#!/usr/bin/python
import HTSeq
import sys
from Bio.Seq import Seq
import os
import argparse
import time
import pickle
import shutil
import multiprocessing
import subprocess

def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True

    return "Not found"


def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        strDNArevC = ''
        for l in strDNA:

        	strDNArevC += basecomplement[l]

        return strDNArevC[::-1]

def translateSeq(DNASeq):
	seq=DNASeq
	tableid=11
	try:
		myseq= Seq(seq)
		protseq=Seq.translate(myseq, table=tableid,cds=True)
	except:
		try:
			seq=reverseComplement(seq)
			myseq= Seq(seq)
			protseq=Seq.translate(myseq, table=tableid,cds=True)
						
		except:
			try:
				seq=seq[::-1]
				myseq= Seq(seq)
				protseq=Seq.translate(myseq, table=tableid,cds=True)
			except:
				try:
					seq=seq[::-1]							
					seq=reverseComplement(seq)
					myseq= Seq(seq)
					protseq=Seq.translate(myseq, table=tableid,cds=True)
				except Exception as e:
					print "translated error"
					print e
					raise
	return protseq
	
def loci_translation (genesList,listOfGenomes2):
	
	gene_fp = open( genesList, 'r')

	genepath=''
	basepath=''
	lGenesFiles = []
	argumentsList = []
	noshortgeneFile=[]
	
	for gene in gene_fp:
		k=0
		gene = gene.rstrip('\n')
		multiple=True
		shortgene= os.path.join(os.path.dirname(gene),"short",os.path.basename(gene))
		shortgene= shortgene.replace(".fasta","_short.fasta")
		if not os.path.isfile(shortgene):
			noshortgeneFile.append(gene)
			break
			
		gene_fp2 = HTSeq.FastaReader(shortgene)
		for allele in gene_fp2: 
			k+=1
			if (len(allele.seq) % 3 != 0):
				multiple=False
				print "allele "+str(k)+" is not multiple of 3: "+str(gene)+" this gene is to be removed"
				break
			else:
				try:
					protseq=translateSeq(allele.seq)
				except:
					multiple=False
					print "allele "+str(k)+" is not translating : "+ str(gene)+" this gene is to be removed"
					break
		if multiple:
			lGenesFiles.append( gene )
			genepath=os.path.dirname(gene)
			basepath=os.path.join(genepath, "temp")
			if not os.path.exists(basepath):
				os.makedirs(basepath)
			filepath=os.path.join(basepath,str(os.path.basename(gene))+"_argList.txt")
			with open(filepath, 'wb') as f:
				var = [gene, listOfGenomes2,genesList]
				pickle.dump(var, f)
			argumentsList.append(filepath)

	

	gene_fp.close()	
	return	genepath,basepath,lGenesFiles,argumentsList,noshortgeneFile

def call_proc(cmd):
	
	p= subprocess.Popen(args=cmd)
	out, err = p.communicate()
	return p

# ================================================ MAIN ================================================ #

def main():
			
	parser = argparse.ArgumentParser(description="This program screens a set of genes in a fasta file.")
	parser.add_argument('-i', nargs='?', type=str, help='List of genome files (list of fasta files)', required=True)
	parser.add_argument('-g', nargs='?', type=str, help='List of genes (fasta)', required=True)
	parser.add_argument('-o', nargs='?', type=str, help="Name of the output files", required=True)
	parser.add_argument('-p', nargs='?', type=str, help="Path to prodigal exec file", required=True)
	parser.add_argument('--cpu', nargs='?', type=int, help="Number of cpus, if over the maximum uses maximum -2", required=True)
	parser.add_argument("-v", "--verbose", help="increase output verbosity",dest='verbose', action="store_true")
	
	
	
	args = parser.parse_args()
	
	genomeFiles = args.i
	genes = args.g
	prodigalPath = args.p
	cpuToUse=args.cpu
	
	if cpuToUse >= multiprocessing.cpu_count()-2:
		cpuToUse=multiprocessing.cpu_count()-2
	
	try:
		verbose=args.verbose
	except:
		verbose=False

	scripts_path=os.path.dirname(os.path.realpath(__file__))
	
	print ("Checking all programs are installed and usable")
	
	print ("Checking Blast... "+str(which('blastp')))
	print ("Checking Prodigal... "+str(which(prodigalPath)))
		
	
	starttime="\nStarting Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y")
	print ("\nStarting Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	
	
	listOfCDSDicts = []
	listOfGenomes = []
	listOfGenomesDict = []

	fp = open(genomeFiles, 'r')

	for genomeFile in fp:

		genomeFile = genomeFile.rstrip('\n')
		genomeFile = genomeFile.rstrip('\r')
		listOfGenomes.append( genomeFile )
		genomeDict = {}

	fp.close()
	
	#check if remnant files from previous run exist, prompt user if exists to know if it's his run and want to continue or start a new one

	with open(genes, 'r') as f:
		first_gene = f.readline()
	
	lGenesFiles = []
	with open(genes, 'r') as gene_fp:
		for gene in gene_fp:
			gene = gene.rstrip('\n')
			gene = gene.rstrip('\r')
			lGenesFiles.append( gene )
	
	first_gene = first_gene.rstrip('\n')
	first_gene = first_gene.rstrip('\r')
	genepath=os.path.dirname(first_gene)
	basepath=os.path.join(genepath, "temp")
	testVar=""
	if os.path.isdir(basepath):
		testVar = raw_input("We found files belonging to a previous run not finished, If they are yours and want to continue were it stopped type Y or yes")
	continueRun=False
	
	if testVar.lower()== "yes" or testVar.lower()== "y":
		continueRun=True
		
	if continueRun:
		argumentsList=[]
		resultsList=[]
		for file in os.listdir(basepath):
			if file.endswith("_argList.txt"):
				argumentsList.append(os.path.join(basepath,str(file)))
			elif file.endswith("_result.txt"):
				resultsList.append((os.path.basename(os.path.join(basepath,str(file)))).replace("_result.txt","_argList.txt"))
		
		#remove unfinished directories
		todelFolders=next(os.walk(genepath))[1]
		for foldertodel in todelFolders:
			if "short" in foldertodel or "temp" in foldertodel:
				pass
			else:
				shutil.rmtree(os.path.join(genepath,str(foldertodel)))
		
	else:
		if os.path.isdir(basepath):
			print ("removing existing files...")
			shutil.rmtree(basepath)
	
		#translate the loci
		genepath,basepath,lGenesFiles,argumentsList, noShort=loci_translation (genes,listOfGenomes)
	
	
	totgenomes= len(listOfGenomes)
		
	joblist =[]
	results = []
	
	if not continueRun:
		if len(argumentsList) ==0:
			shutil.rmtree(basepath)
			raise ValueError('ERROR! AT LEAST ONE GENE FILE MUST CONTAIN AN ALLELE REPRESENTING A CDS')
		
		if len(noShort)>0:
			shutil.rmtree(basepath)
			raise ValueError('ERROR! These loci have no short gene file: '+str(noShort))
			
		# ------------------------------------------------- #
		#           RUN PRODIGAL OVER ALL GENOMES           #
		# ------------------------------------------------- #




		print ("\nStarting Prodigal at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

		
		#Prodigal run on the genomes, one genome per core using n-2 cores (n number of cores)
	
		pool = multiprocessing.Pool(cpuToUse)
		for genome in listOfGenomes:
			pool.apply_async(call_proc, args=([os.path.join(scripts_path,"runProdigal.py"),str(genome),basepath,prodigalPath],))
			
		pool.close()
		pool.join()
	
	
		print "\nChecking all prodigal processes created the necessary files..."
		
		listOfORFCreated=[]
		for orffile in os.listdir(basepath):
			if orffile.endswith("_ORF.txt"):
				listOfORFCreated.append(orffile)
				
		if len(listOfGenomes) >len(listOfORFCreated):
			message="Missing some files from prodigal. "+str((len(listOfGenomes))-(len(listOfORFCreated)))+" missing files out of "+str(len(listOfGenomes))
			shutil.rmtree(basepath)
			raise ValueError(message)
		else:
			print "All prodigal files necessary were created\n"
			
		print ("Finishing Prodigal at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	
		#---CDS to protein---#
		listOFProt=[]
		listAllCDS=[]
		
		i = 0
		j=0
		
		#translate the genome CDSs, load them into dictionaries and fasta files to be used further ahead
		for genomeFile in listOfGenomes:
			listOfCDS={}
			genomeProts=""
			currentCDSDict = {}
			currentGenomeDict = {}
			filepath=os.path.join(basepath,str(os.path.basename(genomeFile))+"_ORF.txt")
			with open(filepath,'rb') as f:
				currentCDSDict = pickle.load(f)
				
			g_fp = HTSeq.FastaReader( genomeFile )
			for contig in g_fp:
				sequence=str(contig.seq)
				genomeDict[ contig.name ] = sequence
			
			currentGenomeDict = genomeDict
			
			i+=1
			for contigTag,value in currentCDSDict.iteritems():

				for protein in value:
					try:
						seq= currentGenomeDict[ contigTag ][ protein[0]:protein[1] ].upper()
						protseq=translateSeq(seq)
						idstr=">"+contigTag+"&protein"+str(j)+"&"+str(protein[0])+"-"+str(protein[1])
						genomeProts+=idstr+"\n"
						listOfCDS[idstr]=seq
						genomeProts+=str(protseq)+"\n"
						
					except Exception as e:
						print str(e)+" "+str(genomeFile)
						pass
					j+=1
			filepath=os.path.join(basepath,str(os.path.basename(genomeFile))+"_ORF_Protein.txt")
			with open(filepath, 'wb') as f:
				var = listOfCDS
				pickle.dump(var, f)

			filepath=os.path.join(basepath,str(os.path.basename(genomeFile))+"_Protein.fasta")
			with open(filepath, 'wb') as f:
				f.write(genomeProts)
				
		
		print ("Starting Genome Blast Db creation at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

		#creation of the Databases for each genome, one genome per core using n-2 cores (n number of cores)
		
		pool = multiprocessing.Pool(cpuToUse)
		for genomeFile in listOfGenomes:
			filepath=os.path.join(basepath,str(os.path.basename(genomeFile))+"_Protein.fasta")
			os.makedirs(os.path.join(basepath,str(os.path.basename(genomeFile)) ))
			
			pool.apply_async(call_proc, args=([os.path.join(scripts_path,"Create_Genome_Blastdb.py"),filepath,os.path.join(basepath,str(os.path.basename(genomeFile)) ),str(os.path.basename(genomeFile))],))
			
			
		pool.close()
		pool.join()
	
	else:
		for argument in argumentsList:
			if argument in resultsList:
				argumentsList.remove(argument)
	
	print
	print ("Starting Allele Calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	

	# Run the allele call, one gene per core using n-2 cores (n number of cores)
	
	
	pool = multiprocessing.Pool(cpuToUse)
	for argList in argumentsList:
		
		pool.apply_async(call_proc, args=([os.path.join(scripts_path,"callAlleles_protein3.py"),str(argList),basepath,str(verbose)],))

	pool.close()
	pool.join()
	

	print ("Finished Allele Calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	output=[]
	for gene in lGenesFiles:
		filepath=os.path.join(basepath, os.path.basename(gene)+"_result.txt")
		with open(filepath,'rb') as f:
			var = pickle.load(f)
			output.append(var)
	
	output2=[]
	for gene in lGenesFiles:
		filepath2=os.path.join(basepath, os.path.basename(gene)+"_result2.txt")
		with open(filepath2,'rb') as f:
			var = pickle.load(f)
			output2.append(var)


	print ("Wrapping up the results")

	#delete all temp files
	shutil.rmtree(basepath)
	
	print "##################################################\n %s genomes used for %s loci" % (len(output[0][0]),len(output) )
	numberexactmatches=0
	for gene in output:
		for gAllele in gene[0]:
			if("EXC:" in gAllele):
				numberexactmatches+=1
					
	
				
	print "\n %s exact matches found out of %s" % (numberexactmatches,(len(output[0][0])*len(output)) )	
	print "\n %s percent of exact matches \n##################################################" % (float((numberexactmatches*100)/(len(output[0][0])*len(output))) )	
		
	print "\nWriting output files\n"
	
	
	#wrapping up the results into a matrix and save in a file
	try:
		phylovout=[]
		phylovout2=[]
		genesnames=[]
		statistics=[]
		
		
		for gene in lGenesFiles:
			
			genename=gene.split("/")
			genename=genename[len(genename)-1]
			genesnames.append(genename)
		for geneOut in output:
			gene=0
			alleleschema=[]
			while gene<len(output[0][0]): 

				genename=(geneOut[1][gene]).split("_")

				if(len(genename)!=1):
					alleleschema.append(genename[1])
				else:
					alleleschema.append(genename[0])
				gene+=1
			phylovout.append(alleleschema)

		for geneOut in output2:
			gene=0
			alleleschema=[]
			while gene<len(output2[0]): 

				genename=(geneOut[gene])
				
				alleleschema.append(genename)
				gene+=1
			phylovout2.append(alleleschema)

		genome=0
		finalphylovinput= "FILE"+ "\t" 
		finalphylovinput2= "FILE"+ "\t" 
		for geneid in genesnames:
			finalphylovinput+= str(geneid)+ "\t"
			finalphylovinput2+= str(geneid)+ "\t"

		
		while genome<len(listOfGenomes):
			currentGenome = os.path.basename(listOfGenomes[genome])
			statsaux=[0]*8 # EXC INF LNF LOT incomplete SAC
			finalphylovinput+= "\n" + currentGenome + "\t"
			for gene in phylovout:
				
				val= str(gene[genome])
				finalphylovinput+= val + "\t"
				if "INF" in val:
					statsaux[1]+=1
				elif "LNF" in val:
					statsaux[2]+=1
				elif "PLOT" in val:
					statsaux[4]+=1
				elif "LOT" in val:
					statsaux[3]+=1 
				elif "NIPL" in val:
					statsaux[5]+=1
				elif "ALM" in val:
					statsaux[6]+=1
				elif "ASM" in val:
					statsaux[7]+=1
				else:
					statsaux[0]+=1
				
			genome+=1
			statistics.append(statsaux)
		
		genome=0	
		while genome<len(listOfGenomes):
			currentGenome = os.path.basename(listOfGenomes[genome])
			finalphylovinput2+= "\n" + currentGenome + "\t"
			for gene in phylovout2:
				
				val= str(gene[genome])
				finalphylovinput2+= val + "\t"
			
			genome+=1
				
			
		
		gOutFile = args.o
		statswrite='Genome\tEXC\tINF\tLNF\tLOT\tPLOT\tNIPL\tALM\tASM'
		i=0
		genome=0
		while genome<len(listOfGenomes):
			currentGenome = os.path.basename(listOfGenomes[genome])
			statsaux=[0]*8 # EXC NA INF LNF LOT incomplete SAC
			statswrite+= "\n" + currentGenome + "\t"
			for k in statistics[i]:
				statswrite+= str(k) + "\t"
			i+=1	
			genome+=1
		
		outputpath=os.path.dirname(gOutFile)
		
		with open(gOutFile, 'wb') as f:
			f.write(finalphylovinput)
			
		print statswrite	
		with open(os.path.join(outputpath,"statistics.txt"), 'wb') as f:
			f.write(str(statswrite))
		
		with open(os.path.join(outputpath,"contigsInfo.txt"), 'wb') as f:
			f.write(str(finalphylovinput2))
			
	except Exception as e:
		print e
		
		exc_type, exc_obj, tb = sys.exc_info()
		f = tb.tb_frame
		lineno = tb.tb_lineno
		print lineno
			
	
	print (starttime)
	print ("Finished Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

if __name__ == "__main__":
    main()
