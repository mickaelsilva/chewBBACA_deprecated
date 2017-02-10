#!/usr/bin/env python
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

def prepGenomes(genomeFile,basepath):
	
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
		currentGenomeDict[ contig.name ] = sequence
	
	j=0
	for contigTag,value in currentCDSDict.iteritems():

		for protein in value:
			try:
				seq= currentGenomeDict[ contigTag ][ protein[0]:protein[1] ].upper()
				protseq=translateSeq(seq)
				j+=1
				idstr=">"+contigTag+"&protein"+str(j)+"&"+str(protein[0])+"-"+str(protein[1])
				genomeProts+=idstr+"\n"
				listOfCDS[idstr]=seq
				genomeProts+=str(protseq)+"\n"
			except Exception as e:
				print (str(e)+" "+str(genomeFile))
				pass

	filepath=os.path.join(basepath,str(os.path.basename(genomeFile))+"_ORF_Protein.txt")
	with open(filepath, 'wb') as f:
		var = listOfCDS
		pickle.dump(var, f)
	listOfCDS=''

	filepath=os.path.join(basepath,str(os.path.basename(genomeFile))+"_Protein.fasta")
	with open(filepath, 'wb') as f:
		f.write(genomeProts)
	genomeProts=''
	var=''
	currentGenomeDict=''
	currentCDSDict=''

	return True

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
					print "translation error"
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
			
	parser = argparse.ArgumentParser(description="This program call alleles for a set of genomes provided a schema")
	parser.add_argument('-i', nargs='?', type=str, help='List of genome files (list of fasta files)', required=True)
	parser.add_argument('-g', nargs='?', type=str, help='List of genes (fasta)', required=True)
	parser.add_argument('-o', nargs='?', type=str, help="Name of the output files", required=True)
	parser.add_argument('--cpu', nargs='?', type=int, help="Number of cpus, if over the maximum uses maximum -2", required=True)
	parser.add_argument("-v", "--verbose", help="increase output verbosity",dest='verbose', action="store_true",default=False)
	parser.add_argument('-b', nargs='?', type=str, help="BLAST full path", required=False,default='blastp')
	parser.add_argument('--bsr', nargs='?', type=float, help="minimum BSR score", required=False,default=0.6)
	parser.add_argument("--so", help="split the output per genome",dest='divideOutput', action="store_true",default=False)
	parser.add_argument('-t', nargs='?', type=str, help="taxon", required=False,default=False)
	
	args = parser.parse_args()
	
	genomeFiles = args.i
	genes = args.g
	cpuToUse=args.cpu
	BSRTresh=args.bsr
	verbose=args.verbose
	BlastpPath=args.b
	divideOutput=args.divideOutput
	gOutFile = args.o
	chosenTaxon=args.t
	
	# avoid user to run the script with all cores available, could impossibilitate any usage when running on a laptop
	if cpuToUse > multiprocessing.cpu_count()-2:
		print "Warning, you are close to use all your cpus, if you are using a laptop you may be uncapable to perform any action"
	
	taxonList={'Campylobacter_Jejuni':'trained_campyJejuni.trn',
				'Acinetobacter_Baumannii':'trained_acinetoBaumannii.trn',
				'Streptococcus_Agalactiae':'trained_strepAgalactiae.trn'
				}	
	if isinstance(chosenTaxon, basestring):
		trainingFolderPAth=os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'TrainingFiles4Prodigal'))
		try:
			chosenTaxon=os.path.join(trainingFolderPAth,taxonList[chosenTaxon])
			
			if os.path.isfile(chosenTaxon):
				print "will use this training file : "+chosenTaxon
			else:
				print "training file don't exist"
				print chosenTaxon
				return "retry"
		except :
			print "Your chosen taxon is not attributed, select one from:"
			for elem in taxonList.keys():
				print elem
			return "retry"

	print BlastpPath
	
	scripts_path=os.path.dirname(os.path.realpath(__file__))
	
	print ("Will use this number of cpus: "+str(cpuToUse))
	print ("Checking all programs are installed")
	
	print ("Checking Blast installed... "+str(which(str(BlastpPath))))
	print ("Checking Prodigal installed... "+str(which('prodigal')))
	
	#check version of Blast
	
	proc = subprocess.Popen([BlastpPath, '-version'], stdout=subprocess.PIPE)
	line = proc.stdout.readline()
	if not "blastp: 2.5." in str(line):
		print "your blast version is " + str(line)
		print "update your blast to 2.5.0, will exit program"
		sys.exit()
	else:
		print "blast version is up to date, the program will continue"
	
	starttime="\nStarting Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y")
	print (starttime)
	
	
	
	listOfCDSDicts = []
	listOfGenomes = []
	listOfGenomesDict = []
	listOfGenomesBasename = []
	
	print "checking if genome files exist.."
	with open(genomeFiles, 'r') as fp:
		for genomeFile in fp:

			genomeFile = genomeFile.rstrip('\n')
			genomeFile = genomeFile.rstrip('\r')
			if os.path.isfile(genomeFile):
				listOfGenomes.append( genomeFile )
				listOfGenomesBasename.append(os.path.basename(genomeFile) )
			else:
				print "File does not exist, will not be used : "+str(genomeFile)
			genomeDict = {}

	
	if len(listOfGenomes) ==0:
			raise ValueError('ERROR! No usable genome files in '+str(genomeFiles))
	
	#check if remnant files from previous run exist, prompt user if exists to know if it's his run and want to continue or start a new one
	
	lGenesFiles = []
	
	print "checking if gene files exist.."
	with open(genes, 'r') as gene_fp:
		for gene in gene_fp:
			gene = gene.rstrip('\n')
			gene = gene.rstrip('\r')
			if os.path.isfile(gene):
				lGenesFiles.append( gene )
			else:
				print "File does not exist, will not be used : "+str(gene)
	
	
	if len(lGenesFiles) ==0:
			raise ValueError('ERROR! No usable gene files in '+str(genes))
	
	#create temp folder inside the folder where the first gene is located
	first_gene = lGenesFiles[0]
	genepath=os.path.dirname(first_gene)
	basepath=os.path.join(genepath, "temp")
	testVar=""
	if os.path.isdir(basepath):
		testVar = raw_input("We found files belonging to a previous run not finished, If they are yours and want to continue were it stopped type Y or yes")
	continueRun=False
	
	if testVar.lower()== "yes" or testVar.lower()== "y":
		continueRun=True
		
	if continueRun:
		print ("You chose to continue the allele call")
		argumentsList=[]
		resultsList=[]
		for file in os.listdir(basepath):
			if file.endswith("_argList.txt"):
				argumentsList.append(os.path.join(basepath,str(file)))
			elif file.endswith("_result.txt"):
				resultsList.append((os.path.join(basepath,str(file))).replace("_result.txt","_argList.txt"))
		
		#remove unfinished directories
		todelFolders=next(os.walk(genepath))[1]
		for foldertodel in todelFolders:
			if "short" in foldertodel or "temp" in foldertodel:
				pass
			else:
				shutil.rmtree(os.path.join(genepath,str(foldertodel)))


	
	if not continueRun:
		
		#user decided not to continue although there was a run that had been stopped before finishing, files will be removed for a new run
		if os.path.isdir(basepath):
			print ("You chose to do a new allele call")
			print ("removing existing files...")
			shutil.rmtree(basepath)
	
		#translate the loci
		genepath,basepath,lGenesFiles,argumentsList, noShort=loci_translation (genes,listOfGenomes)
		
		
		#starting a fresh allele call process, going to run prodigal for all genomes on the list
		
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
			
			pool.apply_async(call_proc, args=([os.path.join(scripts_path,"runProdigal.py"),str(genome),basepath,str(chosenTaxon)],))
			
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
				
		#translate the genome CDSs, load them into dictionaries and fasta files to be used further ahead
		
		print "Translating genomes"
		pool = multiprocessing.Pool(cpuToUse)
		for genomeFile in listOfGenomes:
			
			pool.apply_async(prepGenomes,args=[str(genomeFile),basepath])
		pool.close()
		pool.join()
		
		
		print ("Starting Genome Blast Db creation at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

		#creation of the Databases for each genome, one genome per core using n cores
		
		pool = multiprocessing.Pool(cpuToUse)
		for genomeFile in listOfGenomes:
			filepath=os.path.join(basepath,str(os.path.basename(genomeFile))+"_Protein.fasta")
			os.makedirs(os.path.join(basepath,str(os.path.basename(genomeFile)) ))
			
			pool.apply_async(call_proc, args=([os.path.join(scripts_path,"Create_Genome_Blastdb.py"),filepath,os.path.join(basepath,str(os.path.basename(genomeFile)) ),str(os.path.basename(genomeFile))],))
			
			
		pool.close()
		pool.join()
	
	else:
		
		#user decided to continue the stopped run, finding out which genes were left to run
		argumentsList= list(set(argumentsList) - set(resultsList))
		argumentsList=sorted(argumentsList)

	
	print
	print ("Starting Allele Calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	

	# Run the allele call, one gene per core using n cores
	
	
	pool = multiprocessing.Pool(cpuToUse)
	for argList in argumentsList:
		
		pool.apply_async(call_proc, args=([os.path.join(scripts_path,"callAlleles_protein3.py"),str(argList),basepath,str(BlastpPath),str(verbose),str(BSRTresh)],))

	pool.close()
	pool.join()
	

	print ("Finished Allele Calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	print ("Wrapping up the results")
	
	output=[]
	output2=[]
	for gene in lGenesFiles:
		filepath=os.path.join(basepath, os.path.basename(gene)+"_result.txt")
		filepath2=os.path.join(basepath, os.path.basename(gene)+"_result2.txt")
		with open(filepath,'rb') as f:
			var = pickle.load(f)
			output.append(var)
		with open(filepath2,'rb') as f:
			var = pickle.load(f)
			output2.append(var)
	
		
	#delete all temp files
	shutil.rmtree(basepath)
	
	print "##################################################\n %s genomes used for %s loci" % (len(output[0][1]),len(output) )
	numberexactmatches=0
	for gene in output:
		for gAllele in gene[1]:
			try:
				int(gAllele)
				numberexactmatches+=1
			except:
				pass
					
	
	print "\n used a bsr of : " +str(BSRTresh)	
	print "\n %s exact matches found out of %s" % (numberexactmatches,(len(output[0][1])*len(lGenesFiles)) )	
	print "\n %s percent of exact matches \n##################################################" % (float((numberexactmatches*100)/(len(output[0][1])*len(lGenesFiles))) )	
		
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
		gene=0
		
		containedOutpWrite=''
		for geneOut in output:
			genome=0
			alleleschema=[]
			
			while genome<len(output[0][1]): 
				
				genename=(geneOut[1][genome]).split("_")
				if(len(genename)!=1):
					alleleschema.append(genename[1])
				else:
					alleleschema.append(genename[0])
				genome+=1
			
			phylovout.append(alleleschema)
			if len(geneOut[0])>0:
				containedOutpWrite+= lGenesFiles[gene]+"\n"
			for contained in geneOut[0]:
				containedOutpWrite+=contained[0]+"\t"+contained[1]+"-->"+contained[2]+"\n"
			gene+=1
		
		for geneOut in output2:
			gene=0
			alleleschema=[]
			while gene<len(output2[0]): 

				genename=(geneOut[gene])
				
				alleleschema.append(genename)
				gene+=1
			phylovout2.append(alleleschema)

		genome=0
		
		genesHeader="FILE"+ "\t"+('\t'.join(map(str,genesnames)))
		finalphylovinput=genesHeader
		finalphylovinput2=genesHeader
		
		if divideOutput:
			allelesDict={}
			contigDict={}
			statsDict={}
		
		while genome<len(listOfGenomes):
			auxList=[]
			currentGenome = listOfGenomesBasename[genome]
			statsaux=[0]*7 # EXC INF LNF PLOT NIPL ALM ASM
			finalphylovinput+= "\n" + currentGenome + "\t"
			for gene in phylovout:
				
				val= str(gene[genome])
				auxList.append(val)
				if "INF" in val:
					statsaux[1]+=1
				elif "LNF" in val:
					statsaux[2]+=1
				elif "PLOT" in val:
					statsaux[3]+=1
				elif "NIPL" in val:
					statsaux[4]+=1
				elif "ALM" in val:
					statsaux[5]+=1
				elif "ASM" in val:
					statsaux[6]+=1
				else:
					statsaux[0]+=1
			
			if divideOutput:
				allelesDict[currentGenome]=('\t'.join(map(str,auxList)))	
			finalphylovinput+=('\t'.join(map(str,auxList)))	
			genome+=1
			statistics.append(statsaux)
		
		genome=0	
		while genome<len(listOfGenomes):
			auxList=[]
			currentGenome = listOfGenomesBasename[genome]
			finalphylovinput2+= "\n" + currentGenome + "\t"
			for gene in phylovout2:
				
				val= str(gene[genome])
				auxList.append(val)
			
			if divideOutput:
				contigDict[currentGenome]=('\t'.join(map(str,auxList)))	
			finalphylovinput2+=('\t'.join(map(str,auxList)))	
			genome+=1
				
			
		
		
		statsHeader='Genome\tEXC\tINF\tLNF\tPLOT\tNIPL\tALM\tASM'
		statswrite=statsHeader
		genome=0
		while genome<len(listOfGenomes):
			auxList=[]
			currentGenome = listOfGenomesBasename[genome]
			statsaux=[0]*7 # EXC INF LNF PLOT NIPL ALM ASM
			statswrite+= "\n" + currentGenome + "\t"
			for k in statistics[genome]:
				auxList.append(str(k))
			
			if divideOutput:
				statsDict[currentGenome]=('\t'.join(map(str,auxList)))	
			statswrite+=('\t'.join(map(str,auxList)))	
			genome+=1
		
		if not os.path.exists(gOutFile):
			os.makedirs(gOutFile)
		#outputpath=os.path.dirname(gOutFile)
		outputfolder= os.path.join(gOutFile,"results_"+str(time.strftime("%Y%m%dT%H%M%S")) )
		os.makedirs(outputfolder)
		print statswrite
		print
		print containedOutpWrite
		if not divideOutput:
			with open(os.path.join(outputfolder,"results_alleles.txt"), 'wb') as f:
				f.write(finalphylovinput)
				
				
			with open(os.path.join(outputfolder,"results_statistics.txt"), 'wb') as f:
				f.write(str(statswrite))
				f.write("\n_________________________________________\n")
				f.write(starttime)
				f.write("\nFinished Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
				f.write("\nnumber of genomes: "+str(len(listOfGenomes)))
				f.write("\nnumber of loci: "+str(len(lGenesFiles)))
				f.write ("\nused this number of cpus: "+str(cpuToUse))
				f.write("\nused a bsr of : " +str(BSRTresh))
			
			with open(os.path.join(outputfolder,"results_contigsInfo.txt"), 'wb') as f:
				f.write(str(finalphylovinput2))
			with open(os.path.join(outputfolder,"results_contained.txt"), 'wb') as f:
				f.write(str(containedOutpWrite))
		else:
			for genome in listOfGenomesBasename:
				currentGenome = os.path.splitext(genome)[0]
				perGenomeFolder=os.path.join(outputfolder,currentGenome)
				os.makedirs(perGenomeFolder)
				with open(os.path.join(perGenomeFolder,currentGenome+"_statistics.txt"), 'wb') as f:
					f.write(statsHeader+"\n")
					f.write(genome)
					f.write(statsDict[genome])
				with open(os.path.join(perGenomeFolder,currentGenome+"_contigsInfo.txt"), 'wb') as f:
					f.write(genesHeader+"\n")
					f.write(genome)
					f.write(contigDict[genome])
				with open(os.path.join(perGenomeFolder,currentGenome+"_alleles.txt"), 'wb') as f:
					f.write(genesHeader+"\n")
					f.write(genome)
					f.write(allelesDict[genome])
				
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
