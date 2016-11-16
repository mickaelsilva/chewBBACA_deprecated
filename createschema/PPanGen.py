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
import copy

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
	
	

def checkGeneStrings(genome1,genome2,newName,basepath,cpu,blastp):
	
	pathForTemp=os.path.join(basepath,newName)
	if not os.path.exists(pathForTemp):
		os.makedirs(pathForTemp)
	
	listOfGenomes=[genome1,genome2]


	dictprots={}
	dictprotsLen={}
	dictprotsName={}
	newlistOfCDS={}
	proteinsEqual=0
	smallProteins=0
	protid=0
	genomeProts=""
	genomeProtsTrans=""
	
	try:
		
		for genomeFile in listOfGenomes:
						
			genomename=(os.path.basename(genomeFile)).split(".")
			genomename=genomename[0]
			
			listOfCDS={}
			
			currentCDSDict = {}
			currentGenomeDict = {}
			filepath=os.path.join(basepath,str(os.path.basename(genomeFile))+"_ORF.txt")
			newfilepath=os.path.join(basepath,str(newName))
			
			
				
			g_fp = HTSeq.FastaReader( genomeFile )
			for contig in g_fp:
				sequence=str(contig.seq)
				currentGenomeDict[ contig.name ] = sequence
			
			#after the first iteration, genomes are already defined by their cds and no longer have a cds dictionary pickle file
			try:
				with open(filepath,'rb') as f:
					currentCDSDict = pickle.load(f)
			except:
				for k,v in currentGenomeDict.iteritems():
					currentCDSDict[k]=[[v]]
				#currentCDSDict=copy.deepcopy(currentGenomeDict)
			
			j=0
			
			counter=0
			
			for contigTag,value in currentCDSDict.iteritems():
				
				for protein in value:
					protid+=1
					
					#print protein
					#at first iteration we use the genome file and after a cds only multifasta file
					try:
						seq= currentGenomeDict[ contigTag ][ protein[0]:protein[1] ].upper()
					
					except Exception as e:
						#print contigTag
						#print protein
						seq=str(protein[0])
						
					try:	
						protseq,orderedSeq=translateSeq(seq)
						lengthofProt=len(str(protseq))
					except:
						print str(genome1)+" "+str(genome2)
						pass
					#check if any protein with size on dict
					
					try:
						#print len(protseq)
						if len(str(protseq))<67:
							smallProteins+=1
							pass
						
						elif dictprotsLen[lengthofProt]:
							proteinFound=False
							listproteinsid=dictprotsLen[lengthofProt]
							for elem in listproteinsid:
								
								if protseq==dictprots[elem]:
									proteinFound=True
									proteinsEqual+=1
									break
									
							if not proteinFound:
								dictprotsLen[lengthofProt].append(protid)
								dictprots[protid]=protseq
								newlistOfCDS[protid]=orderedSeq
								try:
									idstr=">"+str(genomename)+"|protein"+str(protid)+"|:"+str(protein[0])+"-"+str(protein[1])+" "
									idstr2=">"+str(genomename)+"|protein"+str(protid)+"|:"+str(protein[0])+"-"+str(protein[1])
								except:
									idstr=">"+str(contigTag)
									idstr2=">"+str((contigTag.split(" "))[0])
								genomeProts+=idstr+"\n"
								genomeProtsTrans+=idstr2+"\n"
								genomeProts+=str(orderedSeq)+"\n"
								genomeProtsTrans+=str(protseq)+"\n"
								dictprotsName[protid]=idstr2
					except Exception as e:
						dictprotsLen[lengthofProt]=[protid]	
						dictprots[protid]=protseq
						newlistOfCDS[protid]=orderedSeq
						try:
							idstr=">"+str(genomename)+"|protein"+str(protid)+"|:"+str(protein[0])+"-"+str(protein[1])+" "
							idstr2=">"+str(genomename)+"|protein"+str(protid)+"|:"+str(protein[0])+"-"+str(protein[1])
							
						except:
							idstr=">"+str(contigTag)
							idstr2=">"+str((contigTag.split(" "))[0])
							
						genomeProts+=idstr+"\n"
						genomeProtsTrans+=idstr2+"\n"
						genomeProts+=str(orderedSeq)+"\n"
						genomeProtsTrans+=str(protseq)+"\n"
						dictprotsName[protid]=idstr2
											
					else :
						pass
			

			listOfCDS=''
			currentGenomeDict=''
			currentCDSDict=''
		
		print "Checked equal proteins for: "+str(genome1)+" "+str(genome2)
		print "Starting with a total of loci: "	+str(protid)
		print "equal proteins : "+str(proteinsEqual)
		print "small proteins : "+str(smallProteins)

		fastaFile=os.path.join(pathForTemp,newName+".fasta")
		with open(fastaFile, 'a') as f:
			f.write(genomeProts)	
			
		newlistOfCDS={}
		genomeProtsTrans=''
		genomeProts=''
		
		
		#check if any protein is substring of a larger, mantaining the larger one
		#ordering the sequences by length, larger ones first
		auxlist=dictprotsLen.keys()
		auxlist=sorted(auxlist, key=int)
		auxlist=auxlist[::-1]
		
		finalProtDict={}
		genomeProtsTrans=''
		auxprotlist=[]
		contained=0
		finalnumber=0
		print "Looking for contained proteins in : "+str(genome1)+" "+str(genome2)
		counter=0
		
		for elem in auxlist:
			counter+=1
			#print "processing "+str(counter)+" out of "+str(len(auxlist))
			#print auxprotlist
			for protid in dictprotsLen[elem]:
				str2=str(dictprots[protid])
				
				try:
					auxprotlist[0]
					for elem2 in auxprotlist:
						isContained=False
						if len(elem2)<len(str2):
							break
						else:
							if str2 in elem2:
								isContained=True
								contained+=1
								break
					if not isContained:
						str1=dictprotsName[protid]
						genomeProtsTrans+=str1+"\n"+str2+"\n"
						finalnumber+=1
						auxprotlist.append(str2)
				except Exception as e:
					#print "first element added"
					str1=dictprotsName[protid]
					genomeProtsTrans+=str1+"\n"+str2+"\n"
					finalnumber+=1
					auxprotlist.append(str2)
		
		print "number of contained proteins : " +str(contained)
		print "total of loci to blast : "+str(finalnumber)
		
		proteinFile=os.path.join(pathForTemp,newName+"_proteins.fasta")
		
		
		with open(proteinFile, 'a') as f:
			f.write(genomeProtsTrans)
		

		#sadasd
		# run createschema
		print "running blast will use this number of cpu: "+str(cpu)
		proc = subprocess.Popen(['/home/msilva/chewBBACA/createschema/CreateSchema.py', '-i', fastaFile,'-l', "200",'--cpu', str(cpu),'-p', proteinFile, '-o', fastaFile, "-b", blastp],stdout=subprocess.PIPE)
		p_status = proc.wait()
		print "finished blast"
		
		os.remove(proteinFile)
		
		
	except Exception as e:
		print e
		return e	
	
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
	
	return protseq,seq
	


def call_proc(cmd):
	
	p= subprocess.Popen(args=cmd)
	out, err = p.communicate()
	return p

# ================================================ MAIN ================================================ #

def main():
			
	parser = argparse.ArgumentParser(description="This program call alleles for a set of genomes provided a schema")
	parser.add_argument('-i', nargs='?', type=str, help='List of genome files (list of fasta files)', required=True)
	parser.add_argument('-o', nargs='?', type=str, help="Name of the output files", required=True)
	parser.add_argument('--cpu', nargs='?', type=int, help="Number of cpus, if over the maximum uses maximum -2", required=True)
	parser.add_argument('-b', nargs='?', type=str, help="BLAST full path", required=False,default='blastp')
	
	
	args = parser.parse_args()
	
	genomeFiles = args.i
	cpuToUse=args.cpu
	outputFile=args.o
	BlastpPath=args.b
	
	# avoid user to run the script with all cores available, could impossibilitate any usage when running on a laptop
	if cpuToUse > multiprocessing.cpu_count()-2:
		print "Warning, you are close to use all your cpus, if you are using a laptop you may be uncapable to perform any action"
		
	
		
	scripts_path=os.path.dirname(os.path.realpath(__file__))
	
	print ("Will use this number of cpus: "+str(cpuToUse))
	print ("Checking all programs are installed")
	
	print ("Checking Prodigal installed... "+str(which('prodigal')))
	
	
	starttime="\nStarting Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y")
	print (starttime)
	
	
	
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
	
	basepath=os.path.join((os.path.dirname(outputFile)), "temp")
	if not os.path.exists(basepath):
				os.makedirs(basepath)
		
	# ------------------------------------------------- #
	#           RUN PRODIGAL OVER ALL GENOMES           #
	# ------------------------------------------------- #


	print ("\nStarting Prodigal at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

	
	#Prodigal run on the genomes, one genome per core using n-2 cores (n number of cores)

	pool = multiprocessing.Pool(cpuToUse)
	for genome in listOfGenomes:
		pool.apply_async(call_proc, args=([os.path.join(scripts_path,"runProdigal.py"),str(genome),basepath],))
		
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
	#listpairs=[]
	toprocess=[]
	processed=[]
	
	pairID=0
	#while len(processed)<len(toprocess):
	while len(listOfGenomes)>0:
		
		pair=[]
		dictPairs={}
		remaining=[]
		
		for genomeFile in listOfGenomes:
			#toprocess.append(listOfGenomes)
			if len(pair)<2:
				pair.append(genomeFile)
			else:
				dictPairs[pairID]=pair
				pairID+=1
				pair=[]
				pair.append(genomeFile)
		
		
		#if total unpair, keep the remainig
		listOfGenomes=[]
		
		if len(pair)==2:
			dictPairs[pairID]=pair
			pairID+=1
		
		elif len(pair)>0:
			listOfGenomes.append(pair[0])
		
		numberOfPairs=len(dictPairs.items())
		extraCpu=0
		if numberOfPairs>=cpuToUse:
			pool = multiprocessing.Pool(cpuToUse)
		else:
			pool = multiprocessing.Pool(numberOfPairs)	
			extraCpu=cpuToUse-numberOfPairs
		
		#print dictPairs	
		for item in dictPairs.items():
			k=item[0]
			v=item[1]
			
			newgGenome="protogenome"+str(k)
			pathFornewgGenome=os.path.join(basepath,newgGenome,newgGenome+".fasta")
			listOfGenomes.append(pathFornewgGenome)
			extraCpuPerProcess=extraCpu/numberOfPairs
			print "running analysis for pair : "+str(v[0])+" "+str(v[1])
			pool.apply_async(checkGeneStrings,args=[v[0],v[1],newgGenome,basepath,extraCpuPerProcess+1,BlastpPath])
			
		
		pool.close()
		pool.join()	
		print "finished running pair analysis"
		
		if len(listOfGenomes)==1:
			print "Creating the schema"
			lastFile=listOfGenomes.pop()
			proc = subprocess.Popen(['/home/msilva/chewBBACA/createschema/CreateSchema.py', '-i', lastFile,'-l', "200",'--cpu', str(cpuToUse),"-b",BlastpPath],stdout=subprocess.PIPE)
			p_status = proc.wait()
			print "Schema Created sucessfully"

	
	#shutil.rmtree(basepath)		
	
	print (starttime)
	print ("Finished Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

if __name__ == "__main__":
    main()
