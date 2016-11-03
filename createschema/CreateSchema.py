#!/usr/bin/python
import HTSeq
import argparse
import os.path
from Bio.Seq import Seq
from CommonFastaFunctions import Create_Blastdb
from CommonFastaFunctions import runBlast
from CommonFastaFunctions import runBlastParser
from Bio.Blast.Applications import NcbiblastpCommandline
import collections
import shutil
from  init_schema_4_bbaca import get_Short
import time

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

def translateSeq(DNASeq,genename):
	seq=DNASeq
	reversedSeq=False
	try:
		myseq= Seq(seq)
		protseq=Seq.translate(myseq, table=11,cds=True)
	except:
		reversedSeq=True
		try:
			seq=reverseComplement(seq)
			myseq= Seq(seq)
			protseq=Seq.translate(myseq, table=11,cds=True)
						
		except:
			try:
				seq=seq[::-1]
				myseq= Seq(seq)
				protseq=Seq.translate(myseq, table=11,cds=True)
			except:
				reversedSeq=False
				try:
					seq=seq[::-1]							
					seq=reverseComplement(seq)
					myseq= Seq(seq)
					protseq=Seq.translate(myseq, table=11,cds=True)
				except Exception as e:
					print "translation error :",e," on locus "+str(genename)

					protseq=""
	return protseq,seq,reversedSeq

def main():

	parser = argparse.ArgumentParser(description="Given an ffn file, recovers the genes that are not paralogs and have a size bigger than the g parameter provided")
	parser.add_argument('-i', nargs='?', type=str, help='ffn file', required=True)
	parser.add_argument('-l', nargs='?', type=int, help='int minimum length', required=True)
	parser.add_argument('--cpu', nargs='?', type=int, help="Number of cpus, if over the maximum uses maximum -2", required=False)
	
	args = parser.parse_args()
	genes = args.i
	sizethresh = args.l
	cpuToUse = args.cpu
	passSteps = False
	
	starttime="\nStarting Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y")
	print ("\nStarting Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	print ("Checking Blast installed... "+str(which('blastp')))
	
	#translate to protein and create new file
	abspath=os.path.abspath(genes)
	filename=os.path.basename(genes)
	abspath=abspath.replace(filename,'')
	proteinfile=os.path.join(abspath,'proteins.fasta') 
	
	geneDict = {}
	protDict={}
	orderedprotDict=collections.OrderedDict()
	alreadyIn=[]
	totalgenes=0
	repeatedgenes=0
	smallgenes=0
	nottranslatable=0
	
	print "Checking translatability of the loci:\n"
	
	if not passSteps:
		#print "not passing steps"
		with open(proteinfile, "wb") as f:
			g_fp = HTSeq.FastaReader( genes )
			
			for gene in g_fp:
				dnaseq=	str(gene.seq)
				protseq,x,y=translateSeq(dnaseq,gene.name)
				totalgenes+=1
				if len(protseq)>1:
					
					if str(protseq) in alreadyIn:
						repeatedgenes+=1
					
					elif len(str(protseq))<67:
						smallgenes+=1
						
					else:	
						alreadyIn.append(str(protseq))
						protname=">"+str(gene.name)+"\n"
														
						f.write(protname+str(protseq)+"\n")
						protDict[protname] = str(protseq)
						geneDict[str(gene.name)] = gene.seq
				else:
					nottranslatable+=1
					continue
			
			print (str(nottranslatable)+" not translatable out of "+ str(totalgenes))
			
			print
			print "Checking if repeated protein sequences:\n"
			
			orderedprotList=[]
			orderedprotList=sorted(protDict.items(), key=lambda x: len(x[1]), reverse=True)
			
			
			i=0
			while i < len(orderedprotList):
				elem=orderedprotList[i]
				orderedprotDict[elem[0]] = elem[1]
				i+=1
		

		print (str(repeatedgenes) + " repeated loci out of "+ str(totalgenes))
		print (str(smallgenes) + " loci out of "+ str(totalgenes)+ " smaller than "+str(sizethresh)+"bp")
		print "\nprotein file created\n"
				
		# first step -  remove genes contained in other genes or 100% equal genes
		
		# list of results - the output of the function
		resultsList = []
		
		auxDict={}
		g_fp = HTSeq.FastaReader( proteinfile )
		g=0
		j=0
		
		print "Checking if protein sequences are contained in others..."
		
		# for each gene from all the annotated genes - starting with an empty dictionary, only add a new gene if the "to be added gene" is not contained or equal to a gene already added to the dictionary
		auxprot=[]

		for elem in orderedprotDict.items():

			contained=False
			
			prot=str(elem[1])
			if any(prot in x for x in auxprot):
				g+=1
				contained=True
			
			else:
				auxDict[elem[1]] = elem[0]
				auxprot.append(str(elem[1]))
			
				
			j+=1
		print "%s loci are contained in other genes\n" %  (g)
		
		#overwrite the original file, obtaining a new file with unique genes
		
		with open(proteinfile, "wb") as f:
			allsequences=''
			for k,v in auxDict.iteritems():
				allsequences+=v+k+"\n"
			f.write(allsequences)
	
	else:
		
		totalgenes=0
		smallgenes=0
		g_fp = HTSeq.FastaReader( genes )
		
		for gene in g_fp:
			dnaseq=	str(gene.seq)
			protseq,x,y=translateSeq(dnaseq)
			if len(protseq)>1:
				
				if str(protseq) in alreadyIn:
					repeatedgenes+=1
				
				elif len(str(protseq))<67:
					smallgenes+=1
					
				else:	
					alreadyIn.append(str(protseq))
					protname=">"+str(gene.name)+"\n"							
					protDict[protname] = str(protseq)
					geneDict[str(gene.name)] = gene.seq
			else:

				print gene.name
	
	
	print "Blasting the total of "+ str(len(auxDict.keys())) + " loci"
	
	geneFile = os.path.abspath( proteinfile )
	Gene_Blast_DB_name = Create_Blastdb( geneFile, 1, True )
	
	geneF = os.path.splitext( geneFile )[0]
	blast_out_file = geneF + '.xml'
					# ------------------------------ RUNNING BLAST ------------------------------ #
	if cpuToUse:
		cline = NcbiblastpCommandline(query=geneFile, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5, num_threads=int(cpuToUse))
	else:
		cline = NcbiblastpCommandline(query=geneFile, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)
	blast_records = runBlastParser(cline, blast_out_file, geneFile)
	toRemove=[]
	genesToKeep=[]
	log=["removed\tcause\texplanation"]
	for blast_record in blast_records:
		
		allelename=blast_record.query
		allelename=allelename.split(" ")
		allelename=allelename[0]
		alleleLength=len(geneDict[allelename])

		try:
			
			#if gene A is not on the toRemove list yet, add to genesToKeep list
			
			if str(blast_record.query) not in toRemove:
				genesToKeep.append(blast_record.query)
				
				i=0
				#if first alignement is not against self, gene B is bigger than gene A and very simillar - remove gene A from genesToKeep and add gene B instead
				if  not str(blast_record.query) == str((blast_record.alignments[0]).hit_def):
					genesToKeep.remove(str(blast_record.query))
					toRemove.append(str(blast_record.query))
					log.append(str(blast_record.query)+"\t"+str((blast_record.alignments[0]).hit_def)+"\t"+"2 is first best match")
					
					#if gene B is not on the toRemove list, add to genesToKeep list
					if str((blast_record.alignments[0]).hit_def) not in toRemove:
						genesToKeep.append(str((blast_record.alignments[0]).hit_def))

					raise
				
				selfblastscore=(((blast_record.alignments[0]).hsps)[0]).score


				
				
				while i<len(blast_record.alignments):
					align=blast_record.alignments[i]
					
					match=(align.hsps)[0]
					scoreRatio=float(match.score)/float(selfblastscore)
					
					alleleLength2=len(geneDict[str(align.hit_def)])
					
					#if good match and gene B not in toremove list
					if(scoreRatio>0.6 and not str(align.hit_def) == str(blast_record.query) and str(align.hit_def) not in toRemove):
						
						#if gene B is bigger than gene A, keep bigger gene B
						if alleleLength2>alleleLength :
							genesToKeep.append(str(align.hit_def))
							genesToKeep.remove(str(blast_record.query))
							toRemove.append(str(blast_record.query))
							log.append(str(blast_record.query)+"\t"+str(align.hit_def)+"\t"+"2 is bigger and bsr >0.6")
							
							raise
						#else add gene B to toremove list
						elif str(align.hit_def) in genesToKeep:
							genesToKeep.remove(str(align.hit_def))
							toRemove.append(str(align.hit_def))
							log.append(str(align.hit_def)+"\t"+str(blast_record.query)+"\t"+"2 is bigger and bsr >0.6")
							
					i+=1
			
			#else gene A is on toRemove list, add all similar genes (not in genesToKeep) list to the toRemove list
			else:		
						
				i=0
				selfblastscore=0
				for align in blast_record.alignments:
					if not (str(align.hit_def) == str(blast_record.query)):
						selfblastscore=((align.hsps)[0]).score
						#print "gene "+str(align.hit_def)+" is larger than gene "+str(blast_record.query)
						raise
				
				while i<len(blast_record.alignments):
					align=blast_record.alignments[i]
					match=(align.hsps)[0]
					scoreRatio=float(match.score)/float(selfblastscore)
					
					if align.hit_def not in genesToKeep and not str(align.hit_def) == str(blast_record.query) and scoreRatio>0.6 :
						toRemove.append(align.hit_def)
						#log.append(str(align.hit_def)+"\t"+str(blast_record.query)+"\t"+"2 was on the removed list and bsr >0.6")
							
					else:
						pass

					i+=1
			

		except Exception as e:
			#print e
			pass
	with open("schemacreation.log", "wb") as f:
		for elem in log:
			
			f.write(str(elem)+"\n")
	

	
	genesToKeep=list(set(genesToKeep))
	toRemove=list(set(toRemove))
	s = set(toRemove)
	notcommonToKeep= [x for x in genesToKeep if x not in s]

	pathfiles=os.path.dirname(geneFile)
	pathfiles=pathfiles+"/"
	listfiles=[]

	g_fp = HTSeq.FastaReader( genes )
	removedparalogs=0
	removedsize=0
	totalgenes=0
	rest=0
	concatenatedFile=''
	schema_folder_path=os.path.join(pathfiles,'schema_seed')
	
	if not os.path.exists(schema_folder_path):
		os.makedirs(schema_folder_path)
	
	for contig in g_fp:
		totalgenes+=1
		name = contig.name+" "+contig.descr
		name2= contig.name
		
		
		if name2 not in toRemove and name2 in genesToKeep:
			if int(len(contig.seq))>sizethresh:
				namefile=contig.name
				namefile=namefile.replace("|","_")
				newFile=os.path.join(schema_folder_path,namefile+".fasta")
				listfiles.append(newFile)
				with open(newFile, "wb") as f:
					f.write(">1\n"+contig.seq+"\n")
				rest+=1	
				concatenatedFile+=">"+namefile+"\n"+contig.seq+"\n"
			else:
				removedsize+=1
		else:

			removedparalogs+=1
		
	print "\nRemoved %s with a high similarity (BSR>0.6)" % str(removedparalogs)
	print "Total of %s loci that constitute the schema" % str(rest)
	
	

	shutil.rmtree(os.path.join(pathfiles,'blastdbs'))
	os.remove(proteinfile)
	os.remove(blast_out_file)
	
	#create short folder
	get_Short(listfiles)
	print (starttime)
	print ("Finished Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
if __name__ == "__main__":
	main()
