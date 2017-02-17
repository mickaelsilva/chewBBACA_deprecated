#!/usr/bin/env python
import HTSeq
import sys
from Bio.Seq import Seq
from Bio.Blast.Applications import NcbiblastpCommandline
from Counter import Counter 
import os
import argparse
from CommonFastaFunctions import Create_Blastdb
from CommonFastaFunctions import runBlastParser
import time
import pickle
import shutil

def getBlastScoreRatios(genefile,basepath,doAll,verbose,blastPath):
	
	if verbose:
		def verboseprint(*args):
			for arg in args:
			   print arg,
			print
	else:   
		verboseprint = lambda *a: None      # do-nothing function
	
	gene_fp = HTSeq.FastaReader(genefile)
	allelescores=[]
	alleleProt=''
	alleleAllProt=''
	alleleList=[]
	alleleI=0
	alleleIlist=[]
	listAllelesNames=[]
	#calculate bsr for each allele
	for allele in gene_fp: 
		
		#usually first allele name is just >1 and after that it has >gene_id_genome
		aux=(allele.name).split("_")
		if len(aux)<2:
			alleleI=int(aux[0])
		else:
			alleleI=int(aux[-1])

		#try to translate the allele
		alleleIlist.append(alleleI)
		alleleList.append(allele.seq)
		listAllelesNames.append(allele.name)
		translatedSequence,x,y=translateSeq(allele.seq)
		
		if translatedSequence =='':
			print "cannot translate allele on bsr calculation"
			pass
		
		#calculate BSR for the allele	
		else:	
			alleleProt=">"+str(alleleI)+"\n"+str(translatedSequence+"\n")
			alleleAllProt+=">"+str(alleleI)+"\n"+str(translatedSequence+"\n")
			proteinfastaPath=os.path.join(basepath,str(os.path.basename(genefile)+'_protein2.fasta'))
			
			#new db for each allele to blast it against himself
			with open(proteinfastaPath, "wb") as f:
				f.write(alleleProt)
			Gene_Blast_DB_name = Create_Blastdb( proteinfastaPath, 1, True )
			
			#if bsr hasn't been calculated, do the BLAST
			if doAll:
				
				blast_out_file = os.path.join(basepath,'blastdbs/temp.xml')
				verboseprint ("Starting Blast alleles at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
				

				# --- get BLAST score ratio --- #
				cline = NcbiblastpCommandline(cmd=blastPath, query=proteinfastaPath, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)
				allelescore=0
			
				blast_records = runBlastParser(cline,blast_out_file, alleleProt)
			
				verboseprint ("Blasted alleles at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
			
				for blast_record in blast_records:

					for alignment in blast_record.alignments:

						for match in alignment.hsps:
								
							allelescores.append(int(match.score))
							
				geneScorePickle=os.path.abspath(genefile)+'_bsr.txt'
				verboseprint ("________")
				#~ var=[alleleI,allelescores]
				var=dict(zip(alleleIlist, allelescores))
				with open(geneScorePickle,'wb') as f:
					pickle.dump(var, f)			
			
			#bsr had already been calculated, load it to memory
			else:
				geneScorePickle=os.path.abspath(genefile)+'_bsr.txt'
				with open(geneScorePickle,'rb') as f:
					var = pickle.load(f)
					#~ allelescores=var[1]
		
	proteinfastaPath=os.path.join(basepath,str(os.path.basename(genefile)+'_protein.fasta'))
	with open(proteinfastaPath, "wb") as f:
			f.write(alleleAllProt)
			
	#returning all allele BSR scores and list of alleles for this gene
	return var,alleleList,listAllelesNames
	
def reDogetBlastScoreRatios(genefile,basepath,alleleI,allelescores2,newGene_Blast_DB_name,alleleList2,picklepath,verbose,blastPath,listAllelesNames):
	
	if verbose:
		def verboseprint(*args):
			for arg in args:
			   print arg,
			print
	else:   
		verboseprint = lambda *a: None      # do-nothing function
	
	gene_fp = HTSeq.FastaReader(genefile)

	alleleProt=''
	
		
	proteinfastaPath=genefile
	
	verboseprint ("Starting Blast of new alleles to calculate BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	blast_out_file2 = os.path.join(basepath,'blastdbs/temp.xml')

	cline = NcbiblastpCommandline(cmd=blastPath, query=proteinfastaPath, db=newGene_Blast_DB_name, evalue=0.001, out=blast_out_file2, outfmt=5)
	allelescore=0
	blast_records = runBlastParser(cline,blast_out_file2, proteinfastaPath)
	
	verboseprint ("Blasted alleles at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	found =False
	matchscore=0
	for blast_record in blast_records:
		
		for alignment in blast_record.alignments:
			
			
			for match in alignment.hsps:
				matchscore=int(match.score)
				

	allelescores2[alleleI]=matchscore
	with open(picklepath,'wb') as f:
		pickle.dump(allelescores2, f)
	
	return allelescores2,alleleList2,listAllelesNames

def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        strDNArevC = ''
        for l in strDNA:

        	strDNArevC += basecomplement[l]

        return strDNArevC[::-1]

def translateSeq(DNASeq):
	seq=DNASeq
	reversedSeq=False
	tableid=11
	try:
		myseq= Seq(seq)
		protseq=Seq.translate(myseq, table=tableid,cds=True)
	except:
		reversedSeq=True
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
				reversedSeq=False
				try:
					seq=seq[::-1]							
					seq=reverseComplement(seq)
					myseq= Seq(seq)
					protseq=Seq.translate(myseq, table=tableid,cds=True)
				except Exception as e:
					print "translation error"
					print e
					protseq=""
	return protseq,seq,reversedSeq


# ======================================================== #
#            Allele calling and classification             #
# ======================================================== #
def main():
	
	try:
		input_file = sys.argv[1]
		temppath = sys.argv[2]
		blastPath= sys.argv[3]
		verbose= sys.argv[4]
		bsrTresh= sys.argv[5]
		
		if verbose == 'True':
			verbose=True
		else:
			verbose=False
		
	except IndexError:
		print "Error starting the callAlleleles_protein3 script. usage: list_pickle_obj"
	
	bsrTresh=float(bsrTresh)
	
	argumentList=[]
	with open(input_file,'rb') as f:
		argumentList = pickle.load(f)
	
	if verbose:
		def verboseprint(*args):

			for arg in args:
			   print arg,
			print
	else:   
		verboseprint = lambda *a: None
	
	
	geneFile = argumentList[0]

	verboseprint( "Using gene: "+str(geneFile))
	shortgeneFile= os.path.join(os.path.dirname(argumentList[0]),"short",os.path.basename(argumentList[0]))
	shortgeneFile= shortgeneFile.replace(".fasta","_short.fasta")
	genomesList = argumentList[1]
	genesList = argumentList[2]
	
	newListgenes=[]
	with open(genesList,'r') as gene_fp:
		for gene in gene_fp:
			gene = gene.rstrip('\n')
			gene = gene.rstrip('\r')
			newListgenes.append(gene)

	statusbar=float(newListgenes.index(str(geneFile)))/len(newListgenes)
	locusnumber=(newListgenes.index(str(geneFile)))
	totalocusnumber=len(newListgenes)
	basepath=os.path.join(temppath,os.path.splitext(geneFile)[0])
	
	print ("Processing "+os.path.basename(geneFile)+". Start "+time.strftime("%H:%M:%S-%d/%m/%Y")+ " Locus "+str(locusnumber)+" of "+str(totalocusnumber)+". Done "+str(int(statusbar*100))+"%." )
	
	
	if not os.path.exists(basepath):
			os.makedirs(basepath)

	gene_fp = HTSeq.FastaReader(geneFile)
	
	fullAlleleList=[]
	fullAlleleNameList=[]
	alleleI=0
	# get full list of alleles from main gene file and last allele number id
	for allele in gene_fp:
		aux=(allele.name).split("_")
		if len(aux)<2:
			alleleI=int(aux[0])
		else:
			alleleI=int(aux[-1])
		fullAlleleList.append(allele.seq)
		fullAlleleNameList.append(allele.name)
	
	resultsList = []
	i = 0
	perfectMatchIdAllele=[]
	perfectMatchIdAllele2=[]
	allelescores=[]
	listShortAllelesNames=[]
	
	verboseprint ("Getting BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

	geneScorePickle=os.path.abspath(shortgeneFile)+'_bsr.txt'
	
	#check if bsr as arealdy been calculated and recalculate it if necessary

	if os.path.isfile(geneScorePickle) :
		
		allelescores,alleleList,listShortAllelesNames=getBlastScoreRatios(shortgeneFile,basepath,False,verbose,blastPath)
		
	else:	
		allelescores,alleleList,listShortAllelesNames=getBlastScoreRatios(shortgeneFile,basepath,True,verbose,blastPath)
		
			
			
	verboseprint ("Finished BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	verboseprint ("starting allele call blast at: "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	for genomeFile in genomesList:
		verboseprint(genomeFile)
		bestmatch=[0,0,False,'',0] #score, score ratio, perfectmatch, key name of the DNA sequence string, allele ID
		currentGenomeDict={}
		currentCDSDict={}
		
		# load the translated CDS from the genome to a dictionary
		filepath=os.path.join(temppath,str(os.path.basename(genomeFile))+"_ORF_Protein.txt")
		
		with open(filepath,'rb') as f:
			currentCDSDict = pickle.load(f)

		#load the contig info of the genome to a dictionary
		g_fp = HTSeq.FastaReader( genomeFile )
		for contig in g_fp:
			sequence=str(contig.seq)
			currentGenomeDict[ contig.name ] = sequence
		sequence=''
		
		
		#create a dictionary with cds length, for future comparing only the alleles with cds of same size 
		dictCDSLen={}
		for contigTag,value in currentCDSDict.iteritems():
			lengthofCDS=len(value)
			try:
				dictCDSLen[lengthofCDS].append(contigTag)
			except:
				dictCDSLen[lengthofCDS]=[contigTag]

		#check if any CDS is completely equal to an allele without blast -FASTER
		try:
			numberExactAlleles=0
			tempExactResult=[]
			equalmatches=False
			reverse=False
			for alleleAux in fullAlleleList:
				lengthofAllele=len(alleleAux)
				#get list of all CDS with same size as the allele
				if lengthofAllele in dictCDSLen.keys():
					tempList=dictCDSLen[lengthofAllele]
					#only compare the CDS with same size instead of all CDS
					for elem in tempList:
						cds=currentCDSDict[elem]
						try:
							reversedcds=reverseComplement(cds)
						except:
							reversedcds=''
						#~ if reversedcds=='':
							#~ pass
						if str(alleleAux) == str(cds):
							equalmatches=True
						elif str(alleleAux) == str(reversedcds):
							equalmatches=True
							reverse=True
						if equalmatches:
							
							numberExactAlleles+=1
							################################################
							# EXACT MATCH --- MATCH == GENE --- GENE FOUND #
							################################################
							contigname=elem.split("&")
							matchLocation=contigname[2]	
							#starting CDS base need to be +1
							matchLocation=matchLocation.split("-")
							matchLocation=[int(matchLocation[0])+1,matchLocation[1]]
							contigname=(contigname[0]).replace(">","")
							alleleName=''
							alleleMatchid=0
							if reverse:
								alleleName=fullAlleleNameList[fullAlleleList.index(reversedcds)]
								alleleMatchid=int((alleleName.split("_"))[-1])
								
								#~ perfectMatchIdAllele.append(str(alleleMatchid))
								tempExactResult.append(str(alleleMatchid))
								#~ perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"-")
								tempExactResult.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"-")
							else:
								alleleName=fullAlleleNameList[fullAlleleList.index(cds)]
								alleleMatchid=int((alleleName.split("_"))[-1])
								
								#~ perfectMatchIdAllele.append(str(alleleMatchid))
								tempExactResult.append(str(alleleMatchid))
								#~ perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"+")
								tempExactResult.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"+")
								
							#check if atributed allele is contained or contains
							containedInfo=(alleleName.split("_"))[1]
							if "CD" in containedInfo:
								#~ resultsList.append([(os.path.basename(genomeFile)),str(alleleMatchid),containedInfo.rstrip()])
								tempExactResult.append([(os.path.basename(genomeFile)),str(alleleMatchid),containedInfo.rstrip()])
							elif "CS" in containedInfo:
								#~ resultsList.append([(os.path.basename(genomeFile)),str(alleleMatchid),containedInfo.rstrip()])
								tempExactResult.append([(os.path.basename(genomeFile)),str(alleleMatchid),containedInfo.rstrip()])
							else:
								pass
							
							equalmatches=False
							reverse=False	
							#~ resultsList.append('EXC:' + str(alleleMatchid) )
			
			if numberExactAlleles>1:
				perfectMatchIdAllele.append('NIPHEM')
				perfectMatchIdAllele2.append('NIPHEM')
				print os.path.basename(genomeFile)+" has "+str(numberExactAlleles)+" multiple exact match : "+os.path.basename(geneFile)+" MULTIPLE ALLELES as EXACT MATCH"
				raise ValueError("MULTIPLE ALLELES as EXACT MATCH")
				
			elif numberExactAlleles==1:
				perfectMatchIdAllele.append(tempExactResult[0])
				perfectMatchIdAllele2.append(tempExactResult[1])
				try:
					resultsList.append(tempExactResult[2])
				except:
					pass
				raise ValueError("EQUAL")
			
							
		
		except Exception, e:
			#~ exc_type, exc_obj, exc_tb = sys.exc_info()
			#~ fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
			#~ print(exc_tb.tb_lineno)
			#~ print e
			continue
		
		else:	
			verboseprint("Blasting alleles on genome at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
			
			blast_out_file = os.path.join(basepath,"blastdbs/"+os.path.basename(geneFile)+ '_List.xml')

			Gene_Blast_DB_name = os.path.join(temppath,str(os.path.basename(genomeFile))+"/"+str(os.path.basename(genomeFile))+"_db")

			proteinfastaPath=os.path.join(basepath,str(os.path.basename(shortgeneFile)+'_protein.fasta'))
			
			
			#blast the genome CDS against the translated locus
			#cline = NcbiblastpCommandline(query=proteinfastaPath, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5,max_target_seqs=10,max_hsps_per_subject=10)
			#2.2.28 up
			cline = NcbiblastpCommandline(cmd=blastPath, query=proteinfastaPath, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5,max_target_seqs=10,max_hsps=10)

				
			blast_records = runBlastParser(cline, blast_out_file, proteinfastaPath)
			verboseprint("Blasted alleles on genome at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
			
			alleleSizes=[]
			for allele in fullAlleleList:
				alleleSizes.append(len(allele))
			
			biggestSizeAllele=max(alleleSizes)
			
			#get mode allele size
			moda=max(set(alleleSizes), key=alleleSizes.count)
			contador= Counter(alleleSizes).most_common()
			
			#if most common allele size appears 1 time, get first allele size
			if (contador[0])[1] ==1:
				moda= alleleSizes[0]

			try:
				
				# iterate through the blast results
				for blast_record in blast_records:
						
					locationcontigs=[]
					
					for alignment in blast_record.alignments:
						
						# select the best match
						for match in alignment.hsps:
							
							#query id comes with query_id, not name of the allele
							alleleMatchid=int((blast_record.query_id.split("_"))[-1])
							
							#~ scoreRatio=float(match.score)/float(allelescores[int(alleleMatchid)-1])
							#query_id starts with 1
							alleleMatchid2=((listShortAllelesNames[alleleMatchid-1]).split("_"))[-1]
							scoreRatio=float(match.score)/float(allelescores[int(alleleMatchid2)])

							cdsStrName=((alignment.title).split(" "))[1]
							
							DNAstr=str(currentCDSDict[">"+cdsStrName])

							AlleleDNAstr=alleleList[int(alleleMatchid)-1]

			
							if scoreRatio>=bsrTresh:
								locationcontigs.append(cdsStrName)
								
							
							#select the best match from BLAST results
							
							if(scoreRatio == 1 and match.score>bestmatch[0]):
								bestmatch=[match.score,scoreRatio,False,cdsStrName,int(alleleMatchid),match,len(AlleleDNAstr)]

							elif(match.score>bestmatch[0] and scoreRatio>=bsrTresh and scoreRatio>bestmatch[1] and bestmatch[2] is False):
								bestmatch=[match.score,scoreRatio,False,cdsStrName,int(alleleMatchid),match,len(AlleleDNAstr)]
								
											
				verboseprint("Classifying the match at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))		
				

				
				#if no best match was found it's a Locus Not Found
				
				#check for ambiguious bases
				if not bestmatch[0]==0:
				
					alleleStr=currentCDSDict[">"+bestmatch[3]]
					listFoundAmbiguities=[]
					listambiguousBases=['K','M','R','Y','S','W','B','V','H','D','X','N','-','.']
					listFoundAmbiguities=[e for e in listambiguousBases if e in AlleleDNAstr]
				
				if bestmatch[0]==0 or len(listFoundAmbiguities)>0 :
							
							###################
							# LOCUS NOT FOUND #
							###################
					if 	bestmatch[0]==0:		
						#~ resultsList.append('LNF3:-1')
						perfectMatchIdAllele.append('LNF')
						perfectMatchIdAllele2.append('LNF')
						verboseprint( "Locus not found, no matches \n")
					else:

						#~ resultsList.append('LNFN:-1')
						perfectMatchIdAllele.append('LNF')
						perfectMatchIdAllele2.append('LNF')
						verboseprint( "Locus has strange base \n")
				
				#if more than one BSR >0.6 in two different CDSs it's a Non Paralog Locus
				elif len(list(set(locationcontigs)))>1:
					verboseprint("NIPH","")
					#~ resultsList.append('NIPH')            
					perfectMatchIdAllele.append('NIPH')
					perfectMatchIdAllele2.append('NIPH')
					for elem in locationcontigs:
						verboseprint(elem)
					
				

				
				# if match with BSR >0.6 and not equal DNA sequences
				else:
					
					match=bestmatch[5]
					geneLen=bestmatch[6]
					alleleStr=currentCDSDict[">"+bestmatch[3]]
					contigname=bestmatch[3]	
					
					contigname=contigname.split("&")
					matchLocation=contigname[2]	
					matchLocation=matchLocation.split("-")
					matchLocation=[int(matchLocation[0])+1,matchLocation[1]]
					contigname=contigname[0]
					
					seq=currentGenomeDict[ contigname ]
					bestMatchContigLen=len(seq)
					seq=''
					
					protSeq,alleleStr,Reversed=translateSeq(alleleStr)
					
					
					rightmatchContig=bestMatchContigLen-int(matchLocation[1])	
					leftmatchContig=int(matchLocation[0])
					
					
					if Reversed:
						aux=rightmatchContig
						rightmatchContig=leftmatchContig
						leftmatchContig=aux
					
									
					
					# get extra space to the right and left between the allele and match and check if it's still inside the contig
					
					rightmatchAllele=geneLen-((int(match.query_end)+1)*3)	
					leftmatchAllele=((int(match.query_start)-1)*3)
					

							###########################
							# LOCUS ON THE CONTIG TIP #
							###########################
					
					
					#check if contig is smaller than the matched allele
					if leftmatchContig<leftmatchAllele and 	rightmatchContig < rightmatchAllele:
					
						#~ resultsList.append('PLOTSC:-1')
						perfectMatchIdAllele.append('LOTSC')
						perfectMatchIdAllele2.append('LOTSC')
						#~ if not Reversed:
							#~ perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"+")
						#~ else:
							#~ perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"-")
						
						verboseprint(match,contigname,geneFile,leftmatchAllele,rightmatchAllele,"Locus is bigger than the contig \n")
						
					
					elif leftmatchContig<leftmatchAllele:
						
						
						#~ resultsList.append('PLOT3:-1')
						perfectMatchIdAllele.append('PLOT3')
						perfectMatchIdAllele2.append('PLOT3')
						#~ if not Reversed:
							#~ perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"+")
						#~ else:
							#~ perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"-")
						
						verboseprint(match,contigname,geneFile,leftmatchAllele,rightmatchAllele,"Locus is on the 3' tip of the contig \n")
						
					
					elif rightmatchContig < rightmatchAllele:
						
						#~ resultsList.append('PLOT5:-1')
						perfectMatchIdAllele.append('PLOT5')
						perfectMatchIdAllele2.append('PLOT5')
						#~ if not Reversed:
							#~ perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"+")
						#~ else:
							#~ perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"-")
						
						verboseprint(match,contigname,geneFile,leftmatchAllele,rightmatchAllele,"Locus is on the 5' tip of the contig \n")
						
				
								
					elif len(alleleStr) > moda+(moda*0.2) :
						
						verboseprint("Locus is larger than mode", moda, alleleStr)
		
						#~ resultsList.append('ALM')
						perfectMatchIdAllele.append('ALM')
						perfectMatchIdAllele2.append('ALM')
						#~ if not Reversed:
							#~ perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"+")
						#~ else:
							#~ perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"-")
					
					elif len(alleleStr) < moda-(moda*0.2):
						
						verboseprint("Locus is smaller than mode", moda, alleleStr)
			
						#~ resultsList.append('ASM')
						perfectMatchIdAllele.append('ASM')
						perfectMatchIdAllele2.append('ASM')
						#~ if not Reversed:
							#~ perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"+")
						#~ else:
							#~ perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"-")
				
						
					else:
								#######################
								# ADD INFERRED ALLELE #		# a new allele 
								#######################
								
						wasContained=False		
						tagAuxC='S'
						for alleleaux in fullAlleleList:
							
							if alleleStr in alleleaux:
								alleleName=fullAlleleNameList[fullAlleleList.index(alleleaux)]
								alleleMatchid=(alleleName.split("_"))[-1]
								tagAuxC='CD'+alleleMatchid.rstrip()
								resultsList.append([(os.path.basename(genomeFile)),str(alleleI+1),tagAuxC])
								break
							elif alleleaux in alleleStr:
								alleleName=fullAlleleNameList[fullAlleleList.index(alleleaux)]
								alleleMatchid=(alleleName.split("_"))[-1]
								tagAuxC='CS'+alleleMatchid.rstrip()
								resultsList.append([(os.path.basename(genomeFile)),str(alleleI+1),tagAuxC])
								break

						
						if not wasContained:
							tagAux='INF'
				
							perfectMatchIdAllele.append( tagAux +"-"+str(alleleI+1))
							
							if not Reversed:
								perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"+")
							else:
								perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"-")
							
							
							verboseprint( "New allele! Adding allele "+ tagAux + str(alleleI+1) +" to the database\n")
																								
							#~ resultsList.append( tagAux + str(alleleI+1) )
							

														# --- add the new allele to the gene fasta --- #
							
							alleleI+=1
							#appendAllele='>'+str((((os.path.basename(geneFile)).split("."))[0]).replace("_","-"))+"_" + str(alleleI) + "_" + str(os.path.basename(genomeFile)) + '\n'
							appendAllele='>'+str((((os.path.basename(geneFile)).split("."))[0]).replace("_","-"))+"_"+tagAuxC+"_"+ (str(os.path.basename(genomeFile))).replace("_","-") + "_" + str(alleleI) + '\n'
							fG = open( geneFile, 'a' )
							fG.write(appendAllele)
							fG.write( alleleStr + '\n')
							fG.close()
							fullAlleleList.append(alleleStr)
							fullAlleleNameList.append(appendAllele)
							
							if bestmatch[1]>=bsrTresh and bestmatch[1]<bsrTresh+0.1:
								fG = open( shortgeneFile, 'a' )
								fG.write(appendAllele)
								fG.write( alleleStr + '\n')
								fG.close()	
								
								geneTransalatedPath2=os.path.join(basepath,str(os.path.basename(shortgeneFile)+'_protein2.fasta'))
								geneTransalatedPath=os.path.join(basepath,str(os.path.basename(shortgeneFile)+'_protein.fasta'))
								
								fG = open( geneTransalatedPath2, 'w' )
								fG.write('>'+str(alleleI)+'\n'+str(protSeq) + '\n')
								fG.close()
								fG = open( geneTransalatedPath, 'a' )
								fG.write('>'+str(alleleI)+'\n'+str(protSeq) + '\n')
								fG.close()	
								
								match=bestmatch[5]
								
								# --- remake blast DB and recalculate the BSR for the locus --- #
								alleleList.append(alleleStr)
								listShortAllelesNames.append(appendAllele)

								genefile2= geneTransalatedPath2
								Gene_Blast_DB_name2 = Create_Blastdb( genefile2, 1, True )
								verboseprint("Re-calculating BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
								allelescores,alleleList,listShortAllelesNames=reDogetBlastScoreRatios(genefile2,basepath,alleleI,allelescores,Gene_Blast_DB_name2,alleleList,geneScorePickle,verbose,blastPath,listShortAllelesNames)
								verboseprint("Done Re-calculating BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
			
							
			
			except Exception as e:
				print "some error occurred"
				print e
				print 'Error on line {}'.format(sys.exc_info()[-1].tb_lineno)
				perfectMatchIdAllele2.append("ERROR")
				perfectMatchIdAllele.append("ERROR")
				#~ resultsList.append('ERROR')  
		
		
	final =	(resultsList,perfectMatchIdAllele)	
	verboseprint("Finished allele calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	filepath=os.path.join(temppath , os.path.basename(geneFile)+"_result.txt")
	filepath2=os.path.join(temppath , os.path.basename(geneFile)+"_result2.txt")
	with open(filepath, 'wb') as f:
		pickle.dump(final, f)
	with open(filepath2, 'wb') as f:
		pickle.dump(perfectMatchIdAllele2, f)
	shutil.rmtree(basepath)
	return True
	
if __name__ == "__main__":
    main()
