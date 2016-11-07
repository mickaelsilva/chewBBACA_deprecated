#!/usr/bin/python
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
import warnings

def getBlastScoreRatios(genefile,basepath,doAll,verbose):
	
	if verbose:
		def verboseprint(*args):
			# Print each argument separately so caller doesn't need to
			# stuff everything to be printed into a single string
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
	
	for allele in gene_fp: #new db for each allele to blast it against himself
		
		aux=(allele.name).split("_")
		if len(aux)<2:
			alleleI=int(aux[0])
		else:
			alleleI=int(aux[1])

		genome=-1
		alleleList.append(allele.seq)
		translatedSequence,x,y=translateSeq(allele.seq)
		
		if translatedSequence =='':
			pass
			
		else:	
			alleleProt=">"+str(alleleI)+"\n"+str(translatedSequence+"\n")
			alleleAllProt+=">"+str(alleleI)+"\n"+str(translatedSequence+"\n")
			proteinfastaPath=os.path.join(basepath,str(os.path.basename(genefile)+'_protein2.fasta'))
			
			with open(proteinfastaPath, "wb") as f:
				f.write(alleleProt)
			Gene_Blast_DB_name = Create_Blastdb( proteinfastaPath, 1, True )
			if doAll:
				
				blast_out_file = os.path.join(basepath,'blastdbs/temp.xml')
				verboseprint ("Starting Blast alleles at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
				

				# --- get BLAST score ratio --- #
				cline = NcbiblastpCommandline(query=proteinfastaPath, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)
				allelescore=0
			
				blast_records = runBlastParser(cline,blast_out_file, alleleProt)
			
				verboseprint ("Blasted alleles at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
			
				for blast_record in blast_records:

					for alignment in blast_record.alignments:

						for match in alignment.hsps:
								
							allelescores.append(int(match.score))
							
				geneScorePickle=os.path.abspath(genefile)+'_bsr.txt'
				verboseprint ("________")
				var=[alleleI,allelescores]
				with open(geneScorePickle,'wb') as f:
					pickle.dump(var, f)			
			
			else:
				geneScorePickle=os.path.abspath(genefile)+'_bsr.txt'
				with open(geneScorePickle,'rb') as f:
					var = pickle.load(f)
					allelescores=var[1]
				
	proteinfastaPath=os.path.join(basepath,str(os.path.basename(genefile)+'_protein.fasta'))
	with open(proteinfastaPath, "wb") as f:
			f.write(alleleAllProt)
			
	
	return allelescores,alleleList
	
def reDogetBlastScoreRatios(genefile,basepath,alleleI,allelescores2,newGene_Blast_DB_name,alleleList2,picklepath,verbose):
	
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
	
	verboseprint ("Re-starting Blast alleles at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	blast_out_file2 = os.path.join(basepath,'blastdbs/temp.xml')

	cline = NcbiblastpCommandline(query=proteinfastaPath, db=newGene_Blast_DB_name, evalue=0.001, out=blast_out_file2, outfmt=5)
	allelescore=0
	blast_records = runBlastParser(cline,blast_out_file2, proteinfastaPath)
	
	verboseprint ("Blasted alleles at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	found =False
	for blast_record in blast_records:
		
		for alignment in blast_record.alignments:
			
			
			for match in alignment.hsps:
				allelescores2.append(int(match.score))
				

	var=[alleleI,allelescores2]
	with open(picklepath,'wb') as f:
		#currentCDSDict = pickle.dump(var, f)
		pickle.dump(var, f)
	
	return allelescores2,alleleList2

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
		verbose= sys.argv[3]
		
		if verbose == 'True':
			verbose=True
		else:
			verbose=False
		
	except IndexError:
		print "Error starting the callAlleleles_protein3 script. usage: list_pickle_obj"
	
	
	
	argumentList=[]
	with open(input_file,'rb') as f:
		argumentList = pickle.load(f)
	warnings.filterwarnings('ignore')
	
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
	alleleI=0
	# get full list of alleles from main gene file and last allele number id
	for allele in gene_fp:
		aux=(allele.name).split("_")
		if len(aux)<2:
			alleleI=int(aux[0])
		else:
			alleleI=int(aux[1])
		fullAlleleList.append(allele.seq)
	
	resultsList = []
	i = 0
	perfectMatchIdAllele=[]
	perfectMatchIdAllele2=[]
	allelescores=[]
	
	verboseprint ("Getting BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

	geneScorePickle=os.path.abspath(shortgeneFile)+'_bsr.txt'
	
	#check if bsr as arealdy been calculated and recalculate it

	if os.path.isfile(geneScorePickle) :
		
		allelescores,alleleList=getBlastScoreRatios(shortgeneFile,basepath,False,verbose)
		
	else:	
		allelescores,alleleList=getBlastScoreRatios(shortgeneFile,basepath,True,verbose)
		
			
			
	verboseprint ("Finished BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	genome=-1	
	
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

		genome+=1
		#listOfCDS=currentCDSDict
		
		#check if any CDS is completely equal to an allele without blast -FASTER
		try:
			for alleleAux in fullAlleleList:
				for k,cds in currentCDSDict.items():
					if str(alleleAux) == str(cds):
						################################################
						# EXACT MATCH --- MATCH == GENE --- GENE FOUND #
						################################################
						contigname=k.split("&")
						matchLocation=contigname[2]	
						
						contigname=(contigname[0]).replace(">","")
						alleleMatchid=(int(fullAlleleList.index(cds)))+1
						
						perfectMatchIdAllele.append(str(alleleMatchid))
						perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation)+"&"+"+")

						resultsList.append('EXC:' + str(alleleMatchid) )
						
						raise ValueError("EQUAL")
		except Exception, e:

			continue
		
		else:	
			verboseprint("Blasting alleles on genome at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
			
			blast_out_file = os.path.join(basepath,"blastdbs/"+os.path.basename(geneFile)+ '_List.xml')

			Gene_Blast_DB_name = os.path.join(temppath,str(os.path.basename(genomeFile))+"/"+str(os.path.basename(genomeFile))+"_db")

			proteinfastaPath=os.path.join(basepath,str(os.path.basename(shortgeneFile)+'_protein.fasta'))
			
			
			#blast the genome CDS against the translated locus
			#	cline = NcbiblastpCommandline(query=proteinfastaPath, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5,max_target_seqs=10,max_hsps_per_subject=10)
			#2.2.28 up
			cline = NcbiblastpCommandline(query=proteinfastaPath, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5,max_target_seqs=10,max_hsps=10)

				
			blast_records = runBlastParser(cline, blast_out_file, proteinfastaPath)
			verboseprint("Blasted alleles on genome at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
			
			alleleSizes=[]
			for allele in fullAlleleList:
				alleleSizes.append(len(allele))
			
			biggestSizeAllele=0
			
			moda=max(set(alleleSizes), key=alleleSizes.count)
			contador= Counter(alleleSizes).most_common()
			
			if (contador[0])[1] ==1:
				moda= alleleSizes[0]

			try:
				
				# iterate through the blast results
				for blast_record in blast_records:
						
					locationcontigs=[]
					
					for alignment in blast_record.alignments:
						
						# select the best match
						for match in alignment.hsps:
							
							alleleMatchid=str(blast_record.query_id).split("_")[1]
							
							scoreRatio=float(match.score)/float(allelescores[int(alleleMatchid)-1])

							cdsStrName=((alignment.title).split(" "))[1]
							
							DNAstr=str(currentCDSDict[">"+cdsStrName])

							AlleleDNAstr=alleleList[int(alleleMatchid)-1]
							if len(AlleleDNAstr)>biggestSizeAllele:
								biggestSizeAllele=len(AlleleDNAstr)
								
							compare=False
							
							#check if DNA sequence is already defined as an allele
							
							
							
							if scoreRatio>0.6:
								locationcontigs.append(cdsStrName)
								
								if DNAstr in fullAlleleList:
									compare=True
									alleleMatchid=(int(fullAlleleList.index(DNAstr)))+1
								
								else:
									try:
										DNAstr=reverseComplement(DNAstr)
										if DNAstr in fullAlleleList:
											compare=True
											alleleMatchid=(int(fullAlleleList.index(DNAstr)))+1
									except:
										#some ambiguous base found cant reversecomplement it
										continue
								
							
							notCDS=False
							try:
								protseq=translateSeq(DNAstr)
							except:
								notCDS=True
							if notCDS:
								pass
							
							elif(compare is True):
								bestmatch=[match.score,1,True,cdsStrName,int(alleleMatchid),match,len(AlleleDNAstr),AlleleDNAstr]
							
							elif(scoreRatio == 1 and bestmatch[2] is False and compare is True):
								bestmatch=[match.score,scoreRatio,True,cdsStrName,int(alleleMatchid),match,len(AlleleDNAstr)]

							elif(scoreRatio == 1 and match.score>bestmatch[0] and compare is True and bestmatch[2] is False):
								bestmatch=[match.score,scoreRatio,True,cdsStrName,int(alleleMatchid),match,len(AlleleDNAstr)]

							elif(scoreRatio == 1 and bestmatch[2] is False and compare is False):
								bestmatch=[match.score,scoreRatio,False,cdsStrName,int(alleleMatchid),match,len(AlleleDNAstr)]
							
							elif(scoreRatio == 1 and match.score>bestmatch[0] and compare is False):
								bestmatch=[match.score,scoreRatio,False,cdsStrName,int(alleleMatchid),match,len(AlleleDNAstr)]

							elif(match.score>bestmatch[0] and scoreRatio>0.6 and scoreRatio>bestmatch[1] and bestmatch[2] is False):
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
						resultsList.append('LNF3:-1')
						perfectMatchIdAllele.append('LNF')
						perfectMatchIdAllele2.append('LNF')
						verboseprint( "Locus not found, no matches \n")
					else:

						resultsList.append('LNFN:-1')
						perfectMatchIdAllele.append('LNF')
						perfectMatchIdAllele2.append('LNF')
						verboseprint( "Locus has strange base \n")
				
				#if more than one BSR >0.6 in two different CDSs it's a Non Paralog Locus
				elif len(list(set(locationcontigs)))>1:
					verboseprint("NIPL","")
					resultsList.append('NIPL')            
					perfectMatchIdAllele.append('NIPL')
					perfectMatchIdAllele2.append('NIPL')
					for elem in locationcontigs:
						verboseprint(elem)
					
				
				#in case the DNA match sequence equal to the DNA sequence of the comparing allele
				elif bestmatch[2] is True:
					contigname=bestmatch[3]	
					
					contigname=contigname.split("&")
					matchLocation=contigname[2]	
					contigname=contigname[0]	
					protSeq,alleleStr,Reversed=translateSeq(alleleStr)
					

					#check for possible locus on tip
					match=bestmatch[5]
					matchLocation2=matchLocation.split("-")			
					seq=currentGenomeDict[ contigname ]
					bestMatchContigLen=len(seq)
					seq=''
					
					rightmatchContig=bestMatchContigLen-int(matchLocation2[1])	
					leftmatchContig=int(matchLocation2[0])
					
					if Reversed:
						aux=rightmatchContig
						rightmatchContig=leftmatchContig
						leftmatchContig=aux
					
					
					
					
					
					# get extra space to the right and left between the allele and match
					
					possibleExtra=int(moda)-((int(match.query_end)*3)-(int(match.query_start)*3))
					
					if possibleExtra<0:
						perfectMatchIdAllele.append(str(bestmatch[4]))
						if not Reversed:
							perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation)+"&"+"+")
						else:
							perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation)+"&"+"-")
						resultsList.append('EXC:' + str(bestmatch[4]) )
					
					else:	
						rightmatchAllele=possibleExtra
						leftmatchAllele=possibleExtra
						
						if leftmatchContig<leftmatchAllele and 	rightmatchContig < rightmatchAllele:
					
							resultsList.append('PLOTSC:-1')
							perfectMatchIdAllele.append('PLOTSC')
							perfectMatchIdAllele2.append('PLOTSC')

							verboseprint(match, "contig extras (l,r)",leftmatchContig,rightmatchContig,"allele extras (l,r)",leftmatchAllele,rightmatchAllele,"Locus is possibly bigger than the contig \n")
							
						
						elif leftmatchContig<leftmatchAllele:
							
							
							resultsList.append('PLOT3:-1')
							perfectMatchIdAllele.append('PLOT3')
							perfectMatchIdAllele2.append('PLOT3')
							
							verboseprint(match, "contig extras (l,r)",leftmatchContig,rightmatchContig,"allele extras (l,r)",leftmatchAllele,rightmatchAllele,"Locus is possibly on the 3' tip of the contig \n")
							
							
						
						elif 	rightmatchContig < rightmatchAllele:
							
							resultsList.append('PLOT5:-1')
							perfectMatchIdAllele.append('PLOT5')
							perfectMatchIdAllele2.append('PLOT5')
							
							verboseprint(match, "contig extras (l,r)",leftmatchContig,rightmatchContig,"allele extras (l,r)",leftmatchAllele,rightmatchAllele,"Locus is possibly on the 5' tip of the contig \n")
			
					
						else:
							#if a perfect match was found
									
							################################################
							# EXACT MATCH --- MATCH == GENE --- GENE FOUND #
							################################################
									
							perfectMatchIdAllele.append(str(bestmatch[4]))
							if not Reversed:
								perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation)+"&"+"+")
							else:
								perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation)+"&"+"-")
							resultsList.append('EXC:' + str(bestmatch[4]) )

				
				# if match with BSR >0.6 and not equal DNA sequences
				else:
					
					match=bestmatch[5]
					geneLen=bestmatch[6]
					alleleStr=currentCDSDict[">"+bestmatch[3]]
					contigname=bestmatch[3]	
					
					contigname=contigname.split("&")
					matchLocation=contigname[2]	
					matchLocation=matchLocation.split("-")
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
					
					
					
					if leftmatchContig<leftmatchAllele and 	rightmatchContig < rightmatchAllele:
					
						resultsList.append('LOTSC:-1')
						perfectMatchIdAllele.append('LOTSC')
						perfectMatchIdAllele2.append('LOTSC')
						
						
						verboseprint(match,contigname,geneFile,leftmatchAllele,rightmatchAllele,"Locus is bigger than the contig \n")
						

					
					elif leftmatchContig<leftmatchAllele:
						
						
						resultsList.append('LOT3:-1')
						perfectMatchIdAllele.append('LOT3')
						perfectMatchIdAllele2.append('LOT3')
						
						verboseprint(match,contigname,geneFile,leftmatchAllele,rightmatchAllele,"Locus is on the 3' tip of the contig \n")
						
					
					elif 	rightmatchContig < rightmatchAllele:
						
						resultsList.append('LOT5:-1')
						perfectMatchIdAllele.append('LOT5')
						perfectMatchIdAllele2.append('LOT5')
						
						verboseprint(match,contigname,geneFile,leftmatchAllele,rightmatchAllele,"Locus is on the 5' tip of the contig \n")
						
				
								
					elif len(alleleStr) > moda+(moda*0.2) :
						
						verboseprint("Locus is larger than mode", moda, alleleStr)
		
						resultsList.append('ALM')
						perfectMatchIdAllele.append('ALM')
						perfectMatchIdAllele2.append('ALM')
					
					elif len(alleleStr) < moda-(moda*0.2):
						
						verboseprint("Locus is smaller than mode", moda, alleleStr)
			
						resultsList.append('ASM')
						perfectMatchIdAllele.append('ASM')
						perfectMatchIdAllele2.append('ASM')
				
						
					else:
								#######################
								# ADD INFERRED ALLELE #		# a new allele 
								#######################
								
														
						tagAux='INF'
						perfectMatchIdAllele.append( tagAux +"-"+str(alleleI+1))
						
						if not Reversed:
							perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"+")
						else:
							perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"-")
						
						
						verboseprint( "New allele! Adding allele "+ tagAux + str(alleleI+1) +" to the database\n")
																							
						resultsList.append( tagAux + str(alleleI+1) )

													# --- add the new allele to the gene fasta --- #
						
						alleleI+=1
						appendAllele='>'+str((((os.path.basename(geneFile)).split("."))[0]).replace("_","-"))+"_" + str(alleleI) + "_" + str(os.path.basename(genomesList[genome])) + '\n'
						fG = open( geneFile, 'a' )
						fG.write(appendAllele)
						fG.write( alleleStr + '\n')
						fG.close()
						fullAlleleList.append(alleleStr)
						
						if bestmatch[1]>0.6 and bestmatch[1]<0.7:
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

							genefile2= geneTransalatedPath2
							Gene_Blast_DB_name2 = Create_Blastdb( genefile2, 1, True )
							verboseprint("Re-calculating BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
							allelescores,alleleList=reDogetBlastScoreRatios(genefile2,basepath,alleleI,allelescores,Gene_Blast_DB_name2,alleleList,geneScorePickle,verbose)
							verboseprint("Done Re-calculating BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
			
			
			
			except Exception as e:
				print "some error occurred"
				print e

				#print 'Error on line {}'.format(sys.exc_info()[-1].tb_lineno)
				perfectMatchIdAllele2.append("ERROR")
				perfectMatchIdAllele.append("ERROR")
				resultsList.append('ERROR')  
		
		
		
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
