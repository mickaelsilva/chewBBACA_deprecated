#!/usr/bin/python
import HTSeq
from Bio.Seq import Seq
import os
import argparse
import mpld3
import matplotlib.pyplot as plt
import json

def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        strDNArevC = ''
        for l in strDNA:

        	strDNArevC += basecomplement[l]

        return strDNArevC[::-1]

def translateSeq(DNASeq,transTable):
	seq=DNASeq
	tableid=transTable
	reversedSeq=False
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

					raise ValueError(e)
	return protseq,seq,reversedSeq

def main():
	
	parser = argparse.ArgumentParser(description="This program analyses cds")
	parser.add_argument('-i', nargs='?', type=str, help='list genes', required=True)
	parser.add_argument('-r', nargs='?', type=bool, help='Return values', required=False)
	
	args=parser.parse_args()

	
	genes = args.i
	
	try:
		ReturnValues=bool(args.r)
	except:
		ReturnValues=False
		pass
		
	analyzeCDS(genes,ReturnValues)

def analyzeCDS(genes,transTable,ReturnValues,outputpath):
	
	gene_fp = open( genes, 'r')
	
	stopc=0
	notStart=0
	notMultiple=0
	totalalleles=0
	
	htmlgenespath=os.path.join(outputpath,"genes_html/")
	
	if not os.path.exists(htmlgenespath):
		os.makedirs(htmlgenespath)
	
	
	statsPerGene={}
	
	for gene in gene_fp:
		
		listStopc=[]
		listnotStart=[]
		listnotMultiple=[]
		
		
		print "####################"
		print str(os.path.basename(gene))
		
		k=0
		gene = gene.rstrip('\n')
		multiple=True
		gene_fp2 = HTSeq.FastaReader(gene)
		
		alleleSizes=[]
		# translate each allele and report the error if unable to translate
		for allele in gene_fp2: 
			
			k+=1
			alleleSizes.append(len(allele.seq))
			# if allele is not multiple of 3 it's useless to try to translate
			if (len(allele.seq) % 3 != 0):
				multiple=False
				listnotMultiple.append(str(k))
				print "allele "+str(k)+" is not multiple of 3"
				pass
			else:
				try:
					protseq,seq,reversedSeq=translateSeq(allele.seq, transTable)
					
				except Exception, err:
					if "Extra in frame stop codon found" in str(err):
						stopc+=1
						listStopc.append(str(k))
					elif "is not a start codon" in str(err):
						notStart+=1
						listnotStart.append(str(k))
					elif "is not a stop codon" in str(err):
						notStart+=1
						listnotStart.append(str(k))
					else:
						print err
					print "allele "+str(k)+" is not translating"
					pass
		
		relpath=os.path.relpath(gene,outputpath)
		statsPerGene[relpath]=listnotMultiple,listStopc,listnotStart,k
		totalalleles+=k
		
		#create html per gene
		#genename=((os.path.basename(gene)).split("."))[0]
		genename=(os.path.basename(gene)).split(".")
		genename.pop()
		genename=".".join(genename)
		with open(htmlgenespath+genename+".html", "wb") as f:
			f.write("<!DOCTYPE html>\n<html>\n<head><script type='text/javascript' src='https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.5/d3.min.js'></script><script type='text/javascript' src='https://mpld3.github.io/js/mpld3.v0.2.js'></script>\n")
			f.write("<script type='text/javascript' src='https://ajax.googleapis.com/ajax/libs/jquery/2.1.3/jquery.min.js'></script>")
			f.write("""<!-- Latest compiled and minified JavaScript -->
	<script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/js/bootstrap.min.js" integrity="sha384-0mSbJDEHialfmuBBQP6A4Qrprq5OVfW37PRR3j5ELqxss1yVqOtnepnHVP9aJ7xS" crossorigin="anonymous"></script>""")
			f.write("""<!-- Latest compiled and minified CSS -->
	<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css" integrity="sha384-1q8mTJOASx8j1Au+a5WDVnPi2lkFfwwEAa8hDDdjZlpLegxhjVME1fgjWPGmkzs7" crossorigin="anonymous">""")
			
			f.write("</head>\n<body><div id='histdiv'></div>\n")
			relpath=os.path.relpath(gene,htmlgenespath)
			#print relpath
			f.write("<input type=button onClick=window.open('"+relpath+"') value='click here to get fasta file'>\n""")
			fig, ax = plt.subplots(figsize=(20,10))
			bp=plt.scatter(list(range(1,len(alleleSizes)+1,1)),alleleSizes)
			plt.ylabel('DNA bp allele length ')
			plt.xlabel('Allele number')
			plt.title('Allele size scatter plot')
			ax.yaxis.labelpad = 40
			#start, end = ax.get_xlim()

			plt.grid(True)
			histplothtml=mpld3.fig_to_dict(fig)
			histplothtml=str(json.dumps(histplothtml))
			
			f.write("<script type='text/javascript'>var hist ="+str(histplothtml)+";mpld3.draw_figure('histdiv', hist);</script></body></html>")
			plt.close('all')
			
	print str(stopc) + " alleles have stop codons inside"
	print str(notStart) + " alleles don't have start codons"
	print "total of alleles : " + str(totalalleles)
	
	if not ReturnValues:
		with open("CheckCDSResults.txt", "wb") as f:
			f.write("Alleles with stop codons inside: \n")
			for item in listStopc:
				f.write(item)
				f.write("\n")
			f.write("\nAlleles without start/stop codon: \n")
			for item in listnotStart:
				f.write(item)
				f.write("\n")
	else:
		return statsPerGene
	
if __name__ == "__main__":
    main()

