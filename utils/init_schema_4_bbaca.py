#!/usr/bin/python
import HTSeq
import sys
from Bio.Seq import Seq
import os
import argparse


def loci_translation (genesList):
	
	gene_fp = open( genesList, 'r')

	genepath=''
	basepath=''
	lGenesFiles = []
	argumentsList = []

	
	for gene in gene_fp:
		gene = gene.rstrip('\n')
		shortgene= os.path.join(os.path.dirname(gene),"short",os.path.basename(gene))
		shortgene= shortgene.replace(".fasta","_short.fasta")

			
		gene_fp2 = HTSeq.FastaReader(gene)
		for allele in gene_fp2: 
			fG = open( shortgene, 'w' )
			fG.write('>'+str(allele.name)+'\n'+str(allele.seq) + '\n')
			fG.close()
			break
		

	

	gene_fp.close()	
	return	True

def main():
			
	parser = argparse.ArgumentParser(description="This program prepares a schema for a chewBBACA allele call, creating a short version of each fast with only the 1st allele")
	parser.add_argument('-i', nargs='?', type=str, help='List of genes files (list of fasta files)', required=True)

	
	
	
	args = parser.parse_args()
	
	geneFiles = args.i
	
	
	loci_translation(geneFiles)

if __name__ == "__main__":
    main()
