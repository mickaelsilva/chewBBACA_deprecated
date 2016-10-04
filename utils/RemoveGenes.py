import csv
import numpy as np
from numpy import array
import argparse
import collections
import operator

def main():

	parser = argparse.ArgumentParser(description="This program removes gens from a matrix")
	parser.add_argument('-i', nargs='?', type=str, help='main matrix file from which to remove', required=True)
	parser.add_argument('-o', nargs='?', type=str, help='list of genes to remove', required=True)
	
	args = parser.parse_args()
	mainListFile = args.i
	toRemoveListFile = args.o
	
	FilesToRemove=[]
	with open(toRemoveListFile) as f:
		for File in f:
			FilesToRemove.append( File.replace("\n",'' ))
	#print FilesToRemove
	
	
	with open(mainListFile,'rb') as tsvin, open("new.tsv", "a") as csvout:
		tsvin = csv.reader(tsvin, delimiter='\t')

		listindextoremove=[]
		for firstline in tsvin:
			for gene in firstline:
				if gene in FilesToRemove:
					listindextoremove.append(firstline.index(gene))
			for elem in reversed(listindextoremove):
				del firstline[elem]
			csvout.write(('\t'.join(firstline))+"\n")
			break
			
		for line in tsvin:
			for elem in reversed(listindextoremove):
				del line[elem]
			csvout.write(('\t'.join(line))+"\n")

	
if __name__ == "__main__":
	main()
