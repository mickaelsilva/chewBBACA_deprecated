#!/usr/bin/python
import sys
import csv
import numpy as np
from numpy import array
import argparse
import collections
from collections import OrderedDict

def presAbs (d3,listgenomesRemove):	

	d2c=np.copy(d3)
	

	geneslist= d2c[:1,:]
	genomeslist= d2c[:,:1]
	genomeslist=(genomeslist.tolist())
	
	#remove genomes
	print "removing the following genomes..."
	for genome in genomeslist:
		if genome[0] in listgenomesRemove:
			print genome[0]
			rowid=genomeslist.index(genome)
			genomeslist.pop(rowid)
			d2c=np.delete(d2c, rowid, 0)
	
	print "all genomes removed"
	
	"building the presence and abscence matrix..."
	row=1
	while row<d2c.shape[0]:
		column=1
		while column<d2c.shape[1]:
			try:

				aux=int(d2c[row,column])
				if aux>0:
					d2c[row,column]=1
				else :
					d2c[row,column]=0
			except:
				try:
					aux=str((d2c[row,column])).replace ("INF-","")
					aux=int(aux)
					d2c[row,column]=1
				except Exception as e:
					d2c[row,column]=0
				
			column+=1
		row+=1
	
	"presence and abscence matrix built"

	d2d=d2c.tolist()
	
	with open ("presence.txt","wb") as f:
		
		f.write('\t'.join(geneslist[0]))

		writer = csv.writer(f,delimiter='	')
		writer.writerows(d2d)
		
	
	return d2c

def clean (inputfile,outputfile,totaldeletedgenes,rangeFloat,toremovegenes,toremovegenomes):
	
	#open the raw file to be clean
	
	with open(inputfile) as f:
		reader = csv.reader(f, delimiter="\t")
		d = list(reader)
	
	originald2 = array(d)
	
	#get presence abscence matrix
	d2=presAbs (originald2,toremovegenomes)
	
	genomeslist= d2[1:,:1]
	
	originald2=originald2.T
	d2=d2.T
	rowid=1
	deleted=0
	abscenceMatrix=True
	lnfdel=0
	balldel=0
	
	
	
	numbergenomes=len(genomeslist)
	pontuationmatrix=[0]*numbergenomes	
	lostgenesList=[]
	
	#clean the original matrix, using the information on the presence/abscence matrix
	
	while rowid< d2.shape[0]:
		columnid=1
		genomeindex=0
		
		
		while columnid< d2.shape[1]:

			if d2[rowid][0] in toremovegenes:
				
				originald2=np.delete(originald2, rowid, 0)
				d2=np.delete(d2, rowid, 0)
				totaldeletedgenes+=1
				deleted+=1
				rowid-=1
				columnid=1
				lnfdel+=1
				break
			
			elif not abscenceMatrix:
				
				if ( "PLOT" in d2[rowid][columnid] or "NIPL" in d2[rowid][columnid] or "LNF" in d2[rowid][columnid] or "LOT" in d2[rowid][columnid]  or "ALM" in d2[rowid][columnid]  or "ASM" in d2[rowid][columnid] or "ERROR" in d2[rowid][columnid]  or "ABM" in d2[rowid][columnid]) or "undefined allele" in d2[rowid][columnid] or "small match" in d2[rowid][columnid] or "allele incomplete" in d2[rowid][columnid]:	
					
					originald2=np.delete(originald2, rowid, 0)
					d2=np.delete(d2, rowid, 0)
					totaldeletedgenes+=1
					deleted+=1
					rowid-=1
					columnid=1
					lnfdel+=1
					break
			else:
				if int(d2[rowid,columnid]) == 0:
					
					originald2=np.delete(originald2, rowid, 0)
					d2=np.delete(d2, rowid, 0)
					totaldeletedgenes+=1
					deleted+=1
					rowid-=1
					columnid=1
					lnfdel+=1
					break


			
			columnid+=1
			genomeindex+=1

		rowid+=1
	
	originald2=originald2.T
	originald2=originald2.tolist()
	
	#write the output file
	
	with open(outputfile, "wb") as f:
		writer = csv.writer(f,delimiter='	')
		writer.writerows(originald2)
	
	#chewbbaca files have INF that needs to be removed
	file = open(outputfile)
	contents = file.read()
	contents = contents.replace('INF-', '')

	
	with open(outputfile, 'w') as f:
			f.write(contents)
	

	print "deleted : %s loci" % totaldeletedgenes
	print "total loci remaining : "+ str(rowid-2)


def main():

	parser = argparse.ArgumentParser(description="This program cleans an output file for phyloviz")
	parser.add_argument('-i', nargs='?', type=str, help='output to clean', required=True)
	parser.add_argument('-o', nargs='?', type=str, help='name of the clean file', required=True)
	parser.add_argument('-r', nargs='?', type=str, help='listgenes to remove', required=False)
	parser.add_argument('-g', nargs='?', type=str, help='listgenomes to remove', required=False)
	
	args = parser.parse_args()

	
	pathOutputfile = args.i
	newfile = args.o
	
	
	
	genesToRemove=[]
	genomesToRemove=[]
	
	try:
		genesToRemoveFile = (args.r)
		fp = open(genesToRemoveFile, 'r')
		
		for geneFile in fp:

			geneFile = geneFile.rstrip('\n')
			geneFile = geneFile.rstrip('\r')
			
			genesToRemove.append( geneFile )
			
	except:		
		pass
	
	try:
		genomesToRemoveFile = (args.g)
		fp = open(genomesToRemoveFile, 'r')
		
		for genomeFile in fp:

			genomeFile = genomeFile.rstrip('\n')
			genomeFile = genomeFile.rstrip('\r')
			
			genomesToRemove.append( genomeFile )
			
	except:		
		pass
		
	clean(pathOutputfile,newfile,0,0.2,genesToRemove,genomesToRemove)

	

	
	
	
	
		
if __name__ == "__main__":
	main()
