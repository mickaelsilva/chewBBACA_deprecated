import csv
import numpy as np
from numpy import array
import argparse
import collections
import operator

def main():

	parser = argparse.ArgumentParser(description="Check if locus are being represented more than once")
	parser.add_argument('-i', nargs='?', type=str, help='contig info file', required=True)
	
	
	args = parser.parse_args()
	contigsfile = args.i
	
	
	with open(contigsfile) as f:
		reader = csv.reader(f, delimiter="\t")
		d = list(reader)
	
	d2 = array(d)
	
	genelist= d2[:1,:]
	genelist=genelist.tolist()[0]
	genomeslist= d2[:,:1]
	genomeslist=genomeslist.tolist()
	genomeslist2=[]

	i=1
	pontuationDict={}
	
	
	while i<len(genomeslist):

		genomeSchema2=d2[i]
		genomeSchema=genomeSchema2.tolist()
		genomeSchema= "\t".join(genomeSchema)
		
		numberLOT=genomeSchema.count("LOT")
		numberLNF=genomeSchema.count("LNF")
		numberNIPL=genomeSchema.count("NIPL")
		numberASM=genomeSchema.count("ASM")
		numberALM=genomeSchema.count("ALM")
		
		numberMissingData=numberLOT+numberLNF+numberNIPL+numberASM+numberALM
		if numberMissingData>0:
			pontuationDict[(genomeslist[i])[0]]=numberMissingData
			genomeslist2.append((genomeslist[i])[0])
		i+=1	
	
	for k,v in pontuationDict.items():
		print k,v

	print genomeslist2
	
	"""with open("RepeatedLoci.txt", "wb") as f:
		f.write("gene\toverrepresented\tproblems\ttotal\n")
		for k,v in ordered:
			try:
				troubledLocus=str(pontuationDict2[k])
			except:
				troubledLocus="0"
				pass	
			f.write( (k+ "\t"+str( v)+"\t"+ troubledLocus+"\t"+str(int(v)+int(troubledLocus))+"\n"))"""
			
	
if __name__ == "__main__":
	main()
