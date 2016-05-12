#!/usr/bin/python
import HTSeq
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
#from Counter import Counter #Counter.py is needed to run this script
#from collections import Counter
import os
import string
from datetime import datetime
from Bio import SearchIO

def extractPrefix(filename):
	
	prefix=''
	l=len(filename)
	while l > 0:
		if filename[l-1] == '.':
			break
		l-=1
	
	if l > 1:
		l-=1
		return filename[0:l]
	else:
		print "File prefix not found. Using full name."
		return filename

def testFile(filename):

	try:
		with open(filename): pass
	except IOError:
		print 'This file does not exist: ' + filename
		sys.exit()

def ensure_dir(f):
	if not os.path.isdir(f):
		#print "No BLAST db dir was found. Creating it..."
		os.makedirs(f)



# overwrite has to be 1 or 0 (True or False)
def Create_Blastdb( questionDB, overwrite, dbtypeProt ):
	base = os.path.basename(questionDB)
	dirname = os.path.dirname( questionDB )
	isProt=dbtypeProt
	
	if len(dirname)==0:
		dirname='.'
	basename = os.path.splitext(base)[0]
	ensure_dir(dirname + "/blastdbs")
	name = dirname+"/blastdbs/" + basename + "_db"

	if not os.path.isfile(name + ".nin") and not os.path.isfile(name + ".nhr") and not os.path.isfile(name + ".nsq"):
		
		if not isProt:
			os.system( "makeblastdb -in " + questionDB + " -out " + name + " -dbtype nucl -logfile " + name + "_blast.log" )
		else:
			os.system( "makeblastdb -in " + questionDB + " -out " + name + " -dbtype prot -logfile " + name + "_blast.log" )

	elif overwrite:
		if not isProt:
			os.system( "makeblastdb -in " + questionDB + " -out " + name + " -dbtype nucl -logfile " + name + "_blast.log" )
		else:
			os.system( "makeblastdb -in " + questionDB + " -out " + name + " -dbtype prot -logfile " + name + "_blast.log" )

	else:
		print "BLAST DB files found. Using existing DBs.."  
	return( name )

# def LoadAlelleFasta( Allele_fasta_file ):
# 	fp = HTSeq.FastaReader( Allele_fasta_file )
# 	#allele_dict = dict( ( s.seq, s.name ) for s in fp )
# 	allele_dict = {}  
# 	allele_sizes = []
# 	for s in fp:
# 		allele_dict[s.seq.rstrip()]=s.name
# 		allele_sizes.append(len(s.seq.rstrip()))
# 		#print "Allele :" + s.name + " size: " + str(len(s.seq.rstrip()))#DEBUG
# 	data = Counter(allele_sizes)
# 	most_common_size = data.most_common(1)
# 	allele_sizes = list(set(allele_sizes))	
# 	# create BLAST db
# 	db_name=Create_Blastdb( Allele_fasta_file )
# 	return (allele_dict, allele_sizes, most_common_size, db_name)
	
def LoadAlellicProfile( AP_file ):
	AP_dict={}
	fp = open( AP_file, "r" )
	for line in fp:
		line = line.rstrip()
		st, g1, g2, g3, g4 , g5, g6, g7 = line.split("\t")
		AP_dict[ g1 + '-' + g2 + '-'+ g3 + '-'+ g4 + '-'+ g5 + '-'+ g6 + '-'+ g7 ] = st
	fp.close()
	return( AP_dict )

# similar to LoadAlellicProfile but outputs the first line of the file, which contains the right order of the MLST genes
def LoadAlellicProfile( AP_file ):
        AP_dict={}
        fp = open( AP_file, "r" )

	firstLineStr = ''
	firstLine = True
        for line in fp:
		if firstLine:
			firstLineStr = line
			firstLine = False
		else:
                	line = line.rstrip()
                	st, g1, g2, g3, g4 , g5, g6, g7 = line.split("\t")
                	AP_dict[ g1 + '-' + g2 + '-'+ g3 + '-'+ g4 + '-'+ g5 + '-'+ g6 + '-'+ g7 ] = st
        fp.close()
        return( AP_dict, firstLineStr )

	
def LoadAlellicProfileGeneric(AP_file):
	AP_dict={}
	AP_list=[]
	fp = open(AP_file, "r")
	for line in fp:
		key=''
		line = line.rstrip()
		AP_list = line.split("\t")
		for i in AP_list:
			key=key+i+'-'
		key=key[:-1]
		AP_dict[key] = st
	fp.close()
	return(AP_dict)
	
def WriteFasta(filename,tag,sequence,mode):
	fp = open (filename, mode) 
	fp.write(tag)
	fp.write(sequence)
	fp.close()

def getSeqsFromAP(AP, filesList):
	
	APList=AP.split('-')
	i=0
	seqList=[]
	
	for f in filesList:
		record_dict = SeqIO.index(f, "fasta")
		iniName=extractPrefix(f)+'_'+APList[i]
		seqList.append(SeqRecord('', iniName, iniName, iniName, None, None, None, None))
		for k in record_dict.keys():
			if k.split('-')[-1]==APList[i]:
				seqList[i]=record_dict[k]
				break
		i+=1
	
	return seqList

def runBlast(cline, bOutFile, locus_sbjct):
	os.system(str(cline))
	rec = open(bOutFile)
	blast_record = NCBIXML.read(rec)
	
	if os.path.isfile(locus_sbjct):
		os.remove(locus_sbjct)
	os.remove(bOutFile)
	
	return blast_record

def runBlastParser(cline, bOutFile, locus_sbjct):


	os.system(str(cline))
	rec = open(bOutFile)
	blast_records = NCBIXML.parse(rec)

#	if os.path.isfile(locus_sbjct):
#		os.remove(locus_sbjct)

	#os.remove(bOutFile)

	return blast_records


