#!/usr/bin/env python
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
import json

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

def call_proc(cmd):
	
	p= subprocess.Popen(args=cmd)
	out, err = p.communicate()
	return p

# ================================================ MAIN ================================================ #

def main():
			
	parser = argparse.ArgumentParser(description="This program call alleles for a set of genomes provided a schema")
	#~ parser.add_argument('-i', nargs='?', type=str, help='List of genome files (list of fasta files)', required=True)
	#~ parser.add_argument('-g', nargs='?', type=str, help='List of genes (fasta)', required=True)
	#~ parser.add_argument('-o', nargs='?', type=str, help="Name of the output files", required=True)
	#~ 
	#~ parser.add_argument("-v", "--verbose", help="increase output verbosity",dest='verbose', action="store_true",default=False)
	#~ 
	#~ 
	#~ parser.add_argument("--so", help="split the output per genome",dest='divideOutput', action="store_true",default=False)
	#~ parser.add_argument('-t', nargs='?', type=str, help="taxon", required=False,default=False)
	#~ parser.add_argument("--fc", help="force continue",required=False, action="store_true",default=False)
	#~ parser.add_argument("--json", help="report in json file",required=False, action="store_true",default=False)
	
	parser.add_argument("--cs",  nargs='?', type=str, help='create schema', required=True)
	parser.add_argument('--bsr', nargs='?', type=float, help="minimum BSR score", required=False,default=0.6)
	parser.add_argument('--cpu', nargs='?', type=int, help="Number of cpus, if over the maximum uses maximum -2", required=True)
	parser.add_argument('-b', nargs='?', type=str, help="BLAST full path", required=False,default='blastp')
	
	args = parser.parse_args()
	
	#~ genomeFiles = args.i
	#~ genes = args.g
	#~ 
	#~ 
	#~ verbose=args.verbose
	#~ 
	#~ divideOutput=args.divideOutput
	#~ gOutFile = args.o
	#~ chosenTaxon=args.t
	#~ forceContinue=args.fc
	#~ jsonReport=args.json
	
	
	genomes2CreateSchema = args.cs
	BSRTresh=args.bsr
	cpuToUse=args.cpu
	BlastpPath=args.b
	
	# avoid user to run the script with all cores available, could impossibilitate any usage when running on a laptop
	if cpuToUse > multiprocessing.cpu_count()-2:
		print "Warning, you are close to use all your cpus, if you are using a laptop you may be uncapable to perform any action"
	
	taxonList={'Campylobacter_Jejuni':'trained_campyJejuni.trn',
				'Acinetobacter_Baumannii':'trained_acinetoBaumannii.trn',
				'Streptococcus_Agalactiae':'trained_strepAgalactiae.trn',
				'Haemophilus_Influenzae':'trained_haemoInfluenzae_A.trn',
				'Yersinia_Enterocolitica':'trained_yersiniaEnterocolitica.trn',
				'Escherichia_Coli':'trained_eColi.trn',
				'Enterococcus_Faecium':'trained_enteroFaecium.trn',
				'Staphylococcus_Haemolyticus':'trained_staphHaemolyticus.trn',
				'Salmonella_Enterica_enteritidis':'trained_salmonellaEnterica_enteritidis.trn'
				}	
	if isinstance(chosenTaxon, basestring):
		trainingFolderPAth=os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'TrainingFiles4Prodigal'))
		try:
			chosenTaxon=os.path.join(trainingFolderPAth,taxonList[chosenTaxon])
			
			if os.path.isfile(chosenTaxon):
				print "will use this training file : "+chosenTaxon
			else:
				print "training file don't exist"
				print chosenTaxon
				return "retry"
		except :
			print "Your chosen taxon is not attributed, select one from:"
			for elem in taxonList.keys():
				print elem
			return "retry"
	
	print BlastpPath
	
	scripts_path=os.path.dirname(os.path.realpath(__file__))
	
	print ("Will use this number of cpus: "+str(cpuToUse))
	print ("Checking all programs are installed")
	
	print ("Checking Blast installed... "+str(which(str(BlastpPath))))
	print ("Checking Prodigal installed... "+str(which('prodigal')))
	
	#check version of Blast
	
	proc = subprocess.Popen([BlastpPath, '-version'], stdout=subprocess.PIPE)
	line = proc.stdout.readline()
	if not "blastp: 2.5." in str(line) and not "blastp: 2.6." in str(line):
		print "your blast version is " + str(line)
		print "update your blast to 2.5.0 or above, will exit program"
		sys.exit()
	else:
		print "blast version is up to date, the program will continue"
	
	starttime="\nStarting Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y")
	print (starttime)
	
	
	# check if given a list of genomes paths or a folder to create schema
	try:
		f=open( genomes2CreateSchema, 'r')
		f.close()
	except IOError:
		listbasename=os.path.basename(os.path.normpath(genomes2CreateSchema))
		
		with open("listGenes"+listbasename+".txt", "wb") as f:
			for gene in os.listdir(genomes2CreateSchema):
				try:
					genepath=os.path.join(genomes2CreateSchema,gene)
					gene_fp2 = HTSeq.FastaReader(genepath)
					for allele in gene_fp2:
						break
					f.write( genepath+"\n")
				except Exception as e:
					print e
					pass
		
		genomes2CreateSchema="listGenes"+listbasename+".txt"
	
	# call ppan
	proc = subprocess.Popen([createSchemaPath, '-i', genomes2CreateSchema,'-l', "200",'--cpu', cpuToUse,"-b",BlastpPath,"-o",outputFile,"--bsr",str(bsr)])
	p_status = proc.wait()
	
	# create list with genes to call and call bbaca
	
	
	# use bbaca outputs, remove paralogs and testqualitygenomes, show plot and ask for threshold to use for phyloviz, set default to 20
	
	
	#send data to phyloviz online ????
	
	
	
	# run shema evaluator	
	
	
	
	
	
	print (starttime)
	print ("Finished Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

if __name__ == "__main__":
    main()
