# chewBBACA: Quick and Dirty

**Important Notes before starting:**

 - For **chewBBACA**, the definition of an allele is that each allele
   must represent a complete Coding Domain Sequence, with starting codon and stop codon according to the [NCBI genetic code table 11](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). It will automatically exclude any allele which DNA sequence that does not contain start or stop codon and that it's length is not multiple of three. The allele identification is performed using [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/). 
 - All the referenced lists of files *must contain full path* for the files.
 - Make sure that your fasta files are UNIX format. If they were created in Linux or MacOS systems they should be in the correct format, but if they were created in Windows systems , you should do a a quick conversion using for example [dos2unix](http://linuxcommand.org/man_pages/dos2unix1.html).

----------
##FAQ

###Q: How to do the quick and dirty?  
A: Just follow the first 4 steps bellow.

###Q: I ran all the steps and my cgMLST loci size is smaller than traditional MLST, does this even work?  
A: Try to run step 5, your analysis may contain some genomes responsible for a considerable loss of loci. Remove some of those genomes and check if the cgMLST loci number rises from the ashes.

###Q: Step 2 is taking hours, will it ever end?  
A: Depending on the variability of the strains you are using to create the schema and the number of cpu you are using. The more variable the strains, the more BLAST comparisons are made.

###Q: Step 3 just crashed at 99% after 2 days running and I smashed my computer, will you refund me?  
A: chewBBACA shoul allow you to continue where you stopped, just re-run the same command and you should be prompted to continue the allele call.

----------

## 1. Setting up the analysis

**Installing chewBBACA**

Git clone the whole repository.

You need to install the following dependencies. Prodigal and BLAST must be added to the PATH variables.

You may use this pip command to install the python dependencies automatically :

```
pip install -r requirements.txt
```


Python dependencies:
* [Biopython 1.66 ](http://biopython.org/wiki/Main_Page)
* [HTSeq 0.6.1p1](http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)


Other dependencies:
* BLAST 2.5.0+ ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.5.0/
* [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/)

----------

## 2. wgMLST schema creation

Create your own wgMLST schema. The command is the following:

`CreateSchema.py -i allffnfile.fasta -l 200 --cpu 4`

**Parameters:**

`-i`  Input file with concatenated gene sequences

`-l` minimum locus lenght (removes any loci with length equal or less the specified value

`--cpu` (optional) number of cpus to use for BLAST, will work for the latest BLAST releases


**Input file :**

 `allffnfile.fasta` : a fasta file resulting from concatenating all the fasta files containing all the genes for the genomes we want to use for the creating the wgMLST schema. In Genebank this files usually have the *.fnn*  extension. A commonly used and highly efficient software for bacterial genome annotation, [Prokka](https://github.com/tseemann/prokka), also provide *.ffn* files with its output. Another source for bacteria .ffn files is NCBI at ftp://ftp.ncbi.nih.gov/genomes/archive/old_genbank/Bacteria/

**Outputs:** 

One fasta file per gene in the `schema_seed/`directory that is created in the dir where the files are found. 

----------

## 3.  Allele call using the wgMLST schema 

Create two list of files with the full paths (one path per line), one list for genomes and another for genes (genes are located on the schema seed created in the last step, ignore the short folder)

Then run is the following:

	% BBACA.py -i listGenomes.txt -g listGenes.txt -o OutPrefix --cpu 3 

**Parameters** 

`-i` file containing the path to the list of genomes. One file path (must be full path) to any fasta/multifasta file containing all the complete or draft genomes you want to call alleles for.

`-g` file containing the path to the list of alleles

`-o` prefix for the output files. ID for the allele call run

`--cpu` Number of cpus to use 

`-v`,`--verbose`  verbose mode(optional). Provides more output of the run.

`-b` Blastp full path(optional). In case of slurm system BLAST version being outdated it may be hard to use a different one, use this option using the full path of the blastp executable


This will use by default, number of CPUs available minus 2 and can be called in a SLURM HPC by  srun  

`srun -c 62 --mem 250G BBACA.py -i listgenomes.txt -g listgenes.txt -o OutPrefix --cpu 62`

**Outputs files**:
```
./< outPrefix >_< datestamp>/< outPrefix >_statistics.txt
./< outPrefix >_< datestamp>/< outPrefix >_contigsInfo.txt
./< outPrefix >_< datestamp>/< outPrefix >_Alleles.txt 
```


----------


### 4. Defining the cgMLST schema

 **Creating a clean allelic profile for PHYLOViZ** 
 
Clean a raw output file from an allele calling to a phyloviz readable file. Use the Alleles.txt output file.


Basic usage:

	% Extract_cgAlleles.py -i rawDataToClean.txt -o cleanedOutput.txt -r removeLocusList.txt -g removeGenomesList.txt
	
`-i` raw output file from an allele calling

`-o` path/name of the clean file

`-r` (optional) list of genes to remove, one per line, advised to use the detected overrepresented genes from ParalogPrunning.py

`-g` (optional) list of genomes to remove, one per line, advised to use genomes that perform worst on statistics or use the result from the testGenomeQuality script

=============

## 5. Evaluate wgMLST call quality per genome


Usage:


	% TestGenomeQuality.py -i out.txt -n 12 -t 250
	
`-i` raw output file from an allele calling

`-n` maximum number of iterations, each iteration removes a set of genomes over the threshold and recalculates all variables

`-t` maximum threshold, will start at 5 increasing in a step of 5 until t

The output consists in a set of plots per iteration and a removedGenomes.txt file where its informed of which genomes are removed per threshold when it reaches a stable point (no more genomes are removed)

Example of an output can be seen [here] (http://i.imgur.com/uQDNNkb.png) . This examples uses an original set of 1042 genomes and a scheme of 5266 loci, using a parameter `-n` of 12 and `-t` of 300.

=============
