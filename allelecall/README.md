# chewBBACA: Quick and Dirty

**Important Notes before starting:**

 - For **chewBBACA**, the definition of an allele is that each allele
   must represent a complete Coding Domain Sequence, with starting codon and stop codon according to the [NCBI genetic code table 11](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). It will automatically exclude any allele which DNA sequence that does not contain start or stop codon and that it's length is not multiple of three. The allele identification is performed using [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/). 
 - All the referenced lists of files *must contain full path* for the files.
 - Make sure that your fasta files are UNIX format. If they were created in Linux or MacOS systems they should be in the correct format, but if they were created in Windows systems , you should do a a quick conversion using for example [dos2unix](http://linuxcommand.org/man_pages/dos2unix1.html).

----------
##FAQ

###Q: How to do the quick and dirty?   
A: Just follow the first 4 steps bellow. (new: use step 1.1)

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
* [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/) or above

----------

## 1.1 Fast run

### Run all analysis in a single command**

**if you have no schema**

`cd` to the folder where you want to do the analysis, create folders for:
folder 1 : genomes fasta files as base for schema creation
folder 2 : genomes fasta files to call the alleles (genomes from folder 1 will already be used)

`/home/user/chewBBACA/fullBBACA.py --cs genomes/cg/ 
--cpu 6 -o myAnalysis -t Streptococcus_Agalactiae --genomes ./genomes/other/`

**if you have a schema**

`cd` to the folder where you want to do the analysis, create folders for:
folder 1 : genomes fasta files to call the alleles
folder 2 : target genes fasta files 


`/home/user/chewBBACA/fullBBACA.py --genes schema_seed/ 
--cpu 6 -o myAnalysis -t Streptococcus_Agalactiae --genomes ./genomes/other/`

`--cs` path to folder with genomes to create schema

`--genes` path to folder with target genes

`--genes` path to folder with target genomes

`--cpu` Number of cpus to use

`-o` folder name with the analysis files

`-t` (Optional) taxon to use for prodigal training input

----------

## 2. wgMLST schema creation

Create your own wgMLST schema based on a set of genomes fasta files. The command is the following:

`PPanGen.py -i listGenomes -o OutputFolderName --cpu 4`

**Parameters**

`-i` file containing the path to the list of genomes. One file path (must be full path) to any fasta/multifasta file containing all the complete or draft genomes you want to call alleles for.

`-o` prefix for the output folder for the schema, chewBBACA ready

`--cpu` Number of cpus to use

`--bsr` (Optional) Minimum BSR for locus similarity. Default at 0.6. 

`-t` (Optional) taxon to use for prodigal training input

**Outputs:** 

One fasta file per gene in the `-o`directory that is created in the dir where the files are found. The fasta file names are the given according the FASTA annotation for each coding sequence. For example the locus with the annotation ` >gi|193804931|gb|AE005672.3|:2864-3112 Streptococcus pneumoniae TIGR4, complete genome` will create the fasta file named  `gi_193804931_gb_AE005672.3_:2864-3112.fasta`. It will also create the necessary files for the allele call, by creating the directory named short. The contents of this dir is already explained in the folder structure subsection.

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

`-t` (Optional) taxon to use for prodigal training input

This will use by default, number of CPUs available minus 2 and can be called in a SLURM HPC by  srun  

`srun -c 62 --mem 250G BBACA.py -i listgenomes.txt -g listgenes.txt -o OutPrefix --cpu 62`

**Outputs files**:
```
./< outPrefix >_< datestamp>/< outPrefix >_statistics.txt
./< outPrefix >_< datestamp>/< outPrefix >_contigsInfo.txt
./< outPrefix >_< datestamp>/< outPrefix >_Alleles.txt 
./< outPrefix >_< datestamp>/< outPrefix >logging_info.txt 
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

Example of an output can be seen [here] (http://i.imgur.com/jlTV2vg.png) . This examples uses an original set of 1042 genomes and a scheme of 5266 loci, using a parameter `-n` of 12 and `-t` of 300.

=============
