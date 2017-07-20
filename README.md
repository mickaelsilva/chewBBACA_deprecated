# chewBBACA: Quick Usage

**chewBBACA** stands for "BSR-Based Allele Calling Algorithm". The "chew" part could be thought of as "Comprehensive and  Highly Efficient Workflow" but at this point still it needs a bit of work to make that claim so we just add "chew" to add extra coolness to the software name. BSR stands for BLAST Score Ratio as proposed by [Rasko DA et al.](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-2) 

chewBBACA is a comprehensive pipeline for the creation and validation of whole genome and core genome MultiLocus Sequence Typing (wg/cgMLST) schemas, providing an allele calling algorithm based on Blast Score Ratio that can be run in multiprocessor settings and a set of functions to visualize and validate allele variation in the loci.

----------
## Check the [wiki pages](https://github.com/mickaelsilva/chewBBACA/wiki) ...
...for a much more thorough chewBBACA walktrough. 
In this readme you can find a quick list of commands to be used by experienced users.

## Use [BBACA gitter](https://gitter.im/BBACA/Lobby)...

... if you have any pressing question. Chat can be faster and better than email for troubleshooting purposes 

**Important Notes before starting:**

 - For **chewBBACA**, the definition of an allele is that each allele
   must represent a complete Coding DNA Sequence, with start and stop codon according to the [NCBI genetic code table 11](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). It will automatically exclude any allele for which the DNA sequence does not contain start or stop codons and for which the length is not multiple of three. The CDS identification is performed using [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/). 
 - All the referenced lists of files *must contain full path* for the files.
 - Make sure that your fasta files are UNIX format. If they were created in Linux or MacOS systems they should be in the correct format, but if they were created in Windows systems, you should do a a quick conversion using for example [dos2unix](http://linuxcommand.org/man_pages/dos2unix1.html).


## 0. Setting up the analysis

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
* BLAST 2.5.0+ ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.5.0/ or above
* [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/) or above

----------

## 1. wgMLST schema creation

Create your own wgMLST schema based on a set of genomes fasta files. The command is the following:

`chewBBACA.py CreateSchema -i ./genomes/cg/ -o OutputFolderName --cpu 4`

**Parameters**

`-i` Folder containing the genomes from which to create the schema. Alternatively a file
 containing the path to the list of genomes. One file path (must be full path) 
 to any fasta/multifasta file containing all the complete or draft genomes you want to call alleles for.

`-o` prefix for the output folder for the schema

`--cpu` Number of cpus to use

`-t` (Optional) taxon to use for prodigal training input

`--bsr` (Optional) Minimum BSR for locus similarity. Default at 0.6. 

**Outputs:** 

One fasta file per gene in the `-o`directory that is created in the dir where the files are found. The fasta file names are the given according the FASTA annotation for each coding sequence. 

----------

## 2.  Allele call using the wgMLST schema 

Create two list of files with the full paths (one path per line), one list for genomes and another for genes (genes are located on the schema seed created in the last step, ignore the short folder)

Then run is the following:

`chewBBACA.py Allelecall -i ./genomes/other/ -g genes/ -o OutPrefix --cpu 3 `

**Parameters** 

`-i` Folder containing the genomes that you want to call the alleles. Alternatively a file
 containing the path to the list of genomes. One file path (must be full path) 
 to any fasta/multifasta file containing all the complete or draft genomes you want to call alleles for.
`-g` Folder containing the genes that you want to call. Alternatively a file
 containing the path to the list of genes. One file path (must be full path) 
 to any fasta/multifasta file containing all the complete or draft genomes you want to call alleles for.

`-o` prefix for the output files. ID for the allele call run

`--cpu` Number of cpus to use 

`-t` (Optional but recommended, contact for new species) taxon to use for prodigal training input

`-b` (optional)Blastp full path. In case of slurm system BLAST version being outdated it may 
be hard to use a different one, use this option using the full path of the updated blastp executable



**Outputs files**:
```
./< outPrefix >_< datestamp>/< outPrefix >/results_statistics.txt
./< outPrefix >_< datestamp>/< outPrefix >/results_contigsInfo.txt
./< outPrefix >_< datestamp>/< outPrefix >/results_Alleles.txt 
./< outPrefix >_< datestamp>/< outPrefix >logging_info.txt 
./< outPrefix >_< datestamp>/< outPrefix >RepeatedLoci.txt
```


----------

## 3. Evaluate wgMLST call quality per genome


Usage:


`chewBBACA.py TestGenomeQuality -i alleles.tsv -n 12 -t 200 -s 5 -o OutFolder`
	
`-i` raw output file from an allele calling

`-n` maximum number of iterations, each iteration removes a set of genomes over the threshold and recalculates all variables

`-t` maximum threshold, will start at 5

`-s` step to add to each threshold (suggested 5)

`-o` Folder for the analysis files

The output consists in a set of plots per iteration and a removedGenomes.txt file where its 
informed of which genomes are removed per threshold when it reaches a stable point (no more genomes are removed)

Example of an output can be seen [here] (http://i.imgur.com/jlTV2vg.png) . This examples uses an 
original set of 714 genomes and a scheme of 3266 loci, using a parameter `-n 12`,`-s 5` and `-t 300`.

----------
## 4. Defining the cgMLST schema

 **Creating a clean allelic profile for PHYLOViZ** 
 
Clean a raw output file from an allele calling to a phyloviz readable file. Use the Alleles.txt output file.


Basic usage:

`chewBBACA.py ExtractCgMLST -i rawDataToClean.tsv -o output_folders`
	
`-i` raw output file from an allele calling

`-o` output folder (created by the script if non existant)

`-r` (optional) list of genes to remove, one per line, advised to use the detected overrepresented genes from ParalogPrunning.py

`-g` (optional) list of genomes to remove, one per line, advised to use genomes that perform worst on statistics or use the result from the testGenomeQuality script

----------
## 5. All in one option 

### Run all chewBBACA analyses in a single command: 

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

`--genomes` path to folder with target genomes

`--cpu` Number of cpus to use (dont use all your cpu if using a laptop)

`-o` folder name with the analysis files

`-t` (Optional) taxon to use for prodigal training input

=============

----------
## FAQ

### Q: What is the fast way to run?  
A: Check step 1.1


### Q: Step 2 is taking hours, will it ever end?  
A: Depending on the variability of the strains used to create the schema and the number of CPUs you are using time used will vary. The more variable the strains, the more BLAST comparisons will be made, taking more time.

### Q: Step 3 just crashed at 99% after 2 days running, do I need to start over :(?  
A: chewBBACA should allow you to continue where you stopped, just re-run the same command and you should be prompted to continue the allele call.

### Q: I ran all the steps and my cgMLST loci size is smaller than traditional MLST, does this even work?  
A: Try to run step 4, your analysis may contain some genomes responsible for a considerable loss of loci. Remove some of those genomes and check if the cgMLST loci number rises.

### Q: Which species already have a training file?  
A: At the moment:
 - *Campylobacter jejuni* (Campylobacter_Jejuni)
 - *Acinetobacter baumannii* (Acinetobacter_Baumannii)
 - *Haemophilus influenzae* (Haemophilus_Influenzae)
 - *Streptococcus agalactiae* (Streptococcus_Agalactiae)
 - *Yersinia enterocolitica* (Yersinia_Enterocolitica)
 - *Enterococcus faecium* (Enterococcus_Faecium)
 - *Staphylococcus haemolyticus* (Staphylococcus_Haemolyticus)
 - *Salmonella enterica enteritidis* (Salmonella_Enterica_enteritidis)
 - *Staphylococcus aureus* (Staphylococcus_aureus)

----------

