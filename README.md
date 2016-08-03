#chewBBACA: BSR-Based Allele Calling Algorithm

**chewBBACA** stands for "BSR-Based Allele Calling Algorithm". The "chew" part could be thought as "Comprehensive and  Highly Efficient Workflow" but at this point still needs a bit of work to make that claim so we just add "chew" to add extra coolness to the software name. BSR stands for BLAST Score Ratio as proposed by  [Rasko DA et al.](http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-6-2) 

chewBBACA is a comprehensive pipeline for the creation and validation of whole genome and core genome MultiLocus Sequence Typing (wg/cgMLST) schemas, and also providing an allele calling algorithm, that can be run in multiprocessor settings.

This page was created with the dual purpose of serving as readme file and tutorial for the software.

##Dependencies:
* [git](https://git-scm.com/)
* [Biopython 1.66 ](http://biopython.org/wiki/Main_Page)
* [HTSeq 0.6.1p1](http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)
* BLAST 2.2.28+ ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.2.28/
* [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/)

##Tutorial contents
 
 1. Setting up the analysis  
 2. wgMLST schema creation  
 3. Allele call using the wgMLST schema
 4. Selecting a cgMLST schema from the wgMLST schema
 5. Validating the cgMLST schema
 6. Allele calling using the cgMLST schema

**Important Notes before starting:**

 - For **chewBBACA**, the definition of an allele is that each allele
   must represent a complete Coding Domain Sequence, with starting codon and stop codon according to the [NCBI genetic code table 11](http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi). It will automatically exclude any allele which DNA sequence that does not contain start or stop codon and that it's length is not multiple of three. The allele identification is performed using [Prodigal 2.6.0 ](https://github.com/hyattpd/prodigal/releases/). 
 - All the referenced lists of files *must contain full path* for the files.
 - Make sure that your fasta files are UNIX format. If they were created in Linux or MacOS systems they should be in the correct format, but if they were created in Windows systems , you should do a a quick conversion using for example [dos2unix](http://linuxcommand.org/man_pages/dos2unix1.html).

### 1. Setting up the analysis

**Installing chewBACCA**

To install chewBACCA, simply select the directory where you want it to be installed and run the following command:

```
 git clone https://github.com/mickaelsilva/chewBBACA.git
``` 

This will create the chewBBACA dir in the present directory with all the needed scripts. 

to update chewBBACA you simply 

```
 git pull
```

You will also need to install all the dependencies indicated above and make sure they are on the PATH.

**PATH settings for chewBBACA**

All the subdirectories in the chewBBACA dir  must be on the PATH if you want to be able to access the scripts from any dir. Replace the `<install_dir>` with the full path to chewBBACA instalation directory

```
export PATH="<install_dir>/chewBBACA/utils/:$PATH"
export PATH="<install_dir>/chewBBACA/allelecall/:$PATH"
export PATH="<install_dir>/chewBBACA/createschema/:$PATH"
export PATH="<install_dir>/SchemaEvaluator/:$PATH"
```
**Folder structure**
*This step is optional and is directed to users without much experience in running scripts.*
We suggest that  for each analysis for a given schema, chewBBACA should be run in a directory (chewbacca_wrkDIR) containing the following directory structure: 
```
	.../chewbacca_wrkDIR/
    .../chewbacca_wrkDIR/genomes 
    .../chewbacca_wrkDIR/schema
```
the `.../chewbacca_wrkDIR/genomes` dir will contain the fasta files with the genomes to be analysed (complete or draft genomes).
In `.../chewbacca_wrkDIR/schema ` will contain a fasta file with the alleles for each loci. This directory will also contain a subdir named ` .../chewbacca_wrkDIR/schema/short ` with the fasta file with all the alleles to be used in the BLAST step of the allele call. This files should have the  name "< gene >_short.fasta" with < gene > matching the filenames in `.../chewbacca_wrkDIR/schema`.
 
### 2. wgMLST schema creation

To create a wgMLST schema for a set of pre-determined genomes, you can run the following command:

`CreateSchema.py -i allffnfile.fasta -g 200`

The parameters for the command are: 

`-i` file with concatenated gene sequences

`-g` minimum locus lenght (removes any loci with length equal or less the specified value


**Input:**  `allffnfile.fasta` : a fasta file resulting from concatenating all the fasta files containing all the genes for the genomes we want to use for the creating the wgMLST schema. In Genebank this files usually have the *.fnn*  extension. A commonly used and highly efficient software for bacterial genome annotation ,[Prokka](https://github.com/tseemann/prokka/blob/master/README.md), also provide *.ffn* files with its output.

**Output:** One fasta file per gene in the `schema_seed/`directory that is created in the dir where the files are found. The fasta file names are the given according the FASTA annotation for each coding sequence. For example the locus with the annotation ` >gi|193804931|gb|AE005672.3|:2864-3112 Streptococcus pneumoniae TIGR4, complete genome` will create the fasta file named  `gi_193804931_gb_AE005672.3_:2864-3112.fasta`.

The script creates a selection of unique loci present in the input file. Firstly, it removes genes that are substring of larger genes (i.e. the CDS are identical but  are annotated with different start codons) and genes with DNA sequences smaller than chosen in the `-g` parameter. 
The second step is grouping all the genes by BLASTIng all the genes against each other. Pairwise gene comparisons with Blast Score Ratio (BSR) greater than 0.6 are considered alleles of the same locus and the allele with larger gene length is kept in the list.

Finally  in the resulting directory `schema_seed` ,  run the `init_bbaca_schema.sh` shellscript. This script will create the necessary files for the allele call, by creating the  directory named `short`. The contents of this dir is already explained in the **folder structure** subsection. It will also contain a  file named `listGenes.txt` , that will contain the list of loci in the schema.

### 3.  Allele call using the wgMLST schema 

Having defined the wgMLST schema with one allele per locus, one can proceed to use it to call alleles.  The command is the following

	% BBACA.py -i listGenomes.txt -g listGenes.txt -o OutPrefix --cpu 3 -p /<PathToProdigal>/prodigal

**Parameters** 
`-i` path to the list of genomes file. One file path (must be full path) to any fasta/multifasta file containing all the complete or draft genomes you want to call alleles for.

`-g` path to the list of alleles file i.e. the `listGenes.txt`in the `.../chewbacca_wrkDIR/schema/short` directory

`-o` prefix for the output files

`-p` Prodigal path to execution file (type prodigal if already on PATH) 

`--cpu` Number of cpus to use (if greater than existing cpus-2 uses cpus-2)

`-v`,`--verbose`  verbose mode(optional). Provides more output of the run.


This will use by defaults cpus-2 and can be called in a SLURM HPC by  srun  (**TODO@Mickael**: @Mickael provide command line example) .
If by some reason the process is interrupted (server crash, etc), running the command line with the same inputs  will display the option to resume the  allele call, avoiding lost computing time.   

**TODO@Mickael**: Change name: alleleCalling_ORFbased_protein_main3_local.py --> BBACA

The outputs files are:
```
< outPrefix >_< datestamp>/< outPrefix >_statistics.txt
< outPrefix >_< datestamp>/< outPrefix >_contigsInfo.txt
< outPrefix >_< datestamp>/< outPrefix >_Alleles.txt 
```
short example `< outPrefix >_statistics.txt` file:

```
Stats:	EXC	INF	LNF	LOT	PLOT	NIPL	ALM	ASM
NC_017162.fna	892	2319	1909	0	0	104	5	37	
NC_011586.fna	1563	1697	1809	0	0	116	6	75	
```

The column headers stand for.
* EXC - allele has exact match (100% identity)
* INF - infered allele with prodigal
* LNF - locus not found
* LOT - locus on the tip of the contig (partial match)
* PLOT - possible locus on the tip of the contig (CDS match is on the tip of the contig and other alleles for the loci are large enough to be considered LOT- to be manualy curated)
* NIPL - Non informative paralogous locus . Two or more BLAST matches with BSR>0.6 for the protein, indicate that the locus could have paralogous and should be removed from the analysis
* ALM - allele 20% larger than locus size mode (match CDS lenght> locus mode length + locus mode length * 0.2)
* ASM - allele 20% smaller than locus size mode (match CDS lenght < locus mode length - locus mode length * 0.2)

The ALM and ASM control for the fact that, for the majority of the loci the allele lengths are quite conserved. However some loci can have larger variation in allele length and those should be manually curated.


An abridged example of an output file for `< outPrefix >_< datestamp>/< outPrefix >_Alleles.txt`  with only two loci:

```
FILE	gi_126640115_ref_NC_009085.1_:1032446-1033294.fasta	gi_126640115_ref_NC_009085.1_:103903-104649.fasta	
NC_017162.fna	INF-2	LNF
NC_011586.fna	INF-3	LNF
NC_011595.fna	3	LNF
```

The first column has the filename for which the allele call was performed. The table headers have the filenames of files where the alleles are stored for each loci. (**TODO @JAC** : explain output result )
 
### 4. Selecting a cgMLST schema from the wgMLST schema 


1. Create a list .txt file containing one draft genome file per line with full paths (similar to 3.)

2. Run the allelecall script (local or cluster version) using the list files created at 3. and 4.

3. Run the `ParalogPrunning.py.` script using `contigsInfo.txt` output as input from step 5.

**TODO@JAC** : Reavaluate paralog prunnig after wgMLST schema is created? 

**TODO@Mickael**: whichRepeatedLoci.py -> ParalogPrunning.py

	% ParalogPrunning.py -i contigsInfo.txt

`-i` contigsInfo.txt file




* PC -Paralog Count - number of times a CDS on this locus as been found in another locus
* NDC - Non Determined Locus count -  Number of times the alleles the target locus was not exact match to an allele  (EXC) or a new inferred allele (INF). (i.e. locus can be a LNF,LOT, PLOT, NIPL, ASM, ALM). It helps to decide towards excluding the locus if this number is also large.  

The output file `asasdasdas`lists only loci where at least one paralogous was found. The file format is the following:

```
gene	PC	NDC
gi_22536185_ref_NC_004116.1_:c2045049-2043157.fasta	1	2
gi_406708523_ref_NC_018646.1_:c1944065-1941807.fasta	1	2

```

In this example the allele call was ran for 3 genomes.
Both locus presented had an exact match or an infered allele for one genome, while 2 genomes had issues or didn't have the locus. The CDS returned for the first locus is present in another locus, while the same happens for the second locus, from which we may clearly infer that a locus is being overrepresented by this two locus, since both are catching the same CDS.
**TODO@JAC: Explain better this section. E**

7. Run the Extract_cgAlleles.py using the outputs from 5. and 6.

**TODO@Mickael**Change name:  XpressGetCleanLoci4Phyloviz.py -> Extract_cgAlleles.py

### 5. Validating the cgMLST schema

### 5.1 Evaluate genome quality


Usefull to determine a core genome and remove genomes that may have technical issues. The algorithm description is the following:

1. For each allelic profile generated for a draft genome , let *nl* be the number of loci that are not present in the allelic profile but are present in 99%  or more of the remaining allelic profiles (97% if total number of genomes under 500 and 95% if under 200);
2. For and exclusion threshold (et) remove all allelic profiles that have nml greater than et. If no allelic profiles are removed, proceed to Step 4;
3. Return to Step 1.
4. The locus present in all the draft genomes for the remaining allelic profiles, are defined as the cgMLST schema for the exclusion threshold (et)

Usage:

testQualityGenomes3 -> TestGenomeQuality

	% testQualityGenomes3.py -i out.txt -n 12 -t 250
	
`-i` raw output file from an allele calling

`-n` maximum number of iterations, each iteration removes a set of genomes over the threshold and recalculates all variables

`-t` maximum threshold, will start at 5 increasing in a step of 5 until t

The output consists in a set of plots per iteration and a removedGenomes.txt file where its informed of which genomes are removed per threshold when it reaches a stable point (no more genomes are removed)

Example of an output can be seen [here] (http://i.imgur.com/uQDNNkb.png) . This examples uses an original set of 1042 genomes and a scheme of 5266 loci, using a parameter `-n` of 12 and `-t` of 300.



9. (optional) Use the testQualityGenomes3.py script to reach/analyze the core genome and re-do step 7


### 6. Allele calling using the cgMLST schema

##  6.1 Creating a clean allelic profile for PHYLOViZ 

The 
Clean a raw output file from an allele calling to a phyloviz readable file. Keep the locus with only Exact matches or new alleles found for all genomes.

Basic usage:

	% XpressGetCleanLoci4Phyloviz.py -i rawDataToClean.txt -g cleanedOutput.txt -r removeLocusList.txt
	
`-i` raw output file from an allele calling

`-g` path/name of the clean file

`-r` (optional) list of genes to remove, one per line, advised to use the detected overrepresented genes from whichRepeatedLoci.py

=============



**Notice :**
Previous versions used the full allele list to make a BLAST search over the genomes, increasing exponentially the time cost necessary to perform the allele call with the growing allele database. 
To address this issue the allele call script has been modified using now 2 files per locus. One file will store all allelic sequence forms found, while a "short" version of the allele will store only
new alleles with a significant difference to the closest allele (0.6 < BSR < 0.7). The short gene form will be used to perform the BLAST search, allowing a better time performance allele call, while not losing the wider diversity range search.
