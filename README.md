# wgMLST
tools developed to perform a bacterial wgMLST

Dependencies:
* [biopython] (http://biopython.org/wiki/Main_Page)
* [HTSeq] (http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)
* BLAST
* [Prodigal] (https://github.com/hyattpd/prodigal/releases/) (tested with v. 2.6.0)


**Notice :**

Previous versions used the full allele list to make a BLAST search over the genomes, increasing exponentially the time cost necessary to perform the allele call with the growing allele database. 
To address this issue the allele call script has been modified using now 2 files per locus. One file will store all allelic sequence forms found, while a "short" version of the allele will store only
new alleles with a significant difference to the closest allele (0.6<BSR<0.7). The short gene form will be used to perform the BLAST search, allowing a better time performance allele call, while not losing
the wider diversity range search.

**Important Notes :**
- lists of files MUST contain FULL PATH!
- be sure your fasta files are formated in UNIX, for quick conversion use [dos2unix] (http://linuxcommand.org/man_pages/dos2unix1.html)
- allele sequence of the gene files must represent a complete Coding Domain Sequence, with starting codon and stop codon according to the [NCBI table 11] (http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi)
- gene folder must have a subfolder cointaning duplicated gene files for the short form

suggested folder structure:

1. main folder - scripts and txt files
 1. sub folder - genomes - all genomes fasta files
 2. sub folder - genes - all genes fasta files
  1. sub folder - short - all genes short form, extension ends with "_short.fasta"


How to perform a complete wgMLST:

1. Concatenate gene sequences in a single fasta file. You can use the .ffn files as source of genes sequences available [here] (http://ftp.ncbi.nih.gov/genomes/archive/old_genbank/Bacteria/)
2. Run CreateSchema.py over the concatenated single fasta file
3. Go to the schema_seed and run the init_bacca_schema.sh. Use the listGenes.txt for the allele call.
4. Create a list .txt file containing one draft genome file per line with full paths (similar to 3.)
5. Run the allelecall script (local or cluster version) using the list files created at 3. and 4.
6. Run the whichRepeatedLoci.py over the contigsInfo.txt output from step 5.
7. Run the XpressGetCleanLoci4Phyloviz.py using the outputs from 5. and 6.
8. Use [phyloviz software] (http://www.phyloviz.net/) to build trees based on the profiles generated or the [online tool] (https://node.phyloviz.net/)
9. (optional) Use the testQualityGenomes3.py script to reach/analyze the core genome and re-do step 7

=============
#Create Schema

dependencies:
* biopython
* HTSeq
* BLAST

Given a concatenated ffn file, removes genes that are substring of bigger genes and genes with DNA sequences smaller than chosen in the -g parameter. Blasts all the genes against each other and saves the bigger genes, removing the smaller genes with a 0.6>BSR (BSR calculated according to the original [paper](https://peerj.com/articles/332/) )

	% CreateSchema.py -i allffnfile.fasta -g 200

`-i` file with concatenated gene sequences

`-g` minimum DNA sequence size


Output:

* x.fasta large set of .fasta files, 1 file per gene inside the folder schema_seed



=============
##Allele Call


`alleleCalling_ORFbased_protein_main3_local.py` - short version of ORF based allele call to be run in a SLURM grid based cluster or a local machine. Uses 2 files per gene.



Performing a **SLURM cluster** allele call short version (uses the maximum number of local cores minus 2, one per gene):

	% srun -c 32 --mem-per-cpu 4G  alleleCalling_ORFbased_protein_main3_local.py -i listGenomes.txt -g listGenes.txt -o outputFileName.txt -p /home/user/prodigal/Prodigal-2.60/prodigal


`-c` number of CPU (for SLURM runs only)

`--mem-per-cpu` memory available per CPU in Gb (for SLURM runs only)
	
`-i` path to the list of genomes file

`-g` path to the list of alleles file

`-o` output file name

`-p` Prodigal path to execution file (included) 

short example statistics file:

* EXC - allele has exact match (100% identity)
* INF - infered allele with prodigal
* LNF - locus not found
* LOT - locus on the tip of the contig (partial match)
* PLOT - locus possibly on the tip of the contig (CDS match is on the tip of the contig - to be manualy curated)
* NIPL - Non informative paralog locus (two or more good blast matches (bsr >0.6) for the protein)
* ALM - allele much larger than gene size mode (match CDS lenght> gene mode length + gene mode length * 0.2)
* ASM - allele much smaller than gene size mode (match CDS lenght < gene mode length - gene mode length * 0.2)

```
Stats:	EXC	INF	LNF	LOT	PLOT	NIPL	ALM	ASM
NC_017162.fna	892	2319	1909	0	0	104	5	37	
NC_011586.fna	1563	1697	1809	0	0	116	6	75	
```

short example file output:

```
FILE	gi_126640115_ref_NC_009085.1_:1032446-1033294.fasta	gi_126640115_ref_NC_009085.1_:103903-104649.fasta	gi_126640115_ref_NC_009085.1_:1056402-1057004.fasta	gi_12664011510_S10_L001.fasta
NC_017162.fna	INF-2	LNF
NC_011586.fna	INF-3	LNF
NC_011595.fna	3	LNF
```
=============
## Evaluate overrepresented loci

Using the contigsInfo.txt output from the allele call, check if the same CDS is being called for different locus

	% whichRepeatedLoci.py -i contigsInfo.txt

`-i` contigsInfo.txt file

short example file output:

* overrepresented - number of times a CDS on this locus as been found in another locus
* problems - neither exact match or infered allele found

```
gene	overrepresented	problems	total
gi_22536185_ref_NC_004116.1_:c2045049-2043157.fasta	1	2	3
gi_406708523_ref_NC_018646.1_:c1944065-1941807.fasta	1	2	3

```
In this example the allele call was ran for 3 genomes.
Both locus presented had an exact match or an infered allele for one genome, while 2 genomes had issues or didn't have the locus. The CDS returned for the first locus is present in another locus, while the same happens for the second locus, from which we may clearly infer that a locus is being overrepresented by this two locus, since both are catching the same CDS.

=============
## Clean raw output profiles from none exact matches/new alleles


Clean a raw output file from an allele calling to a phyloviz readable file. Keep the locus with only Exact matches or new alleles found for all genomes.

Basic usage:

	% XpressGetCleanLoci4Phyloviz.py -i rawDataToClean.txt -g cleanedOutput.txt -r removeLocusList.txt
	
`-i` raw output file from an allele calling

`-g` name of the clean file

`-r` (optional) list of genes to remove, one per line, advised to use the detected overrepresented genes from whichRepeatedLoci.py

=============
## Evaluate genome quality

Dependencies:
* [matplotlib](http://matplotlib.org/)

Usefull to determine a core genome and remove genomes that may have technical issues. The algorithm description is the following:

1. For each allelic profile generated for a draft genome , let nl be the number of loci that are not present in the allelic profile but are present in 99% (97% if total number of genomes under 500 and 95% if under 200) or more of the remaining allelic profiles;
2. For and exclusion threshold (et) remove all allelic profiles that have nml greater than et. If no allelic profiles are removed, proceed to Step 4;
3. Return to Step 1.
4. The locus present in all the draft genomes for the remaining allelic profiles, are defined as the cgMLST schema for the exclusion threshold (et)

Usage:

	% testQualityGenomes3.py -i out.txt -n 12 -t 250
	
`-i` raw output file from an allele calling

`-n` maximum number of iterations, each iteration removes a set of genomes over the threshold and recalculates all variables

`-t` maximum threshold, will start at 5 increasing in a step of 5 until t

The output consists in a set of plots per iteration and a removedGenomes.txt file where its informed of which genomes are removed per threshold when it reaches a stable point (no more genomes are removed)

Example of an output can be seen [here] (http://i.imgur.com/uQDNNkb.png) . This examples uses an original set of 1042 genomes and a scheme of 5266 loci, using a parameter `-n` of 12 and `-t` of 300.
