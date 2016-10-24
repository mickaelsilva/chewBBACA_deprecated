#Validate Schema

Analyze your alleles for a set of parameters, taking special consideration on the allele CDS and the allele sizes per gene/locus

install python dependencies with:

```
pip install -r requirements.txt
```

Python dependencies:
* numpy
* [Biopython 1.66 ](http://biopython.org/wiki/Main_Page)
* [HTSeq 0.6.1p1](http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html)
* [matplotlib](http://matplotlib.org/)
* [mpld3](http://mpld3.github.io/)

Other dependencies:
* [mafft](http://mafft.cbrc.jp/alignment/software/linux.html)
* [clustalw2](http://www.clustal.org/download/2.1/)


#How to run

Main script that calls alleSizeStats.py and CheckCDS.py creating an html file for an easier results reading.

The optional parameters refer to the analysis of the alleles size. This analysis will calculate the mode size per gene and using that value -/+ a threshold (0.05 default) it will consider the allele as being a good allele if it is still considered within the threshold. It is given to the user the choice of the threshold and the choice to consider that a gene is considered "Conserved" if only one of the alleles is outside the threshold (default) or all the alleles must be within the threshold to be considered a "Conserved allele.

	% ValidateSchema.py -i /home/Pneumoniae_Genes/ -ta 11 -l /home/user/results/ratemyschema.html
	
`-i` directory where the genes .fasta files are located or alternatively a .txt file containing the full path for each gene .fasta file per line

`-ta` which translation table to use ( 11 in case of bacteria)

`-l` Location/name of the final html output

`--cpu` number of cpu to use, will be used for mafft and clustal

`--log` (optional) number of alleles per locus plot will be ploted in log10 scale

`-p` (optional) True if all alleles must be within threshold (default=False)

`-t` (optional) threshold used to calculate the range at which the allele is considered good (default=0.05)

`-s` (optional) split the boxplot in subsets of (usefull for >1000 genes scheme)
