#!/usr/bin/env bash
mkdir short
cp -r *.fasta short/
cd short/
for file in *; do mv "$file" "`basename $file .fasta`_short.fasta"; done
cd ..
ls *.fasta | xargs -I {} sh -c 'echo $(pwd)/{}' > listGenes.txt
