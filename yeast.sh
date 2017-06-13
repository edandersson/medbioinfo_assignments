#!/bin/bash
# Extract genes
genes=$(grep -h -e '^>' -r yeast_genes/ | cut -c '2-5' | sort | uniq | tr '\n' ':') # Introduce ":" for later split
# Assess lenght of genes, since not a vector
genes_char=$(echo $genes|wc -c)
genes_char=$(expr $genes_char - 1)
genes_num=$(expr $genes_char / 5)
index=$(seq 1 $genes_num)
for i in $index
do
gene=$(echo $genes | cut -d : -f "$i") # Cut gene name from string
echo $i
touch "$gene.txt"
# Merge all samples and extract data by gene and allign with MAFFT
cat yeast_genes/* | sed -n "/^>$gene*/,/^>*/p" | mafft /dev/stdin > "$gene.txt"
done
