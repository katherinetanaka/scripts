#! /bin/bash

# A short, specialized script to map my contigs to the pAsa5 plasmid
# Required that I lastdb the pAsa5_propre database first
# Gunzip the contig fasta file first

# Pass args, only targetted directory
MYDIR=$1
MYFILE=$2

lastal pAsa5_propre $MYFILE > $MYDIR/rawAL.maf # align
cd $MYDIR
maf-sort rawAl.maf > sort.maf # Sort by pAsa5 positions

# Keep only alignments greater than 2000 bp --> remove some noise caused by rep. regions
# Only print coords in pAsa5 and in contig
awk 'BEGIN {ORS=""} $2 == "pAsa5" && $4 > 2000 {print id, $2, $3, $4, $5, $8;getline;
print $2, $3, $4, $5, $6, "\n"}' id=$MYDIR sort.maf > pAsa5.tab

cd ..

