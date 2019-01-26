#!/bin/bash
#./wrapper/run_$1.R -c B1,B2,B3 -t A1,A2,A3 -o output.txt data/Evers/CRISPRn-RT112.csv
./wrapper/run_$1.R -c pDNA -t RepA,RepB,RepC -o output.txt data/Sanson/CRISPRi-A375.tsv
./scripts/generate_auc_output.R -e data/Sanson/essential-genes.txt -n data/Sanson/non-essential-genes.txt  -o output2.txt output.txt
#./scripts/generate_auc_output.R -e data/Evers/essential-genes.txt -n data/Evers/non-essential-genes.txt  -o output2.txt output.txt
head output2.txt
#rm output.txt
