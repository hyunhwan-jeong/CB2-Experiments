#!/bin/bash
cat ../02_sgRNA-level-analysis/dat/CRISPRn_readcount.txt | sed 's/Gene/gene/g' > input.tsv
#./wrapper/run_CRISPhieRmix.R -c pDNA -t RepA,RepB,RepC -o results/Sanson/CRISPRn-A375/FDR/CRISPhieRmix.csv -p "neg_prefix <-\"NO_CURRENT\"" input.tsv
./wrapper/run_CRISPhieRmix.R -c pDNA -t RepA,RepB,RepC -o results/Sanson/CRISPRn-A375/FDR/CRISPhieRmix.csv -p "neg_prefix <-\"NO_CURRENT\"" data/Sanson/CRISPRn-A375.tsv
cat ../02_sgRNA-level-analysis/dat/CRISPRi_readcount.txt | sed 's/Gene/gene/g' > input.tsv
#./wrapper/run_CRISPhieRmix.R -c pDNA -t RepA,RepB,RepC -o results/Sanson/CRISPRi-A375/FDR/CRISPhieRmix.csv -p "neg_prefix <-\"CONTROL\"" input.tsv
./wrapper/run_CRISPhieRmix.R -c pDNA -t RepA,RepB,RepC -o results/Sanson/CRISPRn-A375/FDR/CRISPhieRmix.csv -p "neg_prefix <-\"NO_CURRENT\"" data/Sanson/CRISPRi-A375.tsv
