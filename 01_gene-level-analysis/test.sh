#!/bin/bash
./wrapper/run_$1.R -c B1,B2,B3 -t A1,A2,A3 -o output.txt data/Evers/CRISPRn-RT112.csv
head output.txt
rm output.txt
