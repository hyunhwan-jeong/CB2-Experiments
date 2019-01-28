library(argparse)
library(tidyverse)
parser <- ArgumentParser()
parser$add_argument("infile", nargs=1, help="File to be analyzed.")
parser$add_argument("-t", "--trt", help="The list of sample names in the treatment group.",  required=T)
parser$add_argument("-c", "--ctl", help="The list of sample names in the control group. ",  required=T)
parser$add_argument("-o", "--outfile", help="File Path to the output",  required=T)
parser$add_argument("-p", "--params", nargs=1, help="Additional parameters can be used in the method.")
args <- parser$parse_args()
infile <<- args$infile
outfile <<- args$outfile
ctl <<- args$ctl %>% strsplit(split = ',') %>% .[[1]]
trt <<- args$trt %>% strsplit(split = ',') %>% .[[1]]
params <<- args$params
if(is.null(params)) params <<- ""

if( file.access(infile) == -1) {
    stop(sprintf("Specified input file ( %s ) does not exist", infile))
    quit(1)
}

if(infile %>% tolower() %>% endsWith(".csv")) {
    df_count <<- read_csv(infile)
} else {
    df_count <<- read_delim(infile, delim="\t")
}


