library(argparse)
library(tidyverse)
parser <- ArgumentParser()
parser$add_argument("-i", "--infile", help="The original file.",  required=T)
parser$add_argument("-o", "--outfile", help="File Path to the output.",  required=T)
parser$add_argument("-m", "--method", help="The analysis method.",  required=T)
args <- parser$parse_args()
infile <- args$infile
outfile <- args$outfile
method <- args$method

if( file.access(infile) == -1) {
    stop(sprintf("Specified input file ( %s ) does not exist", infile))
}

if(method == "HitSelect") {
    read_csv(infile) %>% 
        select(gene = gene, fdr = fdr, stat = rank) %>% 
        write_csv(outfile)
} else if(method == "PinAPL-Py") {
    # read_csv(infile) %>% print()
    read_csv(infile) %>% select(gene = gene, stat = `p-value combined`) %>%
        mutate(fdr = p.adjust(stat, method="fdr")) %>% 
        write_csv(outfile)
} else {
    stop("The method parameter was not properly set.")
}
