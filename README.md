# CB2-Experiments

`CB2-Experiments` is the repository of scripts and data of experiments were performed to benchmark the performance of CB<sup>2</sup>.

## Structure of the repository

`CB2-Experiments` contains four different directories. Each directory contains `README.md` to help running scripts.

| Folder name               | Contents                                                                                                                           |
|---------------------------|------------------------------------------------------------------------------------------------------------------------------------|
| `01_gene-level-analysis`  | Contains scripts and data used for "CB<sup>2</sup> is more sensitive in target gene identification than existing methods" section. |
| `02_sgRNA-level-analysis` | Contains scripts and data used for "CB<sup>2</sup>" is more specific in target gene detection than existing methods" section.      |
| `03_quantification`       | Contains scripts used for "CB <sup>2</sup> provides more accurate alignment without parameter tuning" section.                     |
| `util`                    | A collection of utility scripts.                                                                                                   |

## Requirements

Below bash command line enumerate the list of required/used R packages to use the script in `CB2-Experiments`.

```
grep -e "^library" * -R | sed 's/:/ /g' | awk '{ print $2 }' | sort -u | sed 's/library(//g' | sed 's/)//g'
CB2
DESeq2
PBNPA
RColorBrewer
ScreenBEAM
argparse
cowplot
edgeR
eulerr
gghighlight
ggsci
glue
here
pheatmap
precrec
sgRSEA
tidyverse
```


