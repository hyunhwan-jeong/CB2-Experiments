# CB2-Experiments

`CB2-Experiments` is the repository of scripts and data of experiments were performed to benchmark the performance of CB<sup>2</sup>.

## Structure of the repository

`CB2-Experiments` contains four different directories. Each directory contains `README.md` to help running scripts.

| Folder name               | Contents                                                                                                                           |
|---------------------------|------------------------------------------------------------------------------------------------------------------------------------|
| `01_gene-level-analysis`  | Contains scripts and data used for "CB<sup>2</sup> is more sensitive in target gene identification than existing methods" section. |
| `02_sgRNA-level-analysis` | Contains scripts and data used for "CB<sup>2</sup>" is more specific in target gene detection than existing methods" section.      |
| `03_quantification`       | Contains scripts used for "CB<sup>2</sup> provides more accurate alignment without parameter tuning" section.                     |
| `util`                    | A collection of utility scripts.                                                                                                   |

## Requirements

Below bash command line enumerate the list of required/used R packages to use the script in `CB2-Experiments`.

```
$ grep -e "^library" * -R | sed 's/:/ /g' | awk '{ print $2 }' | sort -u | sed 's/library(//g' | sed 's/)//g'
CB2
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

If you don't have some of the packages than please use below installation code in `R`:

```
# Install packages in `cran`
install.packages(c("PBNPA", "RColorBrewer", "argparse", "cowplot", "eulerr", "gghighlight", "ggsci", "glue", "here", "pheatmap", "precrec", "sgRSEA", "tidyverse"))

# Install packages `devtools`
install.packages("devtools")
devtools::install_github("hyunhwaj/CB2")
devtools::install_github("jyyu/ScreenBEAM")
```

Since many of scripts written in `Snakemake`, so it is necessary to install the `Snakemake` package.

```
pip install snakemake --user
```

If you plan to run wrapper scripts of `RIGER` and `RSA`, you have to download and install the programs and have to specify the location and version in the wrapper scripts.

* `RIGER` is available to download at https://github.com/broadinstitute/rigerj.
* `RSA` is available to download at https://admin-ext.gnf.org/publications/RSA/.

Please use the below command if you haven't install `MAGeCK` and want to run any scripts related to the method:

```
pip install mageck --user
```

In `01_gene-level-analysis/wrapper/run_RIGER.R`, please change line 20 and 21 like below:

```
RIGER_PATH <- "/Users/hyunhwan/Projects/InProgress/rigerj" # path of `RIGER`
RIGER_VERSION <- "2.0.2" # path of the version of `RIGER`
```

In `01_gene-level-analysis/wrapper/run_RSA`, please change line 21 and 22 like below:

```
PY3_PATH <- "/Users/hyunhwan/miniconda3/bin/python3" # path of `python3`, can lookup using `which python3`
RSA_PATH <- "/Users/hyunhwan/Projects/InProgress/RSA" # path of `RSA`
  
  
## Questions

If you have any questions or issues with the respository, please add issue or send an email to `hyunhwaj@bcm.edu`.
