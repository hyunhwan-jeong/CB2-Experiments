## Scripts of the sgRNA-level quantification experiments.

### Data download

* Please create `FASTQ` folder and download `FASTQ` files listed below to the `FASTQ`. If any of `FASTQ` files are gziped, then please `ungzip`.

#### Ever et al.'s data

Evers et al.'s data is available at https://drive.google.com/drive/folders/1x8in4AQh-BcGQpYyldIC14v6ZEnC0pOk

#### Golden et al.'s data
You can execte the following code snippet in your terminal to download Golden's et al.' data:

```
cd FASTQ/
fastq-dump SRR5027845
fastq-dump SRR5027846
fastq-dump SRR5027849
fastq-dump SRR5027850
```
#### Koike-Yusa et al.'s data

You can execte the following code snippet in your terminal to download Golden's et al.' data:

```
cd FASTQ/
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR376/ERR376998/ERR376998.fastq.gz
gunzip ERR376998.fastq.gz
```


### A notebook to create figure 5B.

[worknote/notebook.Rmd](worknote/notebook.Rmd) contains how to create the figure 5A.

### The jupyter notebook to create figure 5B.

[Notebook-for-figure-5B.ipynb](Notebook-for-figure-5B.ipynb) contains codes to generate the figure 5B.
