import glob
import os
DIR = "FASTQ"
LIB = "library/nbt3536-S3.tsv"

FASTQ = [ "RT112_T0_R1"]

rule all:
    input:
        expand("Evers/{sample}_CC2.txt", sample=FASTQ) +
        expand("Evers/{sample}_mageck.txt", sample=FASTQ) +

rule run_CC2:
    input:
        library = LIB,
        fastq = os.path.join(DIR, "{sample}.fastq")
    output:
        "Evers/{sample}_CC2.txt"
    benchmark:
        "Evers/benchmark/{sample}_CC2.txt"
    shell:
        "node js/main_CC2.js {input.library} {input.fastq} > {output}"


rule run_mageck:
    input:
        library = LIB,
        fastq = os.path.join(DIR, "{sample}.fastq")
    output:
        "Evers/{sample}_mageck.txt"
    benchmark:
        "Evers/benchmark/{sample}_mageck.txt"
    shell:
        "node js/main_mageck.js {input.library} {input.fastq} > {output}"



