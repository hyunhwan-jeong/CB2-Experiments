import glob
import os
DIR = "FASTQ"

FASTQ = [ "SRR5027845", "SRR5027846", "SRR5027849", "SRR5027850" ]

rule all:
    input:
        expand("AGO/{sample}_CC2.txt", sample=FASTQ) +
        expand("AGO/{sample}_mageck.txt", sample=FASTQ) 

rule run_CC2:
    input:
        library = "library/human_gecko_v2.tsv",
        fastq = os.path.join(DIR, "{sample}.fastq")
    output:
        "AGO/{sample}_CC2.txt"
    benchmark:
        "AGO/benchmark/{sample}_CC2.txt"
    shell:
        "node main_CC2.js {input.library} {input.fastq} > {output}"


rule run_mageck:
    input:
        library = "library/human_gecko_v2.tsv",
        fastq = os.path.join(DIR, "{sample}.fastq")
    output:
        "AGO/{sample}_mageck.txt"
    benchmark:
        "AGO/benchmark/{sample}_mageck.txt"
    shell:
        "node main_mageck.js {input.library} {input.fastq} > {output}"

