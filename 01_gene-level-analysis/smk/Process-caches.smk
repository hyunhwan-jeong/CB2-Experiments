methods = ["HitSelect", "PinAPL-Py"]
datasets = ["CRISPRn-RT112", "CRISPRi-RT112", "CRISPRn-UMUC3"]

rule all:
    input:
        expand("results/Evers/{dataset}/AUC/{method}.csv", dataset=datasets, method = methods)



rule convert_external:
    input:
        "results_from_gui/{dataset}/{method}.csv"
    output:
        "results/Evers/{dataset}/FDR/{method}.csv"
    params:
        method = "{method}"
    shell:
        "Rscript scripts/convert_external.R -i {input} -o {output} -m {params.method}"
 

rule generate_auc_output:
    input:
        "results/Evers/{dataset}/FDR/{method}.csv"
    output:
        "results/Evers/{dataset}/AUC/{method}.csv"
    shell:
        "./scripts/generate_auc_output.R -e data/Evers/essential-genes.txt -n data/Evers/non-essential-genes.txt -o {output} {input}"

   
