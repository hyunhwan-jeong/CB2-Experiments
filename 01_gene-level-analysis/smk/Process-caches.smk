outputs = expand("results/{study}/{dataset}/AUC/{method}.csv", study = config['study'], dataset=config['datasets'], method = config['methods'])
print(outputs)

rule all:
    input:
        outputs



rule convert_external:
    input:
        "results_from_gui/{study}/{dataset}/{method}.csv"
    output:
        "results/{study}/{dataset}/FDR/{method}.csv"
    params:
        method = "{method}"
    shell:
        "Rscript scripts/convert_external.R -i {input} -o {output} -m {params.method}"
 

rule generate_auc_output:
    input:
        "results/{study}/{dataset}/FDR/{method}.csv"
    output:
        "results/{study}/{dataset}/AUC/{method}.csv"
    params:
        study = "{study}"
    shell:
        "./scripts/generate_auc_output.R -e data/{params.study}/essential-genes.txt -n data/{params.study}/non-essential-genes.txt -o {output} {input}"

  
