methods = ["CB2", "MAGeCK", "sgRSEA", "RIGER", "RSA", "PBNPA"]
datasets = ["CRISPRi-A375", "CRISPRn-A375"]



rule all:
    input:
        expand("results/Sanson/{dataset}/AUC/{method}.csv", dataset=datasets, method = methods) 


rule generate_auc_output:
    input:
        "results/Sanson/{dataset}/FDR/{method}.csv"
    output:
        "results/Sanson/{dataset}/AUC/{method}.csv"
    shell:
        "./scripts/generate_auc_output.R -e data/Sanson/essential-genes.txt -n data/Sanson/non-essential-genes.txt -o {output} {input}"

    
rule run_methods:
    input: 
        "data/Sanson/{dataset}.tsv"
    output:
        "results/Sanson/{dataset}/FDR/{method}.csv"
    params:
        exec = "wrapper/run_{method}.R"
    shell:
        "./{params.exec} -c pDNA -t RepA,RepB,RepC -o {output} {input}"


