methods = ["CB2", "MAGeCK", "PBNPA", "sgRSEA", "RIGER", "RSA", "ScreenBEAM"]
datasets = ["CRISPRn-RT112", "CRISPRi-RT112", "CRISPRn-UMUC3"]

rule all:
    input:
        expand("results/Evers/{dataset}/AUC/{method}.csv", dataset=datasets, method = methods)

rule generate_auc_output:
    input:
        "results/Evers/{dataset}/FDR/{method}.csv"
    output:
        "results/Evers/{dataset}/AUC/{method}.csv"
    shell:
        "./scripts/generate_auc_output.R -e data/Evers/essential-genes.txt -n data/Evers/non-essential-genes.txt -o {output} {input}"

    
rule run_methods:
    input: 
        "data/Evers/{dataset}.csv"
    output:
        "results/Evers/{dataset}/FDR/{method}.csv"
    params:
        exec = "../wrapper/run_{method}.R"
    benchmark:
        "benchmarks_runtime/Evers/{dataset}_{method}.txt"
    shell:
        "./{params.exec} -c B1,B2,B3 -t A1,A2,A3 -o {output} {input}"


