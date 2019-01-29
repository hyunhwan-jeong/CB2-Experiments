input_data = "data/Evers/CRISPRn-RT112.csv"

# setting for mageck
param_mageck = [100, 1000, 10000, 100000]

# setting for ScreenBEAM
param_burnin = [50, 500, 5000, 10000]
param_nitt = [15000]

# setting for PBNPA
param_PBNPA = [10, 50, 100, 500]

# setting for sgRSEA
param_sgRSEA = [10, 50, 100, 500]


# setting for RIGER
param_RIGER = [0.1, 0.5, 1.0, 1.5, 2.0]

# setting for RSA


rule all:
    input:
        expand('results_params_tune/mageck-{param_mageck}.txt', param_mageck = param_mageck) +
        expand('results_params_tune/ScreenBEAM-{param_burnin}-{param_nitt}.txt', param_burnin = param_burnin, param_nitt = param_nitt) +
        expand('results_params_tune/PBNPA-{param_PBNPA}.txt', param_PBNPA = param_PBNPA) +
        expand('results_params_tune/sgRSEA-{param_sgRSEA}.txt', param_sgRSEA = param_sgRSEA) +
        expand('results_params_tune/RIGER-{param_RIGER}.txt', param_RIGER = param_RIGER) 






rule run_mageck:
    input:
        input_data
    output:
        "results_params_tune/mageck-{param_mageck}.txt"
    shell:
        './wrapper/run_MAGeCK.R -c B1,B2,B3 -t A1,A2,A3 -o {output} -p "--adjust-method fdr --norm-method median --additional-rra-parameters \\"--permutation {wildcards.param_mageck}\\"" {input}'



rule run_ScreenBEAM:
    input:
        input_data
    output:
        'results_params_tune/ScreenBEAM-{param_burnin}-{param_nitt}.txt'
    shell:
        './wrapper/run_ScreenBEAM.R -c B1,B2,B3 -t A1,A2,A3 -o {output} -p "burnin <- {wildcards.param_burnin}; nitt <- {wildcards.param_nitt}" {input}'


rule run_PBNPA:
    input:
        input_data
    output:
        'results_params_tune/PBNPA-{param_PBNPA}.txt'
    shell:
        './wrapper/run_PBNPA.R -c B1,B2,B3 -t A1,A2,A3 -o {output} -p "sim.no <- {wildcards.param_PBNPA}" {input}'



rule run_sgRSEA:
    input:
        input_data
    output:
        'results_params_tune/sgRSEA-{param_sgRSEA}.txt'
    shell:
        './wrapper/run_sgRSEA.R -c B1,B2,B3 -t A1,A2,A3 -o {output} -p "multiplier <- {wildcards.param_sgRSEA}" {input}'


rule run_RIGER:
    input:
        input_data
    output:
        'results_params_tune/RIGER-{param_RIGER}.txt'
    shell:
        './wrapper/run_RIGER.R -c B1,B2,B3 -t A1,A2,A3 -o {output} -p "-alpha {wildcards.param_RIGER}" {input}'


