## Scripts and data for the target gene identification benchmarking

Below three lines will be used to generate results of the CRISPR screening analysis methods, including CB<sup>2</sup>.

```
# This will generate results for Evers et al.'s data
snakemake --snakefile smk/Run-Evers.smk 
# This will generate results for Sanson et al.'s data
snakemake --snakefile smk/Run-Sanson.smk
# This will grap some caches from HitSelect and PinAPL-Py
snakemake --snakefile smk/Process-caches.smk --configfile smk/Process-caches.yaml
# This will generate results of the parameter turning experiments.
snakemake --snakefile smk/Run-Evers-param-tune.smk    
```

Once you have results, then you can create figures using scripts in `scripts` folder. You can execute the scripts in Rstudio, or run `Rscript scripts/<file-name-of-the-script>.R`. The following table describes the function of each script in the `scripts` directory.


File name             |  Executable  | Function
----------------------|--------------|-------------------------------------------------------------------------------------
convert_external.R    |  No          | The script to convert output files from HitSelect and PinAPL-Py (We couldn't add them since they only provide GUI to run the analysis)
draw-fig-S1.R         |  Yes         | The script to generate a figure of the ranking comparison (Figure S1)
draw-fig-S16.R        |  Yes         | The script to generate a figure of the parameter tuning (Figure S16)
draw-fig-S2-4-S8-9.R  |  Yes         | The script to generate AUC curves (Figure S2, S3, S4, S8, and S9)
draw-fig-S5.R         |  Yes         | The script to generate accuracy plots for Ever et al.'s data (Figure S5)
draw-fig-S6.R         |  Yes         | The script to generate a heatmap of FDRs of gene statistics (Figure S6)
draw-fig-S7.R         |  Yes         | The script to generate accuracy plots for Sanson's et al.'s data (Figure S5)
draw-fig1.R           |  Yes         | The script to generate the main figure (heatmaps + F1-score plots) (Figure 1)
draw-fig1_Evers.R     |  No          | The script to generate a part of Figure 1 (Ever et al.'s data)
draw-fig1_Sanson.R    |  No          | The script to generate a part of Figure 1 (Sanson et al.'s data)
draw-fig2-S10-S11.R   |  Yes         | The script to generate dot plots and boxplots of RPL5, COPS8, and RPL27 (Figure 2, S10, S11)
draw-upset.R          |  Yes         | The script to visualize how many essential genes were newly found by CB<sup>2</sup>
draw_auc.R            |  No          | The file contains a function to visualize an AUC curve.
generate_auc_output.R |  No          | The script to convert an output from the Snakemake pipeline to another proper format to visualize an AUC curve

Please read [RRA_experiment/README.md](RRA_experiment/README.md) to see details regarding the performance comparison of p-value combined methods (Fisher's method vs. alpha-RRA).
