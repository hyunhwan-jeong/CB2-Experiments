## Scripts and data for "CB<sup>2</sup> is more sensitive in target gene identification than existing methods" section.

Below three lines will generate results of the CRISPR screening analysis methods, including CB<sup>2</sup>.

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

Once you have results, then you can create figures using scripts in `scripts` folder.
