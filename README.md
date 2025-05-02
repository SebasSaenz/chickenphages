# Chickenphages

1.  Use SRA to download samples

```         
prefetch sample_id -O sra_temp --verify yes --max-size 50G

fasterq-dump sra_temp/sample_id --split-files --outdir out_dir --threads 1
```

2.  Quality control

    ```{python}

    # Create index
    bowtie2-build {refer} {index} --threads {THREADS}

    # Read QC and host DNA removal
    trim_galore --paired $WS/$RR1 $WS/$RR2 -j {THREADS}

    ```
