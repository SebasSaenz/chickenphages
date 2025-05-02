# Chickenphages

## Table of Contents

1.  [Download metagenomic samples](#download-metagenomic-samples)
2.  [Quality control and host removal](#quality-control-and-host-removal)
3.  [Assembly and filtering](#assembly-and-filtering)
4.  [Virus identification and quality](#virus-identification-and-quality)
5.  [High-quality viral genomes](#high-quality-viral-genomes)

------------------------------------------------------------------------

## Download metagenomic samples

```         
prefetch sample_id -O sra_temp --verify yes --max-size 50G

fasterq-dump sra_temp/sample_id --split-files --outdir out_dir --threads 1
```

## Quality control and host removal 

```         
# Create a QC report
fastqc -q -t {THREADS} -o $preQC -f fastq $RR1 $RR2

# Create index
bowtie2-build {refer} {index} --threads {THREADS}

# Read QC and host DNA removal
trim_galore --paired $WS/$RR1 $WS/$RR2 -j {THREADS}

# Remove host DNA
bowtie2 -p {THREADS} -x $ind_host -1 $WS/$outdir/*_val_1.f*q.gz -2 $WS/$outdir/*_val_2.f*q.gz \--un-conc-gz $WS/$outdir/no_host > $WS/$outdir/host.sam
```

## Assembly and filtering

```         
#
megahit -1 $CR1 -2 $CR2 -o $MH/$name -t {THREADS} --presets meta-sensitive

# Filter contigs shorter than 5kb
bbduk.sh in=%s out=%s minlen=5000
```

## Virus identification and quality

```         
genomad end-to-end --cleanup --conservative %s %s %s --threads 16

checkv end_to_end %s %s -t 18 -d %s
```

## High-quality viral genomes

```         
seqkit grep -n -f %s %s > %s
```

```         
```
