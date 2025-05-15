# Diversity of phages in the chicken gut

This repository contains the supporting code for the manuscript **"A metagenomic collection of 19,778 viruses reveals the diverse virome of the chicken gut"**

Johan S. Sáenz^1,2,^, Timur Yergaliyev^1,2^,  Bibiana Rios-Galicia^1,2^,  Jana Seifert^1,2^  & Amélia Camarinha-Silva^1,2^

^1^Institute of Animal Science, University of Hohenheim, Emil-Wolff-Str. 6-10, 70593 Stuttgart, Germany

^2^HoLMiR—Hohenheim Center for Livestock Microbiome Research, University of Hohenheim, Leonore-Blosser-Reisen Weg 3, 70593 Stuttgart, Germany

## Data availability

The collection of 47,092 viral genomes and the 19,778 dereplicated species-level vOTUs can be found in Zenodo, <https://zenodo.org/records/15063554>. Metagenomic samples used in this study are publicly available and can be downloaded using the provided Bioproject IDs in data folder.

## Sample summary

## Workflow

### Table of Contents

1.  [Download metagenomic samples](#download-metagenomic-samples)
2.  [Quality control and host removal](#quality-control-and-host-removal)
3.  [Assembly and filtering](#assembly-and-filtering)
4.  [Virus identification and quality](#virus-identification-and-quality)
5.  [High-quality viral genomes and vOTUs](#high-quality-viral-genomes-and-vOTUs)
6.  [Cluster family and genus](#cluster-family-and-genus)
7.  [Taxonomic annotation](#taxonomic-annotation)
8.  [Predict host and lifestyle](#predict-host-and-lifestyle)
9.  [Map reads](#map-reads)
10. [Functional annotation](#functional-annotation)

------------------------------------------------------------------------

## Download metagenomic samples {#download-metagenomic-samples}

```         
prefetch sample_id -O sra_temp --verify yes --max-size 50G

fasterq-dump sra_temp/sample_id --split-files --outdir out_dir --threads 1
```

## Quality control and host removal {#quality-control-and-host-removal}

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

## Assembly and filtering {#assembly-and-filtering}

```         
#
megahit -1 $CR1 -2 $CR2 -o $MH/$name -t {THREADS} --presets meta-sensitive

# Filter contigs shorter than 5kb
bbduk.sh in=%s out=%s minlen=5000
```

## Virus identification and quality {#virus-identification-and-quality}

```         
genomad end-to-end --cleanup --conservative %s %s %s --threads 16

checkv end_to_end %s %s -t 18 -d %s
```

## High-quality viral genomes and vOTUs {#high-quality-viral-genomes-and-votus}

```         
# Filter Medium and High-quality genomes
seqkit grep -n -f %s %s > %s

#Cluster vOTUs
makeblastdb -in %s -dbtype nucl -out %s

blastn -query %s -db %s -outfmt '6 std qlen slen' -max_target_seqs 10000 -out %s -num_threads 20

anicalc.py -i %s -o %s

aniclust.py --fna %s --ani %s --out %s --min_ani 95 --min_tcov 85 --min_qcov 0

seqkit grep -n -f %s %s > %s
```

## Cluster family and genus {#cluster-family-and-genus}

```         
# Create a blast+ databasa using your viral contigs
diamond makedb --in %s --db %s --threads 10

# Compute AAI from BLAST results
diamond blastp --query %s \
               --db %s \
               --out %s.tsv \
               --outfmt 6 \
               --threads 20 \
               --evalue 1e-5 \
               --max-target-seqs 10000 \
               --query-cover 50 \
               --subject-cover 50

/beegfs/work/workspace/ws/ho_kezau83-phage_db-0/amino_acid_identity.py --in_faa %s --in_blast %s --out_tsv %s

# Filter edges and prepare MCL inputs GENUS

python /beegfs/work/workspace/ws/ho_kezau83-phage_db-0/filter_aai.py --in_aai %s --min_percent_shared 20 --min_num_shared 16 --min_aai 40 --out_tsv %s

# Perform MCL-based clustering

mcl %s -te 8 -I 2.0 --abc -o %s
```

## Taxonomic annotation {#taxonomic-annotation}

```         
genomad annotate --cleanup --lenient-taxonomy %s %s %s --threads 16
```

## Predict host and lifestyle {#predict-host-and-lifestyle}

```         
# Split dataset so you can run faster
iphop split --input_file %s --split_dir %s

# Run iPhop
iphop predict -t 20 --fa_file %s --db_dir %s --out_dir %s

bacphlip -i %s --multi_fasta
```

## Map reads {#map-reads}

```         
TMPDIR=. coverm contig --methods rpkm -t 14 -c %s %s -r %s -o %s --min-read-percent-identity 90 --min-read-aligned-percent 75 --min-covered-fraction 75
```

## Functional annotation {#functional-annotation}

```         
pharokka.py -i %s -o %s -d %s -g prodigal-gv -m -t 14

defense-finder run %s -o %s -a
```
