# Diversity of phages in the chicken gut

This repository contains the supporting code for the manuscript **"A metagenomic collection of 19,778 viruses reveals the diverse virome of the chicken gut"**

Johan S. Sáenz^1,2^, Timur Yergaliyev^1,2^, Bibiana Rios-Galicia^1,2^, Jana Seifert^1,2^ & Amélia Camarinha-Silva^1,2^

^1^Institute of Animal Science, University of Hohenheim, Emil-Wolff-Str. 6-10, 70593 Stuttgart, Germany

^2^HoLMiR—Hohenheim Center for Livestock Microbiome Research, University of Hohenheim, Leonore-Blosser-Reisen Weg 3, 70593 Stuttgart, Germany

## Data availability

The collection of 47,092 viral genomes and the 19,778 dereplicated species-level vOTUs can be found in Zenodo, <https://zenodo.org/records/15063554>. Metagenomic samples used in this study are publicly available and can be downloaded using the provided Bioproject IDs in data folder.

## Sample summary

In total we collected 1,514 chicken gut metagenomic samples from 15 different countries, including samples from faeces and 6 gut regions.

![](plots_repo/summary_samples.png)

## Workflow

Data was analyse using the High Performance and Cloud Computing BINAC at the Zentrum für Datenverarbeitung of the University of Tübingen, at the state of Baden-Württemberg, Germany.

Bioinformatic tools were installed using Conda or Mamba.

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

Samples were download using the [SRA-toolkit](https://github.com/ncbi/sra-tools) v3.1.1 and the SRA ID provided in the [metadata data-frame](https://github.com/SebasSaenz/chickenphages/tree/main/data_availability). Requesting more than 1 thread can crash the process.

```         
prefetch \
  SRR12345678 \
  -O sra_temp \
  --verify yes \
  --max-size 50G
  
fasterq-dump \
  sra_temp/SRR12345678 \
  --split-files \
  --outdir fastq_files \
  --threads 4
```

## Quality control and host removal {#quality-control-and-host-removal}

```         
# Create a QC report

fastqc \
  -q \
  -t 4 \
  -o out_dir \
  -f fastq \
  reads/sample_R1.fq.gz \
  reads/sample_R2.fq.gz
  
# Create index

bowtie2-build \
  --threads 8 \
  reference_genome.fasta \
  index/host
  
# Read QC and host DNA removal

trim_galore \
  --paired \
  -j 4 \
  reads/sample_R1.fq.gz \
  reads/sample_R2.fq.gz
  
# Remove host DNA

bowtie2 \
  -p 8 \
  -x index/host \
  -1 out_dir/sample_val_1.fq.gz \
  -2 out_dir/sample_val_2.fq.gz \
  --un-conc-gz out_dir/no_host \
  > out_dir/host.sam
```

## Assembly and filtering {#assembly-and-filtering}

```         
# Create the assemblies

megahit \
  -1 clean_reads/sample_R1.fq.gz \
  -2 clean_reads/sample_R2.fq.gz \
  -o megahit_results/sample_001 \
  -t 8 \
  --presets meta-sensitive
  
# Filter contigs shorter than 5kb

bbduk.sh \
  in=megahit_results/sample_contigs.fasta \
  out=filtered/sample_filtered.fasta \
  minlen=5000
```

## Virus identification and quality {#virus-identification-and-quality}

--cleanup: removes intermediate files to save space. --conservative: uses stricter thresholds for virus identification. genomad-db must be downloaded and the path correctly specified.

`checkv end_to_end` performs quality assessment and completeness estimation of viral genomes. Ensure the CheckV database is downloaded and path is correctly set. Input FASTA files should contain contigs ≥1,500 bp for best results.

```         
genomad end-to-end \
  --cleanup \
  --conservative \
  assemblies/sample_filtered_contigs.fasta \
  genomad_output/sample \
  /path/to/genomad-db \
  --threads 16
  
checkv end_to_end \
  viral_contigs/ \
  checkv_output/ \
  -t 18 \
  -d /path/to/checkv-db
```

## High-quality viral genomes and vOTUs {#high-quality-viral-genomes-and-votus}

Scripts to cluster the viral genomes at 95% ANI and 85% min-coveerage can be found in : <https://github.com/snayfach/MGV>

```         
# Filter Medium and High-quality genomes

seqkit grep \
  -n \
  -f ids.txt \
  contigs.fasta \
  > filtered_contigs.fasta

#Cluster vOTUs

makeblastdb \
  -in contigs.fasta \
  -dbtype nucl \
  -out blast_db/contigs_db
  
blastn \
  -query viral_contigs.fasta \
  -db blast_db/nt \
  -outfmt '6 std qlen slen' \
  -max_target_seqs 10000 \
  -out results/blast_hits.tsv \
  -num_threads 20

anicalc.py \
  -i blast_hits.tsv \
  -o ani_summary.tsv


aniclust.py \
  --fna genomes/fna_files \
  --ani ani_summary.tsv \
  --out clustering_results \
  --min_ani 95 \
  --min_tcov 85 \
  --min_qcov 0
```

## Cluster family and genus {#cluster-family-and-genus}

This script filters AAI-based pairwise comparisons to retain only those with sufficient biological similarity. It is useful for defining genome clusters or taxonomic thresholds based on shared protein content and identity. Script can be found in : <https://github.com/snayfach/MGV>

```         
# Create a blast+ databasa using your viral contigs

diamond makedb \
  --in proteins.faa \
  --db diamond_db/protein_db \
  --threads 10
  
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

python amino_acid_identity.py \
  --in_faa proteins/sample.faa \
  --in_blast blast_results/sample_vs_db.tsv \
  --out_tsv results/aa_identity.tsv
  
# Filter edges and prepare MCL inputs GENUS

python filter_aai.py \
  --in_aai results/aa_identity.tsv \
  --min_percent_shared 20 \
  --min_num_shared 16 \
  --min_aai 40 \
  --out_tsv results/aa_identity_filtered.tsv

# Perform MCL-based clustering

mcl similarity_matrix.txt \
  -te 8 \
  -I 2.0 \
  --abc \
  -o mcl_clusters.txt
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
