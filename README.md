![](https://pathdet.hgc.jp/files/static/images/pathdet_logo.png)

# PATHDET: Command-Line Pipeline for PATHogen DETection (for long-read seq version)
**PATHDET** is a command-line pipeline for *“PATHogen DETection”* from extracellular DNA/RNA using high-throughput sequencing data. Unlike the original web application, this version is designed for local or HPC (High-Performance Computing) environments and can be executed entirely from the command line.

**PATHDET** pipeline (step-by-step workflow shown below) detects pathogenic microorganisms from sequencing data derived from infected hosts (humans, animals, plants, etc.).  It is also applicable to the investigation of microbial contaminants in biological experiments and clinical instruments.

Through the execution of **PATHDET**, users can obtain a comprehensive list of microorganisms present in their samples.  
The pipeline generates two main outputs:

- **Rapid Report**: A concise overview of the detected pathogens, enabling a rapid assessment.  
- **Full Report**: A downloadable file containing taxonomic classification, abundance estimation, and additional QC metrics.

This command-line version has the potential to provide reproducibility, scalability, and flexibility for integration into automated workflows and public health surveillance systems.

![](https://pathdet.hgc.jp/files/static/images/main.png)

## Installation

PATHDET requires a prepared conda environment and several reference databases.  
Follow the steps below to set up the environment before running the pipeline.

### 1. Setup conda environment
```
# Load conda initialization script (example path may vary)
$ source /opt/miniforge/etc/profile.d/conda.sh

# Create environment from provided environment files (if distributed)
$ conda env create -f pathdet_env.yaml

# Activate PATHDET environment
$ conda activate pathdet_env
```

### 2. Prepare reference databases

The following databases are required. Place them under a common directory (e.g. `/path/to/datadir`):

- **Kraken2 indices**  
  - Human (GRCh38)  
  - PlusPF (2024/12/28 version)  

- **BLAST databases**  
  - `nt`  
  - `T2T-CHM13v2.0`  

- **Taxonkit SQL database**  
  - `taxdump`  

- **RefSeq genomes**  
  - Downloaded via `ncbi-genome-download`  

- **Host genome indices (example: Human GRCh38)**  
  - Bowtie2 index  
  - Hisat2 index  
  - *By default, the host genome is set to Human (GRCh38).  
    If analyzing non-human samples, users should prepare and specify the appropriate host genome index.*  

Edit the variable `datadir` in the script `pathdet_short.sh` (around line 67) to specify the parent directory:

```
## Database
datadir="/path/to/datadir"   # <-- change this line
```
- If you set up the environment manually, follow the steps below.


#### Database for both short- and long-read sequencing

```
(pathdet_env)$ kraken2-build --download-taxonomy --db db/kraken2/human --threads 20
(pathdet_env)$ kraken2-build --download-library human --db db/kraken2/human --threads 20
(pathdet_env)$ kraken2-build --build --db db/kraken2/human –-threads 20
(pathdet_env)$ wget \
https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_20241228.tar.gz
(pathdet_env)$ mkdir -p db/kraken2/k2_pluspf_20241228
(pathdet_env)$ tar xvf k2_pluspf_20241228.tar.gz -C db/kraken2/k2_pluspf_20241228
(pathdet_env)$ wget \
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-
CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
(pathdet_env)$ gunzip GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
(pathdet_env)$ makeblastdb -dbtype nucl -in GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
-title t2t -parse_seqids -out db/blast/t2t/GCF_009914755.1_T2T-CHM13v2.0_genomic
(pathdet_env)$ mkdir -p db/blast/nt
(pathdet_env)$ pushd db/blast/nt
(pathdet_env)$ update_blastdb.pl --decompress nt
(pathdet_env)$ popd
(pathdet_env)$ ${CONDA_PREFIX}/opt/krona/updateTaxonomy.sh
(pathdet_env)$ ${CONDA_PREFIX}/opt/krona/updateAccessions.sh
(pathdet_env)$ wget https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
(pathdet_env)$ mkdir -p db/taxonkit/taxdump
(pathdet_env)$ tar xvf taxdump.tar.gz -C db/taxonkit/taxdump
(pathdet_env)$ ncbi-genome-download -p 20 -r 10 \
--format fasta,assembly-report,assembly-stats -l complete -R reference \
-o db/ncbi_genome/reference bacteria,fungi,viral
(pathdet_env)$ ncbi-genome-download -p 20 -r 10 \
--format fasta,assembly-report,assembly-stats -l complete -R all \
-o db/ncbi_genome/all bacteria,fungi,viral
(pathdet_env)$ taxid.pl db/ncbi_genome/reference/refseq db/taxonkit/taxdump/nodes.dmp \
> db/ncbi_genome/reference/refseq_lineage_taxid_fasta.tsv
(pathdet_env)$ taxid.pl db/ncbi_genome/all/refseq db/taxonkit/taxdump/nodes.dmp \
> db/ncbi_genome/all/refseq_lineage_taxid_fasta.tsv
```


#### Database for long-read sequencing data analysis  
```
(pathdet_env)$ wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
(pathdet_env)$ gunzip hg38.fa.gz
(pathdet_env)$ wget https://obj.umiacs.umd.edu/mash/screen/RefSeq88n.msh.gz
(pathdet_env)$ mkdir -p db/mash
(pathdet_env)$ gzip -dc RefSeq88n.msh.gz > db/mash/RefSeq88n.msh
```

### 3. Input files

Prepare sequencing data from Illumina (or other-type short-read sequencers) in **paired-end FASTQ format**.

#### Requirements
- **Read 1 FASTQ file** (e.g., `sample_R1.fastq.gz`)
- **Read 2 FASTQ file** (e.g., `sample_R2.fastq.gz`)
- Both files must correspond to the same sample.

### 4. Example: Running the Analysis Pipeline
Before running, make sure the conda environment is active:

```
$ conda activate pathdet_env
```
####Command
```
(pathdet_env)$ pathdet_long.sh -i sample.fastq.gz -o sample \
-t 16 -s 10000 -l 1000
```
####Option

- `-i <file>`  
  Path to **FASTQ file** (required)

- `-o, --out <string>`  
  Prefix for output files and directories (required)

- `-t, --threads <int>`  
  Number of threads to use (default: 20)

- `-s, --split <int>`  
  Number of lines per chunk for BLAST-nt parallelization  
  (default: 50,000)

- `-l, --length <int>`  
  Minimum read length.
  Reads shorter than the specified value will be filtered out.
  (default: 1000)

- `--debug`  
  Keep intermediate files instead of deleting them after analysis

### Expected Outputs (directory layout)

When running with `-o sample`, the following directories and files are generated:

- `sample_qc/`  
  - `sample_trimmed2_R1.fastq.gz`  
  - `sample_trimmed2_R2.fastq.gz`  
  *→ Cleaned read data after quality control*

- `sample_hgs/`  
  - `sample_nohit.fasta`  
  *→ Reads after host genome subtraction*

- `sample_rapid/`  
  - `sample_pre-report.csv`  
  - `sample.plot.html` (Krona chart)  
  *→ Rapid analysis summary*

- `sample_blast/`  
  - `sample_blast_hit.txt`  
  *→ BLAST analysis results*

- `sample_tbl/`  
  - `sample_report.csv`  
  - `sample_bacteria_report.csv`  
  *→ Full analysis summaries*

  - `sample_fRep/`  
        - `<OUT>_f*.csv`  
  *→ Abundance tables for each taxonomic rank*

  - `sample_bRep/`  
        - `<OUT>_f*.csv`  
  *→ Bacteria-focused abundance tables*

- `sample_map/`  
  - `sample_map.csv`  
  *→ Mapping summary for top 3 detected genomes*
   
  - `spark_*.svg`  
  *→ Coverage plots*
  
#### * Taxonomy-level Output Files

In addition to summary tables, the pipeline generates **taxonomy-resolved abundance files**.  
These are stored under `sample_fRep/` and `sample_bRep/` directories.

- **`sample_fRep/`**  
  Contains abundance tables for all detected microorganisms across different taxonomic ranks.  
  Typical files include:  
  - `sample_fK.csv` : Kingdom level  
  - `sample_fP.csv` : Phylum level  
  - `sample_fF.csv` : Family level  
  - `sample_fG.csv` : Genus level  
  - `sample_fS.csv` : Species level  

- **`sample_bRep/`**  
  Contains bacteria-focused abundance tables, formatted in the same way as above.  
  Typical files include:  
  - `sample_bK.csv` : Kingdom level (bacteria only)  
  - `sample_bP.csv` : Phylum level (bacteria only)  
  - `sample_bF.csv` : Family level (bacteria only)  
  - `sample_bG.csv` : Genus level (bacteria only)  
  - `sample_bS.csv` : Species level (bacteria only)  

### Notes
- The suffixes **K, P, F, G, S** correspond to standard taxonomic ranks:  
  - `K` = Kingdom  
  - `P` = Phylum  
  - `F` = Family  
  - `G` = Genus  
  - `S` = Species  
- These files allow users to summarize microbial composition at different levels of taxonomy depending on their analysis needs.

## Citation
Please cite the following if you are using PATHDET for short-read:
[Horiba K, Torii Y, Okumura T, Takeuchi S, Suzuki T, Kawada JI, Muramatsu H, Takahashi Y, Ogi T, Ito Y. Next-Generation Sequencing to Detect Pathogens in Pediatric Febrile Neutropenia: A Single-Center Retrospective Study of 112 Cases. Open Forum Infect Dis. 2021 May 4;8(11):ofab223. doi: 10.1093/ofid/ofab223.](https://doi.org/10.1093/ofid/ofab223)

Please cite the following if you are using PATHDET for long-read:
[Horiba K, Torii Y, Aizawa Y, Yamaguchi M, Haruta K, Okumura T, Suzuki T, Kawano Y, Kawada JI, Hara S, Saitoh A, Giske CG, Ogi T, Ito Y. Performance of Nanopore and Illumina metagenomic sequencing for pathogen detection and transcriptome analysis in infantile central nervous system infections. Open Forum Infect Dis. 2022 ofab504. doi: 10.1093/ofid/ofac504 .](https://doi.org/10.1093/ofid/ofac504)

## Dependencies

PATHDET requires a conda environment with the following software tools installed.  
These dependencies ensure reproducibility of the pipeline.

| Software / Package       | Version    | Purpose / Notes |
|---------------------------|------------|-----------------|
| **NanoPlot**             | 1.44.1     | QC and visualization of read length distributions |
| **seqkit**               | 2.9.0      | FASTA/FASTQ file utilities |
| **nanofilt**             | 2.8.0      | Read filtering for Nanopore data |
| **kraken2**              | 2.1.3      | Taxonomic classification (k-mer based) |
| **minimap2**             | 2.28       | Long-read alignment |
| **samtools**             | 1.21       | BAM/CRAM/SAM processing |
| **mash**                 | 2.3        | Fast genome distance estimation |
| **blast**                | 2.16.0     | Sequence similarity search |
| **krona**                | 2.8.1      | Interactive hierarchical visualization |
| **parallel**             | 20241222   | Parallelized execution of tasks |
| **taxonkit**             | 0.18.0     | Taxonomic lineage operations |
| **bbmap**                | 39.17      | Read mapping and QC utilities |
| **bcftools**             | 1.21       | Variant calling and filtering |
| **deeptools**            | 3.5.5      | NGS data visualization |
| **ncbi-genome-download** | 0.3.3      | Batch genome downloads from NCBI RefSeq |
| **SparK**                | 2.6.2_1    | Coverage plotting |

## Disclaimer
* The users are responsible for all risks associated with the use of and the data transfer to PATHDET. PATHDET makes no legal liability for any claims or damages arising from the use of PATHDET.
* PATHDET is not responsible for the accuracy, equivalence, and completeness of the analysis results.

## Usage Policy
* PATHDET is provided as a **command-line pipeline** to be executed in local or HPC environments.  
  Users are fully responsible for ensuring appropriate handling of sequencing data.
* Prior to analysis, **human or other host-derived reads containing identifiable personal information** should be removed or masked as required by applicable laws, regulations, and institutional policies.
* PATHDET is free to use for **academic and non-profit purposes**.  
  Commercial use is not permitted without prior approval.
* Since this version runs entirely in local or institutional computing environments, **no sequence data is transmitted externally**.  
  Users are responsible for appropriate data storage, retention, and deletion policies within their institution.

![](https://pathdet.hgc.jp/files/static/images/footer.png)
