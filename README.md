# CircRNA Detection Pipeline for FECD

[![Medium Article](https://img.shields.io/badge/Medium-Read%20Article-black?style=flat&logo=medium)](https://medium.com/@tejaswissh/building-a-robust-computational-pipeline-for-circrna-analysis-in-linux-0445888eb84e)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-Connect-0077B5?style=flat&logo=linkedin)](https://in.linkedin.com/in/tejaswi-velugapally-7451b1222)

A complete computational pipeline for detecting and analyzing circular RNAs (circRNAs) in Fuchs Endothelial Corneal Dystrophy (FECD) using paired-end RNA-Seq data.



## 📋 About This Pipeline

Circular RNAs (circRNAs) are endogenous non-coding RNAs with a closed-loop structure that makes them exceptionally stable. They play crucial roles in gene regulation by acting as microRNA sponges and interacting with RNA-binding proteins. This pipeline enables researchers to identify and analyze circRNAs from paired-end RNA-Seq data, with a specific focus on understanding circRNA-mediated mechanisms in FECD.

### What This Pipeline Does

Detects circRNAs from RNA-Seq data using back-splicing junction analysis
Performs differential expression analysis between disease and control samples
Constructs regulatory networks to understand circRNA-miRNA-mRNA interactions
Provides reproducible workflows for bioinformatics analysis in Linux environments

### Pipeline Capabilities

This pipeline performs:
- ✅ **circRNA Detection** from RNA-Seq data using back-splicing junction identification
- ✅ **Quality Control** of raw sequencing reads
- ✅ **Read Alignment** to the human reference genome (hg38)
- ✅ **Differential Expression Analysis** between patient and control samples
- ✅ **Annotation** of detected circRNAs with gene information

---

## 📁 Repository Structure

```
CircRna-pipeline/
├── scripts/
│   ├── circrna.sh              # Main automated pipeline script
│   ├── circ_rna_pipeline_readme.md  # Detailed technical documentation
│   └── deseq2.R                # Differential expression analysis in R
├── results/
│   └── DRR228774.xlsx          # Sample output file
└── README.md                   # This file
```



## 🛠️ Tools Used

| Tool | Version | Purpose |
|------|---------|---------|
| **Miniconda** | Latest | Environment and package management |
| **FastQC** | 0.11+ | Quality assessment of sequencing reads |
| **Trimmomatic** | 0.39+ | Adapter removal and quality trimming |
| **STAR** | 2.7+ | Splice-aware read alignment |
| **CIRCexplorer2** | 2.3+ | circRNA identification and annotation |
| **DESeq2** | 1.30+ | Differential expression analysis (R) |

**Reference Genome**: GRCh38 (hg38) from GENCODE release 47

---

##  Quick Setup

### Prerequisites

- **Operating System**: Linux (Ubuntu 18.04+ or CentOS 7+ recommended)
- **RAM**: Minimum 32 GB (64 GB recommended)
- **Storage**: 100+ GB free disk space
- **Software**: Bash shell, internet connection

---

## 📖 Step-by-Step Pipeline

### 1. Install Miniconda

Set up Miniconda for managing bioinformatics tools:

```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/miniconda3/etc/profile.d/conda.sh
echo 'source ~/miniconda3/etc/profile.d/conda.sh' >> ~/.bashrc
conda --version
```

---

### 2. Create Environment and Install Tools

Create an isolated environment with all required tools:

```bash
conda create -n circrna_env -c bioconda -c conda-forge \
    fastqc trimmomatic star circexplorer2 ucsc-gtftogenepred perl wget -y
    
conda activate circrna_env
```

---

### 3. Download FASTQ Files

Download your paired-end RNA-Seq data. Example using DDBJ database:

```bash
# Replace DRRXXXXXX with your actual dataset ID
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRRXXX/DRRXXXXXX/DRRXXXXXX_1.fastq.gz
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRRXXX/DRRXXXXXX/DRRXXXXXX_2.fastq.gz
```

Decompress the files:

```bash
gunzip DRRXXXXXX_1.fastq.gz
gunzip DRRXXXXXX_2.fastq.gz
```

---

### 4. Quality Control with FastQC

Assess the quality of your sequencing reads:

```bash
mkdir -p fastqc_reports

fastqc -o fastqc_reports -t 4 \
    DRRXXXXXX_1.fastq \
    DRRXXXXXX_2.fastq
```

**Check**: Open the HTML reports to verify read quality, adapter content, and sequence duplication levels.

---

### 5. Adapter Trimming with Trimmomatic

Remove adapter sequences and low-quality bases:

```bash
mkdir -p trimmed

trimmomatic PE -threads 8 \
    DRRXXXXXX_1.fastq DRRXXXXXX_2.fastq \
    trimmed/DRRXXXXXX_1_paired.fastq trimmed/DRRXXXXXX_1_unpaired.fastq \
    trimmed/DRRXXXXXX_2_paired.fastq trimmed/DRRXXXXXX_2_unpaired.fastq \
    ILLUMINACLIP:/path/to/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

**Parameters explained**:
- `ILLUMINACLIP`: Remove Illumina adapters
- `LEADING:3`: Cut bases with quality < 3 from the start
- `TRAILING:3`: Cut bases with quality < 3 from the end
- `SLIDINGWINDOW:4:15`: Scan with 4-base window, cutting when average quality < 15
- `MINLEN:36`: Drop reads shorter than 36 bases

---

### 6. Download Reference Genome

Download human genome reference (hg38) and annotation files:

```bash
mkdir -p ~/genome/hg38 && cd ~/genome/hg38

# Download genome FASTA
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz

# Download gene annotation GTF
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz

# Decompress
gunzip *.gz
```

---

### 7. Build STAR Genome Index

Create a genome index for fast alignment (only needs to be done once):

```bash
mkdir -p ~/genome/hg38_index

STAR --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir ~/genome/hg38_index \
    --genomeFastaFiles ~/genome/hg38/GRCh38.primary_assembly.genome.fa \
    --sjdbGTFfile ~/genome/hg38/gencode.v47.annotation.gtf \
    --sjdbOverhang 100
```

**Note**: This step requires ~32 GB RAM and takes 1-2 hours. The index can be reused for future analyses.

---

### 8. Align Reads with STAR

Align trimmed reads to the reference genome with chimeric junction detection:

```bash
mkdir -p aligned

STAR --runThreadN 8 \
    --genomeDir ~/genome/hg38_index \
    --readFilesIn trimmed/DRRXXXXXX_1_paired.fastq trimmed/DRRXXXXXX_2_paired.fastq \
    --outFileNamePrefix aligned/DRRXXXXXX_ \
    --outSAMtype BAM SortedByCoordinate \
    --chimSegmentMin 10 \
    --chimOutType Junctions SeparateSAMold
```

**Key parameters for circRNA detection**:
- `--chimSegmentMin 10`: Minimum length of chimeric segment to detect back-splicing
- `--chimOutType Junctions`: Output chimeric junctions separately

**Output files**:
- `DRRXXXXXX_Aligned.sortedByCoord.out.bam` - Aligned reads
- `DRRXXXXXX_Chimeric.out.junction` - Chimeric junctions (used for circRNA detection)

---

### 9. Convert GTF to GenePred Format

Prepare annotation files for CIRCexplorer2:

```bash
cd ~/genome/hg38

# Convert GTF to GenePred format
gtfToGenePred -genePredExt -geneNameAsName2 \
    gencode.v47.annotation.gtf \
    gencode.v47.annotation.genepred

# Convert GenePred to CIRCexplorer2-compatible TXT format
perl -alne '$"="\t";print "@F[11,0..9]"' \
    gencode.v47.annotation.genepred > gencode.v47.annotation.txt
```

---

### 10. Parse Chimeric Junctions

Extract back-spliced junctions from STAR output:

```bash
mkdir -p circexplorer_out

CIRCexplorer2 parse -t STAR \
    aligned/DRRXXXXXX_Chimeric.out.junction | \
    grep -P "^chr" > circexplorer_out/back_spliced_junction.bed
```

**What this does**: Identifies back-splicing events where the downstream exon connects to an upstream exon, forming a circular RNA.

---

### 11. Annotate circRNAs

Annotate detected circRNAs with gene information:

```bash
CIRCexplorer2 annotate \
    -r ~/genome/hg38/gencode.v47.annotation.txt \
    -g ~/genome/hg38/GRCh38.primary_assembly.genome.fa \
    -b circexplorer_out/back_spliced_junction.bed
```

**Parameters**:
- `-r`: Reference annotation (TXT format)
- `-g`: Genome sequence (FASTA format)
- `-b`: Back-spliced junction BED file

**Output**: `circularRNA_known.txt` - Annotated circRNAs with genomic coordinates, host gene names, and junction information.

---

### 12. Differential Expression Analysis (R)

After running the pipeline on all samples, analyze differential expression using DESeq2:

#### Install DESeq2 (one-time setup)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2", ask = FALSE)
```

#### Run the Analysis

Use the provided `deseq2.R` script:

```bash
Rscript scripts/deseq2.R
```

Or run interactively in R (see `scripts/deseq2.R` for full code).

**Requirements**:
- Create a `circRNA_matrix.csv` file with circRNA counts
- Rows = circRNA IDs
- Columns = Sample IDs
- First 7 columns = Control samples
- Next 10 columns = Patient samples

**Output**: `DESeq2_circRNA_results.csv` with differential expression statistics.

---

## 📊 Understanding Your Results

### Key Output Files

| File | Description |
|------|-------------|
| `*_Aligned.sortedByCoord.out.bam` | Aligned reads in BAM format |
| `*_Chimeric.out.junction` | Raw chimeric junctions from STAR |
| `back_spliced_junction.bed` | Filtered back-spliced junctions |
| `circularRNA_known.txt` | Annotated circRNAs with gene names |
| `DESeq2_circRNA_results.csv` | Differential expression results |

### Interpreting DESeq2 Results

The output table contains:

- **baseMean**: Average expression across all samples
- **log2FoldChange**: Expression change (Patient vs Control)
  - Positive = upregulated in patients
  - Negative = downregulated in patients
- **pvalue**: Statistical significance (uncorrected)
- **padj**: Adjusted p-value (FDR-corrected)

**Significance criteria** (commonly used):
- `padj < 0.05` (statistically significant)
- `|log2FoldChange| > 1` (at least 2-fold change)

### Sample Results

Example results are provided in `results/DRR228774.xlsx` showing typical output format and circRNA annotations.

---

## 🔧 Troubleshooting

### Common Issues and Solutions

| Issue | Solution |
|-------|----------|
| **STAR out of memory** | Use a system with ≥32 GB RAM or reduce `--limitGenomeGenerateRAM` |
| **Low alignment rate (<50%)** | Check read quality with FastQC; verify correct genome reference |
| **Few circRNAs detected** | Ensure RNA-Seq library is NOT poly(A)-selected (removes circRNAs) |
| **Trimmomatic adapter file not found** | Locate adapters: `conda list trimmomatic` then check installation directory |
| **"Command not found" errors** | Activate conda environment: `conda activate circrna_env` |

---

## 📚 Additional Resources

### Documentation

- **Detailed Technical Guide**: See `scripts/circ_rna_pipeline_readme.md` for in-depth parameter explanations
- **Medium Article**: [Building a Robust Computational Pipeline for circRNA Analysis in Linux](https://medium.com/@tejaswissh/building-a-robust-computational-pipeline-for-circrna-analysis-in-linux-0445888eb84e)

### External Links

- [GENCODE Database](https://www.gencodegenes.org/)
- [STAR Aligner Manual](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- [CIRCexplorer2 Documentation](https://circexplorer2.readthedocs.io/)
- [DESeq2 Vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

---

## 📖 Citation

If you use this pipeline in your research, please cite:

**Tools:**
- STAR: Dobin et al., *Bioinformatics* 2013
- CIRCexplorer2: Zhang et al., *Cell* 2016  
- DESeq2: Love et al., *Genome Biology* 2014
- Trimmomatic: Bolger et al., *Bioinformatics* 2014


---

## 📝 License

This project is open source and available for academic and research purposes.



## 🤝 Contributing

Contributions and suggestions are welcome! Feel free to:
- Open an issue for bugs or feature requests
- Submit pull requests for improvements
- Share your analysis results and experiences

