# circRNA Detection Pipeline (hg38)

This README provides a complete guide to perform circular RNA (circRNA) detection using paired-end RNA-seq data, from downloading reference files to annotation using **CIRCexplorer2**.

---

## 1. Install Miniconda
```bash
cd ~
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
source ~/miniconda3/etc/profile.d/conda.sh
echo 'source ~/miniconda3/etc/profile.d/conda.sh' >> ~/.bashrc
conda --version
```

---

## 2. Download FASTQ Files
Replace `DRRXXXXXX` with your dataset ID.
```bash
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRRXXX/DRRXXXXXX/DRRXXXXXX_1.fastq.gz
wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/fastq/DRRXXX/DRRXXXXXX/DRRXXXXXX_2.fastq.gz
```

Unzip files:
```bash
gunzip DRRXXXXXX_1.fastq.gz
gunzip DRRXXXXXX_2.fastq.gz
```

---

## 3. Quality Check using FastQC
```bash
fastqc -o /path/to/output_directory -t 4 /path/to/fastq/SAMPLE_ID_1.fastq /path/to/fastq/SAMPLE_ID_2.fastq
```

---

## 4. Install Trimmomatic via Conda
```bash
conda install -c bioconda trimmomatic
```

### Paired-End Trimming Command
```bash
trimmomatic PE -threads 8 \
fastq/SAMPLE_ID_1.fastq fastq/SAMPLE_ID_2.fastq \
trimmed/SAMPLE_ID_1_paired.fastq trimmed/SAMPLE_ID_1_unpaired.fastq \
trimmed/SAMPLE_ID_2_paired.fastq trimmed/SAMPLE_ID_2_unpaired.fastq \
ILLUMINACLIP:/path/to/adapters/TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

---

## 5. Download GENCODE hg38 Reference Files
```bash
mkdir -p ~/genome/hg38 && cd ~/genome/hg38

# Download FASTA and GTF files
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz

# Unzip
gunzip *.gz
```

---

## 6. Build STAR Genome Index
```bash
STAR --runThreadN 8 \
--runMode genomeGenerate \
--genomeDir /path/to/genome_index \
--genomeFastaFiles GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile gencode.v47.annotation.gtf \
--sjdbOverhang 100
```

---

## 7. Align Reads with Reference Genome (for circRNA detection)
```bash
STAR --runThreadN 8 \
--genomeDir /path/to/genome_index \
--readFilesIn trimmed/SAMPLE_1_paired.fastq trimmed/SAMPLE_2_paired.fastq \
--outFileNamePrefix aligned/SAMPLE_ \
--outSAMtype BAM SortedByCoordinate \
--chimSegmentMin 10 \
--chimOutType Junctions SeparateSAMold
```

---

## 8. Convert GTF to GenePred Format
```bash
gtfToGenePred -genePredExt -geneNameAsName2 gencode.v47.annotation.gtf gencode.v47.annotation.genepred
```

Convert GenePred to CIRCexplorer2-compatible TXT:
```bash
perl -alne '$"="\t";print "@F[11,0..9]"' gencode.v47.annotation.genepred > gencode.v47.annotation.txt
```

---

## 9. Run CIRCexplorer2 Parse
```bash
CIRCexplorer2 parse -t STAR ~/hstar_out/DRR415361_Chimeric.out.junction | grep -P "^chr" > ~/hstar_out/back_spliced_junction.bed
```

---

## 10. circRNA Annotation using CIRCexplorer2
```bash
CIRCexplorer2 annotate \
-r ~/genome/hg38/gencode.v47.annotation.txt \
-g ~/genome/hg38/GRCh38.primary_assembly.genome.fa \
-b ~/hstar_out/back_spliced_junction.bed
```

### Parameters
- `-r`: Reference annotation in TXT format  
- `-g`: Genome sequence in FASTA format  
- `-b`: BED file with back-spliced junctions from parse

---

## Output Files
- **Aligned BAM files** – from STAR alignment
- **back_spliced_junction.bed** – from CIRCexplorer2 parse
- **Annotated circRNAs** – from CIRCexplorer2 annotate

---

## Notes
- Replace all placeholder paths (e.g., `/path/to/`) with your actual directory locations.
- Ensure adequate disk space and CPU cores for STAR alignment.
- Use release_47 (hg38) from GENCODE for the latest annotation support.

---
**Pipeline Summary:**  
Miniconda → FASTQ QC (FastQC) → Trimming (Trimmomatic) → Reference Setup (GENCODE hg38) → STAR Index + Alignment → circRNA Detection (CIRCexplorer2).