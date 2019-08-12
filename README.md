# `h3avarcall` - H3ABioNet Variant Calling Pipeline
`h3avarcall` is an H3ABioNet Variant Calling [```Nextflow```](https://www.nextflow.io/) pipeleline developed for blah blah blah...

## 1. Data Pipeline Preparation
Clone `h3avarcall` repository onto you machine:
```
git clone https://github.com/h3abionet/h3avarcall.git
cd h3avarcall
```
Contents of the repository:
```
h3avarcall
  |-containers
  |  |--Singularity.bwa
  |  |--Singularity.fastqc
  |  |--Singularity.gatk
  |  |--Singularity.trimmomatic
  |--gatk-b37-bundle
  |  |--b37_files_minimal.txt
  |--templates
  |  |--download_bundles.sh
  |--LICENSE
  |--README.md
  |--main.config
  |--main.nf
  |--nextflow.config
```

### 1.1. Download test datasets:
Create a data directory within the `h3avarcall` repository:
```
mkdir data
cd data
```
Download the test data from [THE_SITE](http://thesite.com):
#### 1.1.2. Using LFTP (faster)
```
lftp -e "pget -n 20 ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz; bye"
lftp -e "pget -n 20 ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz; bye"
```

#### 1.1.3. Using WGET (slower)
```
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
```

### 1.2. Download the `Singularity` containers required to execute the pipeline:
```
nextflow run main.nf -profile slurm --mode do.GetContainers

```

### 1.3. Download the GATK b37 bundle required to execute the workflow:
This step takes **FOREVER** to run - run it only once!

```
nextflow run main.nf -profile slurm --mode do.GenomeIndexing
```

## 2. QC & READ TRIMMING OF THE DATA

### 2.1. Perform QC on data (optional):
```
nextflow run main.nf -profile slurm --mode do.QC
```

#### 2.3. Perform trimming of the data (optional):
```
nextflow run main.nf -profile slurm --mode do.Trimming
```

## 3. ALIGNMENT OF DATA
### 3.1. Perform alignment
```
nextflow run main.nf -profile slurm --mode do.Alignment
```
