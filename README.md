# `h3avarcall` - H3ABioNet Variant Calling Pipeline
`h3avarcall` is an H3ABioNet Variant Calling [```Nextflow```](https://www.nextflow.io/) pipeleline developed for blah blah blah...

## 1. Obtaining pipeline and preparing Data
Clone `h3avarcall` repository onto you machine:
```bash
git clone https://github.com/h3abionet/h3avarcall.git
cd h3avarcall
```
Contents of the repository:
```bash
h3avarcall
  |-containers                        ## Folder for Singularity images and recipes (in case you want to build yourself). All downloaded images go here!
  |  |--Singularity.bwa               ## Singularity recipe file for BWA and Samtools.
  |  |--Singularity.fastqc            ## Singularity recipe file for FastQC.
  |  |--Singularity.gatk              ## Singularity recipe file for GATK and tabix.
  |  |--Singularity.trimmomatic       ## Singularity recipe file for Trimmimatic.
  |--gatk-b37-bundle                  ## Folder for stoding downloaded GATK-b37-bundle files.
  |  |--b37_files_minimal.txt         ## LList of GATK-b37-bundle files to be downloaded (bundle TOO BIG! Only selected files needed for the workflow). 
  |--templates                        ## Folder for extra scripts for the workflow.
  |  |--download_bundles.sh           ## Script for downloading GATK-b37-bundle.
  |--LICENSE                          ## Duh!
  |--README.md                        ## Duh!
  |--main.config                      ## User configuration file! All inputs, outputs and options GO HERE!! ONLY file that SHOULD be modified by user!
  |--main.nf                          ## Main h3avarcall nextflow scripts.
  |--nextflow.config                  ## Pipeline configuration file! DO NOT EDIT!!!
```

### 1.1. Download test datasets:
Create a data directory within the `h3avarcall` repository:
```bash
mkdir data
cd data
```
Download the test data from [THE_SITE](http://thesite.com):
#### 1.1.2. Using LFTP (faster)
```bash
lftp -e "pget -n 20 ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz; bye"
lftp -e "pget -n 20 ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz; bye"
```

#### 1.1.3. Using WGET (slower)
```bash
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
```

### 1.2. Download the `Singularity` containers required to execute the pipeline:
```bash
nextflow run main.nf -profile slurm --mode do.GetContainers
```

### 1.3. Download the GATK b37 bundle required to execute the workflow:
This step takes **FOREVER** to run - run it only once!

```bash
nextflow run main.nf -profile slurm --mode do.GenomeIndexing
```

## 2. Executing the main `h3avarcall` pipeline

### 2.1. Read QC (optional):
```bash
nextflow run main.nf -profile slurm --mode do.QC
```

### 2.2. Read Trimming (optional):
```bash
nextflow run main.nf -profile slurm --mode do.ReadTrimming
```

### 2.3. Read Alignment
Can be run with `--from do.ReadTrimming` or `--from do.QC` depending on whether these step were run! 
```bash
nextflow run main.nf -profile slurm --mode do.ReadAlignment
```

### 2.4. Variant Calling
```bash
nextflow run main.nf -profile slurm --mode do.VariantCalling --from do.ReadAlignment
```
### 2.5. Variant Filtering
```bash
nextflow run main.nf -profile slurm --mode do.VariantFiltering --from do.VariantCalling 
```
