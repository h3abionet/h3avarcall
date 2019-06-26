# h3avarcall
H3A variant calling pipeline

### 1. Download the test datasets first:
Assuming you have no data to play with, make a directory in the workflow directory and call it ```data``` before continuing, or else provide your own directory where the ```FASTQ``` files are and skip this step. ``cd`` into ``data`` and download the test data using one of the following methods.

#### 1.1. Using LFTP (Faster)
```
lftp -e "pget -n 20 ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz; bye"
lftp -e "pget -n 20 ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz; bye"
```

#### 1.2. Using WGET (Slower)
```
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
```

### 2. Download the ```Singularity``` containers for executing the pipeline:
```
nextflow run main.nf --mode do.GetContainers

```

### 3. Download the GATK bundle needed for execution of the workflow:
```
nextflow run main.nf --mode do.GenomeIndexing
```

### 4. Perform QC on data:
```
nextflow run main.nf --mode do.QC
```

#### 4.1. Perform trimming of the data
```
nextflow run main.nf --mode do.Trimming
```

### 6. Perform alignment
```
nextflow run main.nf --mode do.Alignment
```




