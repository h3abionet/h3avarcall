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
  |--containers                       ## Folder for Singularity images and recipes (in case you want to build yourself). All downloaded images go here!
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
The `main.config` file:
```groovy
/*==================================================================================================
 * THIS FILE IS USED TO SPECIFY INPUT, OUTPUTS AND PARAMETERS. THE FOLLOWING OPTIONS ARE THE ALLOWED:
 * ==================================================================================================
 * data         : Path to where the data is (FASTQ files).
 * out          : Path to store output results.
 * bundle       : GATK-b37-bundle list file.
 * mode         : Worflow step to perform. Can be any of [ do.GetContainers | do.GenomeIndexing | do.QC | do.ReadTrimming | do.ReadAlignment | do.VarianCalling | do.VariantFiltering ].
 * trim         : Trimming options for Trimmomatic.
 * resources    : Location of the GATK-b37-bundle folder.
 * from         : Workflow step to resume workflow from. Can be any of [ do.QC | do.ReadTrimming | do.ReadAlignment | do.VarianCalling | do.VariantFiltering ].
 * params.help  : Print help menu.
 * ==================================================================================================
 * BELOW ARE THE DEFAULT PARAMETERS! YOU'RE MORE THAN WELCOME TO CHANGE AS DESIRED!
 * ==================================================================================================
 */

params {
    data         = "$baseDir/data"
    out          = "$baseDir/results"
    bundle       = "$baseDir/gatk-b37-bundle/b37_files_minimal.txt"
    mode         = "do.QC"
    trim         = "ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:28 MINLEN:40"
    resources    = "$baseDir/gatk-b37-bundle"
    from         = null
    params.help  = null
}

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
If by some miracle you happen to have access to the WITS Cluster, you do not need to download the GATK-b37-bundle! Simply `cd` into the `gatk-b37-bundle` folder of the `h3avarcall` repo and soft-link the GATK-b37-bundle data as follows: 
```
cd gatk-b37-bundle
ln -s /global/blast/gatk-bundle/b37/* .
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
Depends on Read Alignment STEP! **MUST** run STEP 2.3 (`--mode do.ReadAlignment`) before running this step.
```bash
nextflow run main.nf -profile slurm --mode do.VariantCalling 
```
### 2.5. Variant Filtering
Depends on Variant Calling STEP! **MUST** run STEP 2.4 (`--mode do.VariantCalling`) before running this step.
```bash
nextflow run main.nf -profile slurm --mode do.VariantFiltering 
```

## 3. `h3avarcall` results
Assuming the output folder was left as default (in the `main.config` file), the results for running the `h3avarcall` will be found in the `results` folder of the `h3avarcall` repository. The results for each of the main workflow steps (`2.1` - `2.5`) are grouped as follows:
```
- [1] Read QC (optional)         =>    `results/1_QC`
- [2] Read Trimming (optional)   =>    `results/2_Read_Trimming`
- [3] Read Alignment             =>    `results/3_Read_Alignment`
- [4] Variant Calling            =>    `results/4_Variant_Calling`
- [5] Variant Filtering          =>    `results/5_Variant_Filtering`
```

Nested withing these results folders are ouput files from each step of the workflow, including a folder "`workflow_report`" containing `h3avarcall_report.html`, `h3avarcall_timeline.html`, `h3avarcall_workflow.dot` and `h3avarcall_trace.txt` that contain detailed information on the resources (CPU, MEMORY and TIME) usage of each process in the steps. The `results` directory structure within `h3avarcall` repository can be summarized as below:

```bash
h3avarcall
  |--results
  |  |--1_QC
  |  |  |--workflow_report
  |  |  |  |--h3avarcall_report.html
  |  |  |  |--h3avarcall_timeline.html
  |  |  |  |--h3avarcall_workflow.dot
  |  |  |  |--h3avarcall_trace.txt
  |  |  |--<sample_1>_R1.fastqc.html .. <sample_N>_R1.fastqc.html
  |  |  |--<sample_1>_R2.fastqc.html .. <sample_N>_R1.fastqc.html
  |  |--2_Read_Trimming
  |  |  |--workflow_report
  |  |  |  |--h3avarcall_report.html
  |  |  |  |--h3avarcall_timeline.html
  |  |  |  |--h3avarcall_workflow.dot
  |  |  |  |--h3avarcall_trace.txt
  |  |  |--<sample_1>.1P.fastq.gz .. <sample_N>.1P.fastq.gz
  |  |  |--<sample_1>.2P.fastq.gz .. <sample_N>.2P.fastq.gz
  |  |--3_Read_Alignment
  |  |  |--workflow_report
  |  |  |  |--h3avarcall_report.html
  |  |  |  |--h3avarcall_timeline.html
  |  |  |  |--h3avarcall_workflow.dot
  |  |  |  |--h3avarcall_trace.txt
  |  |  |--<sample_1>_md.recal.bam .. <sample_N>_md.recal.bam
  |  |  |--<sample_1>_md.recal.bai .. <sample_N>_md.recal.bai
  |  |--4_Variant_Calling
  |  |  |--workflow_report
  |  |  |  |--h3avarcall_report.html
  |  |  |  |--h3avarcall_timeline.html
  |  |  |  |--h3avarcall_workflow.dot
  |  |  |  |--h3avarcall_trace.txt
  |  |  |--chr_1_genotyped.vcf.gz .. chr_22_genotyped.vcf.gz
  |  |  |--chr_1_genotyped.vcf.gz.tbi .. chr_22_genotyped.vcf.gz.tbi
  |  |--5_Variant_Filtering
  |  |  |--workflow_report
  |  |  |  |--h3avarcall_report.html
  |  |  |  |--h3avarcall_timeline.html
  |  |  |  |--h3avarcall_workflow.dot
  |  |  |  |--h3avarcall_trace.txt
  |  |  |--genome.SNP-recal.vcf.gz
  |  |  |--genome.SNP-recal.vcf.gz.tbi
  |--work
  |  |--<There's a lot of folders here! Lets not worry about them for today!>
```
##MORE TO BE ADDED
