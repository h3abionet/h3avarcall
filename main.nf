#!/usr/bin/env nextflow

// HELP MENU

// PARAMETERS
params.data          = "/home/phelelani/h3abionet/data"
params.out           = "/home/phelelani/h3abionet/results"
params.extension     = "fastq.gz"
params.mode          = "do.alignment"
params.trim          = "ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:28 MINLEN:40"
params.resources     = "/dataA/courses/uom/ref"

// SET PARAMETERS
data_dir             = file(params.data, type: 'dir')
out_dir              = file(params.out, type: 'dir')
ext                  = params.extension
trim_params          = params.trim
resources            = file(params.resources, type: 'dir')
genome               = file("${resources}/ucsc.hg19.fasta", type: 'file')

out_dir.mkdir()

// GET DATA
read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.${ext}", type: 'file')

// START PROCESSING READS

switch (params.mode) {
    case['do.GetContainers']:
        println "\nDownloading Singularity containers."
        
        // process run_GetContainers {
            
        // }
        break
        // --------------------
        
    case['do.GenomeIndexing']:
        println "\nDownloading reference genome and indexing"

        // process run_GenomeIndexind {            

        // }
        break
        // --------------------

    case['do.qc']:
        println "\nPerforming QC for all the samples"
        
        process run_QualityChecks {
            cpus 6
            memory '5 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file(reads) from read_pairs
            
            output:
            set sample, file("${sample}*.html") into qc_results
            
            """
            fastqc ${reads.get(0)} ${reads.get(1)} \
                --threads 5 \
                --noextract
            """
        }     
        break
        // --------------------

    case['do.trimming']:
        println "Performing trimming for all the samples"

        process run_ReadTrimming {
            cpus 6
            memory '20 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file(reads) from read_pairs
            
            output:
            set sample, file("${sample}*") into trim_results
            
            """
            java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
                -threads 5 \
                -trimlog trimlog_${sample}.log \
                -basein ${reads.get(0)} \
                -baseout ${sample}_clean \
                \$(sed 's|ILLUMINACLIP:|ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/|' <<< "${trim_params}")
            """
        }
        break
        // --------------------        
        
    case['do.alignment']:
        println "Performing alignment for all the samples"

        process run_ReadAlignment {
            cpus 6
            memory '20 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file(reads) from read_pairs
            
            output:
            set sample, file("${sample}.bam") into bam_results
            
            """
            bwa mem -R \"@RG\\tID:${sample}\\tLB:LIBA\\tSM:${sample}\\tPL:Illumina\" \
                -t 5 -M ${genome} \
                ${reads.get(0)} ${reads.get(1)} | samtools sort \
                --threads 5 -m 2G - > ${sample}.bam
            """
        }
        break
        // --------------------

    case['do.something']:
        println "Performing something for all the samples"
        
        break
}

bam_results.subscribe { println it }
