#!/usr/bin/env nextflow

// HELP MENU

// PARAMETERS
params.data          = "/home/phelelani/h3abionet/data"
params.out           = "/home/phelelani/h3abionet/results"
params.extension     = "fastq.gz"
params.mode          = "do.alignment"
//params.mode          = "do.qc"
params.trim          = "ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:28 MINLEN:40"
params.resources     = "/global/blast/gatk-bundle/b37"

// SET PARAMETERS
data_dir             = file(params.data, type: 'dir')
out_dir              = file(params.out, type: 'dir')
ext                  = params.extension
trim_params          = params.trim
resources            = file(params.resources, type: 'dir')
//
genome               = file("${resources}/human_g1k_v37_decoy.fasta", type: 'file')

out_dir.mkdir()

// GET DATA
read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.${ext}", type: 'file')

// START PROCESSING READS

switch (params.mode) {
    case['do.GetContainers']:
        println "\nDownloading Singularity containers."
        
        base = "shub://phelelani/nf-rnaSeqMetagen:"
        shub_images = Channel.from( ["${base}gatk", "${base}bwa", "${base}trimmomatic", "${base}fastqc"] )
        
        process run_DownloadContainers {
            cpus 1
            memory '2 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { "Downloading: $link" }
            publishDir "$baseDir/containers", mode: 'copy', overwrite: true, pattern: "*.simg"
            
            input:
            each link from shub_images
            
            output:
            file("*.simg") into containers
        
            """
            singularity pull ${link}
            """
        }  
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
                --threads 5 - > ${sample}.bam
            """
        }

        process run_MarkDuplicates {
            cpus 6
            memory '20 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file(raw_bam) from bam_results

            output:
            set sample, file("${sample}_md.bam") into md_bam_results

            """
            gatk --java-options -Xmx4g MarkDuplicates \
                --MAX_RECORDS_IN_RAM 50000 \
                --INPUT ${raw_bam} \
                --METRICS_FILE ${raw_bam}.metrics \
                --ASSUME_SORT_ORDER coordinate \
                --CREATE_INDEX true \
                --OUTPUT ${sample}_md.bam
            """
        }

        process run_CreateRecalibrationTable {
            cpus 6
            memory '20 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file(md_bam) from md_bam_results

            output:
            set sample, file("${sample}*") into recalibration_results

            """
            gatk --java-options -Xmx10g BaseRecalibrator \
                --input ${md_bam} \
                --output ${sample}_recal.table \
                -R ${genome} \
                --known-sites ${resources}/dbsnp_138.b37.vcf \
                --known-sites ${resources}/1000G_phase1.indels.b37.vcf \
                --known-sites ${resources}/Mills_and_1000G_gold_standard.indels.b37.vcf
            """
        }
        
        // process run_RecalibrateBAM {
        //     cpus 6
        //     memory '20 GB'
        //     time '2h'
        //     scratch '$HOME/tmp'
        //     tag { sample }
        //     publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
        //     input:
        //     set sample, file("${sample}*") from recalibration_results
            
        //     output:
            
            
        //     """
        //     gatk --java-options -Xmx10g ApplyBQSR \
        //          --input ${bam_file} \
        //          --output ${sample_id}.md.recal.bam \
        //          --TMP_DIR ${params.gatk_tmp_dir} \
        //          -R ${genome} \
        //          --create-output-bam-index true \
        //          --bqsr-recal-file ${recal_table_file}
        //     """
        // }

        break
        // --------------------

    case['do.something']:
        println "Performing something for all the samples"
        
        break
}


// lftp -e 'mirror --use-pget-n=10 /bundle/hg19 .' -u gsapubftp-anonymous, ftp.broadinstitute.org
// lftp -e 'pget -n20 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz; bye'
