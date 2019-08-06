#!/usr/bin/env nextflow

// HELP MENU GOES HERE

// PARAMETERS GO HERE! WILL BE IN A SEPARATE CONFIG FILE - JUST TESTINF FOR NOW!
params.data          = "$baseDir/data"
params.out           = "$baseDir/results"
params.extension     = "fastq.gz"
params.bundle        = "$baseDir/templates/b37_files_minimal.txt"
//params.mode          = "do.GetContainers" // DIFFENT MODES: do.GetContainers | do.GenomeIndexing | do.QC | do.Trimming | do.Alignment
params.trim          = "ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:28 MINLEN:40"
params.resources     = "/global/blast/gatk-bundle/b37"

// SET PARAMETERS
data_dir             = file(params.data, type: 'dir')
out_dir              = file(params.out, type: 'dir')
ext                  = params.extension
trim_params          = params.trim
resources            = file(params.resources, type: 'dir')

genome               = file("${resources}/human_g1k_v37_decoy.fasta", type: 'file')
dbsnp_sites          = file("${resources}/dbsnp_138.b37.vcf", type: 'file')
b37_bundle           = file(params.bundle, type: 'file')

out_dir.mkdir()

// GET DATA
read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.${ext}", type: 'file')

// START PROCESSING READS
switch (params.mode) {
    // Download the Singularity images required to execute this workflow! 
    case['do.GetContainers']:
        println "\nDownloading Singularity containers."
        
        base = "shub://h3abionet/h3avarcall:"
        shub_images = Channel.from( ["${base}gatk", "${base}bwa", "${base}trimmomatic", "${base}fastqc"] )
        
        process run_DownloadContainers {
            cpus 1
            memory '2 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { "Downloading: $link" }
            publishDir "$baseDir/containers", mode: 'copy', overwrite: true
            
            input:
            each link from shub_images
            
            output:
            file("*") into containers
            
            """
            singularity pull ${link}
            """
        }  
        break
        // --------------------
        
    case['do.GenomeIndexing']:
        println "\nDownloading reference genome and indexing"

        process run_GenomeIndexing {
            cpus 1
            memory '5 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$baseDir/resources", mode: 'copy', overwrite: false
            
            output:
            file("*.gz") into gatk_bundle

            shell:
            b37_list = "${b37_bundle}"
            template 'download_bundles.sh'
        }
        break
        // --------------------

    case['do.QC']:
        println "\nPerforming QC for all the samples"
        
        process run_QualityChecks {
            cpus 11
            memory '10 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: true
            
            input:
            set sample, file(reads) from read_pairs
            
            output:
            set sample, file("${sample}*.html") into qc_results
            set sample, file(reads) into read_pairs_qcd
            
            """
            fastqc ${reads.get(0)} ${reads.get(1)} \
                --threads 10 \
                --noextract
            """
        }
        break
        // --------------------

    case['do.Trimming']:
        println "Performing trimming for all the samples"

        process run_ReadTrimming {
            cpus 11
            memory '50 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: true
            
            input:
            set sample, file(reads) from read_pairs
            
            output:
            set sample, file("${sample}*{1,2}P*") into read_pairs_trimmed
            
            """
            java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
                -threads 10 \
                -trimlog trimlog_${sample}.log \
                -basein ${reads.get(0)} \
                -baseout ${sample}_trimmed.fastq.gz \
                \$(sed 's|ILLUMINACLIP:|ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/|' <<< "${trim_params}")
            """
        }
        break
        // --------------------        
        
    case['do.Alignment']:
        println "Performing alignment for all the samples"

        process run_ReadAlignment {
            cpus 11
            memory '50 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file(reads) from read_pairs
            
            output:
            set sample, file("${sample}*") into bam_results
            
            """
            bwa mem -R \"@RG\\tID:${sample}\\tLB:LIBA\\tSM:${sample}\\tPL:Illumina\" \
                -t 10 -M ${genome} \
                ${reads.get(0)} ${reads.get(1)} | samtools sort \
                --threads 10 - > ${sample}.bam
            """
        }

        process run_MarkDuplicates {
            cpus 8
            memory '5 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file("*") from bam_results

            output:
            set sample, file("${sample}_md.{bam,bai}") into md_bam_results_crt, md_bam_results_rb

            """
            gatk --java-options \"-Xmx4G\" MarkDuplicates \
                --MAX_RECORDS_IN_RAM 50000 \
                --INPUT ${sample}.bam \
                --METRICS_FILE ${sample}.metrics \
                --ASSUME_SORT_ORDER coordinate \
                --CREATE_INDEX true \
                --OUTPUT ${sample}_md.bam
            """
        }

        process run_CreateRecalibrationTable {
            cpus 8
            memory '5 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file("*") from md_bam_results_crt

            output:
            set sample, file("${sample}_recal.table") into recalibration_results

            """
            gatk --java-options \"-Xmx4G\" BaseRecalibrator \
                --input ${sample}_md.bam \
                --output ${sample}_recal.table \
                -R ${genome} \
                --known-sites ${resources}/dbsnp_138.b37.vcf \
                --known-sites ${resources}/1000G_phase1.indels.b37.vcf \
                --known-sites ${resources}/Mills_and_1000G_gold_standard.indels.b37.vcf
            """
        }
        
        process run_RecalibrateBAM {
            cpus 8
            memory '5 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file("*") from recalibration_results
            set sample, file("*") from md_bam_results_rb

            output:
            set sample, file("${sample}_md.recal.{bam,bai}") into recalibrated_results_stats, recalibrated_results
            
            """
            gatk --java-options \"-Xmx4G\" ApplyBQSR \
                 --input ${sample}_md.bam \
                 --output ${sample}_md.recal.bam \
                 -R ${genome} \
                 --create-output-bam-index true \
                 --bqsr-recal-file ${sample}_recal.table
            """
        }

        process run_SamtoolsStats {
            cpus 11
            memory '50 GB'
            time '2h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file("*") from recalibrated_results_stats

            output:
            set sample, file("${sample}_md.recal.stats")  into recal_stats

            """
            samtools stats \
                --threads 10 \
                ${sample}_md.recal.bam > ${sample}_md.recal.stats
            """
        }

        break
        // --------------------

    case['do.HaplotypeCaller']:
        println "Performing something for all the samples"    
        recalibrated_results = Channel.fromFilePairs("$out_dir/NIST7035_TAAGGCGA_L001/*_md.recal{.bam,.bai}", type: 'file', size: -1)
        chromosomes = Channel.from ( 1..22 )

        process run_HaplotypeCaller {
            cpus 1
            memory '10 GB'
            time '5h'
            scratch '$HOME/tmp'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file("*") from recalibrated_results
            each chrom from chromosomes

            output:
            set val("chr_${chrom}"), file("*{.gz,.gz.tbi}") into HaplotypeCaller_results
            
            """
            gatk --java-options \"-Xmx4G\" HaplotypeCaller \
                -R ${genome} \
                -I ${sample}_md.recal.bam \
                --emit-ref-confidence GVCF \
                --dbsnp ${dbsnp_sites} \
                --L ${chrom} \
                --genotyping-mode DISCOVERY \
                -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
                -stand-call-conf 30 \
                --sample-ploidy 2 \
                -O ${sample}_${chrom}.g.vcf.gz
            """
        }

        // HaplotypeCaller_results.subscribe { println "$it" }

        // HaplotypeCaller_results.collectFile(sort: 'true') {
        //     item -> [ "${item[0]}.txt", "${item[1]}" + '\n' ]
        // }.set { sample_gvcf_list }

       HaplotypeCaller_results
           .groupTuple(by: 0, sort: 'true')
           .set { HaplotypeCaller_per_chrom }
//            .subscribe { println "$it" }

        
        // // sample_gvcf_list.subscribe {
        // //     println "Contents of $it:\n"
        // //     println "$it.text"
        // // }

        // process run_GenotypeGVCF {
        //     cpus 1
        //     memory '10 GB'
        //     time '2h'
        //     scratch '$HOME/tmp'
        //     tag { sample }
        //     publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
        //     input:
        //     set chrom, file(list) from HaplotypeCaller_per_chrom
            
        //     output:
        //     set chrom, file("*") into combined_gvcf

        // """
        // gatk --java-options \"-Xmx4G\" GenotypeGVCFs \
        //     -R ${genome} \
        //     -L ${chrom.substring(2,)} \
        //     -V ${list} \
        //     -stand-call-conf 30 \
        //     -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
        //     -O "${chrom.substring(2,)}.g.vcf.gz" 
        // """
        // }

        break

}

//combined_gvcf.subscribe { println "$it" }

// lftp -e 'mirror --use-pget-n=10 /bundle/hg19 .' -u gsapubftp-anonymous, ftp.broadinstitute.org
// lftp -e 'pget -n20 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg19/dbsnp_138.hg19.vcf.gz; bye'
