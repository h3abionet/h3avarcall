#!/usr/bin/env nextflow

// HELP MENU GOES HERE

// STAGING INPUT FILES AND SETTING PARAMETERS FROM MAIN.CONFIG
data_dir             = file(params.data, type: 'dir')
out_dir              = file(params.out, type: 'dir')
trim_params          = params.trim
resources            = file(params.resources, type: 'dir')
mode                 = params.mode
resume_from          = params.from
//==========
b37_bundle           = file(params.bundle, type: 'file')
genome               = file("${resources}/human_g1k_v37_decoy.fasta", type: 'file')
dbsnp_sites          = file("${resources}/dbsnp_138.b37.vcf", type: 'file')
hapmap               = file("${resources}/hapmap_3.3.b37.vcf", type: 'file')
omni                 = file("${resources}/1000G_omni2.5.b37.vcf", type: 'file')
phase1_indels        = file("${resources}/1000G_phase1.indels.b37.vcf", type: 'file')
phase1_snps          = file("${resources}/1000G_phase1.snps.high_confidence.b37.vcf", type: 'file')
golden_indels        = file("${resources}/Mills_and_1000G_gold_standard.indels.b37.vcf", type: 'file')
ext                  = "fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2"
//bind_dir             = "/" + params.resources.split('/')[1].toString()
bind_dir             = [params.data, params.out, params.resources].collect { it -> "-B ${it}"}.join(" ").toString()

// OUTPUT DIRECTORIES
out_dir.mkdir()
//==========
qc_dir             = file("${out_dir}/1_QC", type: 'dir') 
trim_dir           = file("${out_dir}/2_Read_Trimming", type: 'dir') 
align_dir          = file("${out_dir}/3_Read_Alignment", type: 'dir') 
var_call_dir       = file("${out_dir}/4_Variant_Calling", type: 'dir') 
var_filter_dir     = file("${out_dir}/5_Variant_Filtering", type: 'dir') 
multi_qc_dir       = file("${out_dir}/MultiQC", type: 'dir') 

// INPUT ERROR MESSAGES!
main_data_error  = """
=============================================================================================
Ooops!! Looks like there's an ERROR in your input files! There are no FASTQ files in the directory:
\t${data_dir}
Please ensure that you have give the correct DIRECTORY for you FASTQ input files using the '--data' option!
=============================================================================================
"""
trim_data_error  = """
=============================================================================================
Ooops!! Looks like there's an ERROR in your input files! There are no FASTQ files in the directory:
\t${trim_dir}
Are you sure you ran the READ TRIMMING STEP (--mode do.ReadTrimming) ??
Please ensure that you have ran the READ TRIMMING STEP successfully and try again!
=============================================================================================
"""
align_data_error ="""
=============================================================================================
Ooops!! Looks like there's an ERROR in your input files! There are no BAM files in the directory:
\t${align_dir}
Are you sure you ran the READ ALIGNMENT STEP (--mode do.ReadAlignment) ??
Please ensure that you have ran the READ ALIGNMENT STEP successfully and try again!
=============================================================================================
"""
call_data_error  = """
=============================================================================================
Ooops!! Looks like there's an ERROR in your input files! There are no GVCF files in the directory:
\t${var_call_dir}
Are you sure you ran the VARIANT CALLING STEP (--mode do.VariantCalling) ??
Please ensure that you have ran the VARIANT CALLING STEP successfully and try again!
=============================================================================================
"""
multi_qc_error = """
=============================================================================================
Ooops!! Looks like there's an ERROR in your input files! There are NO FILES in the directory:
\t${out_dir}
You are trying to run MultiQC, but it seems like you have NO OUTPUT FILES, most probably because you haven't run any step of the workflow!
Please ensure that you have ran at least ONE STEP of the workflow successfully and try again!
=============================================================================================
"""

// GET DATA - SPIT OUT ERROR MESSAGES IF THERE IS SOMETHING WRONG WITH THE INPUT OR COMMAND OPTIONS
if(mode == null || mode == 'do.GetContainers' || mode == 'do.GenomeIndexing' ) {
} else if(mode == 'do.QC') {
    if(resume_from == null) {
        read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.{$ext}", type: 'file')
            .ifEmpty {
            error "$main_data_error"
        }
    } else if(resume_from == 'do.Trimming') {
        read_pairs = Channel.fromFilePairs("${trim_dir}/*[1,2]P.{$ext}", type: 'file')
            .ifEmpty {
            error "$trim_data_error"
        }
    } else {
        println "\n============================================================================================="
        println "Ooops!! Looks like there's an ERROR in you command!"
        println "I do not recognise the \'--from ${resume_from}\' option you have given me!"
        println "The allowed options for \'--from\' that can be used with \'--mode ${mode}\' are:"
        println "\t[ do.ReadTrimming ]"
        println "Please use one of the above options, or leave the \'--from\' option out, so I can continue!"
        println "=============================================================================================\n"
        exit 1
    }
} else if(mode == 'do.ReadTrimming') {
    if(resume_from == null) {
        read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.{$ext}", type: 'file')
            .ifEmpty {
            error "$main_data_error"
        }
    } else if(resume_from == 'do.QC') {
        read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.{$ext}", type: 'file')
            .ifEmpty {
            error "$main_data_error"
        }
    } else {
        println "\n============================================================================================="
        println "Ooops!! Looks like there's an ERROR in you command!"
        println "I do not recognise the \'--from ${resume_from}\' option you have given me!"
        println "The allowed options for \'--from\' that can be used with \'--mode ${mode}\' are:"
        println "\t[ do.QC ]"
        println "Please use one of the above options, or leave the \'--from\' option out, so I can continue!"
        println "=============================================================================================\n"
        exit 1
    }
} else if(mode == 'do.ReadAlignment') {
    if(resume_from == null || resume_from == 'do.QC') {
        read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.{$ext}", type: 'file')
            .ifEmpty {
            error "$main_data_error"
        }
    } else if(resume_from == 'do.Trimming') {
        read_pairs = Channel.fromFilePairs("${trim_dir}/*[1,2]P.{$ext}", type: 'file')
            .ifEmpty {
            error "$trim_data_error"
        }
    } else {
        println "\n============================================================================================="
        println "Ooops!! Looks like there's an ERROR in you command!"
        println "I do not recognise the \'--from ${resume_from}\' option you have given me!"
        println "The allowed options for \'--from\' that can be used with \'--mode ${mode}\' are:"
        println "\t[ do.QC | do.Trimming ]"
        println "Please use one of the above options, or leave the \'--from\' option out, so I can continue!"
        println "=============================================================================================\n"
        exit 1
    }
} else if(mode == 'do.VariantCalling') {
    if(resume_from == null || resume_from == 'do.Alignment') {
        bam_recal = Channel.fromFilePairs("${align_dir}/*_md.recal.{bai,bam}")
            .ifEmpty {
            error "$align_data_error"
        }
    } else {
        println "\n============================================================================================="
        println "Ooops!! Looks like there's an ERROR in you command!"
        println "I do not recognise the \'--from ${resume_from}\' option you have given me!"
        println "The allowed options for \'--from\' that can be used with \'--mode ${mode}\' are:"
        println "\t[ do.Alignment ]"
        println "Please us the above option to resume from the ALIGNMENT STEP!"
        println "=============================================================================================\n"
        exit 1
    }
} else if(mode == 'do.VariantFiltering') {
    if(resume_from == null || resume_from == 'do.VariantCalling') {
        genotype_vcf_list = Channel.fromFilePairs("${var_call_dir}/*genotyped.vcf.{gz,gz.tbi}", size: -1) { it -> "${it.name.substring(0,3)}" }
            .ifEmpty {
            error "$call_data_error"
        }
    } else {
        println "\n============================================================================================="
        println "Ooops!! Looks like there's an ERROR in you command!"
        println "I do not recognise the \'--from ${resume_from}\' option you have given me!"
        println "The allowed options for \'--from\' that can be used with \'--mode ${mode}\' are:"
        println "\t[ do.VariantCalling ]"
        println "Please use the above option to resume from the VARIANT CALLING STEP!"
        println "=============================================================================================\n"
        exit 1
    }
} else if(mode == 'do.MultiQC') {
    multi_qc = Channel.fromPath("${out_dir}")
        .ifEmpty {
        error "$multi_qc_error"
    }
} else {
    println "\n============================================================================================="
    println "Ooops!! Looks like there's an ERROR in you command!"
    println "I do not recognise the \'--mode ${mode}\' option you have given me!"
    println "The allowed options for \'--from\' that can be used with \'--mode ${mode}\' are:"
    println "\t[ do.GetContainers | do.GenomeIndexing | do.QC | do.ReadTrimming | do.ReadAlignment | do.VarianCalling | do.VariantFiltering ]"
    println "Please us the above option to resume from the VARIANT CALLING STEP!"
    println "=============================================================================================\n"
    exit 1
}

/*  ======================================================================================================
 *  RUN INFO
 *  ======================================================================================================
 */
log.info "====================================      "
log.info "          H3AVarCall                      "
log.info "===================================="
log.info "Input data                : ${data_dir}"
log.info "Output path               : ${out_dir}"
log.info "GATK-b37-bundle directory : ${resources}"
log.info "Path to bind              : ${bind_dir} "
log.info "====================================\n"

// START PROCESSING READS
switch (mode) {
        // PREPROCESSING - DOWNLOAD THE SINGULARITY IMAGES REQUIRED TO EXECUTE THIS WORKFLOW!
    case['do.GetContainers']:
        base = "shub://h3abionet/h3avarcall:"
        shub_images = Channel.from( ["${base}gatk", "${base}bwa", "${base}trimmomatic", "${base}fastqc"] ) //"${base}multiqc"
        
        process run_DownloadContainers {
            label 'noimage'
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
        
        // PREPROCESSING - DOWNLOAD GATK-B37-BUNDLE (ONLY THE NECESSARY FILES) AND INDEX GENOME WITH BWA. 
    case['do.GenomeIndexing']:
        process run_GenomeIndexing {
            label 'bwa'
            tag { "Downloading b37" }
            publishDir "$baseDir/resources", mode: 'copy', overwrite: true
            
            output:
            file("*.gz") into gatk_bundle

            shell:
            b37_list = "${b37_bundle}"
            template 'download_bundles.sh'
        }
        break
        // --------------------

        // MAIN WORKFLOW - STEP 1 (OPTIONAL): PERFORM QC ON INPUT FASTQ FILES! 
    case['do.QC']:
        process run_QualityChecks {
            label 'fastqc'
            tag { sample }
            publishDir "${qc_dir}", mode: 'copy', overwrite: true
            
            input:
            set sample, file(reads) from read_pairs
            
            output:
            set sample, file("${sample}*.html") into qc_html
            set sample, file("${sample}*.zip") into qc_multiqc

            """
            fastqc ${reads.get(0)} ${reads.get(1)} \
                --threads ${task.cpus} \
                --noextract
            """
        }
        break
        // --------------------

        // MAIN WORKFLOW - STEP 2 (OPTIONAL): TRIMMING OF INPUT FASTQ FILES
    case['do.ReadTrimming']:
        process run_ReadTrimming {
            label 'trimmomatic'
            tag { sample }
            publishDir "${trim_dir}", mode: 'copy', overwrite: true
            
            input:
            set sample, file(reads) from read_pairs
            
            output:
            set sample, file("${sample}*{1,2}P*") into read_pairs_trimmed

            """
            java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
                -threads ${task.cpus} \
                -trimlog trimlog_${sample}.log \
                -basein ${reads.get(0)} \
                -baseout ${sample}_trimmed.fastq.gz \
                \$(sed 's|ILLUMINACLIP:|ILLUMINACLIP:/opt/Trimmomatic-0.39/adapters/|' <<< "${trim_params}")
            """
        }
        break
        // --------------------        
        
        // MAIN WORKFLOW - STEP 3 (COMPULSORY): ALIGNMENT OF FASTQ FILES TO THE REFERENCE GENOME
    case['do.ReadAlignment']:
        println "Performing alignment for all the samples"

        process run_ReadAlignment {
            label 'bwa'
            tag { sample }
            
            input:
            set sample, file(reads) from read_pairs
            
            output:
            set sample, file("${sample}.{bam,bai}") into bam
            
            """
            bwa mem -R \"@RG\\tID:${sample}\\tLB:LIBA\\tSM:${sample}\\tPL:Illumina\" \
                -t ${task.cpus} -M ${genome} \
                ${reads.get(0)} ${reads.get(1)} | samtools sort \
                --threads 10 - > ${sample}.bam
            """
        }

        process run_MarkDuplicates {
            label 'gatk'
            tag { sample }
            
            input:
            set sample, file(list) from bam

            output:
            set sample, file("${sample}_md.{bam,bai}") into bam_md

            """
            gatk --java-options \"-Xmx${task.memory.toGiga()}G\" MarkDuplicates \
                --MAX_RECORDS_IN_RAM 50000 \
                --INPUT ${list.find { it =~ '.bam$' } } \
                --METRICS_FILE ${sample}.metrics \
                --ASSUME_SORT_ORDER coordinate \
                --CREATE_INDEX true \
                --OUTPUT ${sample}_md.bam
            """
        }

        process run_CreateRecalibrationTable {
            label 'gatk'
            tag { sample }
            
            input:
            set sample, file(list_bam) from bam_md

            output:
            set sample, file(list_bam), file("${sample}_recal.table") into recal_table 

            """
            gatk --java-options \"-Xmx${task.memory.toGiga()}G\" BaseRecalibrator \
                --input ${list_bam.find { it =~ '_md.bam$' } } \
                --output ${sample}_recal.table \
                -R ${genome} \
                --known-sites ${dbsnp_sites} \
                --known-sites ${phase1_indels} \
                --known-sites ${golden_indels}
            """
        }
  
        process run_RecalibrateBAM {
            label 'gatk'
            tag { sample }
            publishDir "${align_dir}", mode: 'copy', overwrite: true
            
            input:
            set sample, file(list_bam), file(table) from recal_table

            output:
            set sample, file("${sample}_md.recal.{bam,bai}") into bam_recal, bam_recal_stats

            """
            gatk --java-options \"-Xmx${task.memory.toGiga()}G\" ApplyBQSR \
                 --input ${list_bam.find { it =~ '_md.bam$' } } \
                 --output ${sample}_md.recal.bam \
                 -R ${genome} \
                 --create-output-bam-index true \
                 --bqsr-recal-file ${table}
            """
        }

        process run_SamtoolsStats {
            label 'bwa'
            tag { sample }
            publishDir "$out_dir/3_Read_Alignment", mode: 'copy', overwrite: true
            
            input:
            set sample, file(list_bam) from bam_recal_stats

            output:
            set sample, file("${sample}_md.recal.stats")  into samtools_stats, samtools_multiqc

            """
            samtools stats \
                --threads ${task.cpus} \
                ${list_bam.find { it =~ '_md.recal.bam$' } } > ${sample}_md.recal.stats
            """
        }
        break
        // --------------------

        // MAIN WORKFLOW - STEP 4 (COMPULSORY): VARIANT CALLING ON ALIGNED/MARKDUPLICATE/RECALIBRATED BAM FILES
    case['do.VariantCalling']:
        chromosomes = Channel.from ( 1..22 )
        process run_HaplotypeCaller {
            label 'gatk'
            tag { "${sample}_chr_${chrom}" }
            
            input:
            set sample, file(list_bam) from bam_recal
            each chrom from chromosomes

            output:
            set val("chr_${chrom}"), file("*.g.vcf.{gz,gz.tbi}") into haplotype_calls mode flatten
            
            """
            gatk --java-options \"-Xmx${task.memory.toGiga()}G\" HaplotypeCaller \
                -R ${genome} \
                -I ${list_bam.find { it =~ '_md.recal.bam$' } } \
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

        // ====================== COLLECTION =========================== 
        haplotype_calls
            .groupTuple(by: 0, sort: 'true')
            .set { haplotype_calls_chrom }
        // ====================== COLLECTION =========================== 
                
        process run_CombineGVCF {
            label 'gatk'
            tag { "chr_${chrom}" }
            
            input:
            set chrom, file(list_gvcf) from haplotype_calls_chrom
            
            output:
            set chrom, file("*.g.vcf.{gz,gz.tbi}") into gvcfs_chrom
            
        """
        gatk --java-options \"-Xmx${task.memory.toGiga()}G\" CombineGVCFs \
            -R ${genome} \
            -L ${chrom.substring(4,)} \
            -V ${list_gvcf.findAll { it =~ '.g.vcf.gz$' }.join(' -V ') } \
            -O "${chrom}.g.vcf.gz"
        """            
        }
        
        process run_GenotypeGVCF {
            label 'gatk'
            tag { "chr_${chrom}" }
            publishDir "${var_call_dir}", mode: 'copy', overwrite: true
            
            input:
            set chrom, file(list_gvcf) from gvcfs_chrom
            
            output:
            set val("genome"), file("*.vcf.{gz,gz.tbi}") into genotype_vcfs_chrom mode flatten

        """
        gatk --java-options \"-Xmx${task.memory.toGiga()}G\" GenotypeGVCFs \
            -R ${genome} \
            -L ${chrom.substring(4,)} \
            -V ${list_gvcf.findAll { it =~ '.g.vcf.gz$' }.join() } \
            -stand-call-conf 30 \
            -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
            -O "${chrom}_genotyped.vcf.gz"
        """
        }
        
        // ====================== COLLECTION =========================== 
        genotype_vcfs_chrom
            .groupTuple(by: 0, sort: 'true')
            .set { genotype_vcf_list }
        // ====================== COLLECTION =========================== 

        break
        // --------------------
        
        // MAIN WORKFLOW - STEP 5 (COMPULSORY): VARIANT FILTERING ON gVCF
    case['do.VariantFiltering']:
        println "Performing something for all the samples"        
        
        process run_CombineChromVCFs {
            label 'gatk'
            tag { "Genome" }
            
            input:
            set tuple, file(list_vcf) from genotype_vcf_list
            
            output:
            file("genome.vcf.{gz,gz.tbi}") into genome_genotype_vcf
            
        """
        gatk --java-options \"-Xmx${task.memory.toGiga()}G\" GatherVcfs \
            -R ${genome} \
            -I ${list_vcf.findAll { it =~ '.vcf.gz$' }.collect { (it=~/\d+|\D+/).findAll() }.toSorted().collect{ it.join() }.join(' -I ') } \
            -O "genome.vcf.gz"
        tabix -p vcf genome.vcf.gz
        """
        }

        process run_VQSRonSNPs {
            label 'gatk'
            tag { "Genome" }

            input:
            file(list_vcf) from genome_genotype_vcf

            output:
            set val("genome"), file(list_vcf), file("*.{recal,recal.idx,tranches}") into vqsr_snp_recal

        """
        gatk --java-options \"-Xmx${task.memory.toGiga()}G\" VariantRecalibrator \
            -R ${genome} \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
            -resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni} \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${phase1_snps} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp_sites} \
            -an DP -an FS -an SOR -an MQ -an MQRankSum -an QD -an ReadPosRankSum \
            -mode SNP --max-gaussians 4 \
            -V ${list_vcf.findAll { it =~ '.vcf.gz$' }.join(' -V ') } \
            -O genome.recal-SNP.recal \
            --tranches-file genome.recal-SNP.tranches
        """
        }

        process run_ApplyVQSRonSNPs {
            label 'gatk'
            tag { "Genome" }
            publishDir "${var_filter_dir}", mode: 'copy', overwrite: true

            input:
            set tuple_name, file(list_vcf), file(list_recal) from vqsr_snp_recal

            output:
            file("genome.SNP-recal.vcf.{gz,gz.tbi}") into vqsr_snp_apply

        """
        gatk --java-options \"-Xmx${task.memory.toGiga()}G\" ApplyVQSR \
            -R ${genome} \
            --recal-file ${list_recal.find { it =~ '.recal$' } } \
            --tranches-file ${list_recal.find { it =~ '.tranches$' } } \
            -mode SNP \
            -ts-filter-level 99.5 \
            -V ${list_vcf.find { it =~ '.vcf.gz$' } } \
            -O genome.SNP-recal.vcf.gz
        """
        }

        process run_VQSRonINDELs {
            label 'gatk'
            tag { "Genome" }

            input:
            file(list_vcf) from vqsr_snp_apply

            output:
            set val("genome"), file(list_vcf), file("*.{recal,recal.idx,tranches}") into vqsr_indel_recal

        """
        gatk --java-options \"-Xmx${task.memory.toGiga()}G\" VariantRecalibrator \
            -R ${genome} \
            -resource:mills,known=false,training=true,truth=true,prior=12.0 ${golden_indels} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp_sites} \
            -an DP -an FS -an SOR -an MQ -an MQRankSum -an QD -an ReadPosRankSum \
            -mode INDEL --max-gaussians 4 \
            -V ${list_vcf.find { it =~ 'vcf.gz$' } } \
            -O genome.recal-SNP.recal-INDEL.recal \
            --tranches-file genome.recal-SNP.recal-INDEL.tranches
        """
        }

        process run_ApplyVQSRonINDELs {
            label 'gatk'
            tag { "Genome" }
            publishDir "${var_filter_dir}", mode: 'copy', overwrite: true

            input:
            set tumple_name, file(list_vcf), file(list_recal) from vqsr_indel_recal

            output:
            file("genome.SNP-recal.INDEL-recal.vcf.{gz,gz.tbi}") into vqsr_indel_apply

        """
        gatk --java-options \"-Xmx${task.memory.toGiga()}G\" ApplyVQSR \
            -R ${genome} \
            --recal-file ${list_recal.find { it =~ '.recal$' } } \
            --tranches-file ${list_recal.find { it =~ '.tranches$' } } \
            -mode INDEL \
            -ts-filter-level 99.0 \
            -V ${list_vcf.find { it =~ '.vcf.gz$' } } \
            -O genome.SNP-recal.INDEL-recal.vcf.gz
        """
        }
        break
        // --------------------
        
    case['do.MultiQC']:
        process run_MultiQC {
            tag { sample }
            publishDir "${multi_qc_dir}", mode: 'copy', overwrite: true
            
            input:
            file(folder) from multi_qc
            
            output:
            file("*") into multi_qc_output
            
            """
            multiqc ${folder}
            """
        }
        break
        // --------------------        
}
