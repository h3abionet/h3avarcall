#!/usr/bin/env nextflow

// HELP MENU GOES HERE

// SET PARAMETERS
data_dir             = file(params.data, type: 'dir')
out_dir              = file(params.out, type: 'dir')
trim_params          = params.trim
resources            = file(params.resources, type: 'dir')
mode                 = params.mode
resume_from          = params.from

b37_bundle           = file(params.bundle, type: 'file')
genome               = file("${resources}/human_g1k_v37_decoy.fasta", type: 'file')
dbsnp_sites          = file("${resources}/dbsnp_138.b37.vcf", type: 'file')
hapmap               = file("${resources}/hapmap_3.3.b37.vcf", type: 'file')
omni                 = file("${resources}/1000G_omni2.5.b37.vcf", type: 'file')
phase1_indels        = file("${resources}/1000G_phase1.indels.b37.vcf", type: 'file')
phase1_snps          = file("${resources}/1000G_phase1.snps.high_confidence.b37.vcf", type: 'file')
golden_indels        = file("${resources}/Mills_and_1000G_gold_standard.indels.b37.vcf", type: 'file')

out_dir.mkdir()

// // GET DATA
if(mode == null){
    
} else if(mode == 'do.GetContainers') {
    
} else if(mode == 'do.GenomeIndexing') {
    
} else if(mode == 'do.QC') {
    if(resume_from == null) {
        read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.{fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2}", type: 'file')
    } else if(resume_from == 'do.Trimming') {
        read_pairs = Channel.fromFilePairs("${out_dir}/2_Read_Trimming/*{R,read,P}[1,2]*.{fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2}", type: 'file')
    } else {
        exit 1,"DAMMIT DUDE, DO READ ALIGNMENT!!"
    }
} else if(mode == 'do.Trimming') {
    if(resume_from == null) {
        read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.{fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2}", type: 'file')
    } else if(resume_from == 'do.QC') {
        read_pairs = Channel.fromFilePairs("${out_dir}/2_Read_Trimming/*{R,read,P}[1,2]*.{fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2}", type: 'file')
    } else {
        exit 1,"DAMMIT DUDE, DO READ ALIGNMENT!!"
    }
} else if(mode == 'do.Alignment') {
    if(resume_from == null) {
        read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.{fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2}", type: 'file')
    } else if(resume_from == 'do.QC') {
        read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.{fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2}", type: 'file')
    } else if(resume_from == 'do.Trimming') {
        read_pairs = Channel.fromFilePairs("${out_dir}/2_Read_Trimming/*[1,2]P.{fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2}", type: 'file')
    } else {
        exit 1,"DAMMIT DUDE, DO READ ALIGNMENT!!"
    }
} else if(mode == 'do.VariantCalling') {
    if(resume_from == null) {
        bam_recal = Channel.fromFilePairs("${out_dir}/3_Read_Alignment/*_md.recal.{bai,bam}")
    } else if(resume_from == 'do.Alignment') {
        bam_recal = Channel.fromFilePairs("${out_dir}/3_Read_Alignment/*_md.recal.{bai,bam}")
    } else {
        exit 1, "DAMMIT DUDE, DO READ ALIGNMENT!!"
    }
} else if(mode == 'do.VariantFiltering') {
    if(resume_from == null) {
        genotype_vcf_list = Channel.fromFilePairs("${out_dir}/4_Variant_Calling/*genotyped.vcf.{gz,gz.tbi}", size: -1) { it -> "${it.name.substring(0,3)}" }
    } else if(resume_from == 'do.VariantCalling') {
        genotype_vcf_list = Channel.fromFilePairs("${out_dir}/4_Variant_Calling/*genotyped.vcf.{gz,gz.tbi}", size: -1) { it -> "${it.name.substring(0,3)}" }
    } else {
        exit 1,"DAMMIT DUDE, DO READ ALIGNMENT!!"
    }
} else {
    
}

// START PROCESSING READS
switch (mode) {
        // Download the Singularity images required to execute this workflow! 
    case['do.GetContainers']:
        println "\nDownloading Singularity containers."
        
        base = "shub://h3abionet/h3avarcall:"
        shub_images = Channel.from( ["${base}gatk", "${base}bwa", "${base}trimmomatic", "${base}fastqc"] )
        
        process run_DownloadContainers {
            cpus 1
            memory '2 GB'
            time '2h'
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
            label 'bwa'
            cpus 1
            memory '5 GB'
            time '2h'
            tag { sample }
            publishDir "$baseDir/resources", mode: 'copy', overwrite: true
            
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
            label 'fastqc'
            cpus 11
            memory '10 GB'
            time '2h'
            tag { sample }
            publishDir "$out_dir/1_QC", mode: 'copy', overwrite: true
            
            input:
            set sample, file(reads) from read_pairs
            
            output:
            set sample, file("${sample}*.html") into qc_html
            
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
            label 'trimmomatic'
            cpus 11
            memory '50 GB'
            time '2h'
            tag { sample }
            publishDir "$out_dir/2_Read_Trimming", mode: 'copy', overwrite: true
            
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
            label 'bwa'
            cpus 11
            memory '50 GB'
            time '2h'
            tag { sample }
            
            input:
            set sample, file(reads) from read_pairs
            
            output:
            set sample, file("${sample}.{bam,bai}") into bam
            
            """
            bwa mem -R \"@RG\\tID:${sample}\\tLB:LIBA\\tSM:${sample}\\tPL:Illumina\" \
                -t 10 -M ${genome} \
                ${reads.get(0)} ${reads.get(1)} | samtools sort \
                --threads 10 - > ${sample}.bam
            """
        }

        process run_MarkDuplicates {
            label 'gatk'
            cpus 8
            memory '5 GB'
            time '2h'
            tag { sample }
            
            input:
            set sample, file(list) from bam

            output:
            set sample, file("${sample}_md.{bam,bai}") into bam_md

            """
            gatk --java-options \"-Xmx4G\" MarkDuplicates \
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
            cpus 2
            memory '5 GB'
            time '2h'
            tag { sample }
            
            input:
            set sample, file(list_bam) from bam_md

            output:
            set sample, file(list_bam), file("${sample}_recal.table") into recal_table 

            """
            gatk --java-options \"-Xmx4G\" BaseRecalibrator \
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
            cpus 8
            memory '5 GB'
            time '2h'
            tag { sample }
            publishDir "$out_dir/3_Read_Alignment", mode: 'copy', overwrite: true
            
            input:
            set sample, file(list_bam), file(table) from recal_table

            output:
            set sample, file("${sample}_md.recal.{bam,bai}") into bam_recal, bam_recal_stats
            
            """
            gatk --java-options \"-Xmx4G\" ApplyBQSR \
                 --input ${list_bam.find { it =~ '_md.bam$' } } \
                 --output ${sample}_md.recal.bam \
                 -R ${genome} \
                 --create-output-bam-index true \
                 --bqsr-recal-file ${table}
            """
        }

        process run_SamtoolsStats {
            label 'bwa'
            cpus 11
            memory '50 GB'
            time '2h'
            tag { sample }
            publishDir "$out_dir/3_Read_Alignment", mode: 'copy', overwrite: true
            
            input:
            set sample, file(list_bam) from bam_recal_stats

            output:
            set sample, file("${sample}_md.recal.stats")  into samtools_stats

            """
            samtools stats \
                --threads 10 \
                ${list_bam.find { it =~ '_md.recal.bam$' } } > ${sample}_md.recal.stats
            """
        }
        break
        // --------------------

    case['do.VariantCalling']:
        println "Performing something for all the samples"    

        chromosomes = Channel.from ( 1..22 )

        process run_HaplotypeCaller {
            label 'gatk'
            cpus 1
            memory '10 GB'
            time '5h'
            tag { "${sample}_chr_${chrom}" }
            
            input:
            set sample, file(list_bam) from bam_recal
            each chrom from chromosomes

            output:
            set val("chr_${chrom}"), file("*.g.vcf.{gz,gz.tbi}") into haplotype_calls mode flatten
            
            """
            gatk --java-options \"-Xmx8G\" HaplotypeCaller \
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
            cpus 1
            memory '10 GB'
            time '2h'
            tag { "chr_${chrom}" }
            
            input:
            set chrom, file(list_gvcf) from haplotype_calls_chrom
            
            output:
            set chrom, file("*.g.vcf.{gz,gz.tbi}") into gvcfs_chrom
            
        """
        gatk --java-options \"-Xmx4G\" CombineGVCFs \
            -R ${genome} \
            -L ${chrom.substring(4,)} \
            -V ${list_gvcf.findAll { it =~ '.g.vcf.gz$' }.join(' -V ') } \
            -O "${chrom}.g.vcf.gz"
        """            
        }
        
        process run_GenotypeGVCF {
            label 'gatk'
            cpus 1
            memory '10 GB'
            time '2h'
            tag { "chr_${chrom}" }
            publishDir "$out_dir/4_Variant_Calling", mode: 'copy', overwrite: true
            
            input:
            set chrom, file(list_gvcf) from gvcfs_chrom
            
            output:
            set val("genome"), file("*.vcf.{gz,gz.tbi}") into genotype_vcfs_chrom mode flatten

        """
        gatk --java-options \"-Xmx4G\" GenotypeGVCFs \
            -R ${genome} \
            -L ${chrom.substring(4,)} \
            -V ${list_gvcf.findAll { it =~ '.g.vcf.gz$' }.join() } \
            -stand-call-conf 30 \
            -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
            -O "${chrom}_genotyped.vcf.gz"
        """
        }
        
        genotype_vcfs_chrom
            .groupTuple(by: 0, sort: 'true')
            .set { genotype_vcf_list }

        break
        // --------------------
        
    case['do.VariantFiltering']:
        println "Performing something for all the samples"        
        
        // GATHER!
        process run_CombineChromVCFs {
            label 'gatk'
            cpus 1
            memory '10 GB'
            time '2h'
            tag { "Genome" }
            
            input:
            set tuple, file(list_vcf) from genotype_vcf_list
            
            output:
            file("genome.vcf.{gz,gz.tbi}") into genome_genotype_vcf
            
        """
        gatk --java-options \"-Xmx8G\" GatherVcfs \
            -R ${genome} \
            -I ${list_vcf.findAll { it =~ '.vcf.gz$' }.collect { (it=~/\d+|\D+/).findAll() }.toSorted().collect{ it.join() }.join(' -I ') } \
            -O "genome.vcf.gz"
        tabix -p vcf genome.vcf.gz
        """
        }

        process run_VQSRonSNPs {
            label 'gatk'
            cpus 1
            memory '10 GB'
            time '2h'
            tag { "Genome" }

            input:
            file(list_vcf) from genome_genotype_vcf

            output:
            set val("genome"), file(list_vcf), file("*.{recal,recal.idx,tranches}") into vqsr_snp_recal 

        """
        gatk --java-options \"-Xmx8G\" VariantRecalibrator \
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
            cpus 1
            memory '10 GB'
            time '2h'
            tag { "Genome" }
            publishDir "$out_dir/5_Variant_Filtering", mode: 'copy', overwrite: true

            input:
            set tuple_name, file(list_vcf), file(list_recal) from vqsr_snp_recal

            output:
            file("genome.SNP-recal.vcf.{gz,gz.tbi}") into vqsr_snp_apply

        """
        gatk --java-options \"-Xmx8G\" ApplyVQSR \
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
            cpus 1
            memory '10 GB'
            time '2h'
            tag { "Genome" }

            input:
            file(list_vcf) from vqsr_snp_apply

            output:
            set val("genome"), file(list_vcf), file("*.{recal,recal.idx,tranches}") into vqsr_indel_recal

        """
        gatk --java-options \"-Xmx8G\" VariantRecalibrator \
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
            cpus 1
            memory '10 GB'
            time '2h'
            tag { "Genome" }
            publishDir "$out_dir/5_Variant_Filtering", mode: 'copy', overwrite: true

            input:
            set tumple_name, file(list_vcf), file(list_recal) from vqsr_indel_recal

            output:
            file("genome.SNP-recal.INDEL-recal.vcf.{gz,gz.tbi}") into vqsr_indel_apply

        """
        gatk --java-options \"-Xmx8G\" ApplyVQSR \
            -R ${genome} \
            --recal-file ${list_recal.find { it =~ '.recal$' } } \
            --tranches-file ${list_recal.find { it =~ '.tranches$' } } \
            -mode INDEL \
            -ts-filter-level 99.0 \
            -V ${list_vcf.find { it =~ '.vcf.gz$' } } \
            -O genome.SNP-recal.INDEL-recal.vcf.gz
        """
        }
}
