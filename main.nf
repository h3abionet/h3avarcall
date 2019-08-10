#!/usr/bin/env nextflow

// HELP MENU GOES HERE

// PARAMETERS GO HERE! WILL BE IN A SEPARATE CONFIG FILE - JUST TESTINF FOR NOW!
params.data          = "$baseDir/data/welcome_trust"
params.out           = "$baseDir/results"
params.bundle        = "$baseDir/templates/b37_files_minimal.txt"
//params.mode          = "do.GetContainers" // DIFFENT MODES: do.GetContainers | do.GenomeIndexing | do.QC | do.Trimming | do.Alignment
params.trim          = "ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:8:true TRAILING:28 MINLEN:40"
params.resources     = "/global/blast/gatk-bundle/b37"

// SET PARAMETERS
data_dir             = file(params.data, type: 'dir')
out_dir              = file(params.out, type: 'dir')
trim_params          = params.trim
resources            = file(params.resources, type: 'dir')

genome               = file("${resources}/human_g1k_v37_decoy.fasta", type: 'file')
dbsnp_sites          = file("${resources}/dbsnp_138.b37.vcf", type: 'file')
hapmap               = file("${resources}/hapmap_3.3.b37.vcf", type: 'file')
omni                 = file("${resources}/1000G_omni2.5.b37.vcf", type: 'file')
phase1_snps          = file("${resources}/1000G_phase1.snps.high_confidence.b37.vcf", type: 'file')
golden_indels        = file("${resources}/Mills_and_1000G_gold_standard.indels.b37.vcf", type: 'file')
b37_bundle           = file(params.bundle, type: 'file')

out_dir.mkdir()

// GET DATA
//read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]_*trimmed.${ext}", type: 'file')
read_pairs = Channel.fromFilePairs("${data_dir}/*{R,read}[1,2]*.{fastq,fastq.gz,fastq.bz2,fq,fq.gz,fq.bz2}", type: 'file')

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
            tag { "Downloading: $link" }
            publishDir "$baseDir/containers", mode: 'copy', overwrite: false
            
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
            publishDir "$baseDir/resources", mode: 'copy', overwrite: false
            
            output:
            file("*.gz") into gatk_bundle

            shell:
            b37_list = "${b37_bundle}"
            template 'download_bundles.sh'
        }
        // break
        // --------------------

    case['do.QC']:
        println "\nPerforming QC for all the samples"
        
        process run_QualityChecks {
            label 'fastqc'
            cpus 11
            memory '10 GB'
            time '2h'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
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
        // break
        // --------------------

    case['do.Trimming']:
        println "Performing trimming for all the samples"

        process run_ReadTrimming {
            label 'trimmomatic'
            cpus 11
            memory '50 GB'
            time '2h'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
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
        // break
        // --------------------        
        
    case['do.Alignment']:
        println "Performing alignment for all the samples"

        process run_ReadAlignment {
            label 'bwa'
            cpus 11
            memory '50 GB'
            time '2h'
            tag { sample }
            // publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
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
            label 'gatk'
            cpus 8
            memory '5 GB'
            time '2h'
            tag { sample }
            // publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
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
            label 'gatk'
            cpus 1
            memory '5 GB'
            time '2h'
            tag { sample }
            // publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file("*") from md_bam_results_crt

            output:
            set sample, file("${sample}_recal.table") into recalibration_results mode flatten

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
        
        // ====================== COLLECTION =========================== 
        md_bam_results_rb
            .map { item -> [ item[0], item[1][0], item[1][1] ] }
            .set { md_bam_results_rb_flat }
        
        recalibration_results.join(md_bam_results_rb_flat)
            .map { item -> [ item[0], item[1..3] ] }
            .set { recalibration_results_table }
        // ====================== COLLECTION =========================== 

        process run_RecalibrateBAM {
            label 'gatk'
            cpus 8
            memory '5 GB'
            time '2h'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file("*") from recalibration_results_table

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
            label 'bwa'
            cpus 11
            memory '50 GB'
            time '2h'
            tag { sample }
            publishDir "$out_dir/${sample}", mode: 'copy', overwrite: true
            
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
        
    //     recal_stats.subscribe { println }
    //     break
    //     // --------------------

    // case['do.HaplotypeCaller']:
    //     println "Performing something for all the samples"    
    //     recalibrated_results = Channel.fromFilePairs("$out_dir/**/*_md.recal{.bam,.bai}", type: 'file', size: 2)
        chromosomes = Channel.from ( 1..22 )

        process run_HaplotypeCaller {
            label 'gatk'
            cpus 1
            memory '10 GB'
            time '5h'
            tag { sample }
            // publishDir "$out_dir/${sample}", mode: 'copy', overwrite: false
            
            input:
            set sample, file("*") from recalibrated_results
            each chrom from chromosomes

            output:
            set val("chr_${chrom}"), file("*{.gz,.gz.tbi}") into HaplotypeCaller_results mode flatten
            
            """
            gatk --java-options \"-Xmx8G\" HaplotypeCaller \
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

        // ====================== COLLECTION =========================== 
        HaplotypeCaller_results
            .groupTuple(by: 0, sort: 'true')
            .set { HaplotypeCaller_per_chrom }
        // ====================== COLLECTION =========================== 
                
        process run_CombineGVCF {
            label 'gatk'
            cpus 1
            memory '10 GB'
            time '2h'
            tag { sample }
            // publishDir "$out_dir/GVCFs_per_chrom", mode: 'copy', overwrite: false
            
            input:
            set chrom, file(list) from HaplotypeCaller_per_chrom
            
            output:
            set chrom, file("*{.gz,.gz.tbi}") into combined_gvcf
            
        """
        gatk --java-options \"-Xmx4G\" CombineGVCFs \
            -R ${genome} \
            -L ${chrom.substring(4,)} \
            -V ${list.findAll { it =~ '.g.vcf.gz$' }.join(' -V ') } \
            -O "${chrom}.g.vcf.gz" 
        """            
        }
        
        process run_GenotypeGVCF {
            label 'gatk'
            cpus 1
            memory '10 GB'
            time '2h'
            tag { sample }
            publishDir "$out_dir/GVCFs_per_chrom", mode: 'copy', overwrite: false
            
            input:
            set chrom, file(list) from combined_gvcf
            
            output:
            set val("genome"), file("*") into geno_gvcf mode flatten

        """
        gatk --java-options \"-Xmx4G\" GenotypeGVCFs \
            -R ${genome} \
            -L ${chrom.substring(4,)} \
            -V ${list.findAll { it =~ '.g.vcf.gz$' }.join() } \
            -stand-call-conf 30 \
            -A Coverage -A FisherStrand -A StrandOddsRatio -A MappingQualityRankSumTest -A QualByDepth -A RMSMappingQuality -A ReadPosRankSumTest \
            -O "${chrom}_genotyped.g.vcf.gz" 
        """
        }
        
        geno_gvcf
            .groupTuple(by: 0, sort: 'true')
            .set { genome_chr_genotyped }

    //     break
    //     // --------------------
        
    // case['do.VQSRCaller']:
    //     println "Performing something for all the samples"        
        

        // GATHER!
        process run_CombineChromGVCF {
            label 'gatk'
            cpus 1
            memory '10 GB'
            time '2h'
            tag { sample }
            publishDir "$out_dir/GVCF_genome", mode: 'copy', overwrite: true
            
            input:
            set tuple_name, file(chr_gvcf) from genome_chr_genotyped
            
            output:
            file("*") into geno_gvcf_full
            
        """
        gatk --java-options \"-Xmx8G\" GatherVcfs \
            -R ${genome} \
            -I ${chr_gvcf.findAll { it =~ '.g.vcf.gz$' }.collect { (it=~/\d+|\D+/).findAll() }.toSorted().collect{ it.join() }.join(' -I ') } \
            -O "genome_full.g.vcf.gz"
        /home/phelelani/applications/tabix -p vcf genome_full.g.vcf.gz
        """            
        }

        process run_VQSRonSNPs {
            label 'gatk'
            cpus 1
            memory '10 GB'
            time '2h'
            tag { sample }
            publishDir "$out_dir/GVCF_genome", mode: 'copy', overwrite: false

            input:
            file(list) from geno_gvcf_full

            output:
            set val("geonme"), file("genome*"), file(list) into snps_vqsr_recal 

        """
        gatk --java-options \"-Xmx8G\" VariantRecalibrator \
            -R ${genome} \
            -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
            -resource:omni,known=false,training=true,truth=true,prior=12.0 ${omni} \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 ${phase1_snps} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp_sites} \
            -an DP -an FS -an SOR -an MQ -an MQRankSum -an QD -an ReadPosRankSum \
            -mode SNP --max-gaussians 4 \
            -V ${list.findAll { it =~ '.g.vcf.gz$' }.join(' -V ') } \
            -O genome.recal-SNP.recal \
            --tranches-file genome.recal-SNP.tranches
        """
        }

        process run_ApplyVQSRonSNPs {
            label 'gatk'
            cpus 1
            memory '10 GB'
            time '2h'
            tag { sample }
            publishDir "$out_dir/GVCF_genome", mode: 'copy', overwrite: false

            input:
            set tuple_name, file(list), file(list2) from snps_vqsr_recal

            output:
            file("genome.recal-SNP.vcf.{gz,gz.tbi}") into snps_vqsr_vcf

        """
        gatk --java-options \"-Xmx8G\" ApplyVQSR \
            -R ${genome} \
            --recal-file ${list.find { it =~ '.recal$' } } \
            --tranches-file ${list.find { it =~ '.tranches$' } } \
            -mode SNP \
            -ts-filter-level 99.5 \
            -V ${list2.find { it =~ '.g.vcf.gz$' } } \
            -O genome.recal-SNP.vcf.gz
        """
        }

        process run_VQSRonINDELs {
            label 'gatk'
            cpus 1
            memory '10 GB'
            time '2h'
            tag { sample }
            publishDir "$out_dir/GVCF_genome", mode: 'copy', overwrite: false

            input:
            file(list) from snps_vqsr_vcf

            output:
            set val("genome"), file("genome*"), file(list) into indel_vqsr_recal

        """
        gatk --java-options \"-Xmx8G\" VariantRecalibrator \
            -R ${genome} \
            -resource:mills,known=false,training=true,truth=true,prior=12.0 ${golden_indels} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp_sites} \
            -an DP -an FS -an SOR -an MQ -an MQRankSum -an QD -an ReadPosRankSum \
            -mode INDEL --max-gaussians 4 \
            -V ${list.find { it =~ 'vcf.gz$' } } \
            -O "genome.recal-SNP.vcf.recal-INDEL.recal" \
            --tranches-file "genome.recal-SNP.vcf.recal-INDEL.tranches"
        """
        }

        process run_ApplyVQSRonINDELs {
            label 'gatk'
            cpus 1
            memory '10 GB'
            time '2h'
            tag { sample }
            publishDir "$out_dir/GVCF_genome", mode: 'copy', overwrite: false

            input:
            set tumple_name, file(list), file(list2) from indel_vqsr_recal

            output:
            file("genome*.vcf.{gz,gz.tbi}") into indel_vqsr_vcf

        """
        gatk --java-options \"-Xmx8G\" ApplyVQSR \
            -R ${genome} \
            --recal-file ${list.find { it =~ '.recal$' } } \
            --tranches-file ${list.find { it =~ '.tranches$' } } \
            -mode INDEL \
            -ts-filter-level 99.0 \
            -V ${list2.find { it =~ '.vcf.gz$' } } \
            -O genome.recal-SNP.recal-INDEL.vcf.gz
        """
        }
        
}
