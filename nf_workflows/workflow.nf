echo true

process RunTrim_sequencesTask {	
	
	output:
	stdout into Trim_sequencesOutput	

	shell:
	"""
	nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/trim_sequences.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/workflow.config  
	"""
}

process RunAlignmentTask {

	input:
	val PreAlignmentFlag from Trim_sequencesOutput 	
		
	output:
	stdout into AlignmentOutput

	script:
	if(params.PairedEnd == 'true')
	"""
	nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/alignment.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/workflow.config --TrimmedInputRead1 ${params.Trim_sequencesOutputDirectory}"/"${params.SampleName}".read1.trimmed.fq.gz" --TrimmedInputRead2 ${params.Trim_sequencesOutputDirectory}"/"${params.SampleName}".read2.trimmed.fq.gz"
	"""
	else if(params.PairedEnd == 'false')
        """
        nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/alignment.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/workflow.config --TrimmedInputRead1 ${params.Trim_sequencesOutputDirectory}"/"${params.SampleName}".read1.trimmed.fq.gz"  --TrimmedInputRead2 null
        """
}

process RunDedupTask {

        input:
        val PreDedupFlag from AlignmentOutput

        output:
        stdout into DedupOutput

        shell:
        """
        ls -l
        nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/dedup.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/workflow.config --InputAlignedSortedBam ${params.AlignmentOutputDirectory}"/"${params.SampleName}".aligned.sorted.bam" --InputAlignedSortedBamBai ${params.AlignmentOutputDirectory}"/"${params.SampleName}".aligned.sorted.bam.bai"

        """
}

process RunRealignmentTask {

	input:
	val PreRealignmentFlag from DedupOutput

        output:
        stdout into RealignmentOutput

        shell:
        """
        nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/realignment.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/workflow.config --RefFai ${params.RefFai} --InputAlignedSortedDedupedBam ${params.DedupOutputDirectory}"/"${params.SampleName}".aligned.sorted.deduped.bam" --InputAlignedSortedDedupedBamBai ${params.DedupOutputDirectory}"/"${params.SampleName}".aligned.sorted.deduped.bam.bai"
        """
}

process RunBqsrTask {

        input:
        val PreBqsrFlag from RealignmentOutput

        output:
        stdout into BqsrOutput

        shell:
        """
        nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/bqsr.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/workflow.config --InputAlignedSortedDedupedRealignedBam ${params.RealignmentOutputDirectory}"/"${params.SampleName}".aligned.sorted.deduped.realigned.bam" --InputAlignedSortedDedupedRealignedBamBai ${params.RealignmentOutputDirectory}"/"${params.SampleName}".aligned.sorted.deduped.realigned.bam.bai"
        """                                                                                                              
}

process RunHaplotyperTask {

        input:
        val PreHaplotyperFlag from BqsrOutput

        output:
        stdout into HaplotyperOutput

        shell:
        """
        nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/haplotyper.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/workflow.config --InputAlignedSortedDedupedRealignedBam ${params.RealignmentOutputDirectory}"/"${params.SampleName}".aligned.sorted.deduped.realigned.bam" --InputAlignedSortedDedupedRealignedBamBai ${params.RealignmentOutputDirectory}"/"${params.SampleName}".aligned.sorted.deduped.realigned.bam.bai" --RecalTable ${params.BQSROutputDirectory}"/"${params.SampleName}".recal_data.table"                                                                                                                       
                                                                                                                         
        """                                                                                                              
}

process RunVqsrTask {

        input:
        val PreVqsrFlag from HaplotyperOutput                                                                            
                                                                                                                         
        output:                                                                                                          
        stdout into VqsrOutput                                                                                           
                                                                                                                         
        shell:                                                                                                           
        """                                                                                                              
        nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/vqsr.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/workflow.config --InputVCF ${params.HaplotyperOutputDirectory}"/"${params.SampleName}".vcf" --InputVCFIdx ${params.HaplotyperOutputDirectory}"/"${params.SampleName}".vcf.idx"
        """                                                                                                              
} 
