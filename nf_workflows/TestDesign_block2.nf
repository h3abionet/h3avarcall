echo true

process RunRealignmentTask {	
	
	output:
	stdout into RealignmentOutput	

	shell:
	"""
	nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/realignment.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/design_block2.config
	"""
}

process RunBqsrTask {
	
	input:
	val PreBqsrFlag from RealignmentOutput

        output:
        stdout into BqsrOutput

        shell:
        """
        nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/bqsr.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/design_block2.config --InputAlignedSortedDedupedRealignedBam ${params.RealignmentOutputDirectory}"/"${params.SampleName}".aligned.sorted.deduped.realigned.bam" --InputAlignedSortedDedupedRealignedBamBai ${params.RealignmentOutputDirectory}"/"${params.SampleName}".aligned.sorted.deduped.realigned.bam.bai"
        """
}

process RunHaplotyperTask {

	input:
	val PreHaplotyperFlag from BqsrOutput

	output:
	stdout into HaplotyperOutput

	shell:
        """
        nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/haplotyper.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/design_block2.config --InputAlignedSortedDedupedRealignedBam ${params.RealignmentOutputDirectory}"/"${params.SampleName}".aligned.sorted.deduped.realigned.bam" --InputAlignedSortedDedupedRealignedBamBai ${params.RealignmentOutputDirectory}"/"${params.SampleName}".aligned.sorted.deduped.realigned.bam.bai" --RecalTable ${params.BQSROutputDirectory}"/"${params.SampleName}".recal_data.table"

	"""
}

process RunVqsrTask {

	input:
	val PreVqsrFlag from HaplotyperOutput

        output:
        stdout into VqsrOutput

        shell:
        """
        nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/vqsr.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/design_block2.config --InputVCF ${params.HaplotyperOutputDirectory}"/"${params.SampleName}".vcf" --InputVCFIdx ${params.HaplotyperOutputDirectory}"/"${params.SampleName}".vcf.idx"
        """
}

