echo true

process RunTrim_sequencesTask {	
	
	output:
	stdout into Trim_sequencesOutput	

	shell:
	"""
	nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/trim_sequences.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/design_block1.config  
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
	nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/alignment.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/design_block1.config --TrimmedInputRead1 ${params.Trim_sequencesOutputDirectory}"/"${params.SampleName}".read1.trimmed.fq.gz" --TrimmedInputRead2 ${params.Trim_sequencesOutputDirectory}"/"${params.SampleName}".read2.trimmed.fq.gz"
	"""
	else if(params.PairedEnd == 'false')
        """
        nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/alignment.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/design_block1.config --TrimmedInputRead1 ${params.Trim_sequencesOutputDirectory}"/"${params.SampleName}".read1.trimmed.fq.gz"  --TrimmedInputRead2 null
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
        nextflow run /projects/bioinformatics/PrakruthiWork/Genomics_MGC_VariantCalling_Nextflow/src/nextflow/Tasks/dedup.nf -c /projects/bioinformatics/PrakruthiWork/NextflowConfig/design_block1.config --InputAlignedSortedBam ${params.AlignmentOutputDirectory}"/"${params.SampleName}".aligned.sorted.bam" --InputAlignedSortedBamBai ${params.AlignmentOutputDirectory}"/"${params.SampleName}".aligned.sorted.bam.bai"

        """
}

