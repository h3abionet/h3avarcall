/******************************************************************************/
/*                                                                            */
/*          This Nextflow script is a Germline workflow for 1 sample          */
/*                                                                            */
/******************************************************************************/

/*******        Nextflow option so bash stdout will be displayed       ********/
echo false

/* ******************        Import input variables          **************** */
Platform = params.Platform						// Sequencing platform
Library = params.Library	                    // Library name
PlatformUnit = params.PlatformUnit              // Platform unit/ flowcell ID
CenterName = params.CenterName              	// Name of the sequencing center
PairedEnd = params.PairedEnd					// Is input FASTQ paired ended?
Ref = file(params.Ref)                  		// Reference genome
RefAmb = file(params.RefAmb)                    // Reference indices
RefAnn = file(params.RefAnn)                    // Reference indices
RefBwt = file(params.RefBwt)                    // Reference indices
RefPac = file(params.RefPac)                    // Reference indices
RefSa = file(params.RefSa)                      // Reference indices
RefFai = file(params.RefFai)                    // Reference files- implicit
RefDict = file(params.RefDict)                  // to the GATK tool

BqsrKnownSites = params.BqsrKnownSites          // List of known sites+ dbSNP
BqsrKnownSitesChannel = Channel.from(BqsrKnownSites.tokenize(',')).flatMap{ files(it) }.collect()
BqsrKnownSitesIdxChannel = Channel.from(BqsrKnownSites.tokenize(',')).flatMap{ files(it+'.idx') }.collect()
DBSNP = file(params.DBSNP)                  // DBSNP file
DBSNPIdx = file(params.DBSNPIdx)            // Index file for DBSNP file

BWAExe = file(params.BWAExe)                 // Path to BWA executable
ChunkSizeInBases = params.ChunkSizeInBases 	 // 10000000 to normalize threads
BWAExtraOptionsString = params.BWAExtraOptionsString  // BWA extra options

SamtoolsExe = file(params.SamtoolsExe)       // Path to samtools executable
BwaSamtoolsThreads = params.BwaSamtoolsThreads   // Thread required per run

GATKExe = file(params.GATKExe)                  // GATK executable path
ApplyBQSRExtraOptionsString = params.ApplyBQSRExtraOptionsString
HaplotyperThreads = params.HaplotyperThreads
HaplotyperExtraOptionsString = params.HaplotyperExtraOptionsString
JavaExe = file(params.JavaExe)                  // Java executable path
JavaOptionsString = params.JavaOptionsString    //String of java vm options

BashPreamble = file(params.BashPreamble)        // For zombie processes
BashSharedFunctions = file(params.BashSharedFunctions)	// For variable checks
DebugMode = params.DebugMode					// Flag to enable debug mode

AlignmentScript = file(params.AlignmentScript)	// script running alignment
MergeBamScript = file(params.MergeBamScript)			// script running merge BAMs
DedupScript = file(params.DedupScript)          // script running deduplication
BqsrScript = file(params.BqsrScript)            // script of the bqsr job
HaplotyperScript = file(params.HaplotyperScript) // script of GATK code
MergeGvcfsScript = file(params.MergeGvcfsScript) // script of GATK code


AlignmentMultinode = params.AlignmentMultinode

DeliveryFolder_Alignment = params.DeliveryFolder_Alignment
DeliveryFolder_HaplotyperVC = params.DeliveryFolder_HaplotyperVC


/******************************************************************************/
/*                                                                            */
/*      The following are process definitions and their channels wiring!      */
/*                                                                            */
/******************************************************************************/


/* *********      Alignment: per lane (of a sample) - Required   ************ */
SampleName = params.SampleName					 // Singl sample name
PlatformUnitChannel = Channel
                        .from(PlatformUnit.tokenize(',')) // For each read-pair

InputRead1 = params.InputRead1                              // Input Read File
InputRead2 = params.InputRead2                               // Input Read File
InputRead1Channel = Channel.fromPath(InputRead1.tokenize(','))

InputRead2Channel = ( PairedEnd == 'true'
                      ? Channel.fromPath(InputRead2.tokenize(','))
                      : file('null').fileName )

process Alignment{
  tag "${SampleName}_${PlatformUnit}"

  // publishDir DeliveryFolder_Alignment , mode: 'copy'

	input:
  	val SampleName
  	val Platform
    val Library
    val PlatformUnit from PlatformUnitChannel
    val CenterName
  	val PairedEnd
  	file InputRead1 from InputRead1Channel
  	file InputRead2 from InputRead2Channel
  	file Ref
  	file RefAmb
  	file RefAnn
  	file RefBwt
  	file RefPac
  	file RefSa
    file BWAExe
  	val ChunkSizeInBases
    val BWAExtraOptionsString
    file SamtoolsExe
  	val BwaSamtoolsThreads
    file BashSharedFunctions
   	val DebugMode

    file BashPreamble
  	file AlignmentScript

  output:
      set SampleName, "${SampleName}.${PlatformUnit}.bam",
              "${SampleName}.${PlatformUnit}.bam.bai" into AlignOutput

  script:
     	"""
      
      /bin/bash ${AlignmentScript} -s ${SampleName} -p ${Platform} \
          -L ${Library} -f ${PlatformUnit} -c ${CenterName} -P ${PairedEnd} \
          -l ${InputRead1} -r ${InputRead2} -G ${Ref} -e ${BWAExe} \
          -K ${ChunkSizeInBases} -o \"\'${BWAExtraOptionsString}\'\" \
          -S ${SamtoolsExe} -t ${BwaSamtoolsThreads} -F ${BashSharedFunctions}\
          ${DebugMode}
      """
}

/* ***********    Merge: all lanes (of a sample) - Required      ************ */

process MergeBams {
  tag "${SampleName}_All_lanes"

  publishDir DeliveryFolder_Alignment, mode: 'copy'

	input:
      set SampleName, file(InputBam), file(InputBai) from AlignOutput.groupTuple()    // Link to Alignment
      file SamtoolsExe
      file BashSharedFunctions
      val DebugMode

      file BashPreamble
      file MergeBamScript



  output:
      set SampleName, "${SampleName}.bam",
          "${SampleName}.bam.bai" into MergeBamsOutputToDedup, MergeBamsOutputToBqsr

		"""
      
		/bin/bash ${MergeBamScript} -b ${InputBam.join(',')} -s ${SampleName} \
          -S ${SamtoolsExe} \
          -F ${BashSharedFunctions} ${DebugMode}
		"""
}

/*****************         Dedup: per sample - Optional           *************/
process Deduplication{
 tag "${SampleName}_All_intervals"

 publishDir DeliveryFolder_Alignment, mode: 'copy', overwrite: true

 input:
     set SampleName,
         file(InputBams), file(InputBais) from MergeBamsOutputToDedup.groupTuple()  // Link to MergeBams

     file GATKExe
     val JavaExe
     val JavaOptionsString

     val DebugMode

     file BashPreamble
     file BashSharedFunctions

     file DedupScript

 output:
     set SampleName, "${SampleName}.bam", "${SampleName}.bam.bai" into DedupOutput

 when:
     params.MarkDuplicates == 'true'

 script:
     """
     
     /bin/bash ${DedupScript} -s ${SampleName} -b ${InputBams} -S ${GATKExe} \
          -J ${JavaExe} -e \"\'${JavaOptionsString}\'\" \
          -F ${BashSharedFunctions} {DebugMode}
     """
}

/*********       BQSR: per interval (of a sample) - Required        ***********/

Channel                                         // Chromosome names/intervals
  .from(params.GenomicIntervals.tokenize(','))
  .into{BqsrGenomicIntervals; HCGenomicIntervals; JointCallIntervals}

BQSRInput = (params.MarkDuplicates == 'true'
           ? DedupOutput : MergeBamsOutputToBqsr) // Link MergeBams or Deduplication

process BQSR{
 tag "${SampleName}_${GenomicInterval}"

 input:
      set SampleName, file (InputBams), file (InputBais) from BQSRInput

  	file Ref
    	file RefFai
        file RefDict

  	file BqsrKnownSites from BqsrKnownSitesChannel
  	file BqsrKnownSitesIdx from BqsrKnownSitesIdxChannel

      each GenomicInterval from BqsrGenomicIntervals
      file GATKExe
      val ApplyBQSRExtraOptionsString
      file JavaExe
      val JavaOptionsString

      file BashPreamble
      file BashSharedFunctions
      file BqsrScript

 output:
      set SampleName, "${SampleName}.${GenomicInterval}.bam" ,
           "${SampleName}.${GenomicInterval}.bai" into BqsrOutput

 script:
     """
     
     /bin/bash ${BqsrScript} -s ${SampleName} -b ${InputBams} -G ${Ref} \
       -k ${BqsrKnownSites.join(',')} -I ${GenomicInterval} -S ${GATKExe} \
       -o \"\'${ApplyBQSRExtraOptionsString}\'\" -J ${JavaExe} \
       -e \"\'${JavaOptionsString}\'\" -F ${BashSharedFunctions} ${DebugMode}
     """
}


/********      haplotyper: per interval (of a sample) - Required      *********/

process Haplotyper{
 tag "${SampleName}_${GenomicInterval}"

 input:
  	set SampleName,	file (InputBams), file (InputBais) from BqsrOutput // Link to BQSR

  	file Ref
	    file RefFai
      file RefDict

  	file DBSNP
	    file DBSNPIdx
	    val GenomicInterval from HCGenomicIntervals

      file GATKExe
      val HaplotyperThreads
      val HaplotyperExtraOptionsString
      file JavaExe
      val JavaOptionsString

      file BashPreamble
      file BashSharedFunctions
      file HaplotyperScript

	    val DebugMode

 output:
      set SampleName, "${SampleName}.${GenomicInterval}.g.vcf" ,
          "${SampleName}.${GenomicInterval}.g.vcf.idx" into HCOutput

 script:
     """
     
     /bin/bash ${HaplotyperScript} -s ${SampleName} -b ${InputBams} -G ${Ref} \
          -D ${DBSNP} -I ${GenomicInterval} -S ${GATKExe} \
          -t ${HaplotyperThreads} -o \"\'${HaplotyperExtraOptionsString}\'\" \
          -J ${JavaExe} -e \"\'${JavaOptionsString}\'\" -F ${BashSharedFunctions} \
           ${DebugMode}
     """
}


/********      Merge Gvcfs: all intervals (of a sample) - Required    *********/

process MergeGvcfs {
 tag "${SampleName}_All_intervals"

 publishDir DeliveryFolder_HaplotyperVC, mode: 'copy'

 input:
  	set SampleName, file(InputGvcfs), file(InputIdxs) from HCOutput.groupTuple()  //Link to Haplotyper

      file GATKExe
      file JavaExe
      val JavaOptionsString

      file BashPreamble
      file BashSharedFunctions
      file MergeGvcfsScript

	    val DebugMode

 output:
       file "${SampleName}.g.vcf" into MergeGvcfsOutput
       file "${SampleName}.g.vcf.idx" into MergeGvcfsIdx

 script:
     """
     
     /bin/bash ${MergeGvcfsScript} -s ${SampleName} -b ${InputGvcfs.join(',')}\
      -S ${GATKExe} -J ${JavaExe} -e \"\'${JavaOptionsString}\'\" \
      -F ${BashSharedFunctions} ${DebugMode}
     """
}
