import os
import glob
SAMPLES = glob.glob(os.path.join("/Users/leily/mmk_server/RR/RawDataBackup/RNASeqData/CLLE-ES/EGAF00001719*/*STAR*bam"))

rule all:
    input:
        "done"
# Optional trimming

# minimap2
rule mapping:
	input:
		fastq = "{sample}.fastq"
	params:
        ref =
		gtf =
	output:
		bam = "{sample}.bam"
	shell:
		"""
          minimap2 -ax splice -t 8 -K 1G --junc-bed {params.gtf} {params.ref} {input.fastq} |
          samtools sort -@ 4 -T /tmp/ -o {output.bam} ; samtools index {output.bam}
        """

# deeptools_qc
## coverage
rule bamCoverage:
    input:
        bams = "{sample}.bam"
    params:
        metadata = "metadata.tsv"
    output:
        out = "{sample}.bw"
    shell:
        "bamCoverage -b {input.bam} -o {output.bw} -p 8 --normalizeUsing RPKM "

## bwsummary
rule multiBigWigSummary:
    input:
        bws = expand("{sample}.bw", sample=SAMPLES)
    output:
        matrix = "multiBigWigSummary.npz"
    shell:
        """
        multiBigwigSummary bins -b {input.bws} -o {output.matrix} -p 8
        """

## plotPCA
rule plotPCA:
    input:
    output:
    shell:
# Flair
## bam2bed12
rule bam2bed12:
    
