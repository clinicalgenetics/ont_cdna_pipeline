import os
import glob
import snakemake

rule all:
    input:
        os.path.join(config['outdir'], "pcaplot.png"),
        os.path.join(config['outdir'], "flair_correct", "concat_all_corrected.bed")
# Optional trimming


# minimap2
rule mapping:
    input:
        fastq = os.path.join(config['indir'], "{sample}.fastq.gz")
    output:
        bam = os.path.join(config['outdir'], "bam_files", "{sample}.sorted.bam"),
    params:
        ref = config['ref'],
        gtf = config['gtf']
    shell: """
        minimap2 -ax splice -t 8 -K 1G --junc-bed {params.gtf} {params.ref} {input.fastq} |
        samtools sort -@ 4 -T /tmp/ -o {output.bam} ; samtools index {output.bam}
    """


# deeptools_qc
## coverage
rule bamCoverage:
    input:
        bam = os.path.join(config['outdir'], "bam_files", "{sample}.sorted.bam")
    output:
        bw = os.path.join(config['outdir'], "bamCoverage", "{sample}.bw")
    shell:
        "bamCoverage -b {input.bam} -o {output.bw} -p 8 --normalizeUsing RPKM "

## bwsummary
rule multiBigWigSummary:
    input:
        bws = expand(os.path.join(config['outdir'], "bamCoverage", "{sample}.bw"), sample = config['samples'])
    output:
        matrix = os.path.join(config['outdir'], "multiBigWigSummary.npz")
    shell:
        """
        multiBigwigSummary bins -b {input.bws} -o {output.matrix} -p 8
        """

## plotPCA
rule plotPCA:
    input:
        matrix = os.path.join(config['outdir'], "multiBigWigSummary.npz")
    output:
        plot = os.path.join(config['outdir'], "pcaplot.png")
    shell:
        """
        plotPCA -in {input.matrix} -o {output.plot}
        """
# # Flair
## bam2bed12
rule bam2bed12:
    input:
        bam = os.path.join(config['outdir'], "bam_files", "{sample}.sorted.bam")
    output:
        bed = os.path.join(config['outdir'], "bed12_files", "{sample}.bed")
    shell:
        """
        bam2Bed12 -i {input.bam} > {output.bed}
        """

## flair correct
rule flair_correct:
    input:
        bed = os.path.join(config['outdir'], "bed12_files", "{sample}.bed")
    output:
        correct = os.path.join(config['outdir'], "flair_correct", "{sample}")
    params:
        ref = config['ref'],
        gtf = config['gtf']
    shell:
        """
        flair correct -q  {input.bed} -t 8 -f {params.gtf} -g {params.ref} -o {output.correct}
        """

rule correct_concat:
    input:
        all_correct = expand(os.path.join(config['outdir'], "flair_correct", "{sample}"), sample = config['samples'])
    output:
        concat = os.path.join(config['outdir'], "flair_correct", "concat_all_corrected.bed")
    shell:
        """
        cat {input.all_correct} > {output.concat}
        """


## flair collapse
rule flair_collapse:
    input:
        fastqs = expand(os.path.join(config['indir'], "{sample}.fastq.gz"), sample = config['samples']),
        concat = os.path.join(config['outdir'], "flair_correct", "concat_all_corrected.bed")
    output:
        correct = os.path.join(config['outdir'], "flair_collapse", "concat")
    params:
        ref = config['ref'],
        gtf = config['gtf']
    shell:
        """
        flair collapse -g {params.ref} -f {params.gtf} -q {input.concat} \
        -r {input.fastqs} -o
        """
