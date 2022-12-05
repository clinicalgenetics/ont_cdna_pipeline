import os
import glob
import snakemake

rule all:
    input:
        os.path.join(config['outdir'], "pcaplot.png"),
        os.path.join(config['outdir'], "done")
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
# Flair
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
        all_correct = os.path.join(config['outdir'], "flair_correct", "{sample}_all_corrected.bed")

    params:
        ref = config['ref'],
        gtf = config['gtf'],
        file_prefix = os.path.join(config['outdir'], "flair_correct", "{sample}")
    shell:
        """
        flair correct -q  {input.bed} -t 8 -f {params.gtf} -g {params.ref} -o {params.file_prefix};
        head -n1 {output.all_correct}
        """

rule correct_concat:
    input:
        all_correct = expand(os.path.join(config['outdir'], "flair_correct", "{sample}_all_corrected.bed"), sample = config['samples'])
    output:
        concat = os.path.join(config['outdir'], "flair_correct", "concat_all_corrected.bed")
    shell:
        """
        cat {input.all_correct} > {output.concat}
        """


# flair collapse
rule flair_collapse:
    input:
        fastqs = expand(os.path.join(config['indir'], "{sample}.fastq.gz"), sample = config['samples']),
        concat = os.path.join(config['outdir'], "flair_correct", "concat_all_corrected.bed")
    output:
        isoform_fa = os.path.join(config['outdir'], "flair_collapse", "concat.isoforms.fa"),
        soform_bed = os.path.join(config['outdir'], "flair_collapse", "concat.isoforms.bed")
    params:
        ref = config['ref'],
        gtf = config['gtf'],
        file_prefix = os.path.join(config['outdir'], "flair_collapse", "concat")
    shell:
        """
        flair collapse -g {params.ref} -f {params.gtf} -q {input.concat} -r {input.fastqs} -o {params.file_prefix};
        head -n1 {output.isoform_fa}
        """


# flair quantify
rule flair_quantify:
    input:
        isoform = os.path.join(config['outdir'], "flair_collapse", "concat.isoforms.fa")
    output:
        count_matrix = os.path.join(config['outdir'], "flair_quantify", "flair.quantify.counts.tsv")
    params:
        metadata= config["metadata"],
        file_prefix = os.path.join(config['outdir'], "flair_quantify", "flair.quantify")

    shell:
        """
        flair quantify -r {params.metadata} -i {input.isoform} --tpm -t 6 -o {params.file_prefix};
        """


# flair diffSplice (by default only the first 2 cond in the smaplesheets are taken into account)
rule flair_diffsplice:
    input:
        count_matrix = os.path.join(config['outdir'], "flair_quantify", "flair.quantify.counts.tsv"),
        isoform = os.path.join(config['outdir'], "flair_collapse", "concat.isoforms.bed")
    output:
        path = directory(os.path.join(config['outdir'], "flair_diffsplice")),
        done = os.path.join(config['outdir'], "done")
    shell:
        """
        flair diffSplice -i {input.isoform} -q {input.count_matrix}  --test --threads 4 --batch -o {output.path};
        touch {output.done}
        """
