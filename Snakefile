import os
import glob
SAMPLES = glob.glob(os.path.join("/Users/leily/mmk_server/RR/RawDataBackup/RNASeqData/CLLE-ES/EGAF00001719*/*STAR*bam"))

rule all:
    input:
        "done"

# featurcount
rule count:
	input:
		bam=expand("{sample}", sample=SAMPLES),
	params:
		gtf="/Users/leily/dbgap/ref/gencode.v19.annotation.gtf"
	output:
		"counts0.tsv"
	shell:
		"/Users/leily/.conda/envs/deseq2/bin/featureCounts --extraAttributes gene_name,gene_type,gene_status -O -s 0 -p -t exon -g gene_id -T 14 -a {params.gtf} -o {output} {input}"


# remove extra coulmns before deseq2
rule modify_table:
    input:
        count_matrix = "counts0.tsv"
    params:
        metadata = "metadata.tsv"
    output:
        out = "count.processed.tsv"
    shell:
        "sed -e \"1d\" {input.count_matrix} | cut -f1,10- > {output.out} "
# deseq

rule make_matrix:
	input:
		"count.processed.tsv"
	output:
        path = os.path.join(".")
		flag = temp("done")
	shell:
        """
        Rscript -c "{input} -m {params.metadata} -o {output.path};
        touch {output.flag}
        """
