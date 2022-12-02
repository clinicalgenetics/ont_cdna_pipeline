#!/usr/bin/env python
import snakemake as sm
import argparse
import os
import glob

parser = argparse.ArgumentParser(
                    prog = 'cDNA_pipeline',
                    description = 'ONT cDNA data analysis')

parser.add_argument('-i', '--input',
                    help = "path to fastq files, file names should end with fastq.gz",
                    dest = 'input', required = True)
parser.add_argument('-r', '--ref',
                    help = "reference genome, fasta file, should not be zipped.",
                    dest = 'ref', required = True)
parser.add_argument('-o', '--output',
                    help = "path to the output folder",
                    dest = 'output', required = True)
parser.add_argument('-g', '--gtf',
                    help = "reference genome annotation, gtf file, should not be zipped.",
                    dest = 'gtf', required = False)
parser.add_argument('-v', '--verbose', action='store_true', required = False)


args = parser.parse_args()



# running the snakemake pipeline
wdir = os.path.dirname(os.path.realpath(__file__))
samples = [os.path.basename(i).split(".fastq.gz")[0] for i in glob.glob(os.path.join(args.input, "*.fastq.gz"))]
config = {"samples":samples, "indir": args.input, "outdir": args.output, "ref": args.ref, "gtf": args.gtf}
sm.snakemake(os.path.join(wdir , "Snakefile"), config=config, cores=3, verbose=True, dryrun=False)
