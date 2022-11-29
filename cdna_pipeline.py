#!/usr/bin/env python
import snakemake.snakemake as sm
import argparse


parser = argparse.ArgumentParser(
                    prog = 'cDNA_pipeline',
                    description = 'ONT cDNA data analysis')

parser.add_argument('-i', '--input', nargs='+',
                    help = "path to fastq files",
                    dest = 'input', required = True)
parser.add_argument('-r', '--ref',
                    help = "reference genome, fasta file",
                    dest = 'ref', required = True)
parser.add_argument('-g', '--gtf',
                    help = "reference genome annotation, gtf file",
                    dest = 'gtf', required = False)
parser.add_argument('-v', '--verbose', action='store_true', required = False)


args = parser.parse_args()



# running the snakemake pipeline
sm("path to snakefile")
