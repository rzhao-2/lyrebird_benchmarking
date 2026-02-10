import os
import argparse
import logging

fastq_1 = snakemake.input.r1
fastq_2 = snakemake.input.r2

# fastq_1s = snakemake.input.r1
# fastq_2s = snakemake.input.r2

phanta_config_example = snakemake.params.phanta_config_example
outdir = snakemake.params.outdir
threads = snakemake.params.threads

sample_file = snakemake.output.sample_file
phanta_config = snakemake.output.phanta_config

logging.basicConfig(filename=snakemake.log[0], level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')
logging.info(f"Writing fastq file paths to sample file {sample_file}.")
with open(sample_file, 'w') as f:
    # for fastq_1, fastq_2 in zip(fastq_1s, fastq_2s):
    f.write(f"{fastq_1.split('/')[-1].split('.')[0]}\t{fastq_1}\t{fastq_2}\n")

logging.info(f"Writing phanta config to {phanta_config}.")
with open(phanta_config_example, 'r') as r:
    with open(phanta_config, 'w+') as w:
        for line in r:
            if "sample_file" in line:
                w.write(f"sample_file: {sample_file}\n")
            elif "outdir" in line:
                w.write(f"outdir: {outdir}\n")
            elif "class_threads" in line:
                w.write(f"class_threads: {threads}\n")
            else:
                w.write(line)
logging.info("Done.")

# with open(snakemake.output.done, 'w') as __: pass