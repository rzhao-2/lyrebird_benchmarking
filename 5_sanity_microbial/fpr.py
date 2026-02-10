import os
import logging
from os.path import join
from csv import DictReader
import gzip

viral_reps = set()

viral_rep_paths = snakemake.input.viral_rep_paths
with open(viral_rep_paths) as r:
    for line in r:
        viral_reps.add(line.split('\t')[0])

otu_tables = snakemake.input.otu_tables
genomewise_tables = snakemake.input.genomewise_tables
forward_reads = snakemake.input.r1
reverse_reads = snakemake.input.r2
output_report = snakemake.output.fpr_report
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')
logging.info(f"Processing {len(otu_tables)} OTU tables")

with open(output_report, 'w') as out:
    out.write("OTU_table\tsize_viral\tsize_microbial\tviral_readcount\tmicrobial_readcount\ttrue_viral\ttrue_microbial\tfalse_viral\tfalse_microbial\troot_viral\troot_microbial\n")
    for otu_table, genomewise, r1, r2 in zip(otu_tables, genomewise_tables, forward_reads, reverse_reads):
        filename = os.path.basename(otu_table)
        logging.info(f"Processing {filename} with genomewise table {os.path.basename(genomewise)}")

        true_viral = 0
        true_micro = 0
        false_viral = 0
        false_micro = 0
        root_viral = 0
        root_micro = 0
        
        # Read genomewise table to get the size of viral and microbial genomes
        with open(genomewise) as g:
            genome_sizes = DictReader(g, delimiter='\t')
            for row in genome_sizes:
                if 'd__Viruses' in row['taxonomy']:
                    viral_size = int(row['genome_size'])
                else:
                    microbial_size = int(row['genome_size'])
        
        # Count the number of viral/microbial reads in forward and reverse reads - in fastq.gz format
        viral_reads = 0
        microbial_reads = 0
        with gzip.open(r1) as f1:
            for line in f1:
                if line.startswith(b'@'):
                    header = line.decode('utf-8').strip()
                    taxid = header.split('-')[0][1:]
                    if taxid in viral_reps:
                        viral_reads += 1
                    else:
                        microbial_reads += 1
        with gzip.open(r2) as f2:
            for line in f2:
                if line.startswith(b'@'):
                    header = line.decode('utf-8').strip()
                    taxid = header.split('-')[0][1:]
                    if taxid in viral_reps:
                        viral_reads += 1
                    else:
                        microbial_reads += 1

        with open(otu_table) as r:
            for line in r.readlines()[1:]:
                taxon = line.split('\t')[5]
                readnames = line.split('\t')[6]
                for readname in readnames.split():
                    taxid = readname.split('-')[0]
                    if 'd__Viruses' in taxon:
                        if taxid in viral_reps:
                            true_viral += 1
                        else:
                            false_micro += 1
                    elif 'd__Bacteria' in taxon or 'd__Archaea' in taxon:
                        if taxid in viral_reps:
                            false_viral += 1
                        else:
                            true_micro += 1
                    else:
                        if taxid in viral_reps:
                            root_viral += 1
                        else:
                            root_micro += 1
        out.write(f"{filename}\t{viral_size}\t{microbial_size}\t{viral_reads}\t{microbial_reads}\t{true_viral}\t{true_micro}\t{false_viral}\t{false_micro}\t{root_viral}\t{root_micro}\n")

with open(snakemake.output.done, 'w+') as __:
    pass