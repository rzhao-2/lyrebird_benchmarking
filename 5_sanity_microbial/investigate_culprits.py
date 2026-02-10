import os
import logging
from os.path import join
import csv
from csv import DictReader
import gzip
import numpy as np

csv.field_size_limit(sys.maxsize)

viral_reps = set()

with open("/home/n10927662/rossenzhao/lyrebird-benchmarking/5_sanity_microbial/viral_rep_paths.tsv") as r:
    for line in r:
        viral_reps.add(line.split('\t')[0])

otu_tables = snakemake.input.otu_tables
output_report = snakemake.output.culprits

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s', datefmt='%Y/%m/%d %I:%M:%S %p')
logging.info(f"Processing {len(otu_tables)} OTU tables")

with open(output_report, 'w+') as out:
    out.write("problem_phrog\tcount_otu_tables\tmedian_microbial_readcount\ttotal_microbial_readcount\n")

    bad_phrogs = {}
    for otu_table in otu_tables:
        filename = os.path.basename(otu_table)
        logging.info(f"Processing {filename}")
        phrog_check = {}

        with open(otu_table) as f:
            otu_data = DictReader(f, delimiter='\t')
            for row in otu_data:
                phrog = row['gene']
                if 'd__Viruses' in row['taxonomy']:
                    readnames = row['read_names'].split()
                    for readname in readnames:
                        taxid = readname.split('-')[0][1:]
                        if taxid not in viral_reps:
                            if phrog not in phrog_check:
                                phrog_check[phrog] = 0
                            phrog_check[phrog] += 1
        for phrog in phrog_check:
            if phrog not in bad_phrogs:
                bad_phrogs[phrog] = []
            bad_phrogs[phrog].append(phrog_check[phrog])
    for phrog, counts in bad_phrogs.items():
        median_count = np.median(counts)
        total_count = sum(counts)
        out.write(f"{phrog}\t{len(counts)}\t{median_count}\t{total_count}\n")
        logging.info(f"Phrog {phrog} has {len(counts)} counts, median: {median_count}, total: {total_count}")