import os
from tqdm import tqdm
from Bio import SeqIO
from csv import DictReader

provirus_summary = "/home/n10927662/rossenzhao/genomad_full_gtdb/gtdb_r214_provirus_summary.tsv"
microbial_rep_paths = "/home/n10927662/rossenzhao/lyrebird-benchmarking/5_sanity_microbial/microbial_rep_paths.csv"

proviruses = set()
with open(provirus_summary) as f:
    reader = DictReader(f, delimiter="\t")
    for row in tqdm(reader):
        if row["topology"] == "Provirus":
            proviruses.add(row["seq_name"].split("|")[0])
        else:
            proviruses.add(row["seq_name"])

microbes_without_proviruses = []

with open(microbial_rep_paths) as r:
    for line in tqdm(r.readlines()):
        id, filepath = line.strip().split("\t")
        for record in SeqIO.parse(filepath, "fasta"):
            if record.id in proviruses:
                break
        else:
            microbes_without_proviruses.append(line)

with open("/home/n10927662/rossenzhao/lyrebird-benchmarking/5_sanity_microbial/microbial_rep_paths_no_proviruses.csv", "w+") as w:
    for line in microbes_without_proviruses:
        w.write(line)