#!/usr/bin/env python3

import argparse
import json
from io import StringIO
from pathlib import Path
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML


parser = argparse.ArgumentParser(description="Run BLAST and output a hit matrix")
parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file of top regions")
parser.add_argument("-d", "--blastdb", required=True, help="BLAST database path")
parser.add_argument("-o", "--outdir", required=True, help="Directory to save output files")
parser.add_argument("--blastn", default="blastn", help="Path to blastn binary")
args = parser.parse_args()

fasta_file = Path(args.fasta)
blast_db = Path(args.blastdb)
output_dir = Path(args.outdir)
output_list = output_dir / "seqList.txt"
output_json = output_dir / "blastHitResult.json"

output_dir.mkdir(parents=True, exist_ok=True)

# sequence IDs
master_seq_list = []
final_results = []

# Go through each sequence and BLAST it
for seq_record in SeqIO.parse(str(fasta_file), "fasta"):
    query_name = seq_record.id
    query_seq = seq_record.seq
    master_seq_list.append(query_name)

    query_fasta = f">{query_name}\n{query_seq}"

    # Run BLASTN
    blast_cline = NcbiblastnCommandline(
        cmd=args.blastn,
        db=str(blast_db),
        evalue=0.1,
        outfmt=5
    )
    stdout, stderr = blast_cline(stdin=query_fasta)

    # Parsing BLAST XML output
    blast_results = {}
    blast_xml = StringIO(stdout)
    for blast_record in NCBIXML.parse(blast_xml):
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                target_id = alignment.title.split(' ')[1]
                percent_identity = float(hsp.identities) / float(hsp.align_length) * 100
                normalized_score = percent_identity * float(hsp.align_length) / float(blast_record.query_length)
                blast_results[target_id] = normalized_score

    final_results.append({"Name": query_name, "Hits": blast_results})

# Write sequence ID list
with output_list.open("w") as f:
    for seq_id in master_seq_list:
        f.write(f"{seq_id}\n")

# Write JSON results
with output_json.open("w") as f:
    json.dump(final_results, f)
