#!/usr/bin/env python3

import argparse
from pathlib import Path

# Argument parsing
parser = argparse.ArgumentParser(
    description="Generate probe tiling windows from a FASTA index (.fai) using user-defined probe length and gap size."
)
parser.add_argument("fasta_fai", help="Path to the FASTA index (.fai) file")
parser.add_argument("output_bed", help="Path to output BED file")
parser.add_argument("--gap", "-g", type=int, default=25, help="Gap size between probes (default: 25)")
parser.add_argument("--probe-length", "-p", type=int, default=50, help="Length of each probe. Default is 50. Acceptable range: 25–100 nt")
args = parser.parse_args()

fai_path = Path(args.fasta_fai)
output_path = Path(args.output_bed)
gap = args.gap
probe_length = args.probe_length
with output_path.open("w") as fout, fai_path.open("r") as fin:
    for line in fin:
        chrom, length = line.strip().split('\t')[:2]
        length = int(length)
        probe_num = 1
        i = length
        while i > probe_length:
            start = i - probe_length
            end = i
            fout.write(f"{chrom}\t{start}\t{end}\t{chrom}_Probe{probe_num}\t1\t-\n")
            i -= (probe_length + gap)
            probe_num += 1
