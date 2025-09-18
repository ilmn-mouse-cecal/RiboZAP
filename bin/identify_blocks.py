#!/usr/bin/env python3
import argparse
from pathlib import Path

"""
This script processes a genome-wide coverage file in BED format, to identify continuous high-coverage regions.

Purpose:
To extract blocks of contiguous bases with read coverage above a user-defined threshold (--high)
"""

parser = argparse.ArgumentParser(description="Extract high-coverage blocks from a genomeCoverage.bed file.")
parser.add_argument("-s", "--sample", required=True, help="Sample name (used for file naming)")
parser.add_argument("-c", "--cov", required=True, help="Genome coverage file")
parser.add_argument("--high", type=int, required=True, help="Minimum coverage for high coverage block")

args = parser.parse_args()

input_path = Path(args.cov)
output_path = Path(f"{args.sample}_high_coverage_blocks.tsv")

prev_genome = ""
block_flag = False
block_start = block_end = max_cov = -1

# Write output
with output_path.open("w") as fout:
    fout.write("Genome Name\tBlock Start\tBlock End\tBlock Max Coverage\n")
    print(f"Identifying Blocks for {args.sample}.\nHigh Cov: {args.high}+")

    with input_path.open("r") as fin:
        for line in fin:
            genome, start, end, cov = line.strip().split('\t')
            start, end, cov = int(start), int(end), int(cov)

            if genome == prev_genome:
                if block_flag:
                    if cov >= args.high:
                        block_end = end
                        max_cov = max(max_cov, cov)
                    else:
                        fout.write(f"{prev_genome}\t{block_start}\t{block_end}\t{max_cov}\n")
                        block_flag = False
                else:
                    if cov >= args.high:
                        block_start, block_end = start, end
                        max_cov = cov
                        block_flag = True
            else:
                if block_flag:
                    fout.write(f"{prev_genome}\t{block_start}\t{block_end}\t{max_cov}\n")
                block_flag = False
                if cov >= args.high:
                    block_start, block_end = start, end
                    max_cov = cov
                    block_flag = True

            prev_genome = genome

    if block_flag:
        fout.write(f"{prev_genome}\t{block_start}\t{block_end}\t{max_cov}\n")
