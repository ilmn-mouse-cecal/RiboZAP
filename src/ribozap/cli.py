import argparse
from pathlib import Path
from .rewrite import validate_and_rewrite
from .docker_launcher import run_docker
import logging
import re
import sys
import shutil

MIN_CPUS = 2
MIN_MEM = 10

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s"
)

def main():
    parser = argparse.ArgumentParser(
        prog="ribozap",
        description="Run sample prep and containerized workflow for probe design.",
    )

    # ---------- General options ----------
    opts = parser.add_argument_group('General options')
    opts.add_argument(
        "--analysis-name",
        type=str,
        required=True,
        help="Name of the analysis",
    )

    opts.add_argument(
        "--sample-sheet",
        type=Path,
        required=True,
        help="Path to the input sample sheet (CSV)",
    )

    opts.add_argument(
        "--output-dir",
        type=Path,
        help="Directory to write outputs and mount in Docker",
        required=True
    )

    opts.add_argument(
        "--cpus",
        default=4,
        type=int,
        help="Limit number of CPUs for Docker container (e.g. 4)",
        required=False
    )

    opts.add_argument(
        "--memory", default=16, type=int, help="Limit memory for Docker container (e.g. 16)",
        required=False
    )

    opts.add_argument(
        "--resume",
        action="store_true",
        help="Enable Nextflow -resume mode to continue from previous work directory",
        required=False
    )

    parser.add_argument(
        "--image",
        default="ribozap",
        type=str,
        help="Enter the image name you would like to use. Default is 'ribozap'",
        required=False
    )

    parser.add_argument(
        "--image-tag",
        default='latest',
        type=str,
        help="Enter the docker image tag. Default is 'latest'",
        required=False
    )

    parser.add_argument(
        "--dry-run", 
        action="store_true", 
        help="Show what would be done, but do not run anything",
        required=False
    )

    # ---------- Probe design options ----------
    probes = parser.add_argument_group('Probe design options')
    probes.add_argument(
        "--coverage-threshold",
        default=500,
        type=int,
        help="Minimum read depth required to define a high-coverage block. All consecutive regions with coverage ≥ this value are merged.",
        required=False
    )
    probes.add_argument(
        "--num-cov-regions",
        default=50,
        type=int,
        help="Enter the number of top high-coverage regions from the BED file to keep (default: 50)",
        required=False
    )

    probes.add_argument(
        "--probe-tiling-gap",
        default=25,
        type=int,
        help="Gap size between probes (default: 25)",
        required=False
    )

    # ---------- In-silico testing options ----------
    in_silico = parser.add_argument_group('In-silico testing options')
    in_silico.add_argument(
        "--padding",
        type=int,
        default=50,
        help="Number of bases to extend upstream and downstream of each alignment (default: 50bp). This is useful for in silico probe testing to get coverage beyond the aligned region."
    )

    validate = parser.add_argument_group(
        'Validating designed probes',
        description=(
            "This validation process will test the provided probes against a set of samples. "
            "If either --probes-fasta or --probes-summary is provided, the design workflow "
            "will not be executed."
        )
    )
    validate.add_argument(
        "--probes-fasta",
        type=Path,
        help="For testing designed probes against a set of samples, you can provide FASTA file of probe sequences."
    )
    validate.add_argument(
        "--probes-summary",
        type=Path,
        help="Probes summary sequence CSV file."
    )

    args = parser.parse_args()
    
    if args.cpus < MIN_CPUS:
        logging.error(f"Minimum required CPUs is {MIN_CPUS}. You provided {args.cpus}.")
        sys.exit(1)
    
    if args.memory < MIN_MEM:
        logging.error(f"Minimum required memory is {MIN_MEM} GB. You provided {args.memory} GB.")
        sys.exit(1)

    invalid = re.search(r"[^A-Za-z0-9_.-]", args.analysis_name)
    if invalid:
        logging.error("Invalid analysis name: only letters, numbers, dashes, underscores, and dots are allowed . . .")
        sys.exit(1)

    high_cpus = args.cpus
    high_memory = args.memory
    out_dir = args.output_dir
    analysis_name = args.analysis_name
    out_dir.mkdir(parents=True, exist_ok=True)
    num_coverage_regions = args.num_cov_regions
    padding = args.padding
    cov_threshold = args.coverage_threshold
    probe_tiling_gap = args.probe_tiling_gap
    validation_opt = ""
    other_paths = []
    if args.probes_fasta and args.probes_fasta.exists():
        #fasta_out = out_dir / "probes.fasta"
        #shutil.copy(args.probes_fasta, fasta_out)
        validation_opt += f' --probes_fasta "{args.probes_fasta.resolve()}" '
        other_paths.append(f'"{args.probes_fasta.resolve()}"')

    if args.probes_summary and args.probes_summary.exists():
        #summary_out = out_dir / "probes_summary.csv"
        #shutil.copy(args.probes_summary, summary_out)
        validation_opt += f' --probes_summary "{args.probes_summary.resolve()}" '
        other_paths.append(f'"{args.probes_summary.resolve()}"')

    rewritten_path = out_dir.resolve() / "rewritten_sample_sheet.csv"
    mount_path = out_dir.resolve() / "docker_mounts.txt"

    # Step 1: Rewrite sample sheet
    validate_and_rewrite(args.sample_sheet, rewritten_path, mount_path)

    resume_flag = "-resume" if args.resume else ""

    # Step 2: Run Docker
    res = run_docker(
        image=f"{args.image}:{args.image_tag}",
        mount_file=mount_path,
        out_dir=out_dir,
        analysis_name=analysis_name,
        container_cmd=f"nextflow run main.nf -work-dir {out_dir.resolve()}/{analysis_name}/work/ --analysis_name {analysis_name} --sample_sheet {rewritten_path} --outdir {out_dir.resolve()}/{analysis_name} --trace_dir {out_dir.resolve()}/{analysis_name}/trace_dir --top_coverage_regions {num_coverage_regions} --cpus {high_cpus} --memory '{high_memory} GB' --probe_tiling_gap {probe_tiling_gap} --padding {padding} --coverage_threshold {cov_threshold} {validation_opt} {resume_flag}",
        other_paths = other_paths,
        cpus=args.cpus,
        memory=args.memory,
        dry_run=args.dry_run
    )
    if args.dry_run:
        logging.info(f"""
----------------------------
[DRY RUN]
{res}
----------------------------
        """)


if __name__ == "__main__":
    main()
