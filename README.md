# Create Additional Probes on Unknown Mixed Samples

## Overview

**RiboZAP** is designed to generate custom RNaseH probes for metatranscriptomics samples without prior knowledge of sample compostion. Addition of these probes to RiboZero Plus Microbiome kit will significantly improve depletion of additional ribosomal RNA (rRNA) regions.

Depleting rRNA, which typically constitutes 80–90% of total RNA, before sequencing significantly reduces sequencing costs and increases the coverage of mRNA.


## Installation

**Prerequisites**

- Docker
- Python3

**Clone and install**:

```bash
git clone https://github.com/ilmn-mouse-cecal/RiboZAP.git
cd RiboZAP
pip install .
docker build -t ribozap:latest .
```

**Sample Sheet Format**

Your input must be a CSV with three columns:
```
sample_id,read1,read2
```

An example sample sheet is included [here](https://github.com/ilmn-mouse-cecal/RiboZAP/blob/537c72e1fe69ea05a341eff507c801d5ce370d3a/examples/sample_sheet.csv)

**Run the App**:

Probe design:
```bash
ribozap \
  --sample-sheet samplesheet.csv \
  --output-dir ./my_out/ \
  --analysis-name my_analysis \
  --cpus 4 \
  --memory 16 \
  --image ribozap \
  --image-tag 'latest' \
  --resume
```

Validation:
For validation, `--probes-fasta` and `--probes-summary` are required. These are produced by the probe design step (see outputs under your --output-dir).
Example:
```bash
ribozap \
  --sample-sheet samplesheet.csv \
  --output-dir ./my_out/ \
  --analysis-name my_analysis \
  --probes-fasta ./my_out/reports/probes.fasta \
  --probes-summary ./my_out/probes/probes.summary.csv \
  --cpus 4 \
  --memory 16 \
  --image ribozap \
  --image-tag 'latest' \
  --resume
```

## Code Style

The codebase follows [PEP 8](https://peps.python.org/pep-0008/) conventions and uses [Black](https://github.com/psf/black) for formatting.

## Citation

If you use `ribozap` in your research, please cite our publication:

Bunga, S., Tan, A., Roos, M., & Kuersten, S. (2025). RiboZAP: Custom probe design for rRNA depletion in complex metatranscriptomes. bioRxiv, 2025-11.