# ribozap

**Create Additional Probes on Unknown Mixed Samples**

## Overview

**ribozap** is designed to generate custom RNaseH probes for metatranscriptomics samples without prior knowledge of sample compostion. Addition of these probes to RiboZero Plus Microbiome kit will significantly improve depletion of additional ribosomal RNA (rRNA) regions.

Depleting rRNA, which typically constitutes 80–90% of total RNA, before sequencing significantly reduces sequencing costs and increases the coverage of mRNA.


## Installation

**Prerequisites**

- Docker
- Python3

**Clone the repository**:

```bash
git clone https://github.com/yourusername/ribozap.git
cd RiboZAP
pip install .
docker build -t ribozap:latest .

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


## Code Style

The codebase follows [PEP 8](https://peps.python.org/pep-0008/) conventions and uses [Black](https://github.com/psf/black) for formatting.

## Citation

If you use `ribozap` in your research, please cite our publication:

Roos M, Bunga S, Tan A, Maissy E, Skola D, Richter A, Whittaker DS, Desplats P, Zarrinpar A, Conrad R, Kuersten S.0.Optimizing mouse metatranscriptome profiling by selective removal of redundant nucleic acid sequences. mSystems0:e00167-25.https://doi.org/10.1128/msystems.00167-25