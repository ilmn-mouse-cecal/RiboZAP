#!/bin/bash
set -euo pipefail

# Create resource directories
mkdir -p resources

echo "Downloading rRNA FASTA database..."
curl -L -o resources/rRNA_databases.tar.gz \
  https://github.com/ilmn-mouse-cecal/RiboZAP/releases/download/v0.0.1/rRNA_databases.tar.gz

echo "Extracting databases..."
tar -xzf resources/rRNA_databases.tar.gz -C resources/

echo "Downloading Genomes..."
curl -L -o resources/genomes.tar.gz \
  https://github.com/ilmn-mouse-cecal/RiboZAP/releases/download/v0.0.1/genome.tar.gz

echo "Extracting Genomic Databases . . ."
tar -xzf resources/genomes.tar.gz -C resources/

rm resources/*.tar.gz

echo "✅ Done. Databases ready in resources/"
