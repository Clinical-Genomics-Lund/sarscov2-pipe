#!/usr/bin/env bash
set -e

# Download actual pipeline if it doesn't exist
if [ ! -f "sarscov2.sh" ]; then
    wget https://raw.githubusercontent.com/Clinical-Genomics-Lund/sarscov2-pipe/main/sarscov2.sh
fi

# Download reference files
mkdir -p ref && cd ref
wget https://github.com/connor-lab/ncov2019-artic-nf/raw/master/typing/SARS-CoV-2.types.yaml
wget https://github.com/artic-network/primer-schemes/raw/master/nCoV-2019/V3/nCoV-2019.reference.fasta
wget https://github.com/artic-network/primer-schemes/raw/master/nCoV-2019/V3/nCoV-2019.primer.bed
wget https://github.com/connor-lab/ncov2019-artic-nf/raw/master/typing/MN908947.3.gff
cd -

# Download some additional scripts
mkdir -p bin && cd bin
wget https://github.com/connor-lab/ncov2019-artic-nf/raw/master/bin/type_vcf.py
wget https://github.com/connor-lab/ncov2019-artic-nf/raw/master/bin/qc.py
cd -


## Create local installatin of miniconda
mkdir -p miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3/miniconda.sh
bash miniconda3/miniconda.sh -b -u -p miniconda3
rm -rf miniconda3/miniconda.sh

# Activate the conda installation
source miniconda3/etc/profile.d/conda.sh

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

## Setup iVar/bcftools conda environment
conda create --name ivar ivar=1.3 bcftools=1.10.2 bwa=0.7.17 python=3.9 -q -y
conda activate ivar
pip install bio==0.3.0 pandas==1.2.1 matplotlib==3.3.4 PyVCF==0.6.8 PyYAML==5.4.1

ivar version
bcftools --version
conda deactivate


# Install pangolin
git clone https://github.com/cov-lineages/pangolin.git
cd pangolin
conda env create -f environment.yml -q
conda activate pangolin
python setup.py install
cd -
