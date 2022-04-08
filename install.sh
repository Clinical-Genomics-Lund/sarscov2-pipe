#!/usr/bin/env bash
set -e

# Download actual pipeline if it doesn't exist
if [ ! -f "sarscov2.sh" ]; then
    wget https://raw.githubusercontent.com/Clinical-Genomics-Lund/sarscov2-pipe/main/sarscov2.sh
fi

# Download reference files
mkdir -p ref && cd ref
wget https://raw.githubusercontent.com/connor-lab/ncov2019-artic-nf/master/typing/SARS-CoV-2.types.yaml
wget https://raw.githubusercontent.com/connor-lab/ncov2019-artic-nf/master/typing/MN908947.3.gff
sed -i '/^$/d' MN908947.3.gff
grep -v '^###' MN908947.3.gff > MN908947.3.gff.mod
mv MN908947.3.gff.mod MN908947.3.gff
cd -
mkdir -p ref/v3 && cd ref/v3
wget https://raw.githubusercontent.com/artic-network/primer-schemes/master/nCoV-2019/V3/nCoV-2019.reference.fasta
wget https://raw.githubusercontent.com/artic-network/primer-schemes/master/nCoV-2019/V3/nCoV-2019.primer.bed
cd -
mkdir -p ref/v4.1 && cd ref/v4.1
wget https://raw.githubusercontent.com/artic-network/primer-schemes/master/nCoV-2019/V4.1/SARS-CoV-2.reference.fasta
wget https://raw.githubusercontent.com/artic-network/primer-schemes/master/nCoV-2019/V4.1/SARS-CoV-2.primer.bed
cd -

# Download some additional scripts
mkdir -p bin && cd bin
wget https://raw.githubusercontent.com/connor-lab/ncov2019-artic-nf/master/bin/type_vcf.py
wget https://raw.githubusercontent.com/connor-lab/ncov2019-artic-nf/master/bin/qc.py
# wget https://www.epicov.org/content/gisaid_uploader
cd -

# Create local mamba installation and activate it
mkdir -p mamba && cd mamba
wget -qO- https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj bin/micromamba
export MAMBA_ROOT_PREFIX=$PWD
eval "$(./bin/micromamba shell hook -s posix)"
cd -

# Install environmants with the neccessary tools
micromamba create -c bioconda -c conda-forge -c defaults -n pangolin pangolin -y

micromamba create -c bioconda -c conda-forge -c defaults -n ivar ivar bcftools bwa python sambamba freebayes seqtk pyvcf pandas matp\
lotlib PyYAML -y
micromamba activate ivar
pip install bio
micromamba deactivate

micromamba create -c bioconda -c conda-forge -c defaults -n vep ensembl-vep -y

micromamba create -c bioconda -n nextclade nextclade -y
micromamba activate nextclade
nextclade dataset get --name='sars-cov-2' --reference='MN908947' --output-dir='ref/nextclade/sars-cov-2_MN908947'
micromamba deactivate

tar -xf gisaid_cli2.tgz
cd gisaid_cli2
micromamba create -c conda-forge -c defaults -n cli2_env -f environment.yml -y
cd -
rm -r gisaid_cli2


# # Create local installatin of miniconda
# mkdir -p miniconda3
# wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda3/miniconda.sh
# bash miniconda3/miniconda.sh -b -u -p miniconda3
# rm -rf miniconda3/miniconda.sh
# ## Fix bug in WSL
# find miniconda3/ -type f -exec stat {} + > /dev/null

# # Activate the conda installation
# source miniconda3/etc/profile.d/conda.sh

# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

# conda update -q -y -n base -c defaults conda

# ## Setup iVar/bcftools conda environment
# conda create --name ivar ivar=1.3 bcftools=1.10.2 bwa=0.7.17 python=3.9 sambamba=0.8.0 freebayes=1.3.5 seqtk=1.3 pyvcf=0.6.8 -q -y
# conda activate ivar
# pip install bio pandas==1.2.1 matplotlib==3.3.4 PyYAML==5.4.1
# conda deactivate

# # Setup VEP conda environment
# conda create --name vep ensembl-vep=102.0 -q -y

# # Setup nextclade conda environment
# conda create --name nextclade nextclade=1.3.0 -q -y
# conda activate nextclade
# nextclade dataset get --name='sars-cov-2' --reference='MN908947' --output-dir='ref/nextclade/sars-cov-2_MN908947'
# conda deactivate

# # Install pangolin
# git clone https://github.com/cov-lineages/pangolin.git
# cd pangolin
# conda env create -f environment.yml -q
# conda activate pangolin
# pip install .
# cd -
