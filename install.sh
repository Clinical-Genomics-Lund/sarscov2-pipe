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
