#/usr/bin/env bash
set -e

DIR=`dirname $0`
# source $DIR/miniconda3/etc/profile.d/conda.sh
export MAMBA_ROOT_PREFIX=$DIR/mamba
eval "$(${MAMBA_ROOT_PREFIX}/bin/micromamba shell hook -s posix)"

PRIMER_BED=${5:-"${DIR}/ref/v4.1/SARS-CoV-2.primer.bed"}
REF_FASTA="${DIR}/ref/v4.1/SARS-CoV-2.reference.fasta"

REF_NEXTCLADE="${DIR}/ref/nextclade/sars-cov-2_MN908947"
GFF="${DIR}/ref/MN908947.3.gff"
QC_PY="${DIR}/bin/qc.py"

FQ1=$1
FQ2=$2
ID=$3
NO_PANGOLIN=$4

MIN_DEPTH=10
MIN_FREQ=0.75
MIN_QUAL_VAR=20
MAX_READPAIRS=1000000

micromamba activate ivar

# bgzip gff in needed
if [ ! -f "${GFF}.gz" ]; then
    bgzip -i -f -c $GFF > $GFF.gz
    tabix -f $GFF.gz
fi

# Build BWA index if necessary
if [ ! -f "${REF_FASTA}.bwt" ]; then
    bwa index $REF_FASTA
fi

# Subsample fastqs to at most $MAX_READPAIRS
if [ ! -f "${ID}_subsample_R1_001.fastq.gz" ] && [ ! -f "$ID.trim.sort.bam" ]; then
    seqtk sample -s 1314 $FQ1 $MAX_READPAIRS | gzip -c > ${ID}_subsample_R1_001.fastq.gz &
    seqtk sample -s 1314 $FQ2 $MAX_READPAIRS | gzip -c > ${ID}_subsample_R2_001.fastq.gz &
    wait
fi

# Align and trim primers
if [ ! -f "$ID.trim.sort.bam" ]; then
    bwa mem -t 4 $REF_FASTA ${ID}_subsample_R1_001.fastq.gz ${ID}_subsample_R2_001.fastq.gz | samtools sort | tee $ID.qc.bam | samtools view -F 4 -o $ID.sort.bam
    samtools flagstat $ID.qc.bam > $ID.flagstat
    ivar trim -e -i $ID.sort.bam -b $PRIMER_BED -p $ID.trim -m 30 -q 20
    samtools sort $ID.trim.bam -o $ID.trim.sort.bam
    samtools index $ID.trim.sort.bam
    rm $ID.sort.bam $ID.sort.bam.bai $ID.trim.bam $ID.qc.bam
fi

# Create fastq files for distribution
if [ ! -f "${ID}_R1_001.fastq.gz" ]; then
    samtools fastq -1 ${ID}_R1_001.fastq.gz -2 ${ID}_R2_001.fastq.gz $ID.trim.sort.bam &
fi

# Create detailed depth data
if [ ! -f "$ID.depth" ]; then
    sambamba depth base -c0 $ID.trim.sort.bam -o $ID.depth &
fi

# Create consensus sequence
if [ ! -f "$ID.consensus.fa" ]; then
    samtools mpileup -aa -A -B -d 6000000 -Q 0 --reference $REF_FASTA $ID.trim.sort.bam | ivar consensus -p $ID.consensus -n N -m $MIN_DEPTH -t $MIN_FREQ &
fi

# Call variants with freebayes
if [ ! -f "$ID.freebayes.vcf" ]; then
    freebayes -p 1 --min-coverage $MIN_DEPTH --min-base-quality $MIN_QUAL_VAR -f $REF_FASTA -F $MIN_FREQ -m 60 $ID.trim.sort.bam > $ID.freebayes.raw.vcf
    bcftools norm $ID.freebayes.raw.vcf -f $REF_FASTA -o $ID.freebayes.vcf
    rm $ID.freebayes.raw.vcf
fi

# Annotate variants with VEP
if [ ! -f "$ID.freebayes.vep.vcf" ]; then
    bgzip -i -f -c $ID.freebayes.vcf > $ID.freebayes.vcf.gz
    micromamba activate vep
    vep -i $ID.freebayes.vcf.gz --format vcf --gff $GFF.gz --fasta $REF_FASTA -o $ID.freebayes.vep.vcf --vcf --force_overwrite --no_stats --distance 10 --hgvs
    micromamba activate ivar
    rm $ID.freebayes.vcf.gz*
fi

# Generate QC data
if [ ! -f "$ID.qc.csv" ]; then
    wait
    env -u DISPLAY python $QC_PY --illumina --outfile $ID.qc.csv --sample $ID --ref $REF_FASTA --bam $ID.trim.sort.bam --fasta $ID.consensus.fa
fi

# Run nextclade
if [ ! -f "$ID.nextclade.tsv" ] && [ ! -f "$ID.auspice.json" ]; then
    wait
    micromamba activate nextclade
    nextclade --in-order --input-fasta $ID.consensus.fa --input-dataset $REF_NEXTCLADE --output-tsv $ID.nextclade.tsv --output-tree $ID.auspice.json
fi

# Run pangolin
if [ ! "$NO_PANGOLIN" = "NO_PANGOLIN" ] && [ ! -f "$ID.pangolin.csv" ]; then
    wait
    micromamba activate pangolin
    pangolin $ID.consensus.fa -o $ID.pangolin_tmp
    cp $ID.pangolin_tmp/lineage_report.csv ./$ID.pangolin.csv
    rm -rf $ID.pangolin_tmp
fi

