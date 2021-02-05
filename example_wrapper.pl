#!/usr/bin/perl -w
use strict;
use File::Basename qw(basename dirname);
use File::Spec;

my $SCRIPT_ROOT = File::Spec->rel2abs(dirname($0));
my $INDIR = $ARGV[0];
my $OUTDIR = $ARGV[1];
my $N_PARALLEL = $ARGV[2];

my @files = glob("$INDIR/*_R1_001.fastq.gz");

print "mkdir -p $OUTDIR\n";
print "cd $OUTDIR\n";
my $i = 0;
foreach my $fq1 (@files) {
    my $fq2 = $fq1;
    $fq2 =~ s/_R1_/_R2_/;
    my $id = (split /_/, basename($fq1))[0];
    $i++;
    print "bash $SCRIPT_ROOT/sarscov2.sh $fq1 $fq2 $id NO_PANGOLIN &\n";
    print "wait\n" if $i % $N_PARALLEL==0;
}
print "wait\n";

print "cat *.consensus.fa > all.consensus.fa\n";
print "source $SCRIPT_ROOT/miniconda3/etc/profile.d/conda.sh\n";
print "conda activate pangolin\n";
print "pangolin all.consensus.fa -o pangolin_all\n";
print "cp pangolin_all/lineage_report.csv ./pangolin_all.csv\n";

