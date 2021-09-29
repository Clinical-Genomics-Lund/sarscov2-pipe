#!/usr/bin/perl -w
use strict;
use File::Basename qw(basename dirname);
use File::Spec;
use Cwd qw(abs_path);

my $SCRIPT_ROOT = File::Spec->rel2abs(dirname($0));
my $INDIR = $ARGV[0];
my $OUTDIR = $ARGV[1];
my $N_PARALLEL = $ARGV[2];
my $POST_ANALYSIS_SCRIPT = ($ARGV[3] or 0);

my @files = glob("$INDIR/*_R1_001.fastq.gz");

print "source $SCRIPT_ROOT/miniconda3/etc/profile.d/conda.sh\n";
print "conda activate pangolin\n";
print "pangolin --update\n";


print "mkdir -p $OUTDIR\n";
print "cd $OUTDIR\n";
my $i = 0;
foreach my $fq1 (@files) {
    next if basename($fq1) =~ /^Undetermined_/;
    my $fq2 = $fq1;
    $fq2 =~ s/_R1_/_R2_/;
    my $id = (split /_/, basename($fq1))[0];
    $i++;
    print "bash $SCRIPT_ROOT/sarscov2.sh $fq1 $fq2 $id NO_PANGOLIN &\n";
    print "wait\n" if $i % $N_PARALLEL==0;
}
print "wait\n";

unless( -e $OUTDIR.'/pangolin_all.csv' ) {
    print "cat *.consensus.fa > all.consensus.fa\n";
    print "conda activate pangolin\n";
    print "pangolin all.consensus.fa -o pangolin_all\n";
    print "cp pangolin_all/lineage_report.csv ./pangolin_all.csv\n";
}

if( $POST_ANALYSIS_SCRIPT ) {
    print $ARGV[3]." ".abs_path($OUTDIR)."\n";
}
