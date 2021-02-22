# A boring Sars-Cov-2 pipeline

A pipeline for analysing Illumina Sars-Cov-2 WGS data (ARTIC v3) built using the priciples of Choose Boring Technology and KISS. Mostly based on [iVar](https://github.com/andersen-lab/ivar), [Pangolin](https://github.com/cov-lineages/pangolin) and [freebayes](https://github.com/freebayes/freebayes) with some scripts and files stolen from https://github.com/connor-lab/ncov2019-artic-nf/

Does the following:
* Optionally subsamples the data
* Generates consensus sequence
* Calls majority variants with freebayes
* Annotates variants with VEP
* Determines pangolin type
* Generates QC data
* Generates fastqs without host contamination and primers, for distribution

## Requirements

Any fairly modern and normal Linux distribution should be sufficient. Bash is the only requirements I can think of. Wget and an internet connection is needed for installation.


## Installation

Either clone the repository and run the install.sh script:

```
$ git clone https://github.com/Clinical-Genomics-Lund/sarscov2-pipe
$ cd sarscov2-pipe
$ bash install.sh
```

or do :

```
$ curl https://raw.githubusercontent.com/Clinical-Genomics-Lund/sarscov2-pipe/main/install.sh | bash
```
Both options will install a local copy of miniconda3, setup needed software in conda environements and download required reference data.

## Running

```
$ bash /path/to/sarscov2.sh R1.fastq.gz R2.fastq.gz SAMPLE_ID
```

This will create a number of output files in the current working directory

Typically you would run the script from some sort of wrapper. The repository includes a simple Perl example script to start a whole flow cell. It will use everthing up to the first underscore in the fastq filenames as the sample ID. It will generate a bash file which you need to run somehow.

```
$ /path/to/example_wrapper.pl  FASTQ_DIR  OUTPUT_DIR  #SAMPLES_TO_RUN_IN_PARALLEL | bash
```
