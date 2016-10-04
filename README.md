#  PharmGKB NGS Pipeline

This repository contains a simplified version of our work-in-progress grc38 NGS pipeline for pharmacogenomics.  The primary goal of this pipeline is to generate VCF data suitable for use by [PharmCat](https://github.com/PharmGKB/PharmCAT).


## Quick Overview

The pipeline takes FASTQ files, aligns them using BWA, and then calls variants using GATK.  

It is optimized for WGS and Exome analysis, and currently does not handle trios.  [VQSR](http://gatkforums.broadinstitute.org/gatk/discussion/39/variant-quality-score-recalibration-vqsr) filtering has shown to be effective for WGS data, but you may want to switch to [hard filtering](http://gatkforums.broadinstitute.org/wdl/discussion/2806/howto-apply-hard-filters-to-a-call-set) if using other data sources.

Take a look at `run.sh` to see the steps of our pipeline and the parameters used.  It takes a pair of fastq files from an Illumina NGS run and outputs a complete VCF file containing genotypes for *every* position, useful for pharmacogenomics analysis and a smaller VCF containing only the positions needed for [PharmCAT](https://github.com/PharmGKB/PharmCAT).


## Setting Up

This pipeline uses Bwakit, GATK and Picard.  Follow the [GATK installation instructions](http://gatkforums.broadinstitute.org/wdl/discussion/2899/howto-install-all-software-packages-required-to-follow-the-gatk-best-practices) to install bwa (0.7.13), picard (2.0.1), samtools and gatk (3.5.0).  Please take particular note to use Java 8, and to install the R libraries, as these are needed.  Update the paths in the `run.sh` if necessary.

We recommend placing using the following file structure:

````
pipeline:
|-- run.sh
|-- ext_data:
|   |-- hgbundle
|   '-- genome
'-- utils
    |-- bwa-0.7.12
    |-- picard-2.0.1
    '-- gatk-3.5-0
````    

In `/ext_data`:

* uncompress the __GATK hg38 bundle__ from ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hg38bundle.tar.gz

In `/ext_data/genome`:

* uncompress the __grch38 (no-alt) FASTA data__ from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
* download __grch38 (no alt) FASTA index__ from ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai

You will need to also need to index the genome for BWA. The following commands should automate the download and indexing:

```
> curl -o external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
> curl -o external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
> bwa index external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
> curl -o external_data/grc38.tar.gz ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hg38bundle.tar.gz
> tar -zxvf external_data/grc38.tar.gz -C external_data/
```


## Running

The pipeline can be run using the following command:

```
> ./run.sh fasta1.q fasta2.q @RG\tID:fasta1.gz\tPL:ILLUMINA\tSM:1 working_dir output_prefix
```

Where the _fasta1.q_ and _fasta2.q_ are pair reads from Illumia sequencing, _RG_ is the read group header in the format `@RG\tID:fasta1.gz\tPL:ILLUMINA\tSM:1`, _working_dir_ is the directory where the input fastaq files are and where the output will be put, and _output-prefix_ will be added to beginning of each of the output files.

This script has been tested on Ubuntu, takes about 5 days for 30x WGS data, and requires roughly 200GB of space. 


## Test installation

The following commands can be used to download some test files and check in the pipeline is running:

```
> mkdir test
> curl -o test/fasta1.gz ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
> curl -o test/fasta2.gz ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
> ./run.sh test/fasta1.gz test/fasta2.gz "@RG\tID:fasta1.gz\tPL:ILLUMINA\tSM:1" test run1
```

## Docker version

The pipeline can also be installed using this docker [file](DockerFile).  To run this version check out the project and cd into the main directory then create the docker image:

```
> docker build -t pharmngs .
```
    
Then create container using this image, which a shared volume for the data:

```
> docker run -t -v /data/<folder_with_ngs_data>:/usr/share/docker_data pharmngs /bin/bash
```
 
This gives you command line access to the container.  If you would like to try running the test data type:

```
> nohup ./run.sh test/fasta1.gz test/fasta2.gz "@RG\tID:fasta1.gz\tPL:ILLUMINA\tSM:1" test/ test1 > log.txt &
```
