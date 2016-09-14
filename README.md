# pgkb-ngs-pipeline

##  PharmGKB NGS Pipeline

This repository contains a simplified version of our work-in-progress grc38 NGS pipeline for pharmacogenomics - PGxtract.  run.sh is a simple bash script to show the steps of our pipeline and the parameters used. It has been tested on Ubuntu, takes about 5 days for 30x WGS data, and requires roughly 200GB of space.  It takes a pair of fastq files from an Illumina NGS run, aligns them using BWA, and then calls variants using GATK.  The final vcf is genotyped at *every* position in the file, which is useful for pharmacogenomics analysis but produces large files.  A smaller vcf containing only the positions needed for [PharmCAT](https://github.com/PharmGKB/PharmCAT) is also produced. This pipeline is optimized for WGS and Exome analysis, and currently does not handle trios.  [VQSR](http://gatkforums.broadinstitute.org/gatk/discussion/39/variant-quality-score-recalibration-vqsr) filtering has shown to be effective for WGS data, but if you may want to switch to [hard filtering](http://gatkforums.broadinstitute.org/wdl/discussion/2806/howto-apply-hard-filters-to-a-call-set) if using other data sources. 

## Installation

THis pipeline uses bwakit, GATK and Picard.  Use the [GATK installation instructions](http://gatkforums.broadinstitute.org/wdl/discussion/2899/howto-install-all-software-packages-required-to-follow-the-gatk-best-practices) to install bwa (0.7.13), picard (2.0.1), samtools and gatk (3.5.0), then update the paths in the run.sh if necessary.  We recommend placing using the following file structure:

    
    pipeline:
        run.sh
        -> ext_data:
            ->hgbundle
            ->genome
        -> utils:
            ->bwa-0.7.12
            ->picard-2.0.1
            ->gatk-3.5-0
        -> files:
            fastq1.gz
            fastq2.gz
        -> files2
            second__fastq2.gz
            seconde__fastq2.gz
        -> etc
    
Please take particular note to use Java 8, and to install the R libraries, as these are needed. You will also need to download the [GATK hg38 bundle](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hg38bundle.tar.gz) and place it in the external_data folder, and the [grc38 (no-alt) genome](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz) and [index] (ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai), which should be placed in external_data/genome. You will need to also need to index the genome for bwa. The following commands should automate the download and indexing:

    curl -o external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
    curl -o external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    bwa index external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    curl -o external_data/grc38.tar.gz ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hg38bundle.tar.gz
    tar -zxvf external_data/grc38.tar.gz -C external_data/

## Running
 
 THe pipeline can be run using the following command:
  
  `./run.sh fasta1.q fasta2.q @RG\tID:fasta1.gz\tPL:ILLUMINA\tSM:1 working_dir output_prefix`
  
  Where the fasta1.q and fasta2.q are pair reads from Illumia sequencing, RG is the read group header in the format "@RG\tID:fasta1.gz\tPL:ILLUMINA\tSM:1", working_dir is the directory where the input fastaq files are and where the output will be put, and output-prefix will be added to beginning of each of the output files.
  

## Test installation

The following commands can be used to download some test files and check in the pipeline is running:

    mkdir test
    curl -o test/fasta1.gz ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
    curl -o test/fasta2.gz ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz
     ./run.sh test/fasta1.gz test/fasta2.gz "@RG\tID:fasta1.gz\tPL:ILLUMINA\tSM:1" test run1
     

## TODO

This is an experimental work in progress version of this pipeline, made available only to share and refine the procedure with collaborators.  Later versions will be optimized for parallel operation. Please feel free to take the steps in run.sh and adapt them.



