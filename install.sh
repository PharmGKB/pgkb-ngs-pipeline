#!/usr/bin/env bash


# Check if anaconda is installed
hash conda &> /dev/null
if [ $? -eq 1 ]; then
    echo >&2 "Conda not found. You can download it from https://docs.continuum.io/anaconda/install.  This is how we \
    suggest installing bwakit, picard and gatk.  However you can also install manually if you wish."
    exit
fi

# Use anaconda to install bwakit and gatk
conda config --add channels r
conda config --add channels bioconda
conda install -y picard
conda install -y -c r r-essentials

# Download reference files
mkdir external_data
mkdir external_data/genome
curl -o external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
curl -o external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz

# reindex with current bwa version
echo "reindexing reference - this may take 5 to ten minutes"
bwa index external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz


#download gatk grc38 reference pack
curl -o external_data/grc38.tar.gz ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hg38bundle.tar.gz
tar -zxvf external_data/grc38.tar.gz -C external_data/


#download giab example files
mkdir test
curl -o test/fasta1.gz ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz
curl -o test/fasta2.gz ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz

# last step is gatk, as this requires a manual step you will need to do it yourself:
echo "to finish installation type gatk-register at the command line and follow the instructions"
