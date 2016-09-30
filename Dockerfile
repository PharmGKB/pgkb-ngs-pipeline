FROM continuumio/anaconda


MAINTAINER "Lester Carter" lester@pharmgkb.org

# update
RUN apt-get -y update
RUN apt-get -y install zlib1g-dev vim python-software-properties apt-file software-properties-common build-essential curl git python-setuptools ruby

# Install bioconda packages
RUN conda config --add channels r
RUN conda config --add channels bioconda
RUN conda install -y picard
RUN conda install -y -c r r-essentials

# install gatk
RUN conda install -y gatk
ADD GenomeAnalysisTK.jar GenomeAnalysisTK.jar
RUN gatk-register GenomeAnalysisTK.jar

# Install BWA
RUN curl -L -o bwa.tar.bz2 https://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.12_x64-linux.tar.bz2/download && tar -xvf bwa.tar.bz2
ENV PATH /bwa.kit/:$PATH

# get BWA genome and prepare sequence
RUN mkdir external_data \
&& mkdir external_data/genome \
&& curl -o external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai \
&& curl -o external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
&& gunzip -c external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz  > external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
&& bwa index external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna

# Create Picard dict
RUN picard CreateSequenceDictionary R=external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna O=external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.dict

#download gatk grc38 reference pack
RUN curl -o external_data/grc38.tar.gz ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg38/hg38bundle.tar.gz && tar -zxvf external_data/grc38.tar.gz -C external_data/

#download giab example files
RUN mkdir test \
&& curl -o test/fasta1.gz ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz \
&& curl -o test/fasta2.gz ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz

# Add run script and Rscript
ADD run.sh run.sh
ADD Rprofile .Rprofile
ADD Rsetup/install.R Rsetup/install.R
RUN Rscript Rsetup/install.R
