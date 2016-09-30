#!/usr/bin/env bash

# check input
if test "$#" -ne 5; then
    echo "Use: run.sh fasta1.gz fasta2.gz @RG\tID:fasta1.gz\tPL:ILLUMINA\tSM:1 working_dir output_prefix"
    exit
fi

# Get inputs, fasta pairs, RG and working directory and output name
FASTA1=$1
FASTA2=$2
RG=$3
WORKING_DIR=$4
NAME=$5



GATK="gatk -Xms4g -Xmx4g"
PICARD="picard"
BWA="run-bwamem"

REF='external_data/genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
KNOWN_INDELS_MILLS='external_data/hg38bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz'
KNOWN_INDELS='external_data/hg38bundle/Homo_sapiens_assembly38.known_indels.vcf.gz'
KNOWN_DBSNP='external_data/hg38bundle/dbsnp_144.hg38.vcf.gz' # for now same
CURRENT_DBSNP='external_data/hg38bundle/dbsnp_144.hg38.vcf.gz'
HAPMAP='external_data/hg38bundle/hapmap_3.3.hg38.vcf.gz'
OMNI='external_data/hg38bundle/1000G_omni2.5.hg38.vcf.gz'
HC1000G='external_data/hg38bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz'

echo "Running grc38 pipeline"

echo "Running BWA:"
$BWA  -o $WORKING_DIR/output  -R "$RG"  -H -t 36  -a  -d  -s  -S  $REF  $FASTA1  $FASTA2  | sh

echo "indexing bam:"
$PICARD  BuildBamIndex  TMP_DIR=$WORKING_DIR/tmp  VALIDATION_STRINGENCY=STRICT  MAX_RECORDS_IN_RAM=8000000  INPUT=$WORKING_DIR/output.aln.bam  OUTPUT=$WORKING_DIR/output.aln.bam.bai

echo "Moving files"
mv  $WORKING_DIR/output.aln.bam  $WORKING_DIR/output.moved.bam
mv  $WORKING_DIR/output.aln.bam.bai  $WORKING_DIR/output.moved.bam.bai

echo "Realigning"
$GATK  -T RealignerTargetCreator  -I $WORKING_DIR/output.moved.bam  -ip 100  -R $REF  -nt 36  -o $WORKING_DIR/output.realigner.intervals  -known $KNOWN_INDELS_MILLS -known $KNOWN_INDELS
$GATK  -T IndelRealigner  -I $WORKING_DIR/output.moved.bam  -ip 100  -R $REF  -known $KNOWN_INDELS_MILLS -known $KNOWN_INDELS  -targetIntervals $WORKING_DIR/output.realigner.intervals  -o $WORKING_DIR/output.realigned.bam  -filterNoBases

echo "making dirs"
mkdir $WORKING_DIR/stats
mkdir $WORKING_DIR/tranches

echo "Recalibrating alignment"
$GATK  -T BaseRecalibrator  -I $WORKING_DIR/output.realigned.bam  -ip 100  -R $REF  -nct 36  -knownSites $KNOWN_INDELS_MILLS -knownSites $KNOWN_INDELS -knownSites $KNOWN_DBSNP  -o $WORKING_DIR/output.recal.before.table
$GATK  -T PrintReads  -I $WORKING_DIR/output.realigned.bam  -ip 100  -R $REF  -nct 36  -o $WORKING_DIR/output.bam
$GATK  -T BaseRecalibrator  -I $WORKING_DIR/output.realigned.bam  -ip 100  -R $REF  -BQSR $WORKING_DIR/output.recal.before.table  -nct 36  -knownSites $KNOWN_INDELS_MILLS -knownSites $KNOWN_INDELS -knownSites $KNOWN_DBSNP  -o $WORKING_DIR/output.recal.after.table
$PICARD  CollectMultipleMetrics  TMP_DIR=$WORKING_DIR/tmp  VALIDATION_STRINGENCY=STRICT  MAX_RECORDS_IN_RAM=8000000  INPUT=$WORKING_DIR/output.bam  OUTPUT=$WORKING_DIR/stats/output.metrics  PROGRAM= QualityScoreDistribution PROGRAM= MeanQualityByCycle PROGRAM= CollectAlignmentSummaryMetrics PROGRAM= CollectInsertSizeMetrics PROGRAM= CollectBaseDistributionByCycle

echo "Finishing alignment"


echo "Running calling"
$GATK  -T HaplotypeCaller  -I $WORKING_DIR/output.bam -ip 100 -R $REF  -disable_auto_index_creation_and_locking_when_reading_rods  -nct 16  -o $WORKING_DIR/output.g.vcf.gz  -ERC GVCF  -hets 0.001  -indelHeterozygosity 1.25E-4  -gt_mode DISCOVERY  -mbq 10  -maxReadsInRegionPerSample 1000  -minReadsPerAlignStart 5  -activeRegionExtension 100  -activeRegionMaxSize 300

echo "Genotyping gvcf"
$GATK  -T GenotypeGVCFs  -ip 100  -R $REF  -disable_auto_index_creation_and_locking_when_reading_rods  -nt 16  -V $WORKING_DIR/output.g.vcf.gz  -o $WORKING_DIR/$NAME.vcf.calls.vcf.gz  -allSites  -hets 0.001  -indelHeterozygosity 1.25E-4  -stand_call_conf 30.0  -stand_emit_conf 30.0

echo "Running annotator"
$GATK  -T VariantAnnotator  -I $WORKING_DIR/output.bam  -ip 100  -R $REF  -disable_auto_index_creation_and_locking_when_reading_rods  -nt 16  -V $WORKING_DIR/$NAME.vcf.calls.vcf.gz  -D $CURRENT_DBSNP  -o $WORKING_DIR/$NAME.vcf.annotated.vcf.gz  -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A StrandOddsRatio -A Coverage -A QualByDepth -A Coverage -A FisherStrand -A StrandOddsRatio -A ReadPosRankSumTest -A MappingQualityRankSumTest

echo "Variant recalibration"
$GATK  -T VariantRecalibrator  -ip 100  -R $REF  -disable_auto_index_creation_and_locking_when_reading_rods  -nt 16  -mode SNP  -mG 8  -input $WORKING_DIR/$NAME.vcf.annotated.vcf.gz  -resource:known=false,training=true,truth=true,prior=15.0 $HAPMAP -resource:known=false,training=true,truth=true,prior=12.0 $OMNI -resource:known=false,training=true,truth=false,prior=10.0 $HC1000G -resource:known=true,training=false,truth=false,prior=2.0 $KNOWN_DBSNP  -recalFile $WORKING_DIR/$NAME.vcf.snps.recal  -tranchesFile $WORKING_DIR/tranches/snps.tranches  -an MQ -an FS -an SOR -an DP -an QD  -tranche 100.0 -tranche 99.99 -tranche 99.9 -tranche 99.0 -tranche 90.0  -rscriptFile $WORKING_DIR/tranches/snps.R

echo "Apply variant recalibration"
$GATK  -T ApplyRecalibration  -ip 100  -R $REF  -disable_auto_index_creation_and_locking_when_reading_rods  -nt 16  -input $WORKING_DIR/$NAME.vcf.annotated.vcf.gz  -recalFile $WORKING_DIR/$NAME.vcf.snps.recal  -tranchesFile $WORKING_DIR/tranches/snps.tranches  -o $WORKING_DIR/$NAME.vcf.recal-snps.vcf.gz  -ts_filter_level 90.0  -mode SNP


echo "Indel recalibration"
$GATK  -T VariantRecalibrator  -ip 100  -R $REF  -disable_auto_index_creation_and_locking_when_reading_rods  -nt 16  -mode INDEL  -mG 4  -input $WORKING_DIR/$NAME.vcf.recal-snps.vcf.gz  -resource:known=false,training=true,truth=true,prior=12.0 $KNOWN_INDELS_MILLS -resource:known=true,training=false,truth=false,prior=2.0 $KNOWN_DBSNP  -recalFile $WORKING_DIR/$NAME.vcf.indels.recal  -tranchesFile $WORKING_DIR/tranches/indels.tranches  -an QD -an FS -an SOR -an MQ  -tranche 100.0 -tranche 99.99 -tranche 99.9 -tranche 99.0 -tranche 90.0  -rscriptFile $WORKING_DIR/tranches/indels.R

echo "Apply indel recalibration"
$GATK  -T ApplyRecalibration  -ip 100  -R $REF  -disable_auto_index_creation_and_locking_when_reading_rods  -nt 16  -input $WORKING_DIR/$NAME.vcf.recal-snps.vcf.gz  -recalFile $WORKING_DIR/$NAME.vcf.indels.recal  -tranchesFile $WORKING_DIR/tranches/indels.tranches  -o $WORKING_DIR/$NAME.vcf.recal-indels.vcf.gz  -ts_filter_level 90.0  -mode INDEL

echo "Analyze covariates"
$GATK  -T AnalyzeCovariates  -ip 100  -R $REF  -before $WORKING_DIR/output.recal.before.table  -after $WORKING_DIR/output.recal.after.table  -plots $WORKING_DIR/stats/output.recal.pdf  -csv $WORKING_DIR/stats/output.recal.csv


echo "Extract PGX genes for PharmCAT"
#Select a sample and restrict the output vcf to a set of intervals:

