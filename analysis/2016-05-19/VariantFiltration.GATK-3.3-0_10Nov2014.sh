#!/bin/bash

CORE_PATH="/isilon/sequencing/Seq_Proj/"
GATK_DIR="/isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_3/GenomeAnalysisTK-3.3-0"
REF_GENOME="/isilon/sequencing/GATK_resource_bundle/1.5/b37/human_g1k_v37_decoy.fasta"
JAVA_1_7="/isilon/sequencing/Kurt/Programs/Java/jdk1.7.0_25/bin"
KEY="/isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_2/lee.watkins_jhmi.edu.key"

PROJECT="H_Blakeley_NF1_WGHum-SeqWholeExome_160219_1"
PREFIX="H_Blakeley_NF1_WGHum-SeqWholeExome_160219_1"

$JAVA_1_7/java -jar /isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_3/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $REF_GENOME \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".1.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".2.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".3.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".4.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".5.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".6.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".7.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".8.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".9.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".10.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".11.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".12.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".13.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".14.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".15.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".16.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".17.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".18.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".19.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".20.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".21.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".22.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".X.normal.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".Y.normal.vcf" \
-o $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".raw.HC.vcf"

$JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REF_GENOME \
-selectType SNP \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".raw.HC.vcf" \
-et NO_ET \
-K $KEY \
-o $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".raw.HC.SNV.vcf"

$JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T SelectVariants \
-R $REF_GENOME \
-selectType INDEL \
-selectType MNP \
-selectType MIXED \
-selectType SYMBOLIC \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".raw.HC.vcf" \
--disable_auto_index_creation_and_locking_when_reading_rods \
-et NO_ET \
-K $KEY \
-o $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".raw.HC.INDEL.vcf"

# Not enough variants to perform VQSR, will do hard-cutoffs instead.

$JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $REF_GENOME \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".raw.HC.SNV.vcf" \
--filterExpression 'QD < 2.0' \
--filterName 'QD' \
--filterExpression 'MQ < 40.0' \
--filterName 'MQ' \
--filterExpression 'FS > 60.0' \
--filterName 'FS' \
--filterExpression 'MQRankSum < -12.5' \
--filterName 'MQRankSum' \
--filterExpression 'ReadPosRankSum < -8.0' \
--filterName 'ReadPosRankSum' \
--logging_level ERROR \
-et NO_ET \
-K $KEY \
-o $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.SNV.BestPractices.HardCutOffs.vcf"

$JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R $REF_GENOME \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".raw.HC.INDEL.vcf" \
--filterExpression 'QD < 2.0' \
--filterName 'QD' \
--filterExpression 'FS > 200.0' \
--filterName 'FS_INDEL' \
--filterExpression 'ReadPosRankSum < -20.0' \
--filterName 'ReadPosRankSum_INDEL' \
--logging_level ERROR \
-et NO_ET \
-K $KEY \
-o $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.INDEL.BestPractices.HardCutOffs.vcf"

$JAVA_1_7/java -jar /isilon/sequencing/CIDRSeqSuiteSoftware/gatk/GATK_3/GenomeAnalysisTK-3.1-1/GenomeAnalysisTK.jar \
-T CombineVariants \
-R $REF_GENOME \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.SNV.BestPractices.HardCutOffs.vcf" \
--variant $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.INDEL.BestPractices.HardCutOffs.vcf" \
-et NO_ET \
-K $KEY \
-o $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".BEDsuperset.BestPractices.HardCutoffs.vcf"

bgzip-0.2.6 -c $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".BEDsuperset.BestPractices.HardCutoffs.vcf" >| $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".BEDsuperset.BestPractices.HardCutoffs.vcf.gz"

tabix-0.2.6 -p vcf -f $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".BEDsuperset.BestPractices.HardCutoffs.vcf.gz"

# $JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
# -T VariantRecalibrator \
# -R $REF_GENOME \
# --input:VCF $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".raw.HC.SNV.vcf" \
# -resource:hapmap,known=false,training=true,truth=true,prior=15.0 /isilon/sequencing/GATK_resource_bundle/2.5/b37/hapmap_3.3.b37.vcf \
# -resource:omni,known=false,training=true,truth=true,prior=12.0 /isilon/sequencing/GATK_resource_bundle/2.5/b37/1000G_omni2.5.b37.vcf \
# -resource:1000G,known=false,training=true,truth=false,prior=10.0 /isilon/sequencing/GATK_resource_bundle/2.5/b37/1000G_phase1.snps.high_confidence.b37.vcf \
# -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /isilon/sequencing/GATK_resource_bundle/2.8/b37/dbsnp_138.b37.vcf \
# -mode SNP \
# -an QD \
# -an MQRankSum \
# -an ReadPosRankSum \
# -an FS \
# -tranche 100.0 \
# -tranche 99.9 \
# -tranche 99.8 \
# -tranche 99.7 \
# -tranche 99.6 \
# -tranche 99.5 \
# -tranche 99.4 \
# -tranche 99.3 \
# -tranche 99.2 \
# -tranche 99.1 \
# -tranche 99.0 \
# -tranche 98.0 \
# -tranche 97.0 \
# -tranche 96.0 \
# -tranche 95.0 \
# -tranche 90.0 \
# -recalFile $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.SNV.recal" \
# -tranchesFile $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.SNV.tranches" \
# -rscriptFile $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.SNV.R"
# 
# $JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
# -T VariantRecalibrator \
# -R $REF_GENOME \
# --input:VCF $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".raw.HC.INDEL.vcf" \
# -resource:mills,known=true,training=true,truth=true,prior=12.0 /isilon/sequencing/GATK_resource_bundle/2.2/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
# --maxGaussians 4 \
# -an MQRankSum \
# -an FS \
# -an ReadPosRankSum \
# -an QD \
# -mode INDEL \
# -tranche 100.0 \
# -tranche 99.9 \
# -tranche 99.8 \
# -tranche 99.7 \
# -tranche 99.6 \
# -tranche 99.5 \
# -tranche 99.4 \
# -tranche 99.3 \
# -tranche 99.2 \
# -tranche 99.1 \
# -tranche 99.0 \
# -tranche 98.0 \
# -tranche 97.0 \
# -tranche 96.0 \
# -tranche 95.0 \
# -tranche 90.0 \
# -recalFile $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.INDEL.recal" \
# -tranchesFile $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.INDEL.tranches" \
# -rscriptFile $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.INDEL.R"
# 
# $JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
# -T ApplyRecalibration \
# -R $REF_GENOME \
# --input:VCF $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".raw.HC.SNV.vcf" \
# --ts_filter_level 99.9 \
# -recalFile $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.SNV.recal" \
# -tranchesFile $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.SNV.tranches" \
# -mode SNP \
# -o $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.SNV.VQSR.vcf"
# 
# $JAVA_1_7/java -jar $GATK_DIR/GenomeAnalysisTK.jar \
# -T ApplyRecalibration \
# -R $REF_GENOME \
# --input:VCF $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".raw.HC.INDEL.vcf" \
# --ts_filter_level 99.9 \
# -recalFile $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.INDEL.recal" \
# -tranchesFile $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.INDEL.tranches" \
# -mode INDEL \
# -o $CORE_PATH/$PROJECT/MULTI_SAMPLE/$PREFIX".HC.INDEL.VQSR.vcf"
