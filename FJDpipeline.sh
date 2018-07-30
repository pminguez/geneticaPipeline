#!/bin/bash
#Implemented pipepline for trios.
#Set the pathways to work for you.
#MDAP=main directory analysis patwhay (all the data will be generated in this file).
#FQ/   FOLDER     where fastq files rest. (only for pair end).
#HG19: where all the references and indexed files will be + bundle data from UCSC.
#mapped_data obtaining the .SAM files after alignment using BWA.
#sorted_data sorted SAM with PICARD.
#dedupped_data DUPLICATES MARKED data using PICARD.
#recalibrated_base_quality_score_read_data, data after quality check, using GATK.
#applied_base_quality_score_read_data, data after selected the file of qualities and applied the MARK_scores file, using GATK.
#haplotype_caller_data_gvcf, data after calling for variants, using GATK.
#OPTIONAL:combined_data_gvcf IF MORE THAN ONE INPUT.
#genotyped_data_vcf, ONE INPUT or after using CombineGVCF, performs the joint genotyping on GVCFS produced by Hcaller, using GATK, obtaining the VCF.
#selectVariantsData_vcfs, selecting the qualities and the desired FILTERING options.


#software, where all the programs are located.

MDAP='/home/marius/testFJDRP/DataAnalisis3'
FQ='/home/marius/testFJDRP/DataAnalisis3/fastq' #Seleccionar la carpeta
HG19='/home/marius/testFJDRP/hg19'
MD='/home/marius/testFJDRP/DataAnalisis3/mapped_data'
SD='/home/marius/testFJDRP/DataAnalisis3/sorted_data'
DD='/home/marius/testFJDRP/DataAnalisis3/dedupped_data'
RBQSRD='/home/marius/testFJDRP/DataAnalisis3/recalibrated_bqsr_data'
ABQSRD='/home/marius/testFJDRP/DataAnalisis3/applied_bqsr_data'
HCDGVCF='/home/marius/testFJDRP/DataAnalisis3/haplotypeCaller_data_gvcf'
CGVCF='/home/marius/testFJDRP/DataAnalisis3/combined_gvcf'
GDVCF='/home/marius/testFJDRP/DataAnalisis3/genotyped_data_vcf'
VFDVCF='/home/marius/testFJDRP/DataAnalisis3/variant_filtration_data_vcf'
SVDVCF='/home/marius/testFJDRP/DataAnalisis3/selecVariants_data_vcf'
VEPVCFA='/home/marius/testFJDRP/DataAnalisis3/vep_vcf_annotated'
VTTVCF='/home/marius/testFJDRP/DataAnalisis3/variantstotable_vcf'
SFT='/home/marius/software'

#create a register logfile
touch registerFile

echo ············································································································

echo -e                                           "\n \tINDEXING REFERENCE FILES (BWA)\n"

echo ············································································································

#Start BWA INDEX.
#echo "Starts BWA INDEX"
#echo "Starts BWA INDEX" >>  registerFile
#$SFT/bwa/./bwa index $HG19/ucsc.hg19.fasta
#echo -e "\nINDEXADO COMPLETADO" ; paplay /usr/share/sounds/freedesktop/stereo/complete.oga
#echo -e "\nINDEXADO COMPLETADO" >> registerFile
#Creating .FAI in HG19.
#echo "Create .FAI file, using samtools faidx"
#echo "Create .FAI file, using samtools faidx" >> registerFile
#$SFT/samtools/./samtools faidx $HG19/ucsc.hg19.fasta -o $HG19/ucsc.hg19.fai
#echo -e "\n.FAI COMPLETADO"
#echo -e "\n.FAI COMPLETADO" >> registerFile

#Creating .DICT in HG19.
#echo "Create .DICT file, using picardtools CreateSequnceDictionary"
#echo "Create .DICT file, using picardtools CreateSequnceDictionary" >> registerFile
#java -jar $SFT/picard/build/libs/picard.jar CreateSequenceDictionary \
#R=ucsc.hg19.fasta \
#O=$HG19/ucsc.hg19.dict
#echo -e "\.DICT COMPLETADO"
#echo -e "\.DICT COMPLETADO">> registerFile


echo ············································································································

echo -e                                              "\n \tMAPPING (BWA)\n"

echo ············································································································

#mapping data to Reference after BWA INDEX using UCSC.HG19.FASTA
mkdir mapped_data
echo "mkdir mapped_data" >> registerFile
for i in $@
do
	#Run BWA mem -t 12
	#-t threads
	#-P search for Pair mate if not mapped properly, if it found a better hit, skips it.
	echo "Start BWA MEM '$1'"
	echo "Start BWA MEM '$1'">> registerFile
	$SFT/bwa/./bwa mem -t 12 -R '@RG\tID:$1\tPL:illumina\tSM:Analysis' $HG19/ucsc.hg19.fasta \
	$FQ/$i*_1.fastq.gz \
	$FQ/$i*_2.fastq.gz | gzip -3 > $MD/mapped$i.sam.gz
	echo -e "\nBWA MEM '$i'  COMPLETADO" ; paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	echo -e "\nBWA MEM '$i'  COMPLETADO"  >> registerFile
	#Unzip all the files before continuing the process.
	echo "Unzip mapped sam's."
	echo "Unzip mapped sam's." >> registerFile
	gunzip -k $MD/*.sam.gz
	echo "Gunzip completed."
	echo "Gunzip completed." >> registerFile
done



echo ···········································································································

echo -e                                         "\n \tSORTING SAM (PICARD)\n"

echo ···········································································································



mkdir sorted_data
echo "mkdir sorted_data">> registerFile
#Sorting the mapped data.
for i in $@
do
	#SORTING THE SAM FILE.
	echo "Run picard SortSam '$i'"
	echo "Run picard SortSam '$i'">> registerFile
	java -jar $SFT/picard/build/libs/picard.jar SortSam I=$MD/mapped$i.sam \
	O=$SD/sorted$i.bam \
	SORT_ORDER=coordinate
	echo -e "\nPicard SortSam '$i'  COMPLETADO" ; paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	echo -e "\nPicard SortSam '$i'  COMPLETADO" >> registerFile
done


echo ············································································································

echo -e                                         "\n \tMARKING DUPLICATES (PICARD)\n"

echo ············································································································



#Selecting the duplicates reads from the mapped and sorted reads.
mkdir dedupped_data #mkdir deduppe_data_$i
echo "mkdir dedupped_data">> registerFile
for i in $@
do
	#Mark duplicates PICARD
	echo "Start picard MarkDuplicates '$i' "
	echo "Start picard MarkDuplicates '$i' ">>registerFile
	java -jar $SFT/picard/build/libs/picard.jar MarkDuplicates \
	I=$SD/sorted$i.bam \
	O=$DD/dedupped$i.bam \
	M=$DD/marked_dup_metrics$i.txt \
	REMOVE_DUPLICATES=true \
	AS=SortOrder
	echo -e "\n PICARD MarkDuplicates '$i' COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	echo -e "\n PICARD MarkDuplicates '$i' COMPLETADO"  >> registerFile
done


echo ············································································································
echo -e                                                  	"\n \tBQSR (GATK)\n"
echo -e                                                        "\t-.1 Recalibrate Data\"
echo -e								"\t-.2 Apply Recalibration\n"
echo ············································································································


#Create a .BAI file to compare original vs removed duplicates reads.
for i in $@
do
	#Indexing the BAM files.
	#Generating the .BAI files from the DEDUPED (markedDuplicates from the original SAM/BAM file.
	echo "Indexing '$i' BAM files"
	echo "Indexing '$i' BAM files" >> registerFile
	java -jar $SFT/picard/build/libs/picard.jar BuildBamIndex \
	I=$DD/dedupped$i.bam \
	O=$DD/$i.dedupped.bai
	echo -e "\n PICARD BuildBamIndex '$i' COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	echo -e "\n PICARD BuildBamIndex '$i' COMPLETADO" >> registerFile
done

#Recalibrrating  the reads using base quality score reads.
mkdir recalibrated_bqsr_data
echo "mkdir recalibrated_bqsr_data" >> registerFile
for i in $@
do
	#GATK RECALIBRADO
	#BaseRecalibration + table
	echo "Starts GATK '$1' Recalibrator"
	echo "Starts GATK '$1' Recalibrator" >> registerFile
	java -jar $SFT/gatk/build/libs/gatk.jar BaseRecalibrator \
	-I $DD/dedupped$i.bam \
	-R $HG19/ucsc.hg19.fasta \
	--known-sites $HG19/dbsnp_138.hg19.vcf \
	--known-sites $HG19/1000G_phase1.indels.hg19.sites.vcf \
	--known-sites $HG19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
	-O $RBQSRD/recalibrated_bqsr_data$i.table.bam  #--bqsr 1st_racalibrationNIST7035.table
	echo -e "\n GATK BaseRecalibrator '$i' COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	echo -e "\n GATK BaseRecalibrator '$i' COMPLETADO" >>registerFile
done

mkdir applied_bqsr_data
echo "mkdir applied_bqsr_data" >> registerFile
for i in $@
do
	#ApplyBQSR
	echo "Starts picard  '$i' ApplyBQSR"
	echo "Starts picard  '$i' ApplyBQSR" >> registerFile
	java -jar $SFT/gatk/build/libs/gatk.jar ApplyBQSR \
	-R $HG19/ucsc.hg19.fasta \
	-I $DD/dedupped$i.bam \
	--bqsr $RBQSRD/recalibrated_bqsr_data$i.table.bam  \
	-O $ABQSRD/applied_bqsr_data$i.bam

	echo -e "\nGATK ApplyBQSR '$i' COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	echo -e "\nGATK ApplyBQSR '$i' COMPLETADO" >>registerFile
done


echo ············································································································

echo -e                                           "\n \tHAPLOTYPE CALLER (GATK)\n"

echo ············································································································

#Ready to call for Variants.
mkdir haplotypeCaller_data_gvcf
echo "mkdir haplotypeCaller_data_gvcf" >> registerFile
for i in $@
do
	#HaplotypeCaller for each sample for later joint genotyping.
	echo -e "\nGATK HaplotypeCallerGVCF for '$i' STARTS"
	echo -e "\nGATK HaplotypeCallerGVCF for '$i' STARTS">> registerFile
	java -jar $SFT/gatk/build/libs/gatk.jar HaplotypeCaller \
	-R $HG19/ucsc.hg19.fasta \
	-I $ABQSRD/applied_bqsr_data$i.bam \
	-ERC GVCF \
	-bamout $HCDGVCF/HCbamout$i.bam \
        -O $HCDGVCF/HCdata$i.g.vcf \
	-G StandardAnnotation \
	-G AS_StandardAnnotation \
	-G StandardHCAnnotation \
	-A QualByDepth \
	-A FisherStrand \
	-A StrandOddsRatio \
	-A RMSMappingQuality \
	-A MappingQualityRankSumTest \
	-A ReadPosRankSumTest \
	-A DepthPerSampleHC \
	-A BaseQualityRankSumTest \
	-A ExcessHet \
	-A StrandArtifact \
	--annotate-with-num-discovered-alleles=true

	echo -e "\nGATK HaplotypeCallerGVCF ERC GVCF '$i' COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	echo -e "\nGATK HaplotypeCallerGVCF ERC GVCF '$i' COMPLETADO" >> registerFile
done

echo ············································································································

echo -e                                            "\n \tJOINT GENOTYPING (GATK)\n"

echo ············································································································



#If more than one input the pipeline will continue the Joint analysis.
for i in $@
do
	if [ $# > 2 ]
		then
		mkdir combined_gvcf
			echo "mkdir combined_gvcf" >> registerFile
			#Obtenido en GCVF pasamos al Joint Genotyping on one or more samples called with HC.
			echo -e "\nUsing GATK COMBINEGVCFs for merging GVCFs"
			echo -e "\nUsing GATK COMBINEGVCFs for merging GVCFs">> registerFile

			java -jar $SFT/gatk/build/libs/gatk.jar CombineGVCFs \
			-R $HG19/ucsc.hg19.fasta \
			--variant $HCDGVCF/HCaller$1.g.vcf \
			--variant $HCDGVCF/HCaller$2.g.vcf \
			--variant $HCDGVCF/HCaller$3.g.vcf \
			--variant $HCDGVCF/HCaller$4.g.vcf \
			--variant $HCDGVCF/HCaller$5.g.vcf \
			-O $CGVCF/combined.g.vcf

			echo -e "\nGATK COMBINEGVCFs COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
			echo -e "\nGATK COMBINEGVCFs COMPLETADO"  >> registerFile
			mkdir genotyped_data_vcf
	       	        echo "mkdir genotyped_data_vcf">> registerFile
	                #GenotypeGVCFs into final VCF
	                echo -e "\nUsing GATK GenotypeGVCFs for final VCF"
	                java -jar $SFT/gatk/build/libs/gatk.jar GenotypeGVCFs \
	                -R $HG19/ucsc.hg19.fasta \
	                -V $HCDGVCF/combined.g.vcf \
       		        -G StandardAnnotation \
	                -G AS_StandardAnnotation \
	                -O $GDVCF/genotyped_data$i.vcf
	                echo -e "\nGATK GenotypeGVCFs COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	                echo -e "\nGATK GenotypeGVCFs COMPLETADO" >> registerFile
#if only one INPUT continue to joint analysis.
		else
			mkdir genotyped_data_vcf
			echo "mkdir genotyped_data_vcf">> registerFile
			#GenotypeGVCFs into final VCF
			echo -e "\nUsing GATK GenotypeGVCFs for final VCF"
			java -jar $SFT/gatk/build/libs/gatk.jar GenotypeGVCFs \
			-R $HG19/ucsc.hg19.fasta \
			-V $HCDGVCF/HCdata$i.g.vcf \
			-G StandardAnnotation \
			-G AS_StandardAnnotation \
			-O $GDVCF/genotyped_data$i.vcf
			echo -e "\nGATK GenotypeGVCFs COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
			echo -e "\nGATK GenotypeGVCFs COMPLETADO" >> registerFile
	fi
done
#echo ············································································································

#echo -e        "\n \tIf more than 30 samples running VQSR for variants detection using machine learning with the smaples"
#echo -e            "\tConsidering using the BAM from 1000G for VCF calling and using them to train the detector"
#echo -e	     		"\tand being able to search for OUR VARIANTS in OUR SAMPLE wiht VQSR\n"

#echo ············································································································



#echo ············································································································

#echo -e  "\n \tHard filtering if less than 30 samples and doing it in the classical way, selecting variants\n"

#echo ············································································································



echo ············································································································

echo -e               			 "\n \tHard filtering (GATK)"
echo -e                           "\tDATA PREFILTERING  = (vcf_processing_step2.R)\n"

echo ············································································································


#to be added the data for allele and variants frequencies gnomad, exac.


#HARD FILTERING

#Data Heterozygosity and Pedigree specification.
#mkdir variant_filtration_data_vcf
#echo "variantfiltration_data_vcf ">> registerFile
#Label the heterozygous genotypes and select the Heterozygous=1.

#java -jar $SFT/gatk/build/libs/gatk.jar VariantFiltration \
#-R $HG19/ucsc.hg19.fasta \
#-V $GDVCF/genotyped_data$i.vcf \
#-O $VFDVCF/variantFiltered$i.vcf\
#--genotype-filter-expression "isHet == 1" \
#--genotype-filter-name "isHetFilter" \
#--set-filtered-genotype-to-no-call=true 
#--filterName /home/marius/... path to file with myfilters.

#Data HARD-FILTERING.
#mkdir selectVariants_data_vcf
#echo "selectVariants_data_vcf" >> registerFile
#Select filtering values, and set filtered GT  with  isHet=1 to nocall ./. 0/1 and 1/0--> ./.  Only 1/1.
#java -jar $SFT/gatk/build/libs/gatk SelectVariants \
#-R $HG19/ucsc.hg19.fasta \
#-V $GDVCF/variantFiltered$i.vcf \
#-select "QD > 2.0" \
#-select "QUAL > 100" \
#-select "MQ > 40.0" \
#-select "FS < 65.0" \
#-select "ReadPosRankSum > -8.0" \
#--pedigree TOBEADDED
#-O $SVDVCF/selectVariants_data$i.vcf

echo ············································································································

echo -e                  		 "\n \tVARIANT ANNOTATION (VEP ENSEMBL)\n"

echo ············································································································




#Annotate the VCF using VEP(variant effect predictor).

#$SFT/vep/.vep \
#-i selectVariants_data$i.vcf \
#--per_gene \
#--coding_only \
#--offline \
#-o $VEPVCFA/vepAnnotatedVCF$@.vcf

echo ············································································································

echo -e                        	"\n \tVQSR (more than 30 samples (1000genomes + samples)\n"

echo ············································································································


#VQSR

#Building the SNP recalibration model.
#java -jar $SFT/gatk/build/libs/gatk.jar VariantRecalibrator \
# -R $HG19/ucsc.hg19.fasta \
# -V myvariants.vcf \
# --resource hapmap,known=false,training=true,truth=true,prior=15.0:$HG19/hapmap_3.3.hg19.sites.vcf \
# --resource omni,known=false,training=true,truth=false,prior=12.0:$HG19/1000G_omni2.5.hg19.sites.vcf \
# --resource 1000G,known=false,training=true,truth=false,prior=10.0:$HG19/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
# --resource dbsnp,known=true,training=false,truth=false,prior=2.0:$HG19/dbsnp_138.hg19.vcf \
# -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an StrandArtifact -an ExcessHet \
# -mode SNP \
# --output recalibrated.SNP.recal \
# --tranches-file recal.SNP.tranches \
# --rscript-file recal.SNP.R

#Building the INDEL recalibration model.
#java -jar $SFT/gatk/build/libs/gatk.jar VariantRecalibrator \
#-R $HG19/ucsc.hg19.fasta \
#-V myvariants.vcf \
#--resource mills,known=false,training=true,truth=true,prior=12.0:$HG19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
# --resource dbsnp,known=true,training=false,truth=false,prior=2.0:$HG19/dbsnp_138.hg19.vcf \
# -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an StrandArtifact -an ExcessHet \
# -mode INDEL \
# --output recalibrated.INDEL.recal \
# --tranches-file recal.INDEL.tranches \
# --rscript-file recal.INDEL.R

#ApplyRecalibration step, where the output will be input file + the SNPS anotated with their recalibrated quality scores.
# either with PASS or FILTER depending on wheter or not they are included in the selected tranche.
#the biggher is the tsfilter, less variants leaving outside being considerated is inverse..

#java -jar $SFT/gatk/build/libs/gatk.jar ApplyVQSR \
#-R $HG19/ucsc.hg19.fasta \
#-
#TOBEADDED 		TOBEADDED				TOBEADDED

echo ············································································································

echo -e                                  	"\n \tVARIANTS TO TABLE (GATK)"
echo -e 	                         "\tPrefiltering = vcf_processing_step1_v2.py\n"

echo ············································································································

#USING GATK, SELECT DATA FR0M THE VCF VILE TO SEPARATE IN DIFFERENT COLUMNS.
#-F VCF DATA FROM (INFO)
#-GF VCF DATA FROM (FORMAT)

mkdir variantstotable_vcf
echo "mkdir variantstotable_vcf">> registerFile
echo -e "\nUsing GATK VariantsToTable for final VCF sep fileds"
#VariantsToTable avoinding extra scripts.
java -jar $SFT/gatk/build/libs/gatk.jar VariantsToTable \
-V $VEPVCFA/vepAnnotatedVCF$@.vcf \ #VCF file
#select all the FIELDS  wanted to be separated.
-F CHROM -F POS -F TYPE -F REF -F ALT -F QUAL -F QD -F FS -F SOR -F MQ -F RMSMappingQuality -F MappingQualityRankSumTest -F ReadPosRankSumTest -F ExcessHet -F MLEAF -F MLEAC -F HET -F HOM-REF -F HOM-VAR -GF GT -GF AD -GF DP -GF GQ -GF PL  \
-O $VTTVCF/vairantTotablegenotypedvep4$i.table.tsv \
--error-if-missing-data #data NA if missing.
echo -e "\nGATK VariantsToTable COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
echo -e "\nGATK VariantsToTable COMPLETADO" >> registerFile


#info
echo ·····················································································

echo -e \n \tAnalisis COMPLETADO para: '$@', el dia, "$(date +%Y/%m/%d) > READin.txt"
echo -e \n \tAnalisis COMPLETADO para: '$@', el dia, "$(date +%Y/%m/%d) >> registerFile"

echo ·····················································································
