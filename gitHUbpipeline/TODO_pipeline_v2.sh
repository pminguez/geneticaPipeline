#!/bin/bash

##Implemented pipepline for comparing WES and obtainning gVCF for later trios comparissons.

##Set the pathways to work for you.

##MDAP=main directory analysis patwhay (all the data will be generated in this file).
##FQ/   FOLDER     where fastq files rest. (only for pair end).
##HG19: where all the references and indexed files will be + bundle data from UCSC.
##mapped_data obtaining the .SAM files after alignment using BWA.
##sorted_data sorted SAM with PICARD.
##dedupped_data DUPLICATES MARKED data using PICARD.
##recalibrated_base_quality_score_read_data, data after quality check, using GATK.
##applied_base_quality_score_read_data, data after selected the file of qualities and applied the MARK_scores file, using GATK.
##haplotype_caller_data_gvcf, data after calling for variants, using GATK.
##OPTIONAL:combined_data_gvcf IF MORE THAN ONE INPUT.
##genotyped_data_vcf, ONE INPUT or after using CombineGVCF, performs the joint genotyping on GVCFS produced by Hcaller, using GATK, obtaining the VCF.
##selectVariantsData_vcfs, selecting the qualities and the desired FILTERING options.


##software, where all the programs are located.

MDAP='/home/marius/testFJDRP/DataAnalisis8'
FQ='/home/marius/genetica3/marius/fastq/RP-1773/fastq' #Seleccionar la carpeta
HG19='/home/marius/testFJDRP/hg19'
MD='/home/marius/testFJDRP/DataAnalisis8/mapped_data'
SD='/home/marius/testFJDRP/DataAnalisis8/sorted_data'
DD='/home/marius/testFJDRP/DataAnalisis8/dedupped_data'
RBQSRD='/home/marius/testFJDRP/DataAnalisis8/recalibrated_bqsr_data'
PRD='/home/marius/testFJDRP/DataAnalisis8/plotRecalibration_data'
ABQSRD='/home/marius/testFJDRP/DataAnalisis8/applied_bqsr_data'
HCDGVCF='/home/marius/testFJDRP/DataAnalisis8/haplotypeCaller_data_gvcf'
CGVCF='/home/marius/testFJDRP/DataAnalisis8/combined_gvcf'
GDVCF='/home/marius/testFJDRP/DataAnalisis8/genotyped_data_vcf'
VFDVCF='/home/marius/testFJDRP/DataAnalisis8/variant_filtration_data_vcf'
SVDVCF='/home/marius/testFJDRP/DataAnalisis8/selecVariants_data_vcf'
VEPVCFA='/home/marius/testFJDRP/DataAnalisis8/vep_vcf_annotated'
VTTVCF='/home/marius/testFJDRP/DataAnalisis8/variantstotable_vcf'
SFT='/home/marius/software'
###
##SAMPLE_NAME=$ 	todo





#echo ············································································································

#echo -e                                           "\n \tINDEXING REFERENCE FILES (BWA)\n"

#echo ············································································································





##Start BWA INDEX.
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

###Creating .DICT in HG19.
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





##mapping data to Reference after BWA INDEX using UCSC.HG19.FASTA
mkdir mapped_data
echo "mkdir mapped_data" >> registerFile
for i in $@
do
#Run BWA mem -t 12
	#-t threads
	#-P search for Pair mate if not mapped properly, if it found a better hit, skips it.
	echo "Start BWA MEM '$1'"
	echo "Start BWA MEM '$1'">> registerFile
	$SFT/bwa/./bwa mem -t 12 -R '@RG\tID:$@\tPL:illumina\tSM:Analysis' $HG19/ucsc.hg19.fasta \
	$FQ/$i*1.fastq.gz \
	$FQ/$i*2.fastq.gz | gzip -3 > $MD/mapped$i.sam.gz
	echo -e "\nBWA MEM '$i'  COMPLETADO" ; paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	echo -e "\nBWA MEM '$i'  COMPLETADO"  >> registerFile
	#Unzip all the files before continuing the process.
	echo "Unzip mapped sam's."
	echo "Unzip mapped sam's." >> registerFile
	gunzip -k $MD/mapped$i.sam.gz
	echo "Gunzip completed."
	echo "Gunzip completed." >> registerFile
done





echo ···········································································································

echo -e                                         "\n \tSORTING SAM (PICARD)\n"

echo ···········································································································






mkdir sorted_data
echo "mkdir sorted_data">> registerFile
##Sorting the mapped data.
for i in $@
do
##SORTING THE SAM FILE.
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







##Selecting the duplicates reads from the mapped and sorted reads.
mkdir dedupped_data #mkdir deduppe_data_$i
echo "mkdir dedupped_data">> registerFile
for i in $@
do
##Mark duplicates PICARD
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

echo -e                                         "\n \t Dedupped BAM index (BAI) file. \n"

echo ············································································································






##Create a .BAI file to compare original vs removed duplicates reads.
for i in $@
do
	##Indexing the BAM files.
	##Generating the .BAI files from the DEDUPED (markedDuplicates from the original SAM/BAM file.
	echo "Indexing '$i' BAM files"
	echo "Indexing '$i' BAM files" >> registerFile
	java -jar $SFT/picard/build/libs/picard.jar BuildBamIndex \
	I=$DD/dedupped$i.bam \
	O=$DD/$i.dedupped.bai
	echo -e "\n PICARD BuildBamIndex '$i' COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	echo -e "\n PICARD BuildBamIndex '$i' COMPLETADO" >> registerFile
done







echo ············································································································
echo -e                                     '\n \tBQSR (GATK)\n'
echo -e                                      '\t -.1 Recalibration data table\n'
echo -e                                      '\t -.2 Recalibration data table \n'
echo ············································································································





##Recalibrrating  the reads using base quality score reads.
mkdir recalibrated_bqsr_data
echo "mkdir recalibrated_bqsr_data" >> registerFile
for i in $@
do
	##GATK BaseRecalibration first table
	##BaseRecalibration + table
	echo "Starts GATK '$1' Recalibrator"
	echo "Starts GATK '$1' Recalibrator" >> registerFile
	java -jar $SFT/gatk/build/libs/gatk.jar BaseRecalibrator \
     	-R $HG19/ucsc.hg19.fasta \
	-I $DD/dedupped$i.bam \
	--known-sites $HG19/dbsnp_138.hg19.vcf.gz \
	--known-sites $HG19/1000G_phase1.indels.hg19.sites.vcf.gz \
	--known-sites $HG19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz \
	-O $RBQSRD/before_recalibrated_bqsr_data$i.recal.table  #--bqsr 1st_racalibration.table
	echo -e "\n GATK BaseRecalibrator '$i' COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	echo -e "\n GATK BaseRecalibrator '$i' COMPLETADO" >>registerFile
done



echo ············································································································

echo -e                                       '\n \tApplying BQSR GATK\n'

echo ············································································································




##Applying the recalibration table to the bam file to continue the analysis.
mkdir applied_bqsr_data
echo "mkdir applied_bqsr_data" >> registerFile
for i in $@
do
	##ApplyBQSR
	echo "Starts picard  '$i' ApplyBQSR"
	echo "Starts picard  '$i' ApplyBQSR" >> registerFile
	java -jar $SFT/gatk/build/libs/gatk.jar ApplyBQSR \
	-R $HG19/ucsc.hg19.fasta \
	-I $DD/dedupped$i.bam \
	--bqsr $RBQSRD/before_recalibrated_bqsr_data$i.recal.table \
	-O $ABQSRD/applied_bqsr_data$i.bam
	echo -e "\nGATK ApplyBQSR '$i' COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	echo -e "\nGATK ApplyBQSR '$i' COMPLETADO" >>registerFile
done




echo ············································································································

echo -e             '\n \tAnalyze Covariates for PDF plots comparation of the Recalibration GATK\n'

echo ············································································································




mkdir plotRecalibration_data
echo "mkdir plotRecalibration_data" >> registerFile
for i in $@
do
##GATK BaseRecalibration second table for next step AnalyzeCovariates.
##Generates the second pass table.
##instead of second table, we use the BAM created by ApplyBQSR to regenrate a new TABLE for plot.
        echo "Starts GATK '$1' Second Recalibration"
        echo "Starts GATK '$1' Second Recalibration" >> registerFile
        java -jar $SFT/gatk/build/libs/gatk.jar BaseRecalibrator \
        -I $ABQSRD/applied_bqsr_data$i.bam \
        -R $HG19/ucsc.hg19.fasta \
        --known-sites $HG19/dbsnp_138.hg19.vcf.gz \
        --known-sites $HG19/1000G_phase1.indels.hg19.sites.vcf.gz \
        --known-sites $HG19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz \
        -O $PRD/after_recalibrated_bqsr_data$i.recal.table            #second bqsr table for comparation.
done


##Full generating the Plots of the recalibration tables using AnalyzeCovariates and saving a csv copy.
##Analyze  the tables.
for i in $@
do
	java -jar $SFT/gatk/build/libs/gatk.jar AnalyzeCovariates \
	-before $RBQSRD/before_recalibrated_bqsr_data$i.recal.table \
       	-after $PRD/after_recalibrated_bqsr_data$i.recal.table \
        -csv $PRD/BQSR_$i.csv \
        -plots $PRD/AnalyzeCovariates_bqsr_$i.pdf
done
	##Obtaining an CSV and PDF file of the comparrisson between first and second pass of the recalibration applied to the bam file.






echo ············································································································

echo -e                                       '\n \tHAPLOTYPE CALLER GATK\n'

echo ············································································································






##Ready to call for Variants.
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
	-A FisherStrand -A StrandOddsRatio -A RMSMappingQuality -A MappingQualityRankSumTest -A ReadPosRankSumTest -A DepthPerSampleHC -A BaseQualityRankSumTest -A ExcessHet -A StrandArtifact \
	--annotate-with-num-discovered-alleles=true
	echo -e "\nGATK HaplotypeCallerGVCF ERC GVCF '$i' COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	echo -e "\nGATK HaplotypeCallerGVCF ERC GVCF '$i' COMPLETADO" >> registerFile
done





echo ············································································································

echo -e                                            "\n \tJOINT GENOTYPING (GATK)\n"

echo ············································································································





##If more than one input the pipeline will continue the Joint analysis.

if [ $# -ne 1 ];then
	mkdir combined_gvcf
	echo "mkdir combined_gvcf" >> registerFile
	##Obtenido en GCVF pasamos al Joint Genotyping on one or more samples called with HC.
	echo -e "\nUsing GATK COMBINEGVCFs for merging GVCFs"
	echo -e "\nUsing GATK COMBINEGVCFs for merging GVCFs">> registerFile
	java -jar $SFT/gatk/build/libs/gatk.jar CombineGVCFs \
	-R $HG19/ucsc.hg19.fasta \
	--variant $HCDGVCF/HCdata$1.g.vcf \
	--variant $HCDGVCF/HCdata$2.g.vcf \
	-O $CGVCF/combined.g.vcf
	echo -e "\nGATK COMBINEGVCFs COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	echo -e "\nGATK COMBINEGVCFs COMPLETADO"  >> registerFile



	mkdir genotyped_data_vcf
        echo "mkdir genotyped_data_vcf">> registerFile
       	##GenotypeGVCFs into final VCF
        echo -e "\nUsing GATK GenotypeGVCFs for final VCF"
        java -jar $SFT/gatk/build/libs/gatk.jar GenotypeGVCFs \
        -R $HG19/ucsc.hg19.fasta \
        -V $CGVCF/combined.g.vcf \
        -G StandardAnnotation \
  	-O $GDVCF/genotyped_data.vcf
        echo -e "\nGATK GenotypeGVCFs COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
        echo -e "\nGATK GenotypeGVCFs COMPLETADO" >> registerFile
#if only one INPUT continue to joint analysis.
else
	mkdir genotyped_data_vcf
	echo "mkdir genotyped_data_vcf">> registerFile
	##GenotypeGVCFs into final VCF
	echo -e "\nUsing GATK GenotypeGVCFs for final VCF"
	java -jar $SFT/gatk/build/libs/gatk.jar GenotypeGVCFs \
	-R $HG19/ucsc.hg19.fasta \
	-V $HCDGVCF/HCdata$i.g.vcf \
	-G StandardAnnotation \
	-O $GDVCF/genotyped_data$@.vcf
	echo -e "\nGATK GenotypeGVCFs COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
	echo -e "\nGATK GenotypeGVCFs COMPLETADO" >> registerFile

fi


echo ············································································································

echo -e        "\n \tIf more than 30 samples running VQSR for variants detection using machine learning with the smaples"
echo -e            "\tConsidering using the BAM from 1000G for VCF calling and using them to train the detector"
echo -e	     		"\tand being able to search for OUR VARIANTS in OUR SAMPLE wiht VQSR\n"

echo ············································································································







echo ············································································································

echo -e  "\n \tHard filtering if less than 30 samples and doing it in the classical way, selecting variants\n"

echo ············································································································







echo ············································································································

echo -e               			 "\n \tHard filtering (GATK)"
echo -e                           "\tDATA PREFILTERING  = (vcf_processing_step2.R)\n"

echo ············································································································


##HARD FILTERING
##First step extacting the SNP's
##Second step extracting the INDEL's


mkdir variant_filtration_data_vcf
echo "variantfiltration_data_vcf ">> registerFile

##1.Extract the SNP's from the call set.
echo "Extract the SNP's from the call set."

java -jar $SFT/gatk/build/libs/gatk.jar SelectVariants \
-R $HG19/ucsc.hg19.fasta \
-V $GDVCF/genotyped_data.vcf \
--select-type-to-include SNP \
-O $VFDVCF/selected_raw_snp.vcf
Creates the selected_raw_snps vcf containing just the SNP's from the original  callset.


##2.Apply the filters to the SNP's callset.

java -jar $SFT/gatk/build/libs/gatk.jar VariantFiltration \
-R $HG19/ucsc.hg19.fasta \
-V $VFDVCF/selected_raw_snp.vcf \
--filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
--filter-name "my_SNP_filter" \
-O $VFDVCF/filtered_SNP_data.vcf


##3. Extract the INDELS from the ORIGINAL call set.
java -jar $SFT/gatk/build/libs/gatk.jar SelectVariants \
-R $HG19/ucsc.hg19.fasta \
-V $GDVCF/genotyped_data.vcf \
--select-type-to-include INDEL \
-O $VFDVCF/selected_raw_indels.vcf
Creates the selected_raw_indels vcf containing just the INDEL's from the original  callset.


##4.Apply the filters to the INDEL's callset.
java -jar $SFT/gatk/build/libs/gatk.jar VariantFiltration \
-R $HG19/ucsc.hg19.fasta \
-V $VFDVCF/selected_raw_indels.vcf \
--filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
--filter-name "my_INDEL_filter" \
-O $VFDVCF/filtered_INDEL_data.vcf

##Filtered the INDELS and the SNPS of the original file.

##Combine Variants after using SNPS and INDELS filtering.

java -jar $SFT/gatk/build/libs/gatk.jar MergeVcfs \
-R $HG19/ucsc.hg19.fasta \
-I $VFDVCF/filtered_SNP_data.vcf \
-I $VFDVCF/filtered_INDEL_data.vcf \
-O $VFDVCF/filtered_INDEL_SNP_data.vcf


##GENOTYPEPRIORIS
##GENOTYPE following to do.
##FILTERING THE HETERO OR HOMOZYGOTE VARIANTS
##java -jar $SFT/gatk/build/libs/gatk.jar VariantFiltration \
##-R $HG19/ucsc.hg19.fasta \
##--varaint $VFDVCF/filtered_INDEL_SNP_data.vcf \
##--genotypeFilterExpression "isHet == 1" \
##--genotypeFilterName "GT" \
##-o $VFDVCF/variant_HET_raw.vcf

##Will write  a PASS in homozygotes and a GT in heterozygotes.

##java -jar $SFT/gatk/build/libs/gatk.jar SelectVariants \
##-R $HG19/ucsc.hg19.fasta \
##-V $GDVCF/filtered_INDEL_data.vcf \
##--excludeFiltered \
##--setFilteredGtToNoCall \
##-O $VFDVCF/filtered_SNP_INDEL_HET_excluded.vcf

##Exclude the GT/Heterozygotes in this case.
##Select filtering values, and set filtered GT  with  isHet=1 to nocall ./. 0/1 and 1/0--> ./.  Only 1/1.


echo ············································································································

echo -e                  		 "\n \tVARIANT ANNOTATION (VEP ENSEMBL)\n"

echo ············································································································


###BEFORE ANNOTATION SELECTVARIOANTS MUST BE DONE.

####4.Apply the filters to the INDEL's callset.
##java -jar $SFT/gatk/build/libs/gatk.jar VariantFiltration \
##-R $HG19/ucsc.hg19.fasta \
##-V $VFDVCF/filtered_INDEL_SNP_data.vcf \
##--filter-expression "QUAL > 100 | QD > 2.0 || MQ < 40.0 || FS > 65.0 || ReadPosRankSum < -8.0" \
##--filter-name "my_preprocessing_filter" \
##-O $VFDVCF/filtered_INDEL_SNP_data_BeforeAnnotation.vcf


##Annotate the VCF using VEP(variant effect predictor).
mkdir vep_vcf_annotated
echo " mkdir vep_vcf_annotated " >>registerFile

/home/marius/src/ensembl-vep/./vep -i $VFDVCF/filtered_INDEL_SNP_data.vcf \
--offline \
--cache \
--force_overwrite \
--fork 12 \
--symbol \
--ccds \
--protein \
--uniprot \
--hgvs \
--pubmed \
--af \
--af_1kg \
--af_gnomad \
--af_esp \
--biotype \
--canonical \
--sift b \
--polyphen b \
--regulatory \
--ccds \
--protein \
--fasta /home/marius/.vep/homo_sapiens/93_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz \
--numbers \
--domains \
--uniprot \
--gene_phenotype \
--max_af \
--variant_class \
--filter_common \
--plugin dbNSFP,/home/marius/src/ensembl-vep/dbNSFP_hg19.gz,ALL \
--plugin LoF,/home/marius/src/ensembl-vep/loftee,ALL \
--tab \
--fields "VARIANT_CLASS,Location,Uploaded_variation,Allele,Gene,Feature,Feature_type,BIOTYPE,SYMBOL,Consequence,CANONICAL,PolyPhen,SIFT,gnomAD_AF,ExAC_NFE_AF,EUR_AF,cDNA_position,gnomAD_exomes_AC,gnomAD_exomes_AF,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_AF,PUBMED,SWISSPROT,HGVSc,HGVSp,CLIN_SIG" \
-o $VEPVCFA/vep_anotated_$@.tsv

#--plugin dbNSFP,$SFT/ensembl-vep/dbNSFP_hg19.gz,ALL \
#--plugin LoF \
#--vcf \
#--fields Location,Uploaded_variation,Allele,Gene,Feature,Feature_type,BIOTYPE,SYMBOL,Consequence,CANONICAL,SIFT,PolyPhen,gnomAD_AF,ExAC_NFE_AF,EUR_AF,SWISSPROT,gnomAD_exomes_AC,gnomAD_exomes_AF,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_AF,PUBMED,HGVSc,HGVSp,CLIN_SIG \
#-o $VEPVCFA/vepAnnotatedVCF.vcf
##--coding_only
##--fork 8 threads, assembly GRCh37 --species homo_sapines --check_svs see offline donde
##To apply vep filter do the following unhashtag.
##-o STDOUT | \
##filter_vep -filter CANONICAL is YES and SIFT is deleterios  ##grep -E  rs22323|##"  |cut | head -
##--vcf \
##--fields "Uploaded_variation,Location,Gene,Consequence,SYMBOL,CANONICAL,SWISSPROT,SIFT,PolyPhen,gnomAD_AF,ExAC_NFE_AF,EUR_AF,Allele$
##--fields "Location,Uploaded_variation,Allele,Gene,Feature,Feature_type,SYMBOL,Consequence,CANONICAL,SIFT,PolyPhen,gnomAD_AF,ExAC_NFE$
##SWISSPROT,gnomAD_exomes_AC,gnomAD_exomes_AF,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_AF,PUBMED,HGVSc,HGVSp,CLIN_SIG"$



echo ············································································································

echo -e                        	"\n \tVQSR (more than 30 samples (1000genomes + samples)\n"

echo ············································································································









##VQSR

##Building the SNP recalibration model.
##java -jar $SFT/gatk/build/libs/gatk.jar VariantRecalibrator \
##-R $HG19/ucsc.hg19.fasta \
## -V myvariants.vcf \
## --resource hapmap,known=false,training=true,truth=true,prior=15.0:$HG19/hapmap_3.3.hg19.sites.vcf \
## --resource omni,known=false,training=true,truth=false,prior=12.0:$HG19/1000G_omni2.5.hg19.sites.vcf \
## --resource 1000G,known=false,training=true,truth=false,prior=10.0:$HG19/1000G_phase1.snps.high_confidence.hg19.sites.vcf \
## --resource dbsnp,known=true,training=false,truth=false,prior=2.0:$HG19/dbsnp_138.hg19.vcf \
## -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an StrandArtifact -an ExcessHet \
## -mode SNP \
## --output recalibrated.SNP.recal \
##--tranches-file recal.SNP.tranches \
##--rscript-file recal.SNP.R

##Building the INDEL recalibration model.
##java -jar $SFT/gatk/build/libs/gatk.jar VariantRecalibrator \
##-R $HG19/ucsc.hg19.fasta \
##-V myvariants.vcf \
##--resource mills,known=false,training=true,truth=true,prior=12.0:$HG19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
## --resource dbsnp,known=true,training=false,truth=false,prior=2.0:$HG19/dbsnp_138.hg19.vcf \
## -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an StrandArtifact -an ExcessHet \
## -mode INDEL \
## --output recalibrated.INDEL.recal \
## --tranches-file recal.INDEL.tranches \
## --rscript-file recal.INDEL.R

##ApplyRecalibration step, where the output will be input file + the SNPS anotated with their recalibrated quality scores.
##either with PASS or FILTER depending on wheter or not they are included in the selected tranche.
##the biggher is the tsfilter, less variants leaving outside being considerated is inverse..

##java -jar $SFT/gatk/build/libs/gatk.jar ApplyVQSR \
##-R $HG19/ucsc.hg19.fasta \


##TOBEADDED    VQSR		TOBEADDED	VQSR			TOBEADDED VQSR










#echo ············································································································
#echo -e                                  	"\n \tVARIANTS TO TABLE (GATK)"
#echo -e 	                         "\tPrefiltering = vcf_processing_step1_v2.py\n"
#echo ············································································································


##USING GATK, SELECT DATA FR0M THE VCF VILE TO SEPARATE IN DIFFERENT COLUMNS.
##-F VCF DATA FROM (INFO)
##-GF VCF DATA FROM (FORMAT)

#mkdir variantstotable_vcf
#echo "mkdir variantstotable_vcf">> registerFile
#echo -e "\nUsing GATK VariantsToTable for final VCF sep fileds"

##VariantsToTable avoinding extra scripts.

#java -jar $SFT/gatk/build/libs/gatk.jar VariantsToTable \
#-V $VEPVCFA/vepAnnotatedVCF$@.vcf \
#-F CHROM -F Uploaded_variation -F Location -F POS -F TYPE -F Allele -F REF -F ALT -F Gene -F Feature -F Feature_type -F TYPE -F Consequence -F SYMBOL -F CANONICAL -F SWISSPROT -F SIFT -F PolyPhen -F gnomAD_AF -F ExAC_NFE_AF -F EUR_AF -F gnomAD_exomes_AC -F gnomAD_exomes_AF -F gnomAD_exomes_NFE_AC -F gnomAD_exomes_NFE_AF -F PUBMED -F HGVSc -F HGVSp -F CLIN_SIG -F QUAL -F QD -F FS -F SOR -F MQ -F RMSMappingQuality -F MappingQualityRankSumTest -F ReadPosRankSumTest -F ExcessHet -F MLEAF -F MLEAC -F HET -F HOM-REF -F HOM-VAR -GF GT -GF AD -GF DP -GF GQ -GF PL \
#--output $VTTVCF/variants2Table$@.table.tsv





##--error-if-missing-data=true \
##data NA if missing.
##VCF file or change it to  .txt
##select all the FIELDS  wanted to be separated.

#echo -e "\nGATK VariantsToTable COMPLETADO" ;paplay /usr/share/sounds/freedesktop/stereo/complete.oga
#echo -e "\nGATK VariantsToTable COMPLETADO" >> registerFile


##info

#FILTERING the TSV in order to accurately get a readable TSV file.





echo ·····················································································

echo -e "\n \tAnalisis COMPLETADO para:" $@ ", el dia, "$(date +%Y/%m/%d)
echo -e "\n \tAnalisis COMPLETADO para:" $@ ", el dia, "$(date +%Y/%m/%d) > READin.txt
echo -e "\n \tAnalisis COMPLETADO para:" $@ ", el dia, "$(date +%Y/%m/%d) >> registerFile
echo -e "\n \tDONE"
