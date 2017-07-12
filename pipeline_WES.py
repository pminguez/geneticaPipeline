
import sys
from glob import glob
from subprocess import call 

analysis = raw_input("Choose Amplicon(A) or Captured(C):")
runner = raw_input("Local(L) or Server(S):")

#Importing samples
forward_paths = sorted(glob(sys.argv[1] + '*R1.fastq.gz'))
reverse_paths = sorted(glob(sys.argv[1] + '*R2.fastq.gz'))


#Importing bedfile
#bedfile = sys.argv[2]

#Setting paths to work in local or in the server
if runner == 'L':
	genome_ref = "/home/daguilera/Documents/genome_data/hg19/ucsc.hg19.fasta"
	picardtools = "/mnt/datos1/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar"
	gatk = '/mnt/datos1/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar'
	hg19_path = "/home/daguilera/Documents/genome_data/hg19/"
	annovar = "/mnt/datos1/GeneticaPipeDB/software/annovar/table_annovar.pl"
	annovarDB = "/mnt/datos1/GeneticaPipeDB/software/annovar/humandb"
	genome_fai = '/mnt/datos2/d.aguilera/CNVs_analysis/CoNVaDING-1.1.6/ucsc.hg19_convading.fasta.fai'

elif runner == 'S':
	genome_ref = "/home/domin/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta"
	picardtools = "/home/domin/genetica/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar"
	gatk = '/home/domin/genetica/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar'
	hg19_path = "/home/domin/genetica/GeneticaPipeDB/genome_data/hg19/"
	annovar = "/home/domin/genetica/GeneticaPipeDB/software/annovar/table_annovar.pl"
	annovarDB = "/home/domin/genetica/GeneticaPipeDB/software/annovar/humandb"
	genome_fai = '/home/domin/genetica2/d.aguilera/CNVs_analysis/CoNVaDING-1.1.6/ucsc.hg19_convading.fasta.fai'

#Aligning with BWA, sorting and indexing with a pipe process
for i in range(0,len(forward_paths)):
	sample_path = forward_paths[i][:forward_paths[i].rfind('/')+1]
	sample_name = forward_paths[i][forward_paths[i].rfind('/')+1:forward_paths[i].rfind('_R')]
	

	#call('bwa mem -t8 -R "@RG\tID:sample\tLB:library\tPL:illumina\tPU:library\tSM:sample" ' + genome_ref + ' ' + forward_paths[i] + ' ' + reverse_paths[i] + ' > ' + sample_path + '/' +  sample_name + '_bwa.sam',shell =True)

	#call('samtools view -Shu ' + sample_path + '/' + sample_name +  '_bwa.sam | \
	#	samtools sort -@ 8 -m 3G - ' + sample_path + '/' + sample_name + '_sorted && samtools index ' + sample_path + '/' + sample_name + '_sorted.bam',shell =True)

	#Setting input name for indelrealigner just for pcr-based data

	indelrealigner_input = '*_sorted.bam'

#Marking Duplicates by parallelizing samples (only for non pcr-based data)
if analysis == 'C' or analysis == 'c':
	call("find " + sample_path + "*_sorted.bam | parallel -j3 'java -Xmx9g -jar " + picardtools + " \
		MarkDuplicates \
		I= {} \
		O= {}_dedupped.bam \
		CREATE_INDEX=true \
		VALIDATION_STRINGENCY=SILENT \
		TMP_DIR= " + sample_path + "working_temp \
		M= {}_duplicate_metrics.txt'",shell = True)

		#Setting input name for indelrealigner just for non pcr-based data
	indelrealigner_input = '*_dedupped.bam'

#Indelrealigning by parallelizing samples 
call("find " + sample_path +  indelrealigner_input + " | parallel -j3 'java -Xmx9g -jar " + gatk + " \
	-T IndelRealigner \
	-R " + genome_ref + " \
	-I {} \
	-o {}_indelrealigned.bam \
	-targetIntervals " + hg19_path + "hg19_indels_output.intervals \
	-known " + hg19_path + "1000G_phase1.indels.hg19.sites.vcf \
	-known " + hg19_path + "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
	-LOD 0.4 \
	--consensusDeterminationModel KNOWNS_ONLY'",shell=True)
	#--fix_misencoded_quality_scores \

#BQSR (step1) by parallelizing samples 
call("find " + sample_path + "*_indelrealigned.bam | parallel -j3 'java -Xmx9g -jar " + gatk + " \
	-T BaseRecalibrator \
	-R " + genome_ref + " \
	-I {} \
	-knownSites " + hg19_path + "1000G_phase1.indels.hg19.sites.vcf \
	-knownSites " + hg19_path + "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
	-knownSites " + hg19_path + "dbsnp_138.hg19.vcf \
	-o {}_recal_data.table'",shell=True)

#BQSR (step2) by parallelizing samples  
call("find " + sample_path + "*_indelrealigned.bam | parallel -j3 'java -Xmx9g -jar " + gatk + " \
	-T PrintReads \
	-R " + genome_ref + " \
	-I {} \
	-BQSR {}_recal_data.table \
	-o {}_bqsr.bam'",shell=True)

#Variant calling by parallelizing samples
call("find " + sample_path + "*_bqsr.bam | parallel -j3 'java -Xmx9g -jar " + gatk + " \
	-T HaplotypeCaller \
	-R " + genome_ref + " \
	-I {} \
	--emitRefConfidence GVCF \
	--variant_index_type LINEAR \
 	--variant_index_parameter 128000 \
	-o {}.g.vcf'",shell=True)

#Genotyping in single-mode
call("find " + sample_path + "*.g.vcf | parallel -j1 'java -Xmx28g -jar " + gatk + " \
	-T GenotypeGVCFs \
	-nt 8 \
	-R " + genome_ref + " \
	-V {} \
	-o {}_singleGT_raw.vcf'",shell = True) 

#Annovar variant annotation
variants = sorted(glob(sys.argv[1] + '*raw.vcf'))
for vcffile in variants:
	sample_path = vcffile[:vcffile.rfind('/')+1]
	sample_name = vcffile[vcffile.rfind('/')+1:vcffile.rfind('_sorted.bam')]
	output = sample_path + sample_name

	call(annovar + ' ' + vcffile + ' ' + annovarDB + ' -buildver hg19 \
		-out ' + output + ' \
		--remove \
		--otherinfo \
		--protocol refGene,dbnsfp31a_interpro,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_eur,exac03,hrcr1,kaviar_20150923,popfreq_max_20150413,avsnp147,mcap,revel,cadd13,caddindel,dann,eigen,gwava,spidex,dbscsnv11,dbnsfp30a,clinvar_20160302,phastConsElements46way,tfbsConsSites,wgRna,targetScanS,gwasCatalog \
		--operation g,f,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r \
		--nastring . \
		--vcfinput \
		--thread 8',shell = True)

	#Converting from vcf to table
	call('python /mnt/genetica/domin/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step1.py ' + sample_path + '/' + sample_name + '.hg19_multianno.txt' ,shell = True)
	call('Rscript /mnt/genetica/domin/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step2.R ' + sample_path + '/' + sample_name + '_annotated_formatted.txt',shell = True)







'''


call("find " + sample_path + "*_bqsr.bam | parallel -j4 'bedtools coverage \
	-a " + bedfile + " \
	-b {} \
	-sorted \
	-g " + genome_fai + " \
	-mean| sort -k1,1 -k2,2n \
	> {}_regions_mean_coverage.txt",shell = True)


call("find " + sample_path + "*_bqsr.bam | parallel -j4 'java -Xmx6g -jar " + picardtools + " \
	CollectHsMetrics \
	BI= " + interval_list + " \
	TI= " + interval_list + " \
	I= {} \
	O= {}_HSmetrics.txt \
	R= " + genome_ref,shell = True)


call("find " + sample_path + "*_rg_sorted.bam | parallel 'java -jar " + picardtools + " \
	CollectAlignmentSummaryMetrics \
	R= " + genome_ref + " \
	I= {} \
	O= {}_alignment_metrics.txt '",shell = True) 


#Converting bedfile to intervalist for using picard
call('java -jar ' + picardtools + '\
	BedToIntervalList \
	I= ' + bedfile + ' \
	O= ' + bedfile + '_interval_list.txt \
	SD= ' + hg19_path + 'ucsc.hg19.dict',shell = True)

interval_list = bedfile + '_interval_list.txt'



##### FASTQ QUALITY CONTROL #####
call('fastqc ' + forward_paths[i] + ' ' + reverse_paths[i] + ' -t 8 --noextract -o ' + sample_path + '/',shell = True)
call('mkdir ' + sample_path + '/working_temp',shell = True) #Creating temporary working folder needed in following step

call("find " + sample_path + "*_bqsr.bam | parallel 'java -jar " + gatk + " \
	-T HaplotypeCaller \
	-R " + genome_ref + " \
	-I {} \
	-stand_emit_conf 10 \
	-stand_call_conf 30 \
	-o {}_raw.vcf'",shell=True)


'''



'''

'''


'''



			call('bedtools coverage -a ' + sys.argv[1] + '/*_targets.bed -b ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned.bam -mean| sort -k1,1 -k2,2n > ' + sample_path + '/' + sample_name + '_coverage.txt',shell = True)
			call('java -jar /mnt/datos1/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar CollectTargetedPcrMetrics I= ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned.bam O= '+ sample_path + '/' + sample_name + '_pcr_metrics.txt R=/home/daguilera/Documents/hg19/ucsc.hg19.fasta AMPLICON_INTERVALS= ' + sys.argv[1] + '/*_amplicons.interval_list TARGET_INTERVALS= '+ sys.argv[1] + '/*_targets.interval_list',shell = True)
			

		if analysis == "C":
			
			call('rm -r ' + sample_path + '/working_temp',shell = True) #Removing the temporary working folder (not needed anymore)
			call('bedtools coverage -a ' + sys.argv[1] + '/*_targets.bed -b ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned_bqsr.bam -mean| sort -k1,1 -k2,2n > ' + sample_path + '/' + sample_name + '_coverage.txt',shell = True)
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=6 -Xmx4g -jar /mnt/datos1/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar CollectHsMetrics BI= '+ sys.argv[1] + '/*_targets.interval_list TI= '+ sys.argv[1] + '/*_targets.interval_list I= ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned_bqsr.bam O= ' + sample_path + '/' + sample_name + '_coverage.txt R=/home/daguilera/Documents/hg19/ucsc.hg19.fasta',shell = True)
			##### VARIANT CALLING #####
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=6 -Xmx4g -jar /mnt/datos1/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R /home/daguilera/Documents/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned_bqsr.bam -stand_emit_conf 10 -stand_call_conf 30 -o ' + sample_path + '/' + sample_name + '_raw_unformatted.vcf',shell = True)
		
		##### VARIANT ANNOTATION #####
		call('/mnt/datos1/GeneticaPipeDB/software/annovar/table_annovar.pl ' + sample_path + '/' + sample_name + '_raw_unformatted.vcf /mnt/datos1/GeneticaPipeDB/software/annovar/humandb -buildver hg19 -out ' + sample_path + '/' + sample_name + ' --remove --otherinfo --protocol refGene,dbnsfp31a_interpro,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_eur,exac03,hrcr1,kaviar_20150923,popfreq_max_20150413,avsnp147,cadd,dann,eigen,gwava,dbscsnv11,dbnsfp30a,clinvar_20160302,phastConsElements46way,tfbsConsSites,wgRna,targetScanS,gwasCatalog --operation g,f,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r --nastring . --vcfinput',shell = True)
		call('python /mnt/datos1/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step1.py ' + sample_path + '/' + sample_name + '.hg19_multianno.txt' ,shell = True)
		call('Rscript /mnt/datos1/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step2.R ' + sample_path + '/' + sample_name + '_annotated_formatted.txt',shell = True)
		
	if runner == "S":

		##### MAPPING TO REF GENOME #####
		call('bwa mem -t4 /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta ' + fastq1_path + ' ' + fastq2_path + ' > ' + sample_path + '/' + sample_name + '_bwa.sam',shell = True)
		call('java -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx4g -jar /home/domin/genetica/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar AddOrReplaceReadGroups I= ' + sample_path + '/' + sample_name + '_bwa.sam O= ' + sample_path + '/' + sample_name + '_bwa_rg_sorted.bam SO=coordinate RGID=' + sample_name + ' RGLB=library RGPL=illumina RGPU=library RGSM=' + sample_name + ' VALIDATION_STRINGENCY=SILENT TMP_DIR= '+ sample_path + '/working_temp',shell = True)
		call('java -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx4g -jar /home/domin/genetica/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar CollectAlignmentSummaryMetrics R=/mnt/genetica/domin/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta I= ' + sample_path + '/' + sample_name + '_bwa_rg_sorted.bam O= ' + sample_path + '/' + sample_name + '_alignment_metrics.txt',shell = True)
		if analysis == "A":
			call('samtools index ' + sample_path + '/' + sample_name + '_bwa_rg_sorted.bam ' + sample_path + '/' + sample_name + '_bwa_rg_sorted.bai',shell = True)
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx4g -jar /home/domin/genetica/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T IndelRealigner -R /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted.bam -o ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned.bam -targetIntervals /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/hg19_indels_output.intervals -known /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/1000G_phase1.indels.hg19.sites.vcf -known /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -LOD 0.4 --consensusDeterminationModel KNOWNS_ONLY',shell=True)
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx4g -jar /home/domin/genetica/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T BaseRecalibrator -R /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned.bam -knownSites /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/1000G_phase1.indels.hg19.sites.vcf -knownSites /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/dbsnp_138.hg19.vcf -o ' + sample_path + '/' + sample_name + '_recal_data.table',shell =True)
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx4g -jar /home/domin/genetica/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T PrintReads -R /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned.bam -BQSR ' + sample_path + '/' + sample_name + '_recal_data.table -o ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned_bqsr.bam',shell = True)
			call('bedtools coverage -a ' + sys.argv[1] + '/*_targets.bed -b ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned.bam -mean| sort -k1,1 -k2,2n > ' + sample_path + '/' + sample_name + '_coverage.txt',shell = True)
			call('java -jar /home/domin/genetica/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar CollectTargetedPcrMetrics I= ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned.bam O= '+ sample_path + '/' + sample_name + '_pcr_metrics.txt R=/home/domin/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta AMPLICON_INTERVALS= ' + sys.argv[1] + '/*_amplicons.interval_list TARGET_INTERVALS= '+ sys.argv[1] + '/*_targets.interval_list',shell = True)
			
			##### VARIANT CALLING #####
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx4g -jar /home/domin/genetica/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned_bqsr.bam -stand_emit_conf 10 -stand_call_conf 30 -o ' + sample_path + '/' + sample_name + '_raw_unformatted.vcf',shell = True)

		if analysis == "C":
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx4g -jar /home/domin/genetica/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar MarkDuplicates I= ' + sample_path + '/' + sample_name + '_bwa_rg_sorted.bam O= ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR= ' + sample_path + '/working_temp M= ' + sample_path + '/' + sample_name + '_duplicates.metrics',shell = True)
			call('rm -r ' + sample_path + '/working_temp',shell = True) #Removing the temporary working folder (not needed anymore)
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx4g -jar /home/domin/genetica/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T IndelRealigner -R /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped.bam -o ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned.bam -targetIntervals /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/hg19_indels_output.intervals -known /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/1000G_phase1.indels.hg19.sites.vcf -known /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -LOD 0.4 --consensusDeterminationModel KNOWNS_ONLY',shell=True)
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx4g -jar /home/domin/genetica/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T BaseRecalibrator -R /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned.bam -knownSites /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/1000G_phase1.indels.hg19.sites.vcf -knownSites /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/dbsnp_138.hg19.vcf -o ' + sample_path + '/' + sample_name + '_recal_data.table',shell =True)
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx4g -jar /home/domin/genetica/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T PrintReads -R /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned.bam -BQSR ' + sample_path + '/' + sample_name + '_recal_data.table -o ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned_bqsr.bam',shell = True)
		
			#call('java -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx4g -jar /home/domin/genetica/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar CollectHsMetrics BI= '+ sys.argv[1] + '/*_targets.interval_list TI= '+ sys.argv[1] + '/*_targets.interval_list I= ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned_bqsr.bam O= ' + sample_path + '/' + sample_name + '_HS_metrics.txt R=/home/domin/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta',shell = True)
			##### VARIANT CALLING #####
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx4g -jar /home/domin/genetica/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned_bqsr.bam -stand_emit_conf 10 -stand_call_conf 30 -o ' + sample_path + '/' + sample_name + '_raw_unformatted.vcf',shell = True)
	
	##### VARIANT ANNOTATION #####

		
		call('/mnt/genetica/domin/GeneticaPipeDB/software/annovar/table_annovar.pl ' + sample_path + '/' + sample_name + '_raw_unformatted.vcf /mnt/genetica/domin/GeneticaPipeDB/software/annovar/humandb -buildver hg19 -out ' + sample_path + '/' + sample_name + ' --remove --otherinfo --protocol refGene,dbnsfp31a_interpro,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_eur,exac03,hrcr1,kaviar_20150923,popfreq_max_20150413,avsnp147,mcap,revel,cadd13,caddindel,dann,eigen,gwava,spidex,dbscsnv11,dbnsfp30a,clinvar_20160302,phastConsElements46way,tfbsConsSites,wgRna,targetScanS,gwasCatalog --operation g,f,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r --nastring . --vcfinput --thread 6',shell = True)
		
		
	##### SORTING FILES #####
	call('rm ' + sample_path + '/' + sample_name + '.avinput',shell = True) #Removing file generated by annovar for input

	call('mkdir ' + sample_path + '/fastq_files',shell = True) #Creating folder for storing fastq files
	call('mv ' + sample_path + '/*fastq* ' + sample_path + '/fastq_files',shell = True) #Storing fastq files

	call('mkdir ' + sample_path + '/bam_files',shell = True) #Creating folder for storing bam/sam related files
	call('mv ' + sample_path + '/*.ba* ' + sample_path + '/bam_files',shell = True) #Storing bam and bai files
	call('mv ' + sample_path + '/*.sam ' + sample_path + '/bam_files',shell = True) #Storing sam files
	call('mv ' + sample_path + '/*.table ' + sample_path + '/bam_files',shell = True) #Storing bqsr table file
	call('mv ' + sample_path + '/*.sai ' + sample_path + '/bam_files',shell = True) # Storing sai files 
	call('mv ' + sample_path + '/*.metrics ' + sample_path + '/bam_files',shell = True) #Storing duplicate metrics

	call('mkdir ' + sample_path + '/variant_files',shell = True) #Creating folder for storing variant related files
	call('mv ' + sample_path + '/*.vcf* ' + sample_path + '/variant_files',shell = True) #Storing vcf files
'''