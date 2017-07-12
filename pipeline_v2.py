#v2 includes option for selecting Amplicon or Captured approach 

import sys
from glob import glob
from subprocess import call 

analysis = raw_input("Choose Amplicon(A) or Captured(C):")
runner = raw_input("Local(L) or Server(S):")
#Importing samples
samples_path = sys.argv[1] + '*'
samples_path = glob(samples_path)

for sample_path in samples_path:

	#Obtaining sample_name
	sample_name = sample_path[sample_path.rfind('/')+1:]
	


	#Obtaining fastq files paths
	#fastq1_path = sample_path + '/' + sample_name + ('_1.fastq.gz')
	#fastq2_path = sample_path + '/' + sample_name + ('_2.fastq.gz')

	###################################################################################

	fastq1_path = glob(sample_path + '/*R1_merged.fastq.gz')
	fastq2_path = glob(sample_path + '/*R2_merged.fastq.gz')
	
	if fastq1_path == []:
		continue
	fastq1_path = fastq1_path[0]
	fastq2_path = fastq2_path[0]

	###################################################################################
	#Checking if fastq files are gunzipped or not
	test = glob(fastq1_path)
	if test == []:
		fastq1_path = sample_path + '/' + sample_name + ('_1.fastq')
		fastq2_path = sample_path + '/' + sample_name + ('_2.fastq')


	##### FASTQ QUALITY CONTROL #####
	call('fastqc ' + fastq1_path + ' ' + fastq2_path + ' -t 5 --noextract -o ' + sample_path + '/',shell = True)

	call('mkdir ' + sample_path + '/working_temp',shell = True) #Creating temporary working folder needed in following steps



	if runner == "L":
			
		##### MAPPING TO REF GENOME #####
		call('bwa mem -t6 /home/daguilera/Documents/hg19/ucsc.hg19.fasta ' + fastq1_path + ' ' + fastq2_path + ' > ' + sample_path + '/' + sample_name + '_bwa.sam',shell = True)
		call('java -XX:+UseParallelGC -XX:ParallelGCThreads=6 -Xmx4g -jar /mnt/datos1/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar AddOrReplaceReadGroups I= ' + sample_path + '/' + sample_name + '_bwa.sam O= ' + sample_path + '/' + sample_name + '_bwa_rg_sorted.bam SO=coordinate RGID=' + sample_name + ' RGLB=library RGPL=illumina RGPU=library RGSM=' + sample_name + ' VALIDATION_STRINGENCY=SILENT TMP_DIR= '+ sample_path + '/working_temp',shell = True)
		call('java -XX:+UseParallelGC -XX:ParallelGCThreads=6 -Xmx4g -jar /mnt/datos1/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar CollectAlignmentSummaryMetrics R=/home/daguilera/Documents/hg19/ucsc.hg19.fasta I= ' + sample_path + '/' + sample_name + '_bwa_rg_sorted.bam O= ' + sample_path + '/' + sample_name + '_alignment_metrics.txt',shell = True)

		if analysis == "A":
			call('samtools index ' + sample_path + '/' + sample_name + '_bwa_rg_sorted.bam ' + sample_path + '/' + sample_name + '_bwa_rg_sorted.bai',shell = True)
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=6 -Xmx4g -jar /mnt/datos1/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T IndelRealigner -R /home/daguilera/Documents/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted.bam -o ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned.bam -targetIntervals /home/daguilera/Documents/hg19/hg19_indels_output.intervals -known /home/daguilera/Documents/hg19/1000G_phase1.indels.hg19.sites.vcf -known /home/daguilera/Documents/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -LOD 0.4 --consensusDeterminationModel KNOWNS_ONLY',shell=True)
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=6 -Xmx4g -jar /mnt/datos1/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T BaseRecalibrator -R /home/daguilera/Documents/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned.bam -knownSites /home/daguilera/Documents/hg19/1000G_phase1.indels.hg19.sites.vcf -knownSites /home/daguilera/Documents/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites /home/daguilera/Documents/hg19/dbsnp_138.hg19.vcf -o ' + sample_path + '/' + sample_name + '_recal_data.table',shell =True)
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=6 -Xmx4g -jar /mnt/datos1/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T PrintReads -R /home/daguilera/Documents/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned.bam -BQSR ' + sample_path + '/' + sample_name + '_recal_data.table -o ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned_bqsr.bam',shell = True)
			call('bedtools coverage -a ' + sys.argv[1] + '/*_targets.bed -b ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned.bam -mean| sort -k1,1 -k2,2n > ' + sample_path + '/' + sample_name + '_coverage.txt',shell = True)
			call('java -jar /mnt/datos1/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar CollectTargetedPcrMetrics I= ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned.bam O= '+ sample_path + '/' + sample_name + '_pcr_metrics.txt R=/home/daguilera/Documents/hg19/ucsc.hg19.fasta AMPLICON_INTERVALS= ' + sys.argv[1] + '/*_amplicons.interval_list TARGET_INTERVALS= '+ sys.argv[1] + '/*_targets.interval_list',shell = True)
			
			##### VARIANT CALLING #####
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=6 -Xmx4g -jar /mnt/datos1/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R /home/daguilera/Documents/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_notdedupped_realigned_bqsr.bam -stand_emit_conf 10 -stand_call_conf 30 -o ' + sample_path + '/' + sample_name + '_raw_unformatted.vcf',shell = True)

		if analysis == "C":
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=6 -Xmx4g -jar /mnt/datos1/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar MarkDuplicates I= ' + sample_path + '/' + sample_name + '_bwa_rg_sorted.bam O= ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT TMP_DIR= ' + sample_path + '/working_temp M= ' + sample_path + '/' + sample_name + '_duplicates.metrics',shell = True)
			call('rm -r ' + sample_path + '/working_temp',shell = True) #Removing the temporary working folder (not needed anymore)
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=6 -Xmx4g -jar /mnt/datos1/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T IndelRealigner -R /home/daguilera/Documents/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped.bam -o ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned.bam -targetIntervals /home/daguilera/Documents/hg19/hg19_indels_output.intervals -known /home/daguilera/Documents/hg19/1000G_phase1.indels.hg19.sites.vcf -known /home/daguilera/Documents/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -LOD 0.4 --consensusDeterminationModel KNOWNS_ONLY',shell=True)
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=6 -Xmx4g -jar /mnt/datos1/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T BaseRecalibrator -R /home/daguilera/Documents/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned.bam -knownSites /home/daguilera/Documents/hg19/1000G_phase1.indels.hg19.sites.vcf -knownSites /home/daguilera/Documents/hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites /home/daguilera/Documents/hg19/dbsnp_138.hg19.vcf -o ' + sample_path + '/' + sample_name + '_recal_data.table',shell =True)
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=6 -Xmx4g -jar /mnt/datos1/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T PrintReads -R /home/daguilera/Documents/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned.bam -BQSR ' + sample_path + '/' + sample_name + '_recal_data.table -o ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned_bqsr.bam',shell = True)
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
		'''	
			#call('java -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx4g -jar /home/domin/genetica/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar CollectHsMetrics BI= '+ sys.argv[1] + '/*_targets.interval_list TI= '+ sys.argv[1] + '/*_targets.interval_list I= ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned_bqsr.bam O= ' + sample_path + '/' + sample_name + '_HS_metrics.txt R=/home/domin/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta',shell = True)
			##### VARIANT CALLING #####
			call('java -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Xmx4g -jar /home/domin/genetica/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R /home/domin/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta -I ' + sample_path + '/' + sample_name + '_bwa_rg_sorted_dedupped_realigned_bqsr.bam -stand_emit_conf 10 -stand_call_conf 30 -o ' + sample_path + '/' + sample_name + '_raw_unformatted.vcf',shell = True)
	
	##### VARIANT ANNOTATION #####

		
		call('/mnt/genetica/domin/GeneticaPipeDB/software/annovar/table_annovar.pl ' + sample_path + '/' + sample_name + '_raw_unformatted.vcf /mnt/genetica/domin/GeneticaPipeDB/software/annovar/humandb -buildver hg19 -out ' + sample_path + '/' + sample_name + ' --remove --otherinfo --protocol refGene,dbnsfp31a_interpro,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_eur,exac03,hrcr1,kaviar_20150923,popfreq_max_20150413,avsnp147,mcap,revel,cadd13,caddindel,dann,eigen,gwava,spidex,dbscsnv11,dbnsfp30a,clinvar_20160302,phastConsElements46way,tfbsConsSites,wgRna,targetScanS,gwasCatalog --operation g,f,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r --nastring . --vcfinput --thread 6',shell = True)
		call('python /mnt/genetica/domin/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step1.py ' + sample_path + '/' + sample_name + '.hg19_multianno.txt' ,shell = True)
		call('Rscript /mnt/genetica/domin/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step2.R ' + sample_path + '/' + sample_name + '_annotated_formatted.txt',shell = True)
		'''
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
