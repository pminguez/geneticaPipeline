#############################################################################
## Pipeline is now on github: https://github.com/pminguez/geneticaPipeline ##
#############################################################################

import sys
from glob import glob
from subprocess import call 
import argparse
import time
import datetime
import os
import re

def countdown(t):
    for t in range(t,-1,-1):
        mins, secs = divmod(t, 60)
        timeformat = '{:02d}'.format(secs)
        sys.stdout.write('\rRunning in ' + timeformat + ' secs')
        sys.stdout.flush()
        time.sleep(1)

parser = argparse.ArgumentParser(description="Process Fastq files for getting variants")

parser.add_argument("-u", action="store", dest='user',
					help="user name to look/store in the correct dir")

parser.add_argument("-I", action="store", dest='input',
                    help="path to input folder")

parser.add_argument("-o", action="store", dest='output',
					help="path to output results")

parser.add_argument("-T", action="store", dest='threads', type = int, default = 16,
                    help="specify number of threads to use")

parser.add_argument("-J", action="store", dest='parallelization', type = int, default = 5,
                    help="specify number of samples to run in parallel")

parser.add_argument("-duplicates", action="store_true",
                    help="set this flag to markduplicates with picardtools")

parser.add_argument("-local", action="store_true",
                    help="set this flag to run the pipeline using local paths")

parser.add_argument("-gvcf", action="store_true",
					help="set this flag to keep gvcf files")



args = parser.parse_args()


if args.input == None:
	print ''
	print 'ERROR: An input folder containing fastq files is needed'
	print ''
	parser.print_help()
	exit()

if args.output == None:
	print ''
	print 'ERROR: An output folder containing fastq files is needed'
	print ''

	parser.print_help()
	exit()

#Importing samples
forward_paths = sorted(glob(args.input + '*_R1.fastq.gz'))
reverse_paths = sorted(glob(args.input + '*_R2.fastq.gz'))


if forward_paths == []:
	forward_paths = sorted(glob(args.input + '*_1.fastq.gz'))
	reverse_paths = sorted(glob(args.input + '*_2.fastq.gz'))
	

if forward_paths == []:
	print ''
	print 'ERROR: No fastq files detected in ' + args.input + '.\nFastq files names should be named: name_R1.fastq.gz and name_R2.fastq.gz or name_1.fastq.gz and name_2.fastq.gz'
	print ''
	exit()

if len(forward_paths) != len(reverse_paths):
	print ''
	print 'ERROR: Different number of forward and reverse fastq files detected. PLEASE CHECK.'
	print ''
	exit()


print '---------------------------------------------------------------------------------------------'
print '                            *****Running FJD Pipeline*****                             '
print '---------------------------------------------------------------------------------------------'
print ''
print 'Number of samples to analyze: ' + str(len(forward_paths))
print ''
print 'ARGUMENTS:'
print ''
print ' -User: ' + str(args.user)
print '	-Input: ' + str(args.input)
print '	-Output: ' + str(args.output)
print '	-Threads: ' + str(args.threads)
print '	-Sample to parallelizate: ' + str(args.parallelization)
print '	-MarkDuplicates: ' + str(args.duplicates)
print '	-Running local: ' + str(args.local)
print '	-Keep gVCF files: ' + str(args.gvcf)
print ''
print '---------------------------------------------------------------------------------------------'
print 'Please review the arguments and number of samples to process...'
print ''
countdown(5)
print ''
print '---------------------------------------------------------------------------------------------'


if args.local:
	genome_ref = "/mnt/genetica/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta"
	picardtools = "/mnt/genetica/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar"
	gatk = "/mnt/genetica/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar"
	hg19_path = "/mnt/genetica/GeneticaPipeDB/genome_data/hg19/"
	annovar = "/mnt/genetica/GeneticaPipeDB/software/annovar/table_annovar.pl"
	annovarDB = "/mnt/genetica/GeneticaPipeDB/software/annovar/humandb"
	output_path = args.output

else:
	genome_ref = "/mnt/genetica/" + str(args.user) + "/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta"
	picardtools = "/mnt/genetica3/GeneticaPipeDB_updated/picard/build/libs/picard.jar"
	gatk = "/mnt/genetica3/GeneticaPipeDB_updated/gatk-4.0.5.1/gatk-package-4.0.5.1-local.jar"
	hg19_path = "/mnt/genetica/" + str(args.user) + "/GeneticaPipeDB/genome_data/hg19/"
	annovar = "/mnt/genetica3/GeneticaPipeDB_updated/annovar/table_annovar.pl"
	annovarDB = "/mnt/genetica3/GeneticaPipeDB_updated/annovar/humandb"
	output_path = args.output



print '                               Mapping fastq files (BWA)                                      '
print '----------------------------------------------------------------------------------------------'

#Load genome reference to memory

call('/mnt/genetica3/GeneticaPipeDB_updated/bwa/bwa shm ' + genome_ref,shell = True)


#Loop samples for BWA
for i in range(0,len(forward_paths)):
	sample_path = forward_paths[i][:forward_paths[i].rfind('/')+1]
	sample_name = forward_paths[i][forward_paths[i].rfind('/')+1:forward_paths[i].rfind('_R')]

	call('/mnt/genetica3/GeneticaPipeDB_updated/bwa/bwa mem -t' + str(args.threads) + ' -R "@RG\\tID:' + sample_name + '\\tLB:library\\tPL:illumina\\tPU:library\\tSM:' + sample_name + '" ' + genome_ref + ' ' + forward_paths[i] + ' ' + reverse_paths[i] + ' > ' + output_path + '/' + sample_name + '_bwa.sam',shell = True)

	
#Unload genome reference	

call('/mnt/genetica3/GeneticaPipeDB_updated/bwa/bwa shm -d',shell = True)
print '----------------------------------------------------------------------------------------------'

print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Sorting and creating bam and bai files of samples...'
print '----------------------------------------------------------------------------------------------'
call("find " + output_path +  "*.sam | parallel --no-notice -j" + str(args.parallelization) + " '/mnt/genetica3/GeneticaPipeDB_updated/samtools-1.8/samtools sort {} -O BAM -@ " + str(args.threads / 2) + " -o {}_sorted.bam && /mnt/genetica3/GeneticaPipeDB_updated/samtools-1.8/samtools index {}_sorted.bam'", shell = True)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] OK! '
print '----------------------------------------------------------------------------------------------'
baserecalibrator_input = '*_sorted.bam'


#Remove sam files
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing sam files...'
print '----------------------------------------------------------------------------------------------'
for i in glob(output_path + '*.sam'):
	os.remove(i)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing sam files...OK'
print '----------------------------------------------------------------------------------------------'


#MarkDuplicates with picardtools
if args.duplicates:
	print '----------------------------------------------------------------------------------------------'
	print '[FJD_Pipeline] Marking Duplicates... '
	print '----------------------------------------------------------------------------------------------'
	call("find " + output_path + "*_sorted.bam | parallel --no-notice -j" + str(args.parallelization) + " 'java -Xmx9g -jar " + picardtools + " \
		MarkDuplicates \
		I= {} \
		O= {}_dedupped.bam \
		CREATE_INDEX=true \
		VALIDATION_STRINGENCY=SILENT \
		TMP_DIR= " + sample_path + "working_temp \
		M= {}_duplicate_metrics.txt'",shell = True)

	print '----------------------------------------------------------------------------------------------'
	print '[FJD_Pipeline] Marking Duplicates...OK! '
	print '----------------------------------------------------------------------------------------------'

	baserecalibrator_input = '*_dedupped.bam'
	

#Empieza GATK

#Remove intermediary files
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing intermediary files...'
print '----------------------------------------------------------------------------------------------'
intermediary_files = glob(output_path + '*_sorted.bam') + glob(output_path + '*_sorted.bam.bai')
for i in intermediary_files:
	os.remove(i)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing intermediary files...OK'
print '----------------------------------------------------------------------------------------------'


#Quality
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Doing Base Quality Score Recalibration (Step1)...'
print '----------------------------------------------------------------------------------------------'
call("find " + output_path + baserecalibrator_input + " | parallel -j" + str(args.parallelization) + " 'java -Xmx9g -jar " + gatk + " \
	BaseRecalibrator \
	-R " + genome_ref + " \
	-I {} \
	--known-sites " + hg19_path + "1000G_phase1.indels.hg19.sites.vcf \
	--known-sites " + hg19_path + "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
	--known-sites " + hg19_path + "dbsnp_138.hg19.vcf \
	-O {}_recal_data.table'",shell=True)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Doing Base Quality Score Recalibration...(Step1)OK!'
print '----------------------------------------------------------------------------------------------'


print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Doing Base Quality Score Recalibration (Step2)...'
print '----------------------------------------------------------------------------------------------'
call("find " + output_path + baserecalibrator_input + " | parallel -j" + str(args.parallelization) + " 'java -Xmx9g -jar " + gatk + " \
	ApplyBQSR \
	-R " + genome_ref + " \
	-I {} \
	-bqsr {}_recal_data.table \
	-O {}_bqsr.bam'",shell=True)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Doing Base Quality Score Recalibration (Step2)...OK!'
print '----------------------------------------------------------------------------------------------'


#Crea archivo gVCF desde BAM
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Calling the variants...'
print '----------------------------------------------------------------------------------------------'
call("find " + output_path + "*_bqsr.bam | parallel -j" + str(args.parallelization) + " 'java -Xmx9g -jar " + gatk + " \
	HaplotypeCaller \
	-R " + genome_ref + " \
	-I {} \
	-ERC GVCF \
	-OVI \
	-O {}.g.vcf'",shell=True)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Calling the variants...OK!'
print '----------------------------------------------------------------------------------------------'


#Crea VCF desde gVCF
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Genotyping in single mode...'
print '----------------------------------------------------------------------------------------------'
call("find " + output_path + "*.g.vcf | parallel -j1 'java -Xmx28g -jar " + gatk + " \
	GenotypeGVCFs \
	-R " + genome_ref + " \
	-V {} \
	-O {}_singleGT_raw.vcf'",shell = True)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Genotyping in single mode...OK!'
print '----------------------------------------------------------------------------------------------'


#Remove gVCF files
if args.gvcf == None:
	print '----------------------------------------------------------------------------------------------'
	print '[FJD_Pipeline] Removing gVCF files...'
	print '----------------------------------------------------------------------------------------------'
	intermediary_files3 = glob(output_path + '*.g.vcf') + glob(output_path + '*.g.vcf.idx')
	for i in intermediary_files3:
		os.remove(i)
	print '----------------------------------------------------------------------------------------------'
	print '[FJD_Pipeline] Removing gVCF files...OK'
	print '----------------------------------------------------------------------------------------------'


#Empieza annovar
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Annotating Variants...'
print '----------------------------------------------------------------------------------------------'

variants = sorted(glob(args.output + '*singleGT_raw.vcf'))

for vcffile in variants:
	sample_path = vcffile[:vcffile.rfind('/')+1]
	sample_name = vcffile[vcffile.rfind('/')+1:vcffile.rfind('singleGT_raw.vcf')]
	output = sample_path + sample_name
	
	
	call(annovar + ' ' + vcffile + ' ' + annovarDB + ' -buildver hg19 \
		-out ' + output + ' \
		--remove \
		--otherinfo \
		--protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_eur,exac03,gnomad_exome,gnomad_genome,hrcr1,kaviar_20150923,popfreq_max_20150413,avsnp147,intervar_20170202,spidex,dbscsnv11,dbnsfp33a,revel,gwava,clinvar_20170130\
		--operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f \
		--nastring . \
		--vcfinput \
		--thread ' + str(args.threads),shell = True)


	print '----------------------------------------------------------------------------------------------'
	print '[FJD_Pipeline] Annotating Variants...OK!'
	print '----------------------------------------------------------------------------------------------'

	print '----------------------------------------------------------------------------------------------'
	print '[FJD_Pipeline] Formating and Filtering Variants...'
	print '----------------------------------------------------------------------------------------------'	

	if args.local:
		
		call('python /mnt/genetica/' + str(args.user) + '/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step1_v2.py ' + sample_path + '/' + sample_name + '.hg19_multianno.txt' ,shell = True)
		call('Rscript /mnt/genetica/' + str(args.user) + '/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step2.R ' + sample_path + '/' + sample_name + '_annotated_formatted.txt',shell = True)
		
	else:

		call('python /mnt/genetica/ionut/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step1_v2.py ' + sample_path + '/' + sample_name + '.hg19_multianno.txt' ,shell = True)
		call('Rscript /mnt/genetica/ionut/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step2_server.R ' + sample_path + '/' + sample_name + '_annotated_formatted.txt',shell = True)

	print '----------------------------------------------------------------------------------------------'
	print '[FJD_Pipeline] Formating and Filtering Variants...OK!'
	print '----------------------------------------------------------------------------------------------'

#Remove intermediary files
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing intermediary files...'
print '----------------------------------------------------------------------------------------------'
intermediary_files4 = glob(output_path + '*_recal_data.table') + glob(output_path + '*_.avinput') + glob(output_path +'*__annotated_formatted.txt') + glob(output_path + '*_dedupped.bam') + glob(output_path + '*_dedupped.bai')
for i in intermediary_files4:
	os.remove(i)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing intermediary files...OK'
print '----------------------------------------------------------------------------------------------'


print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Renaming files ...'
print '----------------------------------------------------------------------------------------------'

today = str(datetime.date.today())
test_date = re.sub('-', '', today)

#sample_path = glob('/mnt/genetica/ionut/AnalisisPipeline/Muestras_Prueba/fastq/*')
#sample_path = glob('/mnt/genetica/ionut/AnalisisPipeline/Muestras_Prueba/fastq/slice1/results/*')
sample_path = glob(output_path + '*')

while len(sample_path) > 0:

	for sample in sample_path:
		
		name = '_' + test_date + '_v5'

		if sample.endswith("_bqsr.bam"):
			new_name1 = re.sub('.fastq.g_bwa.sam_sorted.bam_dedupped.bam_bqsr.bam', name + '.bam', sample)
			os.rename(sample, new_name1)

		elif sample.endswith("_bqsr.bai"):
			new_name2 = re.sub('.fastq.g_bwa.sam_sorted.bam_dedupped.bam_bqsr.bai', name + '.bai', sample)
			os.rename(sample, new_name2)

		elif sample.endswith(".bam.g.vcf"):
			new_name3 = re.sub('.fastq.g_bwa.sam_sorted.bam_dedupped.bam_bqsr.bam.g.vcf', name + '.g.vcf', sample)
			os.rename(sample, new_name3)

		elif sample.endswith(".g.vcf.idx"):
			new_name4 = re.sub('.fastq.g_bwa.sam_sorted.bam_dedupped.bam_bqsr.bam.g.vcf.idx', name + '.g.vcf.idx', sample)
			os.rename(sample, new_name4)

		elif sample.endswith("_duplicate_metrics.txt"):
			new_name5 = re.sub('.fastq.g_bwa.sam_sorted.bam_duplicate_metrics.txt', name + '_duplicate_metrics.txt', sample)
			os.rename(sample, new_name5)

		elif sample.endswith("_singleGT_raw.vcf"):
			new_name6 = re.sub('.fastq.g_bwa.sam_sorted.bam_dedupped.bam_bqsr.bam.g.vcf_singleGT_raw.vcf', name + '_singleGT_raw.vcf', sample)
			os.rename(sample, new_name6)

		elif sample.endswith("_singleGT_raw.vcf.idx"):
			new_name7 = re.sub('.fastq.g_bwa.sam_sorted.bam_dedupped.bam_bqsr.bam.g.vcf_singleGT_raw.vcf.idx', name + '_singleGT_raw.vcf.idx', sample)
			os.rename(sample, new_name7)

		elif sample.endswith("_multianno.txt"):
			new_name8 = re.sub('.fastq.g_bwa.sam_sorted.bam_dedupped.bam_bqsr.bam.g.vcf_.hg19_multianno.txt', name + '_multianno.txt', sample)
			os.rename(sample, new_name8)

		elif sample.endswith("_multianno.vcf"):
			new_name9 = re.sub('.fastq.g_bwa.sam_sorted.bam_dedupped.bam_bqsr.bam.g.vcf_.hg19_multianno.vcf', name + '_multianno.vcf', sample)
			os.rename(sample, new_name9)

		elif sample.endswith("_raw_variants.txt"):
			new_name10 = re.sub('.fastq.g_bwa.sam_sorted.bam_dedupped.bam_bqsr.bam.g.vcf__annotated_formatted.txt_raw_variants.txt', name + '_raw_variants.txt', sample)
			os.rename(sample, new_name10)

		elif sample.endswith("_prefiltered.txt"):
			new_name11 = re.sub('.fastq.g_bwa.sam_sorted.bam_dedupped.bam_bqsr.bam.g.vcf__annotated_formatted.txt_prefiltered.txt', name + '_prefiltered.txt', sample)
			os.rename(sample, new_name11)


print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Renaming files...OK'
print '----------------------------------------------------------------------------------------------'

print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] #VARIANTS READY FOR ANALYSIS#'
print '----------------------------------------------------------------------------------------------'	
