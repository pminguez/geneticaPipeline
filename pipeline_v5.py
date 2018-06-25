#############################################################################
## Pipeline is now on github: https://github.com/pminguez/geneticaPipeline ##
#############################################################################

import sys
from glob import glob
from subprocess import call 
import argparse
import time
import os

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

parser.add_argument("-I", action="store",dest='input',
                    help="path to input folder")

parser.add_argument("-T", action="store",dest='threads', type = int, default = 16,
                    help="specify number of threads to use")

parser.add_argument("-J", action="store",dest='parallelization', type = int, default = 5,
                    help="specify number of samples to run in parallel")

parser.add_argument("-duplicates", action="store_true",
                    help="set this flag to markduplicates with picardtools")

parser.add_argument("-local", action="store_true",
                    help="set this flag to run the pipeline using local paths")

args = parser.parse_args()


if args.input == None:
	print ''
	print 'ERROR: An input folder containing fastq files is needed'
	print ''
	parser.print_help()
	exit()

#Importing samples
forward_paths = sorted(glob(args.input + '*_R1.fastq.gz'))
reverse_paths = sorted(glob(args.input + '*_R2.fastq.gz'))


if forward_paths == []:
	print ''
	print 'ERROR: No fastq files detected in ' + args.input + '.\nFastq files names should be named: name_R1.fastq.gz and name_R2.fastq.gz'
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
print '	-Threads: ' + str(args.threads)
print '	-Sample to parallelizate: ' + str(args.parallelization)
print '	-MarkDuplicates: ' + str(args.duplicates)
print '	-Running local: ' + str(args.local)
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

else:
	genome_ref = "/mnt/genetica/" + str(args.user) + "/GeneticaPipeDB/genome_data/hg19/ucsc.hg19.fasta"
	picardtools = "/mnt/genetica/"+ str(args.user) + "/GeneticaPipeDB/software/picard-tools-2.1.1/picard.jar"
	gatk = "/mnt/genetica/"+ str(args.user) + "/GeneticaPipeDB/software/GenomeAnalysisTK-3.5/GenomeAnalysisTK.jar"
	hg19_path = "/mnt/genetica/"+ str(args.user) + "/GeneticaPipeDB/genome_data/hg19/"
	annovar = "/mnt/genetica/"+ str(args.user) + "/GeneticaPipeDB/software/annovar/table_annovar.pl"
	annovarDB = "/mnt/genetica/"+ str(args.user) + "/GeneticaPipeDB/software/annovar/humandb"



print '                               Mapping fastq files (BWA)                                      '
print '----------------------------------------------------------------------------------------------'

#Load genome reference to memory

call('bwa shm ' + genome_ref,shell = True)


#Loop samples for BWA
for i in range(0,len(forward_paths)):
	sample_path = forward_paths[i][:forward_paths[i].rfind('/')+1]
	sample_name = forward_paths[i][forward_paths[i].rfind('/')+1:forward_paths[i].rfind('_R')]

	call('bwa mem -t' + str(args.threads) + ' -R "@RG\tID:' + sample_name + '\tLB:library\tPL:illumina\tPU:library\tSM:' + sample_name + '" ' + genome_ref + ' ' + forward_paths[i] + ' ' + reverse_paths[i] + ' > ' + sample_path + '/' + sample_name + '_bwa.sam',shell = True)



#Unload genome reference	
call('bwa shm -d',shell = True)
print '----------------------------------------------------------------------------------------------'

print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Sorting and creating bam and bai files of samples...'
print '----------------------------------------------------------------------------------------------'
call("find " + sample_path +  "*.sam | parallel --no-notice -j" + str(args.parallelization) + " 'samtools sort {} -O BAM -@ " + str(args.threads / 2) + " -o {}_sorted.bam && samtools index {}_sorted.bam'", shell = True)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] OK! '
print '----------------------------------------------------------------------------------------------'

indelrealigner_input = '*_sorted.bam'


#Remove sam files
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing sam files...'
print '----------------------------------------------------------------------------------------------'
for i in glob(sample_path + '*.sam'):
	os.remove(i)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing sam files...OK'
print '----------------------------------------------------------------------------------------------'


#MarkDuplicates with picardtools
if args.duplicates:
	print '----------------------------------------------------------------------------------------------'
	print '[FJD_Pipeline] Marking Duplicates... '
	print '----------------------------------------------------------------------------------------------'
	call("find " + sample_path + "*_sorted.bam | parallel --no-notice -j" + str(args.parallelization) + " 'java -Xmx9g -jar " + picardtools + " \
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

	indelrealigner_input = '*_dedupped.bam'


#Empieza GATK
#Indel Realignment with GATK
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Doing IndelRealigment... '
print '---------------------------------------------------------------------------------------------'
call("find " + sample_path +  indelrealigner_input + " | parallel --no-notice -j" + str(args.parallelization) + " 'java -Xmx9g -jar " + gatk + " \
	-T IndelRealigner \
	-R " + genome_ref + " \
	-I {} \
	-o {}_indelrealigned.bam \
	-targetIntervals " + hg19_path + "hg19_indels_output.intervals \
	-known " + hg19_path + "1000G_phase1.indels.hg19.sites.vcf \
	-known " + hg19_path + "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
	-LOD 0.4 \
	--consensusDeterminationModel KNOWNS_ONLY'",shell=True)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Doing IndelRealigment...OK!'
print '----------------------------------------------------------------------------------------------'


#Remove intermediary files
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing intermediary files...'
print '----------------------------------------------------------------------------------------------'
intermediary_files = glob(sample_path + '*_dedupped.bam') + glob(sample_path + '*_dedupped.bai') + glob(sample_path + '*_sorted.bam') + glob(sample_path + '*_sorted.bam.bai')
for i in intermediary_files:
	os.remove(i)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing intermediary files...OK'
print '----------------------------------------------------------------------------------------------'


#Quality
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Doing Base Quality Score Recalibration (Step1)...'
print '----------------------------------------------------------------------------------------------'
call("find " + sample_path + "*_indelrealigned.bam | parallel -j" + str(args.parallelization) + " 'java -Xmx9g -jar " + gatk + " \
	-T BaseRecalibrator \
	-R " + genome_ref + " \
	-I {} \
	-knownSites " + hg19_path + "1000G_phase1.indels.hg19.sites.vcf \
	-knownSites " + hg19_path + "Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
	-knownSites " + hg19_path + "dbsnp_138.hg19.vcf \
	-o {}_recal_data.table'",shell=True)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Doing Base Quality Score Recalibration...(Step1)OK!'
print '----------------------------------------------------------------------------------------------'


print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Doing Base Quality Score Recalibration (Step2)...'
print '----------------------------------------------------------------------------------------------'
call("find " + sample_path + "*_indelrealigned.bam | parallel -j" + str(args.parallelization) + " 'java -Xmx9g -jar " + gatk + " \
	-T PrintReads \
	-R " + genome_ref + " \
	-I {} \
	-BQSR {}_recal_data.table \
	-o {}_bqsr.bam'",shell=True)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Doing Base Quality Score Recalibration (Step2)...OK!'
print '----------------------------------------------------------------------------------------------'


#Remove intermediary files
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing intermediary files...'
print '----------------------------------------------------------------------------------------------'
intermediary_files2 = glob(sample_path + '*_indelrealigned.bam') + glob(sample_path + '*_indelrealigned.bai')
for i in intermediary_files2:
	os.remove(i)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing intermediary files...OK'
print '----------------------------------------------------------------------------------------------'



#Crea archivo gVCF desde BAM
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Calling the variants...'
print '----------------------------------------------------------------------------------------------'
call("find " + sample_path + "*_bqsr.bam | parallel -j" + str(args.parallelization) + " 'java -Xmx9g -jar " + gatk + " \
	-T HaplotypeCaller \
	-R " + genome_ref + " \
	-I {} \
	--emitRefConfidence GVCF \
	--variant_index_type LINEAR \
 	--variant_index_parameter 128000 \
	-o {}.g.vcf'",shell=True)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Calling the variants...OK!'
print '----------------------------------------------------------------------------------------------'


#PARA ANALIZAR TRIOS.
'''
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Merging the gVCFs...'
print '----------------------------------------------------------------------------------------------'
call("find " + sample_path + "*.g.vcf | parallel -j" + str(args.parallelization) + " 'java -Xmx9g -jar " + gatk + " \
	-T CombineGVCFs \
	-R " + genome_ref + " \
	-I {} \
	-o {}_trio.g.vcf'",shell=True)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Merging the gVCFs...OK'
print '----------------------------------------------------------------------------------------------'

#Si se corre la pipeline para trios usar el siguiente paso


print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Genotyping trio...'
print '----------------------------------------------------------------------------------------------'
call("find " + sample_path + "*_trio.g.vcf | parallel -j1 'java -Xmx28g -jar " + gatk + " \
	-T GenotypeGVCFs \
	-nt " + str(args.threads) + " \
	-R " + genome_ref + " \
	-V {} \
	-o {}_singleGT_trio_raw.vcf'",shell = True)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Genotyping trio...OK!'
print '----------------------------------------------------------------------------------------------'

'''


#Crea VCF desde gVCF
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Genotyping in single mode...'
print '----------------------------------------------------------------------------------------------'
call("find " + sample_path + "*.g.vcf | parallel -j1 'java -Xmx28g -jar " + gatk + " \
	-T GenotypeGVCFs \
	-nt " + str(args.threads) + " \
	-R " + genome_ref + " \
	-V {} \
	-o {}_singleGT_raw.vcf'",shell = True)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Genotyping in single mode...OK!'
print '----------------------------------------------------------------------------------------------'

'''
#Remove gVCF files
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing gVCF files...'
print '----------------------------------------------------------------------------------------------'
intermediary_files3 = glob(sample_path + '*.g.vcf') + glob(sample_path + '*.g.vcf.idx')
for i in intermediary_files3:
	os.remove(i)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing gVCF files...OK'
print '----------------------------------------------------------------------------------------------'
'''

#Empieza annovar
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Annotating Variants...'
print '----------------------------------------------------------------------------------------------'

variants = sorted(glob(args.input + '*singleGT_raw.vcf'))

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

		call('python /mnt/genetica/' + str(args.user) + '/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step1_v2.py ' + sample_path + '/' + sample_name + '.hg19_multianno.txt' ,shell = True)
		call('Rscript /mnt/genetica/' + str(args.user) + '/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step2_server.R ' + sample_path + '/' + sample_name + '_annotated_formatted.txt',shell = True)

	print '----------------------------------------------------------------------------------------------'
	print '[FJD_Pipeline] Formating and Filtering Variants...OK!'
	print '----------------------------------------------------------------------------------------------'

'''
#Remove intermediary files
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing intermediary files...'
print '----------------------------------------------------------------------------------------------'
intermediary_files4 = glob(sample_path + '*.hg19_multianno.vcf') + glob(sample_path + '*_annotated_formatted.txt')
for i in intermediary_files4:
	os.remove(i)
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] Removing intermediary files...OK'
print '----------------------------------------------------------------------------------------------'
'''
print '----------------------------------------------------------------------------------------------'
print '[FJD_Pipeline] #VARIANTS READY FOR ANALYSIS#'
print '----------------------------------------------------------------------------------------------'	
