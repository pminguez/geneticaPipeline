import sys
from glob import glob
from subprocess import call 
import argparse
import time


parser = argparse.ArgumentParser(description="Process Fastq files for getting variants")

parser.add_argument("-I", action="store",dest='input',
                    help="path to input folder")

parser.add_argument("-T", action="store",dest='threads', type = int, default = 12,
                    help="specify number of threads to use")

args = parser.parse_args()

annovar = "/home/domin/genetica/GeneticaPipeDB/software/annovar/table_annovar.pl"
annovarDB = "/home/domin/genetica/GeneticaPipeDB/software/annovar/humandb"

variants = sorted(glob(args.input + '*singleGT_raw.vcf'))


for vcffile in variants:
	sample_path = vcffile[:vcffile.rfind('/')+1]
	#sample_name = vcffile[vcffile.rfind('/')+1:vcffile.rfind('singleGT_raw.vcf')]
	sample_name = vcffile[vcffile.rfind('/')+1:vcffile.rfind('.hg19_multianno.vcf')]

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
	'''

	call(annovar + ' ' + vcffile + ' ' + annovarDB + ' -buildver hg19 \
		-out ' + output + ' \
		--remove \
		--otherinfo \
		--protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_eur,exac03,gnomad_exome,gnomad_genome,hrcr1,kaviar_20150923,popfreq_max_20150413,avsnp147,intervar_20170202,cosmic70,icgc21,spidex,dbscsnv11,dbnsfp33a,revel,gwava,clinvar_20170130 \
		--operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f \
		--nastring . \
		--vcfinput \
		--thread ' + str(args.threads),shell = True)
	'''
	print '----------------------------------------------------------------------------------------------'
	print '[DAguilera_Pipeline] Annotating Variants...OK!'
	print '----------------------------------------------------------------------------------------------'

	print '----------------------------------------------------------------------------------------------'
	print '[DAguilera_Pipeline] Formating and Filtering Variants...'
	print '----------------------------------------------------------------------------------------------'	


	call('python /mnt/genetica/domin/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step1_v2.py ' + sample_path + sample_name + '.hg19_multianno.txt' ,shell = True)
	call('Rscript /mnt/genetica/domin/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step2_server.R ' + sample_path + '/' + sample_name + '_annotated_formatted.txt',shell = True)
	
	print '----------------------------------------------------------------------------------------------'
	print '[DAguilera_Pipeline] Formating and Filtering Variants...OK!'
	print '----------------------------------------------------------------------------------------------'

print '----------------------------------------------------------------------------------------------'
print '[DAguilera_Pipeline] #VARIANTS READY FOR ANALYSIS#'
print '----------------------------------------------------------------------------------------------'	
