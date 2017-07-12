from subprocess import call 
from glob import glob
import sys
import os

bam_folder = glob('/mnt/datos2/Malformaciones_Oculares/Haloplex_MOC/resultados/MOC_2_NextSeq/FJD_pipeline_Analysis_kit2_NextSeq/bam_files/*bam')
#Calculating mean coverage of each target region (bedfile) for every sample and storing it in the working folder

bedfile = '/mnt/datos2/Malformaciones_Oculares/Haloplex_MOC/resultados/MOC_1_MiSeq/FJD_pipeline_Analysis_kit1_MiSeq/MOC_MiSeq_pool1/MOC_Haloplex_targets.bed'

raw_coverages_folder = '/mnt/datos2/d.aguilera/Halo_MOC/Nextseq/'

#Creating raw_coverages folder if doesn't exist
if not os.path.exists(raw_coverages_folder):
    os.makedirs(raw_coverages_folder)


def get_sample_name(sample_path):
	sample_name = sample_path[sample_path.rfind('/') +1:]
	sample_name = sample_name[:sample_name.find('_sorted')]
	return sample_name


for bam_file in bam_folder:
	sample_name = get_sample_name(bam_file)
	print '-----'
	print 'Calculating coverage of sample: ' + sample_name + ' ...'
	#call('bedtools coverage -a ' + bedfile + ' -b ' + bam_file + ' -sorted -g /mnt/datos2/d.aguilera/CNVs_analysis/CoNVaDING-1.1.6/ucsc.hg19_convading.fasta.fai -mean| sort -k1,1 -k2,2n  > ' + raw_coverages_folder + '/' + sample_name + '_mean_coverage.txt',shell = True)
	#call('bedtools coverage -a ' + bedfile + ' -b ' + bam_file + ' -mean| sort -k1,1 -k2,2n > ' + raw_coverages_folder + '/' + sample_name + '_mean_coverage.txt',shell = True)
	call('bedtools coverage -hist -a ' + bedfile + ' -b ' + bam_file + ' | grep ^all > ' + bam_file + '.hist.all.txt',shell = True)
	print 'Saved in ' + raw_coverages_folder



