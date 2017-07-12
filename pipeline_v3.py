
import sys
from glob import glob
from subprocess import call 

analysis = raw_input("Choose Amplicon(A) or Captured(C):")
runner = raw_input("Local(L) or Server(S):")

#Importing samples
forward_paths = sorted(glob(sys.argv[1] + '*R1_001.fastq.gz'))
reverse_paths = sorted(glob(sys.argv[1] + '*R2_001.fastq.gz'))


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


	call('bwa mem -t8 -R "@RG\tID:' + sample_name + '\tLB:library\tPL:illumina\tPU:library\tSM:sample" ' + genome_ref + ' ' + forward_paths[i] + ' ' + reverse_paths[i] + ' | \
		samtools view -Shu - | \
		samtools sort -@ 8 -m 3G - ' + sample_path + '/' + sample_name + '_sorted && samtools index ' + sample_path + '/' + sample_name + '_sorted.bam',shell =True)
	
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
variants = sorted(glob(sys.argv[1] + '*singleGT_raw.vcf'))
for vcffile in variants:
	sample_path = vcffile[:vcffile.rfind('/')+1]
	sample_name = vcffile[vcffile.rfind('/')+1:vcffile.rfind('singleGT_raw.vcf')]
	output = sample_path + sample_name
	
	call(annovar + ' ' + vcffile + ' ' + annovarDB + ' -buildver hg19 \
		-out ' + output + ' \
		--remove \
		--otherinfo \
		--protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2015aug_eur,exac03,gnomad_exome,gnomad_genome,hrcr1,kaviar_20150923,popfreq_max_20150413,avsnp147,intervar_20170202,spidex,dbscsnv11,dbnsfp33a,revel,gwava,clinvar_20170130,phastConsElements46way,tfbsConsSites,wgRna,targetScanS,gwasCatalog \
		--operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,r,r,r,r,r \
		--nastring . \
		--vcfinput \
		--thread 8',shell = True)

	
	#Converting from vcf to table
	call('python /mnt/genetica/domin/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step1.py ' + sample_path + '/' + sample_name + '.hg19_multianno.txt' ,shell = True)
	#call('python /mnt/genetica/domin/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step1_trio.py ' + sample_path + '/' + sample_name + '.hg19_multianno.txt' ,shell = True)
	call('Rscript /mnt/genetica/domin/GeneticaPipeDB/pipeline/annotation_scripts/vcf_processing_step2.R ' + sample_path + '/' + sample_name + '_annotated_formatted.txt',shell = True)




