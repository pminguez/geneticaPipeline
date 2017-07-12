from subprocess import call 
import sys
import os

#Arguments
bam_files_folder = sys.argv[1]
bedfile = sys.argv[2]
output_folder = sys.argv[3]

#Creating Controls folder if doesn't exist
Controls = output_folder + 'Controls'
if not os.path.exists(Controls):
    os.makedirs(Controls)

#Creating StartWithBestScore folder if doesn't exist
StartWithBestScore = output_folder + 'StartWithBestScore'
if not os.path.exists(StartWithBestScore):
    os.makedirs(StartWithBestScore)

#Creating CNVs_calling_results folder if doesn't exist
CNVs_calling_results = output_folder + 'CNVs_calling_results'
if not os.path.exists(CNVs_calling_results):
    os.makedirs(CNVs_calling_results)

#Creating TargetQcList folder if doesn't exist
TargetQcList = output_folder + 'TargetQcList'
if not os.path.exists(TargetQcList):
    os.makedirs(TargetQcList)

#Creating FinalList folder if doesn't exist
FinalList = output_folder + 'FinalList'
if not os.path.exists(FinalList):
    os.makedirs(FinalList)
'''
#Calculating coverages for all samples
call('python /mnt/datos2/d.aguilera/CNVs_analysis/CoNVaDING-1.1.6/CoNVaDING.py \
	' + bam_files_folder + ' \
	' + bedfile + ' \
	' + output_folder,shell = True)

#Copying coverages folder to Controls folder
call('cp ' + output_folder + 'StartWithMatchScore/* ' + Controls, shell = True)


#Calculating 100 Best Controls from StartWithMatchScore folder

call('perl /mnt/datos2/d.aguilera/CNVs_analysis/CoNVaDING-1.1.6/CoNVaDING.pl \
	-mode StartWithMatchScore \
	-inputDir ' + output_folder + 'StartWithMatchScore \
	-controlsDir ' + Controls + ' \
	-outputDir ' + StartWithBestScore + ' \
	-controlSamples 100 \
	-sexChr ',shell = True)
'''
#Making the CNVs calling
call('perl /mnt/datos2/d.aguilera/CNVs_analysis/CoNVaDING-1.1.6/CoNVaDING.pl \
	-mode StartWithBestScore \
	-inputDir ' + StartWithBestScore + ' \
	-controlsDir ' + Controls + ' \
	-outputDir ' + CNVs_calling_results,shell = True)

#Creating TargetQcList
call('perl /mnt/datos2/d.aguilera/CNVs_analysis/CoNVaDING-1.1.6/CoNVaDING.pl \
	-mode GenerateTargetQcList \
	-inputDir ' + Controls + ' \
	-controlsDir ' + Controls + ' \
	-outputDir ' + TargetQcList,shell = True)

#Creating FinalList
call('perl /mnt/datos2/d.aguilera/CNVs_analysis/CoNVaDING-1.1.6/CoNVaDING.pl \
	-mode CreateFinalList \
	-inputDir ' + StartWithBestScore + ' \
	-targetQcList ' + TargetQcList + '/targetQcList.txt \
	-outputDir ' + FinalList,shell = True)


