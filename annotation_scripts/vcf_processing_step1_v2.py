import sys
vcfFile = open(sys.argv[1], 'r')
vcfFile = vcfFile.readlines()
header = vcfFile[0]
vcfFile = vcfFile[1:]
header = header.split('\t')
header = header[0:len(header)-1]
header_newvalues = ['QUAL','DP','AC','AF','AN','ExcessHet','FS','MLEAC','MLEAF','MQ','QD','BaseQRankSum','ClippingRankSum','MQRankSum','ReadPosRankSum','SOR','GT','ALTD','GQ','PL']
for string in header_newvalues:
	header.append(string)
header = str(header).replace("',", "'\t")
header = str(header).replace("'", "")		
result_line = list()
result_line.append(header)

for i in range(0,len(vcfFile)):
	vcfFile_lines = vcfFile[i].split('\t')

	field2split = vcfFile_lines[168]
	field2split_2 = vcfFile_lines[170]
	#field2split = vcfFile_lines[97]
	#field2split_2 = vcfFile_lines[99]
	field2split_2 = field2split_2.split(':')
	QUAL = vcfFile_lines[166]
	DP = vcfFile_lines[160]
	#QUAL = vcfFile_lines[88]
	#DP = vcfFile_lines[89]
	
	AC = field2split[field2split.find('AC=')+len('AC='):field2split.find(';')]
	field2split_temp = field2split[field2split.find(';AF=')+1:]
	AF = field2split_temp[field2split_temp.find('AF=')+len('AF='):field2split_temp.find(';')]
	field2split_temp = field2split[field2split.find(';AN=')+1:]
	AN = field2split_temp[field2split_temp.find('AN=')+len('AN='):field2split_temp.find(';')]
	field2split_temp = field2split[field2split.find(';ExcessHet=')+1:]
	ExcessHet = field2split_temp[field2split_temp.find('ExcessHet=')+len('ExcessHet='):field2split_temp.find(';')]
	field2split_temp = field2split[field2split.find(';FS=')+1:]
	FS = field2split_temp[field2split_temp.find('FS=')+len('FS='):field2split_temp.find(';')]
	field2split_temp = field2split[field2split.find(';MLEAC=')+1:]
	MLEAC = field2split_temp[field2split_temp.find('MLEAC=')+len('MLEAC='):field2split_temp.find(';')]
	field2split_temp = field2split[field2split.find(';MLEAF=')+1:]
	MLEAF = field2split_temp[field2split_temp.find('MLEAF=')+len('MLEAF='):field2split_temp.find(';')]
	field2split_temp = field2split[field2split.find(';MQ=')+1:]
	MQ = field2split_temp[field2split_temp.find('MQ=')+len('MQ='):field2split_temp.find(';')]
	field2split_temp = field2split[field2split.find(';QD=')+1:]
	QD = field2split_temp[field2split_temp.find('QD=')+len('QD='):field2split_temp.find(';')]
	field2split_temp = field2split[field2split.find(';BaseQRankSum=')+1:]
	BaseQRankSum = field2split_temp[field2split_temp.find('BaseQRankSum=')+len('BaseQRankSum='):field2split_temp.find(';')]
	field2split_temp = field2split[field2split.find(';ClippingRankSum=')+1:]
	ClippingRankSum = field2split_temp[field2split_temp.find('ClippingRankSum=')+len('ClippingRankSum='):field2split_temp.find(';')]
	field2split_temp = field2split[field2split.find(';MQRankSum=')+1:]
	MQRankSum = field2split_temp[field2split_temp.find('MQRankSum=')+len('MQRankSum='):field2split_temp.find(';')]
	field2split_temp = field2split[field2split.find(';ReadPosRankSum=')+1:]
	ReadPosRankSum = field2split_temp[field2split_temp.find('ReadPosRankSum=')+len('ReadPosRankSum='):field2split_temp.find(';')]
	field2split_temp = field2split[field2split.find(';SOR=')+1:]
	SOR = field2split_temp[field2split_temp.find('SOR=')+len('SOR='):]

	GT = field2split_2[0]
	
	if len(field2split_2) < 5:
		ALTD = ''
		GQ = field2split_2[1]
		PL = field2split_2[2]

	else:
		ALTD = field2split_2[1][field2split_2[1].rfind(',')+1:]
		GQ = field2split_2[3]
		PL = field2split_2[4][:field2split_2[4].rfind('\n')]

	PL = PL.replace(',',';')

	#vcfFile_lines = list(vcfFile_lines[0:87])
	vcfFile_lines = list(vcfFile_lines[0:158])
	
	
	for element in [QUAL,DP,AC,AF,AN,ExcessHet,FS,MLEAC,MLEAF,MQ,QD,BaseQRankSum,ClippingRankSum,MQRankSum,ReadPosRankSum,SOR,GT,ALTD,GQ,PL]:
		if element == '':
			element = '.'
		vcfFile_lines.append(element)
		
	vcfFile_lines = str(vcfFile_lines).replace("',", "'\t")
	vcfFile_lines = str(vcfFile_lines).replace("'", "")
	result_line.append(vcfFile_lines)

output= open(sys.argv[1][:sys.argv[1].rfind('.hg19')] + '_annotated_formatted.txt','w')
for item in result_line:
  output.write("%s\n" % item)
output.close()
