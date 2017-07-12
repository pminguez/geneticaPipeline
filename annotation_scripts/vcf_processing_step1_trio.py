import sys
vcfFile = open(sys.argv[1], 'r')
vcfFile = vcfFile.readlines()
header = vcfFile[0]
vcfFile = vcfFile[1:]
header = header.split('\t')
header = header[0:len(header)-1]
#header_newvalues = ['GT_father','GT_mother','GT_child','QUAL','DP','AC','AF','AN','ExcessHet','FS','MLEAC','MLEAF','MQ','QD','BaseQRankSum','ClippingRankSum','MQRankSum','ReadPosRankSum','SOR','ALTD_father','GQ_father','PGT_father','PID_father','PL_father','TP_father','ALTD_mother','GQ_mother','PGT_mother','PID_mother','PL_mother','TP_mother','ALTD_child','GQ_child','PGT_child','PID_child','PL_child','TP_child','']
header_newvalues = ['GT_father','GT_mother','GT_sister','GT_child','QUAL','DP','AC','AF','AN','ExcessHet','FS','MLEAC','MLEAF','MQ','QD','BaseQRankSum','ClippingRankSum','MQRankSum','ReadPosRankSum','SOR','DP_father','ALTD_father','GQ_father','PGT_father','PID_father','PL_father','TP_father','DP_mother','ALTD_mother','GQ_mother','PGT_mother','PID_mother','PL_mother','TP_mother','DP_sister','ALTD_sister','GQ_sister','PGT_sister','PID_sister','PL_sister','TP_sister','DP_child','ALTD_child','GQ_child','PGT_child','PID_child','PL_child','TP_child','']
for string in header_newvalues:
	header.append(string)
header = str(header).replace("',", "'\t")
header = str(header).replace("'", "")		
result_line = list()
result_line.append(header)

for i in range(0,len(vcfFile)):
	vcfFile_lines = vcfFile[i].split('\t')	
	field2split = vcfFile_lines[156]
	field2split_2 = vcfFile_lines[160]
	field2split_3 = vcfFile_lines[158]
	field2split_4 = vcfFile_lines[159]
	field2split_5 = vcfFile_lines[161]

	#field2split = vcfFile_lines[97]
	#field2split_2 = vcfFile_lines[99]
	field2split_2 = field2split_2.split(':')
	field2split_3 = field2split_3.split(':')
	field2split_4 = field2split_4.split(':')
	field2split_5 = field2split_5.split(':')

	QUAL = vcfFile_lines[154]
	DP = vcfFile_lines[148]
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

	GT_father = field2split_2[0]
	GT_mother = field2split_3[0]
	GT_child = field2split_4[0]
	GT_sister = field2split_5[0]

	if GT_father == '1/1':
		GT_father = 'Hom'

	elif GT_father == '0/1':
		GT_father = 'Het'

	elif GT_father == '0/0':
		GT_father = 'Wt'

	elif GT_father == './.':
		GT_father = 'Not_covered'

	if GT_mother == '1/1':
		GT_mother = 'Hom'

	elif GT_mother == '0/1':
		GT_mother = 'Het'

	elif GT_mother == '0/0':
		GT_mother = 'Wt'

	elif GT_mother == './.':
		GT_mother = 'Not_covered'	

	if GT_child == '1/1':
		GT_child = 'Hom'

	elif GT_child == '0/1':
		GT_child = 'Het'

	elif GT_child == '0/0':
		GT_child = 'Wt'

	elif GT_child == './.':
		GT_child = 'Not_covered'

	
	if len(field2split_2) == 3:
		ALTD_father = ''
		DP_father = ''
		GQ_father = ''
		PGT_father = ''
		PID_father = ''
		PL_father= ''
		TP_father = ''

	elif len(field2split_2) == 5:
		if field2split_2[1].rfind(',') == field2split_2[1].find(','): 
			ALTD_father = field2split_2[1][field2split_2[1].rfind(',')+1:]
		else:
			ALTD_father = field2split_2[1]

		DP_father = field2split_2[2]
		GQ_father = field2split_2[3]
		PGT_father = ''
		PID_father = ''
		PL_father= field2split_2[4][:field2split_2[4].rfind('\n')]
		TP_father = ''

	elif len(field2split_2) == 6:
		if field2split_2[1].rfind(',') == field2split_2[1].find(','): 
			ALTD_father = field2split_2[1][field2split_2[1].rfind(',')+1:]
		else:
			ALTD_father = field2split_2[1]

		DP_father = field2split_2[2]
		GQ_father = field2split_2[3]
		PGT_father = ''
		PID_father = ''
		PL_father= field2split_2[4]
		TP_father = field2split_2[5][:field2split_2[5].rfind('\n')]
	
	elif len(field2split_2) == 7:
		if field2split_2[1].rfind(',') == field2split_2[1].find(','): 
			ALTD_father = field2split_2[1][field2split_2[1].rfind(',')+1:]
		else:
			ALTD_father = field2split_2[1]

		DP_father = field2split_2[2]
		GQ_father = field2split_2[3]
		PGT_father = field2split_2[4]
		PID_father = field2split_2[5]
		PL_father= field2split_2[6][:field2split_2[6].rfind('\n')]
		TP_father = ''

	elif len(field2split_2) == 8:
		if field2split_2[1].rfind(',') == field2split_2[1].find(','): 
			ALTD_father = field2split_2[1][field2split_2[1].rfind(',')+1:]
		else:
			ALTD_father = field2split_2[1]

		DP_father = field2split_2[2]
		GQ_father = field2split_2[3]
		PGT_father = field2split_2[4]
		PID_father = field2split_2[5]
		PL_father= field2split_2[6]
		TP_father = field2split_2[7][:field2split_2[7].rfind('\n')]



	if len(field2split_3) == 3:
		ALTD_mother = ''
		DP_mother = ''
		GQ_mother = ''
		PGT_mother = ''
		PID_mother = ''
		PL_mother = ''
		TP_mother = ''

	elif len(field2split_3) == 5:
		if field2split_3[1].rfind(',') == field2split_3[1].find(','): 
			ALTD_mother = field2split_3[1][field2split_3[1].rfind(',')+1:]
		else:
			ALTD_mother = field2split_3[1]

		DP_mother = field2split_3[2]
		GQ_mother = field2split_3[3]
		PGT_mother = ''
		PID_mother = ''
		PL_mother = field2split_3[4][:field2split_3[4].rfind('\n')]
		TP_mother = ''

	elif len(field2split_3) == 6:
		if field2split_3[1].rfind(',') == field2split_3[1].find(','): 
			ALTD_mother = field2split_3[1][field2split_3[1].rfind(',')+1:]
		else:
			ALTD_mother = field2split_3[1]

		DP_mother = field2split_3[2]
		GQ_mother = field2split_3[3]
		PGT_mother = ''
		PID_mother = ''
		PL_mother = field2split_3[4]
		TP_mother = field2split_3[5][:field2split_3[5].rfind('\n')]
	
	elif len(field2split_3) == 7:
		if field2split_3[1].rfind(',') == field2split_3[1].find(','): 
			ALTD_mother = field2split_3[1][field2split_3[1].rfind(',')+1:]
		else:
			ALTD_mother = field2split_3[1]

		DP_mother = field2split_3[2]
		GQ_mother = field2split_3[3]
		PGT_mother = field2split_3[4]
		PID_mother = field2split_3[5]
		PL_mother = field2split_3[6][:field2split_3[6].rfind('\n')]
		TP_mother = ''

	elif len(field2split_3) == 8:
		if field2split_3[1].rfind(',') == field2split_3[1].find(','): 
			ALTD_mother = field2split_3[1][field2split_3[1].rfind(',')+1:]
		else:
			ALTD_mother = field2split_3[1]

		DP_mother = field2split_3[2]
		GQ_mother = field2split_3[3]
		PGT_mother = field2split_3[4]
		PID_mother = field2split_3[5]
		PL_mother = field2split_3[6]
		TP_mother = field2split_3[7][:field2split_3[7].rfind('\n')]


	if len(field2split_4) == 3:
		ALTD_child = ''
		DP_child = ''
		GQ_child = ''
		PGT_child = ''
		PID_child = ''
		PL_child = ''
		TP_child = ''

	elif len(field2split_4) == 5:
		if field2split_4[1].rfind(',') == field2split_4[1].find(','): 
			ALTD_child = field2split_4[1][field2split_4[1].rfind(',')+1:]
		else:
			ALTD_child = field2split_4[1]

		DP_child = field2split_4[2]
		GQ_child = field2split_4[3]
		PGT_child = ''
		PID_child = ''
		PL_child = field2split_4[4][:field2split_4[4].rfind('\n')]
		TP_child = ''

	elif len(field2split_4) == 6:
		if field2split_4[1].rfind(',') == field2split_4[1].find(','): 
			ALTD_child = field2split_4[1][field2split_4[1].rfind(',')+1:]
		else:
			ALTD_child = field2split_4[1]

		DP_child = field2split_4[2]
		GQ_child = field2split_4[3]
		PGT_child = ''
		PID_child = ''
		PL_child = field2split_4[4]
		TP_child = field2split_4[5][:field2split_4[5].rfind('\n')]
	
	elif len(field2split_4) == 7:
		if field2split_4[1].rfind(',') == field2split_4[1].find(','): 
			ALTD_child = field2split_4[1][field2split_4[1].rfind(',')+1:]
		else:
			ALTD_child = field2split_4[1]

		DP_child = field2split_4[2]
		GQ_child = field2split_4[3]
		PGT_child = field2split_4[4]
		PID_child = field2split_4[5]
		PL_child = field2split_4[6][:field2split_4[6].rfind('\n')]
		TP_child = ''

	elif len(field2split_4) == 8:
		if field2split_4[1].rfind(',') == field2split_4[1].find(','): 
			ALTD_child = field2split_4[1][field2split_4[1].rfind(',')+1:]
		else:
			ALTD_child = field2split_4[1]

		DP_child = field2split_4[2]
		GQ_child = field2split_4[3]
		PGT_child = field2split_4[4]
		PID_child = field2split_4[5]
		PL_child = field2split_4[6]
		TP_child = field2split_4[7][:field2split_4[7].rfind('\n')]


	if len(field2split_5) == 3:
		ALTD_sister = ''
		DP_sister = ''
		GQ_sister = ''
		PGT_sister = ''
		PID_sister = ''
		PL_sister= ''
		TP_sister = ''

	elif len(field2split_5) == 5:
		if field2split_5[1].rfind(',') == field2split_5[1].find(','): 
			ALTD_sister = field2split_5[1][field2split_5[1].rfind(',')+1:]
		else:
			ALTD_sister = field2split_5[1]

		DP_sister = field2split_5[2]
		GQ_sister = field2split_5[3]
		PGT_sister = ''
		PID_sister = ''
		PL_sister= field2split_5[4][:field2split_5[4].rfind('\n')]
		TP_sister = ''

	elif len(field2split_5) == 6:
		if field2split_5[1].rfind(',') == field2split_5[1].find(','): 
			ALTD_sister = field2split_5[1][field2split_5[1].rfind(',')+1:]
		else:
			ALTD_sister = field2split_5[1]

		DP_sister = field2split_5[2]
		GQ_sister = field2split_5[3]
		PGT_sister = ''
		PID_sister = ''
		PL_sister= field2split_5[4]
		TP_sister = field2split_5[5][:field2split_5[5].rfind('\n')]
	
	elif len(field2split_5) == 7:
		if field2split_5[1].rfind(',') == field2split_5[1].find(','): 
			ALTD_sister = field2split_5[1][field2split_5[1].rfind(',')+1:]
		else:
			ALTD_sister = field2split_5[1]

		DP_sister = field2split_5[2]
		GQ_sister = field2split_5[3]
		PGT_sister = field2split_5[4]
		PID_sister = field2split_5[5]
		PL_sister= field2split_5[6][:field2split_5[6].rfind('\n')]
		TP_sister = ''

	elif len(field2split_5) == 8:
		if field2split_5[1].rfind(',') == field2split_5[1].find(','): 
			ALTD_sister = field2split_5[1][field2split_5[1].rfind(',')+1:]
		else:
			ALTD_sister = field2split_5[1]

		DP_sister = field2split_5[2]
		GQ_sister = field2split_5[3]
		PGT_sister = field2split_5[4]
		PID_sister = field2split_5[5]
		PL_sister= field2split_5[6]
		TP_sister = field2split_5[7][:field2split_5[7].rfind('\n')]


	PL_father = PL_father.replace(',',';')
	PL_mother = PL_mother.replace(',',';')
	PL_child = PL_child.replace(',',';')
	PL_sister = PL_sister.replace(',',';')
	#vcfFile_lines = list(vcfFile_lines[0:87])
	vcfFile_lines = list(vcfFile_lines[0:146])

	
	for element in [GT_father,GT_mother,GT_sister,GT_child,QUAL,DP,AC,AF,AN,ExcessHet,FS,MLEAC,MLEAF,MQ,QD,BaseQRankSum,ClippingRankSum,MQRankSum,ReadPosRankSum,SOR,DP_father,ALTD_father,GQ_father,PGT_father,PID_father,PL_father,TP_father,DP_mother,ALTD_mother,GQ_mother,PGT_mother,PID_mother,PL_mother,TP_mother,DP_sister,ALTD_sister,GQ_sister,PGT_sister,PID_sister,PL_sister,TP_sister,DP_child,ALTD_child,GQ_child,PGT_child,PID_child,PL_child,TP_child]:
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
