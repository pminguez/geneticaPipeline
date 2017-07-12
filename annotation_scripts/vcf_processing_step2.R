args=commandArgs(TRUE)
vcfFile = read.csv(args[1],header = T, sep = '\t',stringsAsFactors = F)
vcfFile[,1] = gsub("\\[|\\]",'',vcfFile[,1])
colnames(vcfFile)[1]= 'Chr' 


vcfFile$PL. = gsub("\\[|\\]",'',vcfFile$PL.)
vcfFile$X. = NULL

for (i in c(1:ncol(vcfFile))){
  vcfFile[,i] = trimws(vcfFile[,i])
}

vcfFile = subset(vcfFile, vcfFile$Chr != 'chrM')
vcfFile = subset(vcfFile, vcfFile$Func.refGene == 'exonic' | 
                   vcfFile$Func.refGene == 'intronic' |
                   vcfFile$Func.refGene == 'UTR3' |
                   vcfFile$Func.refGene == 'UTR5'|
                   vcfFile$Func.refGene == 'splicing')

vcfFile = subset(vcfFile, as.numeric(vcfFile$DP) > 10)

GDI = read.table('/mnt/datos1/GeneticaPipeDB/genome_data/GDI_full_10282015.txt', sep = '\t', header = T)
GDI = GDI[,c(1:7)]
colnames(GDI)[1] <- "Gene.refGene"
vcfFile = merge(vcfFile,GDI,by.x = "Gene.refGene",all.x = TRUE)

LoFtool = read.table('/mnt/datos1/GeneticaPipeDB/genome_data/LoFtool_scores.txt', sep = '\t', header = T)
colnames(LoFtool)[1] <- "Gene.refGene"
vcfFile = merge(vcfFile,LoFtool,by.x = "Gene.refGene",all.x = TRUE)

RVIS = read.table('/mnt/datos1/GeneticaPipeDB/genome_data/RVIS_ExAC_4KW.txt', sep = '\t', header = T)
RVIS = RVIS[,c(5:7)]
colnames(RVIS) = c("Gene.refGene","RVIS_score","RVIS_percentile")
vcfFile = merge(vcfFile,RVIS,by.x = "Gene.refGene",all.x = TRUE)

Phenotype = read.table('/mnt/datos1/GeneticaPipeDB/genome_data/hgnc_genes_phenotypes.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
colnames(Phenotype) = c("Gene.refGene","Phenotype_description")
vcfFile = merge(vcfFile,Phenotype,by.x = "Gene.refGene",all.x = TRUE)

#Introducing NA values 
write.table(vcfFile, paste(args[1],'_raw_variants.txt',sep = ''), sep="\t",quote = F, row.names = F)
vcfFile[vcfFile == '.'] <- NA

#Getting rid of low QUAL variants
vcfFile = subset(vcfFile, as.numeric(vcfFile$QUAL) > 100)
vcfFile = subset(vcfFile, as.numeric(vcfFile$QD) > 2.0)
vcfFile = subset(vcfFile, as.numeric(vcfFile$MQ) > 40 | is.na(vcfFile$MQ))

#Getting rid of variants with strand bias (FS)
vcfFile = subset(vcfFile, as.numeric(vcfFile$FS) < 65 | is.na(vcfFile$FS))

#Getting rid of variants near the end of reads
vcfFile = subset(vcfFile, as.numeric(vcfFile$ReadPosRankSum) > as.numeric(-8.0) | is.na(vcfFile$ReadPosRankSum))

#Getting rid of high frequency variants
vcfFile = subset(vcfFile, as.numeric(vcfFile$ExAC_EAS) < 0.01 | is.na(vcfFile$ExAC_EAS))
vcfFile = subset(vcfFile, as.numeric(vcfFile$Kaviar_AF) < 0.01 | is.na(vcfFile$Kaviar_AF))
vcfFile = subset(vcfFile, as.numeric(vcfFile$X1000g2015aug_eur) < 0.01 | is.na(vcfFile$X1000g2015aug_eur))
vcfFile = subset(vcfFile, as.numeric(vcfFile$gnomAD_exome_ALL) < 0.01 | is.na(vcfFile$gnomAD_exome_ALL))
vcfFile = subset(vcfFile, as.numeric(vcfFile$gnomAD_genome_ALL) < 0.01 | is.na(vcfFile$gnomAD_genome_ALL))
vcfFile$AltFreq = as.numeric(vcfFile$ALTD) / as.numeric(vcfFile$DP)*100
vcfFile = subset(vcfFile, as.numeric(vcfFile$AltFreq) > 10 | is.na(vcfFile$AltFreq))


write.table(vcfFile, paste(args[1],'_prefiltered.txt',sep = ''), sep="\t",quote = F, row.names = F)
