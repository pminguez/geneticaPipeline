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
                   vcfFile$Func.refGene == 'splicing')

vcfFile = subset(vcfFile, as.numeric(vcfFile$DP) > 10)

GDI = read.table('/mnt/genetica/domin/GeneticaPipeDB/genome_data/GDI_full_10282015.txt', sep = '\t', header = T)
GDI = GDI[,c(1:7)]
colnames(GDI)[1] <- "Gene.refGene"
vcfFile = merge(vcfFile,GDI,by.x = "Gene.refGene",all.x = TRUE)

LoFtool = read.table('/mnt/genetica/domin/GeneticaPipeDB/genome_data/LoFtool_scores.txt', sep = '\t', header = T)
colnames(LoFtool)[1] <- "Gene.refGene"
vcfFile = merge(vcfFile,LoFtool,by.x = "Gene.refGene",all.x = TRUE)

RVIS = read.table('/mnt/genetica/domin/GeneticaPipeDB/genome_data/RVIS_ExAC_4KW.txt', sep = '\t', header = T)
RVIS = RVIS[,c(5:7)]
colnames(RVIS) = c("Gene.refGene","RVIS_score","RVIS_percentile")
vcfFile = merge(vcfFile,RVIS,by.x = "Gene.refGene",all.x = TRUE)

Phenotype = read.table('/mnt/genetica/domin/GeneticaPipeDB/genome_data/hgnc_genes_phenotypes.txt', sep = '\t', header = T, stringsAsFactors = F, quote = '')
colnames(Phenotype) = c("Gene.refGene","Phenotype_description")
vcfFile = merge(vcfFile,Phenotype,by.x = "Gene.refGene",all.x = TRUE)


write.table(vcfFile, paste(args[1],'_raw_variants.txt',sep = ''), sep="\t",quote = F, row.names = F)

#Introducing NA values 


