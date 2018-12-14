############################################## Read in the TCGA expression data. ##############################################################
#load('/Share/home/JiFansen/JiFansen/TCGA-Project/differentialAnalysis/working.RData')
library(clusterProfiler)
load('/Share/home/JiFansen/JiFansen/TCGA-Project/Expdata/TCGA_pancancer.noNA.logqq.RData')
#****************************************************************************************************************************
purity_leukocyte <- read.table(file = '/fshare2/JiFansen/Immune_LandScape/Purity_Leukocyte.txt', header = T, sep = '\t')
purity_leukocyte <- purity_leukocyte[,c(1,2,4,11,13)]
purity_leukocyte$Sum <- purity_leukocyte$purity+purity_leukocyte$LeukocyteRatio
purity_leukocyte <- purity_leukocyte[-which(purity_leukocyte$Sum>=1),]
purity_leukocyte$other <- 1 - purity_leukocyte$Sum
sampleName <- colnames(expression.table)
sampleName <- gsub("\\.","-",sampleName)
rawName <- sampleName
sampleName <- substr(sampleName, 1, 15)
match.result <- match(as.character(purity_leukocyte$array),sampleName)
big.index <- which(is.na(match.result)==FALSE)
small.index <- match.result[big.index]
purity_leukocyte <- purity_leukocyte[big.index,]
expression.table <- expression.table[,small.index]
rawName <- rawName[small.index]
cancertype <- as.character(purity_leukocyte$cancertype)
cancername <- unique(cancertype)
colnames(expression.table) <- rawName
#**********************************************************************************************************************************************
# Read the candidate genes.
candidate.genes <- read.csv(file = '/Share/home/JiFansen/pathways/TBCP_candidgenes.csv', header = T)
candidate.genes <- as.character(candidate.genes[,2])
#**********************************************************************************************************************************************
cluster_number <- 3
top.genes <- 40
candidate.genes.top20 <- candidate.genes[1:top.genes]
expression.candidate <- expression.table[match(candidate.genes.top20, rownames(expression.table)), ]
hc<-hclust(dist(t(expression.candidate)))
groups <- cutree(hc,k=cluster_number)
Groups <- groups
tmp <- as.data.frame(table(Groups))
tmp <- tmp[order(tmp$Freq), ]
library(NMI)
hc.pattern <- data.frame(c(1:length(Groups)), Groups)
library(SIMLR)
nmi.values <- c()
for(i in 1:10){
  example_large_scale = SIMLR_Large_Scale(X = expression.candidate, c = cluster_number, k = 5, kk=5)
  newtmp <- example_large_scale$y$cluster
  newtmp <- as.data.frame(table(newtmp))
  newtmp <- newtmp[order(newtmp$Freq), ]
  simlr.pattern <- data.frame(c(1:length(Groups)), example_large_scale$y$cluster)
  one.index <- which(simlr.pattern$example_large_scale.y.cluster==newtmp$newtmp[1])
  two.index <- which(simlr.pattern$example_large_scale.y.cluster==newtmp$newtmp[2])
  three.index <- which(simlr.pattern$example_large_scale.y.cluster==newtmp$newtmp[3])
  simlr.pattern$example_large_scale.y.cluster[one.index] <- tmp$Groups[1]
  simlr.pattern$example_large_scale.y.cluster[two.index] <- tmp$Groups[2]
  simlr.pattern$example_large_scale.y.cluster[three.index] <- tmp$Groups[3]
  nmi.values <- c(nmi.values, NMI(hc.pattern, simlr.pattern)$value)
}
NMI(hc.pattern, simlr.pattern)$value
################################################## Read the 1 year survival of TCGA. #######################################################
clinical <- read.csv(file = '/fshare2/JiFansen/Clinical/three_year_data.csv', header = T)
clinical <- clinical[,c(1,2,3,10,12,13)]
clinical$X.1 <- gsub('\\.','-',as.character(clinical$X.1))
index <- match(clinical$X.1, colnames(expression.table))
immuno.clinical <- clinical[which(is.na(index)==FALSE), ]
index <- index[which(is.na(index)==FALSE)]
immuno.expression <- expression.table[, index]
immuno.genes <- bitr(rownames(immuno.expression), fromType = 'SYMBOL', toType ="ENTREZID", OrgDb="org.Hs.eg.db")
immuno.expression <- immuno.expression[immuno.genes$SYMBOL,]
immuno.pvalues <- apply(immuno.expression,1,FUN = function(x){
  return(wilcox.test(as.numeric(x)~as.factor(immuno.clinical$vital_status), alternative = "less")$p.value)
})
fold.change <- apply(immuno.expression, 1, FUN = function(x){
  zero.samples <- median(x[which(immuno.clinical$vital_status==0)])
  one.samples <- median(x[which(immuno.clinical$vital_status==1)])
  return(zero.samples/one.samples)
})
geneList <- data.frame(gene = immuno.genes$SYMBOL, entrez = immuno.genes$ENTREZID, p.value = immuno.pvalues, fold.change = fold.change)
order.index <- order(geneList$p.value, decreasing = FALSE)
geneList <- geneList[order.index, ]
geneList$rank <- order(c(1:dim(geneList)[1]), decreasing = T)
d <- geneList$rank
names(d) <- as.character(geneList$entrez)
kk2 <- gseKEGG(geneList = d, 
               organism = 'hsa',
               nPerm = 1000,
               minGSSize = 50,
               pvalueCutoff = 1,
               verbose = FALSE)
write.table(kk2[,1:9], file = '/Share/home/JiFansen/JiFansen/TCGA-Project/differentialAnalysis/GSEA-lymphocyte-correction-less.txt', row.names = F, sep = '\t', quote = F)
################################################################################################################################################################################
nonimmuno.clinical <- read.csv(file = '/fshare2/JiFansen/Clinical/40candidatesgpredict_nonimmune_lable.csv', header = T)
nonimmuno.clinical$X <- gsub('\\.', '-', as.character(nonimmuno.clinical$X))
index <- match(nonimmuno.clinical$X, colnames(expression.table))
nonimmuno.clinical <- nonimmuno.clinical[which(is.na(index)==FALSE), ]
index <- index[which(is.na(index)==FALSE)]
nonimmuno.expression <- expression.table[, index]
nonimmuno.genes <- bitr(rownames(nonimmuno.expression), fromType = 'SYMBOL', toType ="ENTREZID", OrgDb="org.Hs.eg.db")
nonimmuno.expression <- nonimmuno.expression[immuno.genes$SYMBOL,]
nonimmuno.pvalues <- apply(nonimmuno.expression,1,FUN = function(x){
  return(wilcox.test(as.numeric(x)~as.factor(nonimmuno.clinical$lab), alternative = "less")$p.value)
})
fold.change <- apply(nonimmuno.expression, 1, FUN = function(x){
  zero.samples <- median(x[which(nonimmuno.clinical$lab==0)])
  one.samples <- median(x[which(nonimmuno.clinical$lab==1)])
  return(zero.samples/one.samples)
})
geneList <- data.frame(gene = nonimmuno.genes$SYMBOL, entrez = nonimmuno.genes$ENTREZID, p.value = nonimmuno.pvalues, fold.change = fold.change)
order.index <- order(geneList$p.value, decreasing = FALSE)
geneList <- geneList[order.index, ]
geneList$rank <- order(c(1:dim(geneList)[1]), decreasing = T)
d <- geneList$rank
names(d) <- as.character(geneList$entrez)
kk2 <- gseKEGG(geneList = d, 
               organism = 'hsa',
               nPerm = 1000,
               minGSSize = 50,
               pvalueCutoff = 1,
               verbose = FALSE)
write.table(kk2[,1:9], file = '/Share/home/JiFansen/JiFansen/TCGA-Project/differentialAnalysis/no-label-GSEA-mutation-correction-less.txt', row.names = F, sep = '\t', quote = F)

















####################################################################################################################################################################################
nonimmuno.palues <- apply(nonimmuno.expression, 1, FUN = function(x){
  return(t.test(as.numeric(x)~as.factor(nonimmuno.clinical$lab))$p.value)
})
nonimmuno.palues.adjust <- p.adjust(nonimmuno.palues, method = 'BH')
nonimmuno.genes.different <- names(nonimmuno.palues.adjust[which(nonimmuno.palues.adjust<0.05)])
nonimmuno.expression <- nonimmuno.expression[nonimmuno.genes.different, ]
fold.change <- apply(nonimmuno.expression, 1, FUN = function(x){
  zero.samples <- mean(x[which(nonimmuno.clinical$lab==0)])
  one.samples <- mean(x[which(nonimmuno.clinical$lab==1)])
  return(zero.samples/one.samples)
})
nonimmuno.genes.different <- nonimmuno.genes.different[-which(fold.change>0.5&fold.change<2)]
write.table(nonimmuno.genes.different, file = '/Share/home/JiFansen/JiFansen/TCGA-Project/differentialAnalysis/nonimmuno.different.genes', sep = '\t', quote = FALSE, col.names = F, row.names = F)
#########################################################################################################################################
#immuno.pvalues.adjust <- p.adjust(immuno.pvalues, method = 'BH')
#immuno.genes.different <- names(immuno.pvalues.adjust[which(immuno.pvalues.adjust<=0.1)])
#geneList <- geneList[which(abs(log2(geneList$fold.change))>1),]
gseaplot(kk2, geneSetID = "hsa03008")
enriched.genes <- as.character(kk2[1,-1])[10]
enriched.genes <- unlist(strsplit(enriched.genes,'/'))
enriched.genes <- as.numeric(enriched.genes)
enriched.genes <- bitr(enriched.genes, fromType = 'ENTREZID', toType = c("ENSEMBL", "SYMBOL"), OrgDb="org.Hs.eg.db")
candidate.genes <- read.csv(file = '/Share/home/JiFansen/pathways/TBCP_candidgenes.csv', header = T)
candidate.genes <- as.character(candidate.genes[,2])
length(intersect(enriched.genes$SYMBOL, candidate.genes))

length(which(immuno.pvalues.adjust[match(candidate.genes, names(immuno.pvalues.adjust))]<0.05))
test <- immuno.pvalues.adjust[match(candidate.genes, names(immuno.pvalues.adjust))]
test[test<0.05]

immuno.genes.different <- immuno.genes.different[which(abs(log2(fold.change))>log2(2))]
write.table(immuno.genes.different, file = '/Share/home/JiFansen/JiFansen/TCGA-Project/differentialAnalysis/immuno.different.genes', sep = '\t', quote = FALSE, col.names = F, row.names = F)
############################################################################################################################################################################################

diff.genes <- read.table(file = '/Share/home/JiFansen/JiFansen/TCGA-Project/differentialAnalysis/nonimmuno.different.genes', sep = '\t', header = F)
diff.genes <- as.character(diff.genes[,1])
expression.table <- read.table(file = '/fshare2/sharedData/TCGA/RNA/TCGA_pancancer.noNA.logqq.exp', header = T, sep = '\t',row.names = 1)
purity_leukocyte <- read.table(file = '/fshare2/JiFansen/Immune_LandScape/Purity_Leukocyte.txt', header = T, sep = '\t')
purity_leukocyte <- purity_leukocyte[,c(1,2,4,11,13)]
purity_leukocyte$Sum <- purity_leukocyte$purity+purity_leukocyte$LeukocyteRatio
purity_leukocyte <- purity_leukocyte[-which(purity_leukocyte$Sum>=1),]
purity_leukocyte$other <- 1 - purity_leukocyte$Sum
sampleName <- colnames(expression.table)
sampleName <- gsub("\\.","-",sampleName)
rawName <- sampleName
sampleName <- substr(sampleName, 1, 15)
match.result <- match(as.character(purity_leukocyte$array),sampleName)
big.index <- which(is.na(match.result)==FALSE)
small.index <- match.result[big.index]
purity_leukocyte <- purity_leukocyte[big.index,]
expression.table <- expression.table[,small.index]
rawName <- rawName[small.index]
cancertype <- as.character(purity_leukocyte$cancertype)
cancername <- unique(cancertype)
colnames(expression.table) <- rawName
gene.index <- match(diff.genes, rownames(expression.table))
expression.diff <- expression.table[gene.index, ]
library('destiny')
library('readxl')
library('Biobase')
library('ggplot2')
expression.diff <- as.data.frame(t(expression.diff))
expression.diff <- as.ExpressionSet(expression.diff)
dm <- DiffusionMap(expression.diff, k = find_dm_k(nrow(expression.diff)-1))
myData <- as.data.frame(dm)
for(i in 1:3){
  test <- c()
  for(j in 21:(dim(myData)[2])){
    test <- c(test, cor(myData[,i], myData[,j], method = 'spearman'))
  }
  test <- as.matrix(test)
  if(i == 1){
    cor.matrix <- test
  }else{
    cor.matrix <- cbind(cor.matrix, test)
  }
}
rownames(cor.matrix) <- rownames(expression.diff)
colnames(cor.matrix) <- c('DC1', 'DC2', 'DC3')
write.table(cor.matrix, file = 'cor.matrix', sep = '\t', quote = FALSE)
DC1.positive <- rownames(cor.matrix)[which(cor.matrix[,1]>0.2)]
DC1.negative <- rownames(cor.matrix)[which(cor.matrix[,1]<(-0.2))]
DC2.positive <- rownames(cor.matrix)[which(cor.matrix[,2]>0.2)]
DC2.negative <- rownames(cor.matrix)[which(cor.matrix[,2]<(-0.2))]
DC3.positive <- rownames(cor.matrix)[which(cor.matrix[,3]>0.2)]
DC3.negative <- rownames(cor.matrix)[which(cor.matrix[,3]<(-0.2))]
write.table(DC1.positive, file = 'DC1.positive', sep = '\t', quote = FALSE, row.names = F, col.names = F)



theme <- theme_set(theme_classic())
qplot(DC1, DC2, data = dm, colour = factor(cancertype))+theme(panel.border = element_blank())
palette(cube_helix(6))
for(i in 21:75){
  pdf(paste0(names(dm)[i], '.pdf'))
  plot(dm, pch = 20, col_by = names(dm)[i], legend_main = names(dm)[i])
  dev.off()
}

palette(colors()[index])
plot(dm, pch = 20, col_by = 'cancertype', legend_main = 'Cancer Type')

library(org.Hs.eg.db)
hs <- org.Hs.eg.db
hc <- select(hs, 
       keys = DC1.positive,
       columns = c("ENTREZID", "SYMBOL"),
       keytype = "SYMBOL")
hc <- hc[-which(is.na(hc$ENTREZID)), ]
ego <- enrichGO(gene = hc$ENTREZID, organism = 'human', pvalueCutoff = 0.05, readable = TRUE)
ego <- enrichKEGG(gene = hc$ENTREZID, organism = 'human', pvalueCutoff = 0.05, readable = TRUE)
#https://zhuanlan.zhihu.com/p/35510434
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
immuno.genes <- read.table(file = 'nonimmuno.different.genes', header = F, sep = '\t')
immuno.genes <- as.character(immuno.genes[,1])
immuno.genes <- immuno.genes.different
immuno.genes <- bitr(immuno.genes, fromType = 'SYMBOL', toType = c("ENSEMBL", "ENTREZID"), OrgDb="org.Hs.eg.db")
ego_ALL <- enrichGO(gene = immuno.genes$ENTREZID, 
                    #universe = names(geneList), #背景基因集
                    OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                    #keytype = 'ENSEMBL',
                    ont = "ALL", #也可以是 CC  BP  MF中的一种
                    pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                    pvalueCutoff = 0.05, #P值会过滤掉很多，可以全部输出
                    qvalueCutoff = 0.05,
                    readable = TRUE) #Gene ID 转成gene Symbol ，易读
head(ego_ALL)
#其中：ONTOLOGY：CC  BP  MF 
#GO ID: Gene Ontology数据库中唯一的标号信息
#Description ：Gene Ontology功能的描述信息
#GeneRatio：差异基因中与该Term相关的基因数与整个差异基因总数的比值
#BgRation：所有（ bg）基因中与该Term相关的基因数与所有（ bg）基因的比值
#pvalue: 富集分析统计学显著水平，一般情况下， P-value < 0.05 该功能为富集项
#p.adjust 矫正后的P-Value
#qvalue：对p值进行统计学检验的q值
#geneID：与该Term相关的基因
#Count：与该Term相关的基因数
dotplot(ego_ALL,title="EnrichmentGO_ALL")#点图，按富集的数从大到小的
kk <- enrichKEGG(gene = immuno.genes$ENTREZID,
                 organism = 'hsa', #KEGG可以用organism = 'hsa'
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH")
dotplot(kk,title="Enrichment_KEGG_dot")#点图，按富集的数从大到小的

##可视化--条形图
barplot(kk, showCategory=20,title="EnrichmentGO_MF")#条状图，按p从小到大排，绘制前20个Term
############################################## Build a cutoff function. #######################################################################
# data: the survival table.
# colTime: the colnames index of survival time.
# colStatus: the colnames index of survival status. The survival status should be 0(1) values for which 0 represents alive and 1 represents dead.
# cutoff: survival cutoff.
cutoffIndex <- function(data = data, colTime = colTime, colStatus = colStatus, cutoff=365){
  test <- data
  censor.index <- which(test[,colTime]<cutoff&test[,colStatus]==0)
  if(length(censor.index)>0){
    return(censor.index)
  }
  
}
survivalCutoff <- function(data = data, colTime = colTime, colStatus = colStatus, cutoff=365){
  test <- data
  censor.index <- which(test[,colTime]<cutoff&test[,colStatus]==0)
  test <- test[-censor.index, ]
  test[,colStatus][which(test[,colTime]>=cutoff)] <- 0
  test[,colTime][which(test[, colTime]>=cutoff)] <- cutoff
  return(test)
}

############################################## Read in the TCGA expression data. ##############################################################
library('preprocessCore')
################################## Read the independent test expression and clinical data. #######################################################
dataset <- read.csv(file = '/Share/home/JiFansen/JiFansen/Clinical/training-TCGA-GSE-bms-Exp.csv', row.names = 1)
survival <- read.csv(file = '/Share/home/JiFansen/JiFansen/Clinical/TCGA-GSE-bms-Clinical.csv')
dataset_42 <- read.table(file = '/Share/home/JiFansen/Machine_learning/MEL_RPKM.txt', header = T, sep = '\t')
survival_42 <- read.table(file = '/Share/home/JiFansen/Machine_learning/sample42_Clinical.txt', header = T, sep = '\t')
name1 <- rownames(dataset_42)
name2 <- colnames(dataset_42)
dataset_42 <- normalize.quantiles(as.matrix(dataset_42), copy = T)
rownames(dataset_42) <- name1
colnames(dataset_42) <- name2
##################################################################################################################################################
############################################### Split the independent test data. #################################################################
dataset.gse <- dataset[, 165:190]
dataset.bms <- dataset[, 191:241]
survival.gse <- survival[165:190, ]
survival.bms <- survival[191:241, ]
name1 <- rownames(dataset.bms)
name2 <- colnames(dataset.bms)
dataset.bms <- normalize.quantiles(as.matrix(dataset.bms), copy = T)
rownames(dataset.bms) <- name1
colnames(dataset.bms) <- name2
name1 <- rownames(dataset.gse)
name2 <- colnames(dataset.gse)
dataset.gse <- normalize.quantiles(as.matrix(dataset.gse), copy = T)
rownames(dataset.gse) <- name1
colnames(dataset.gse) <- name2
dataset.gse <- dataset.gse[, -cutoffIndex(data = survival.gse, colTime = 4, colStatus = 5, cutoff = 365)]
survival.gse <- survivalCutoff(data = survival.gse, colTime = 4, colStatus = 5, cutoff = 365)
dataset.bms <- dataset.bms[, -cutoffIndex(data = survival.bms, colTime = 4, colStatus = 5, cutoff = 365)]
survival.bms <- survivalCutoff(data = survival.bms, colTime = 4, colStatus = 5, cutoff = 365)
if(length(cutoffIndex(data = survival_42, colTime = 2, colStatus = 4, cutoff = 365))>0){
  dataset_42 <- dataset_42[, -cutoffIndex(data = survival_42, colTime = 2, colStatus = 4, cutoff = 365)]
  survival_42 <- survivalCutoff(data = survival_42, colTime = 2, colStatus = 4, cutoff = 365)
}
##################################################################################################################
immuno.pvalues <- apply(dataset_42,1,FUN = function(x){
  return(t.test(as.numeric(x)~as.factor(survival_42$Status))$p.value)
})
immuno.pvalues.adjust <- p.adjust(immuno.pvalues, method = 'BH')
immuno.genes.different <- names(immuno.pvalues.adjust[which(immuno.pvalues.adjust<=0.05)])
immuno.expression <- immuno.expression[immuno.genes.different, ]
fold.change <- apply(immuno.expression, 1, FUN = function(x){
  zero.samples <- mean(x[which(immuno.clinical$vital_status==0)])
  one.samples <- mean(x[which(immuno.clinical$vital_status==1)])
  return(zero.samples/one.samples)
})
immuno.genes.different <- immuno.genes.different[-which(fold.change>0.5&fold.change<2)]
write.table(immuno.genes.different, file = '/Share/home/JiFansen/JiFansen/TCGA-Project/differentialAnalysis/immuno.different.genes', sep = '\t', quote = FALSE, col.names = F, row.names = F)