differential.p <- t(apply(immuno.expression, 1, FUN = function(x){
  test <- cor.test(as.numeric(x), as.numeric(immuno.mutation), method = 'spearman', exact = FALSE)
  return(c(test$p.value, test$estimate))
}))
differential.p <- as.data.frame(differential.p)
differential.p <- differential.p[which(abs(differential.p$rho)>0.3), ]
differential.p <- differential.p[which(differential.p$V1<0.05), ]
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
immuno.genes <- rownames(differential.p)
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
#######################################################################################################################







data.summary.new <- data.summary[which(data.summary$mutation.count<=100), ]
wilcox.test(as.numeric(data.summary.new$mutation.count)~as.factor(data.summary.new$group), exact = FALSE)$p.value
differential.p <- t(apply(nonimmuno.expression, 1, FUN = function(x){
  test <- cor.test(as.numeric(x), as.numeric(nonimmuno.mutation), method = 'spearman', exact = FALSE)
  return(c(test$p.value, test$estimate))
}))
differential.p <- as.data.frame(differential.p)
differential.p <- differential.p[which(abs(differential.p$rho)>0.3), ]
differential.p <- differential.p[which(differential.p$V1<0.05), ]
#############################################################################################################################
mutation.load <- read.table(file = '/Share/home/JiFansen/JiFansen/TCGA-Project/mutation_load/mutation-load_updated.txt', header = T, sep = '\t')
clinical <- read.csv(file = '/fshare2/JiFansen/Clinical/three_year_data.csv', header = T)
clinical <- clinical[,c(1,2,3,10,12,13)]
clinical$X.1 <- gsub('\\.','-',as.character(clinical$X.1))
index <- match(substr(clinical$X.1, 1,15), mutation.load$Tumor_Sample_ID)
immuno.clinical <- clinical[which(is.na(index)==FALSE), ]
index <- index[which(is.na(index)==FALSE)]
immuno.mutation <- mutation.load$Non.silent.per.Mb[index]
data.summary <- data.frame(mutation.count = as.numeric(immuno.mutation), group = as.factor(immuno.clinical$vital_status))
p <- ggplot(data.summary, aes(x=group, y=mutation.count, color=group)) + geom_violin(trim=FALSE)
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.3)+ theme_classic()+ylim(c(0,50))
aggregate(x=data.summary[,1], by = list(data.summary$group), FUN = median)
wilcox.test(as.numeric(data.summary$mutation.count)~as.factor(data.summary$group), exact = FALSE)$p.value
##################################################################################################################################
nonimmuno.clinical <- read.csv(file = '/fshare2/JiFansen/Clinical/40candidatesgpredict_nonimmune_lable.csv', header = T)
nonimmuno.clinical$X <- gsub('\\.', '-', as.character(nonimmuno.clinical$X))
index <- match(substr(nonimmuno.clinical$X, 1,15), mutation.load$Tumor_Sample_ID)
nonimmuno.clinical <- nonimmuno.clinical[which(is.na(index)==FALSE), ]
index <- index[which(is.na(index)==FALSE)]
nonimmuno.mutation <- mutation.load$Non.silent.per.Mb[index]
data.summary <- data.frame(mutation.count = as.numeric(nonimmuno.mutation), group = as.factor(nonimmuno.clinical$lab))
p <- ggplot(data.summary, aes(x=group, y=mutation.count, color=group)) + geom_violin(trim=FALSE)
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1)+ theme_classic()
theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) +
  theme(legend.position="none")
aggregate(x=data.summary[,1], by = list(data.summary$group), FUN = median)














# Convert the variable dose from a numeric to a factor variable
ToothGrowth$dose <- as.factor(ToothGrowth$dose)
head(ToothGrowth)
library(ggplot2)
p <- ggplot(ToothGrowth, aes(x=dose, y=len)) + geom_violin()
# Rotate the violin plot
p + coord_flip()
# Set trim argument to FALSE
ggplot(ToothGrowth, aes(x=dose, y=len)) + geom_violin(trim=FALSE)
p + scale_x_discrete(limits=c("0.5", "2"))
# violin plot with mean points
p + stat_summary(fun.y=mean, geom="point", shape=23, size=2)
# violin plot with median points
p + stat_summary(fun.y=median, geom="point", size=2, color="red")
p + geom_boxplot(width=0.1)
p <- ggplot(ToothGrowth, aes(x=dose, y=len)) + geom_violin(trim=FALSE)
p + stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.2 )
p + stat_summary(fun.data=mean_sdl, mult=1, geom="pointrange", color="red")
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
p + stat_summary(fun.data=data_summary)
# violin plot with dot plot
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
# violin plot with jittered points
# 0.2 : degree of jitter in x direction
p + geom_jitter(shape=16, position=position_jitter(0.2))
# Change violin plot line colors by groups
p<-ggplot(ToothGrowth, aes(x=dose, y=len, color=dose)) +
  geom_violin(trim=FALSE)
p
# Use custom color palettes
p+scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))
# Use brewer color palettes
p+scale_color_brewer(palette="Dark2")
# Use grey scale
p + scale_color_grey() + theme_classic()
# Use single color
ggplot(ToothGrowth, aes(x=dose, y=len)) +
  geom_violin(trim=FALSE, fill='#A4A4A4', color="darkred")+
  geom_boxplot(width=0.1) + theme_minimal()
# Change violin plot colors by groups
p<-ggplot(ToothGrowth, aes(x=dose, y=len, fill=dose)) +
  geom_violin(trim=FALSE)
p
# Use custom color palettes
p+scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))
# Use brewer color palettes
p+scale_fill_brewer(palette="Dark2")
# Use grey scale
p + scale_fill_grey() + theme_classic()




