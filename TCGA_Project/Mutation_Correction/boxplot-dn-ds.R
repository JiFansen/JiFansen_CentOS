library(ggplot2)
library(ggthemes)
m=1
#for(File in list.files('/Share/home/JiFansen/JiFansen/Mutation-Correction/non-sys/', pattern = '*-final.results.txt')){
for(File in list.files('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/dN-dS/')){
  cancertype <- unlist(strsplit(File,'-'))[1]
  print(cancertype)
  tmp <- read.table(file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/dN-dS/', File), header = T, sep = '\t')
  #tmp <- read.table(file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/non-sys/', File),header = T, sep = '\t')
  tmp$cancertype <- cancertype
  if(m==1){
    total.genes <- tmp
    m <- m+1
  }else{
    total.genes <- rbind.data.frame(total.genes, tmp)
    m <- m+1
  }
}
sig.pathway.fraction <- read.csv(file = '/Share/home/JiFansen/pathways/Lymphocyte-Correction-sig-pathway.csv', 
                                 header = F, sep = '\t')
kegg.pathway <- read.csv(file = '/Share/home/JiFansen/Desktop/new-confidence-interval/kegg_gene_symbol.csv', header = F,
                         fill = T)
rownames(kegg.pathway) <- as.character(kegg.pathway[,1])
kegg.pathway <- kegg.pathway[,-1]
match.result <- match(rownames(kegg.pathway),as.character(sig.pathway.fraction[,1]))
big.index1 <- which(is.na(match.result)==FALSE)
d=1
for(i in big.index1){
  print(i)
  genelist <- unique(as.character(as.matrix(kegg.pathway[i,])))
  genelist <- genelist[-length(genelist)]
  common.genes <- intersect(as.character(total.genes$gene),genelist)
  newmatch.result <- match(as.character(total.genes$gene),common.genes)
  big.index <- which(is.na(newmatch.result)==FALSE)
  x <- total.genes[big.index,]
  x$pathways <- as.factor(rownames(kegg.pathway)[i])
  if(d==1){
    figure.data.frame <- x
    d=d+1
  }else{
    figure.data.frame <- rbind.data.frame(figure.data.frame,x)
    d=d+1
  }
}
figure.data.frame$log <- log10(figure.data.frame$dN.dS)
whole.background <- total.genes
whole.background$pathways <- 'wholegenome'
whole.background$pathways <- as.factor(whole.background$pathways)
whole.background$log <- log10(whole.background$dN.dS)
figure.data.frame <- rbind.data.frame(figure.data.frame, whole.background)
title <- "Lymphocyte Correction dN/dS Boxplot"
theme <- theme_set(theme_classic())
#theme <- theme_update(legend.position="top",legend.title=element_blank(), panel.grid.major.x=element_blank())
p <- ggplot(data=figure.data.frame, aes(x=pathways,y=log))
p <- p+geom_boxplot(aes(fill=pathways))+labs(title = title, x = '', y = 'log-transformed dN/dS values')
#p <- p+scale_fill_manual(values=c('deepskyblue3','darkmagenta','springgreen3','violetred2'))
p <- p+theme(panel.border = element_blank())
p <- p+theme(axis.text.x = element_blank())
p <- p+theme(plot.title = element_text(hjust = 0.5, size=20)) 
p <- p+theme(axis.title.y = element_text(size = 16))
p <- p+theme(legend.title = element_text(size=12, face="bold"))
p <- p+theme(legend.text = element_text(color="azure4", size = 10, face = "bold"))
p <- p+ylim(-0.1,0.1)
p

title <- "Lymphocyte Correction dN/dS Boxplot"
theme <- theme_set(theme_tufte())
p <- ggplot(data=figure.data.frame, aes(x=pathways,y=log))+scale_size(range = c(2,8))
p <- p+geom_tufteboxplot(outlier.colour="transparent",width = 16,colour='deepskyblue3')+labs(title = title, x = '', y = 'log-transformed dN/dS values')
p <- p+theme(panel.border = element_blank())
p <- p+theme(axis.text.x = element_text(angle=90, vjust=0.6, size=10))
p <- p+theme(axis.text.y = element_text(size = 15))
p <- p+theme(plot.title = element_text(hjust = 0.5, size=20)) 
p <- p+theme(axis.title.y = element_text(size = 16))
p <- p+theme(legend.text = element_blank())
p <- p+ylim(-0.1,0.1)
p <- p+geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")
p <- p+geom_hline(aes(yintercept=median(whole.background$log)), linetype="dashed")
p
pdf('Lym-dnds.pdf', width=20, height = 8.3)
print(p)
dev.off()
p <- p+geom_boxplot(outlier.colour = NULL, aes_string(colour="group", fill="group")) 
p <- p+geom_hline(aes(yintercept=0))
p <- p+geom_hline(aes(yintercept=median(whole.background$log),col = 'blue'), size = 1)


############################################################# According to cancer type.########################################################################

sig.pathway.fraction <- read.csv(file = '/Share/home/JiFansen/Desktop/new-confidence-interval/SigPathway_fraction_cutoff.csv', 
                                 header = T, sep = ',', fill = T)
sig.pathway.fraction <- sig.pathway.fraction[,-1]
sig.index <- which(as.character(sig.pathway.fraction[,4])=='T')
sig.pathway.fraction <- sig.pathway.fraction[sig.index,]
kegg.pathway <- read.csv(file = '/Share/home/JiFansen/Desktop/new-confidence-interval/kegg_gene_symbol.csv', header = F,
                         fill = T)
rownames(kegg.pathway) <- as.character(kegg.pathway[,1])
kegg.pathway <- kegg.pathway[,-1]
match.result <- match(rownames(kegg.pathway),as.character(sig.pathway.fraction[,1]))
big.index1 <- which(is.na(match.result)==FALSE)
for(File in list.files('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/dN-dS/')){
  cancertype <- unlist(strsplit(File,'-'))[1]
  tmp <- read.table(file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/dN-dS/', File), header = T, sep = '\t')
  tmp$cancertype <- cancertype
  d=1
  t <- 61
  l <- 65
  for(i in big.index1[t:l]){
    print(i)
    genelist <- unique(as.character(as.matrix(kegg.pathway[i,])))
    genelist <- genelist[-length(genelist)]
    common.genes <- intersect(as.character(tmp$gene),genelist)
    if(length(common.genes)>0){
    newmatch.result <- match(as.character(tmp$gene),common.genes)
    big.index <- which(is.na(newmatch.result)==FALSE)
    x <- tmp[big.index,]
    x$group <- as.factor(rownames(kegg.pathway)[i])
    if(d==1){
      figure.data.frame <- x
      d=d+1
    }else{
      figure.data.frame <- rbind.data.frame(figure.data.frame,x)
      d=d+1
    }
    }
  }
  figure.data.frame$log <- log10(figure.data.frame$dN.dS)
  
  #whole.background <- data.frame(gene = unique(as.character(total.genes$gene)), non.count = 1, sys.count = 1, non.background = 1, sys.background = 1,
  #                               cancertype = cancertype, dN.dS = 1)
  #match.result <- match(as.character(whole.background$gene), as.character(tmp$gene))
  #one.index <- which(is.na(match.result)==FALSE)
  #two.index <- match.result[one.index]
  #whole.background[one.index,] <- tmp[two.index,]
  
  whole.background <- tmp
  whole.background$group <- 'wholegenome'
  whole.background$group <- as.factor(whole.background$group)
  whole.background$log <- log10(whole.background$dN.dS)
  figure.data.frame <- rbind.data.frame(figure.data.frame, whole.background)
  p<-ggplot(data=figure.data.frame, aes(x=group,y=log))+geom_boxplot(aes(fill=group))
  p <- p+geom_hline(aes(yintercept=0))
  pdf(file = paste0('/Share/home/JiFansen/Desktop/test/', cancertype, t,'-',l, '.pdf'), width = 28, height = 13)
  print(p)
  dev.off()
}
save.image(file = '/Share/home/JiFansen/JiFansen/working-image/boxplot-dn-ds.RData')
