library(ggplot2)
library(plotrix)
pathway_name <- 'lym_correction_c_type.txt'
# Read the t-cell receptor pathway genes.
t_cell_genes <- read.table(file = paste0('/fshare2/JiFansen/Mutation-Correction/ThreePathways/', pathway_name),
                           header = T, sep = '\t')
t_cell_genes <- as.character(t_cell_genes$gene)
# For each t-cell genes, find the dnds values in each cancer type.
t_cell_genes <- read.csv(file = '/fshare2/JiFansen/Mutation-Correction/ThreePathways/Phosphatidylinositol_signaling_system_candidgenes.csv', header = T)
t_cell_genes <- as.character(t_cell_genes$x)
w <- 1
for(genes in t_cell_genes){
  dn_ds_vector <- c()
  cancertype <- c()
  for(file in list.files('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/GeneCount-Background')){
    cancername <- unlist(strsplit(file,"-"))[1]
    value_table <- read.table(file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/GeneCount-Background/', file), header = T, sep = '\t')
    genes_index <- which(value_table$gene==genes)
    if(length(genes_index)>0){
      dn_ds_vector <- c(dn_ds_vector, ((value_table$non.count[genes_index]/value_table$non.background[genes_index])+1)/((value_table$sys.count[genes_index]/value_table$sys.background[genes_index])+1))
      cancertype <- c(cancertype, cancername)
    }
  }
  gene_summary <- data.frame(gene = genes, dN.dS = dn_ds_vector, cancertype = cancertype)
  if(w==1){
    final_table <- gene_summary
    w <- w+1
  }else{
    final_table <- rbind.data.frame(final_table, gene_summary)
    w <- w+1
  }
}
#df <- list()
#for(genename in unique(final_table$gene)){
#  df[[genename]] <- final_table$dN.dS[which(final_table$gene==genename)]
#}
#boxplot(df)
theme <- theme_set(theme_classic())
p <- ggplot()
p <- p+geom_boxplot(data = final_table, aes(x = gene,y = log10(dN.dS)))+labs(
  title = 'Phosphatidylinositol signalling system pathway candidate genes', x = 'Gene', y = expression(log[10](dN.dS)))
#p <- p+geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")
#p <- p+geom_hline(aes(yintercept=log10(1.015027)), colour="blue", linetype="dashed")
p <- p+theme(panel.border = element_blank())
p <- p+theme(axis.text.x = element_text(angle=60, vjust=0.6, size=20))
p <- p+theme(axis.text.y = element_text(size = 20))
p <- p+theme(plot.title = element_text(hjust = 0.6, size=30)) 
p <- p+theme(axis.title.y = element_text(size = 16))
p <- p+theme(axis.title.x = element_text(size = 20))
p <- p+theme(legend.text = element_text(size = 20))
p <- p+theme(legend.title = element_blank())
#p <- p+theme(legend.text = element_blank())
a <- data.frame(value = c(0,log10(1.015027)), group = c('dN/dS = 1', 'WholeGenome dN/dS'))
p <- p+geom_hline(aes(yintercept=value, col = group), data = a, size = 0.8)
p <- p+scale_color_manual(values=c("#F8766D", "#00BA38"))
p <- p+coord_cartesian(ylim=c(-0.2,0.4))

#pdf(file = '/Share/home/JiFansen/JiFansen/Mutation-Correction/ThreePathways/b_cells_genes.pdf',width = 28, height = 13)
pdf(file = '/Share/home/JiFansen/PaperFigure/Phosphatidylinositol_genes.pdf',width = 16, height = 15)
print(p)
dev.off()
