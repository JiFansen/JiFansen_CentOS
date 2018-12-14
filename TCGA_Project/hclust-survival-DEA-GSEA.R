library("igraph");library("ggplot2");library("Rtsne");library("ggpubr");library("metap");library('NbClust');library('gridExtra');
library("factoextra");library('preprocessCore');library('GSVA');library('doMC');library('survival');library('survminer');
library('sfsmisc');library('SIMLR');library(gplots);library(scales);library(reshape);library('caret');library(pheatmap)
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
load('/Share/home/JiFansen/JiFansen/TCGA-Project/Expdata/TCGA_pancancer.noNA.logqq.RData')
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
##################################################################################################################################################
# Read the candidate genes.
candidate.genes <- read.csv(file = '/Share/home/JiFansen/pathways/fin_TBCP_RFselection.csv', header = T)
candidate.genes <- as.character(candidate.genes[,1])
###################################################################################################################################################
###################################################### All non-immunotherapy patients' survival difference. ###########################################################
for(cluster_number in c(2,3)){
  for(top.genes in c(40)){
    candidate.genes.top20 <- candidate.genes[1:top.genes]
    expression.candidate <- expression.table[match(candidate.genes.top20, rownames(expression.table)), ]
    hc<-hclust(dist(t(expression.candidate)))
    groups <- cutree(hc,k=cluster_number)
    Groups <- groups
    all.non.immune.clinical <- read.csv(file = '/fshare2/JiFansen/Clinical/ALLnon-immuneSurvivalthree_year_data_matchexp.csv', header = T, sep = ',')
    all.non.immune.clinical <- all.non.immune.clinical[,-c(2,3,4,5,6,7,9)]
    all.non.immune.clinical$sample_barcode <- gsub('\\.', '-', all.non.immune.clinical$sample_barcode)
    match.result <- match(all.non.immune.clinical$sample_barcode, names(Groups))
    all.non.immune.clinical$group <- Groups[match.result]
    all.non.immune.clinical <- all.non.immune.clinical[-which(is.na(all.non.immune.clinical$group)), ]
    surv_object <- Surv(time = all.non.immune.clinical$times, event = all.non.immune.clinical$vital_status)
    fit1 <- survdiff(surv_object ~ group, data = all.non.immune.clinical)
    p.fit1 <- 1-pchisq(fit1$chisq, length(fit1$n)-1)
    fit1 <- survfit(surv_object ~ group, data = all.non.immune.clinical)
    p.survival <- ggsurvplot(fit1, data = all.non.immune.clinical, pval = TRUE)
    pdf(paste0('/Share/home/JiFansen/TCGA-', cluster_number, 
               '-', top.genes,'-non-immunotherapy-survival.pdf'))
    print(p.survival$plot)
    dev.off()
  }
}
################################################### Clustering Analysis. #########################################################################
for(cluster_number in c(2, 3)){
  for(top.genes in c(40)){
    SampleRecord <- list()
    candidate.genes.top20 <- candidate.genes[1:top.genes]
    ############################################# Plot TCGA data. ###############################################################
    expression.candidate <- expression.table[match(candidate.genes.top20, rownames(expression.table)), ]
    hc<-hclust(dist(t(expression.candidate)))
    groups <- cutree(hc,k=cluster_number)
    Groups <- groups
    ############################################ Calculate centroid. ############################################################
    for(h in 1:cluster_number){
      expression.subset <- expression.candidate[, which(groups==h)]
      expression.subset <- t(t(apply(expression.subset, 1, mean)))
      if(h==1){
        centriod <- expression.subset
      }else{
        centriod <- cbind(centriod, expression.subset)
      }
    }
    colnames(centriod) <- as.character(c(1:cluster_number))
    ############################################################### TCGA survival analysis. ##########################################################
    clinical <- read.csv(file = '/fshare2/JiFansen/Clinical/three_year_data.csv', header = T, sep = ',')
    clinical <- clinical[,c(1,2,3,10,12,81)]
    clinical$X.1 <- gsub('\\.','-',as.character(clinical$X.1))
    common.clinical <- intersect(as.character(clinical$bcr_patient_barcode),substr(as.character(purity_leukocyte$sample),1,12))
    match.result <- match(clinical$bcr_patient_barcode, common.clinical)
    big.index <- which(is.na(match.result)==FALSE)
    clinical <- clinical[big.index,]
    tumor.type <- as.data.frame(table(as.character(clinical$cancer)))
    tumor.type.number <- as.numeric(tumor.type[,2])
    tumor.type.name <- as.character(tumor.type[,1])
    all.non.immune.clinical <- read.csv(file = '/fshare2/JiFansen/Clinical/ALLnon-immuneSurvivalthree_year_data_matchexp.csv', header = T, sep = ',')
    all.non.immune.clinical <- all.non.immune.clinical[,-c(2,3,4,5,6,7,9)]
    common.clinical <- intersect(as.character(all.non.immune.clinical$X),substr(as.character(purity_leukocyte$sample),1,12))
    match.result <- match(all.non.immune.clinical$X, common.clinical)
    big.index <- which(is.na(match.result)==FALSE)
    all.non.immune.clinical <- all.non.immune.clinical[big.index,]
    group = rep('No-clinical',ncol(expression.table))
    match.result <- match(substr(colnames(expression.table),1,12), as.character(clinical$bcr_patient_barcode))
    immuno.big.index <- which(is.na(match.result)==FALSE)
    immuno.small.index <- match.result[immuno.big.index]
    group[immuno.big.index] <- 'Clinical'
    survival.time <- rep(0,ncol(expression.table))
    survival.time[immuno.big.index] <- clinical$times[immuno.small.index]
    stat <- as.data.frame(table(colnames(expression.table)))
    stat <- stat[stat$Freq>1,]
    titleName <- substr(colnames(expression.table),1,12)
    for(i in 1:dim(stat)[1]){
      titleName[which(titleName==as.character(stat[i,1]))[2]] <- paste(as.character(stat[i,1]),'-1', sep='')
    }
    immune.pvalue <- c()
    nonimmune.pvalue <- c()
    for(s in 1:1000){
      non.immune.clinical <- 0
      for(w in 1:length(tumor.type.name)){
        cancer.index <- which(all.non.immune.clinical$cancer_type==tumor.type.name[w])
        if(w==1){
          non.immune.clinical <- all.non.immune.clinical[cancer.index[sample(1:length(cancer.index),tumor.type.number[w])],]
        }else{
          non.immune.clinical <- rbind.data.frame(non.immune.clinical,all.non.immune.clinical[cancer.index[sample(1:length(cancer.index),tumor.type.number[w])],])
        }
      }
      common.clinical <- intersect(as.character(non.immune.clinical$X),substr(as.character(purity_leukocyte$sample),1,12))
      match.result <- match(non.immune.clinical$X, common.clinical)
      non.immuno.big.index <- which(is.na(match.result)==FALSE)
      non.immune.clinical <- non.immune.clinical[non.immuno.big.index,]
      other.therapy = rep('no-therapy', ncol(expression.table))
      match.result <- match(titleName, as.character(non.immune.clinical$X))
      other.big.index <- which(is.na(match.result)==FALSE)
      other.small.index <- match.result[other.big.index]
      
      other.therapy[other.big.index] <- 'Other-therapy'
      other.thrapy.survival.time <- rep(0,ncol(expression.table))
      other.thrapy.survival.time[other.big.index] <- non.immune.clinical$times[other.small.index]
      dataset <- data.frame(group = group, survival = survival.time, OtherTerapy = other.therapy, 
                            OtherSurvival = other.thrapy.survival.time,
                            shape = as.factor(groups))
      survival.table <- dataset[immuno.big.index,]
      survival.table$status <- clinical$vital_status[immuno.small.index]
      survival.table$shape <- as.factor(survival.table$shape)
      surv_object <- Surv(time = survival.table$survival, event = survival.table$status)
      fit1 <- survdiff(surv_object ~ shape, data = survival.table)
      p.fit1 <- 1-pchisq(fit1$chisq, length(fit1$n)-1)
      #fit1 <- pairwise_survdiff(surv_object ~ shape, data = survival.table)
      #p.fit1 <- as.numeric(fit1$p.value)
      survival.table.new <- dataset[other.big.index, ]
      survival.table.new$status <- as.numeric(non.immune.clinical$vital_status[other.small.index])
      survival.table.new$group <- as.factor(survival.table.new$shape)
      SampleRecord[[s]] <- survival.table.new
      surv_object.new <- Surv(time = survival.table.new$OtherSurvival, event = survival.table.new$status)
      #fit2 <- survfit(surv_object.new ~ shape, data = survival.table.new)
      #p.survival.new <- ggsurvplot(fit2, data = survival.table.new, pval = TRUE)
      #print(p.survival.new)
      #fit2 <- pairwise_survdiff(surv_object.new ~ shape, data = survival.table.new)
      fit2 <- survdiff(surv_object.new ~ shape, data = survival.table.new)
      #p.fit2 <- as.numeric(fit2$p.value)
      p.fit2 <- 1-pchisq(fit2$chisq, length(fit2$n)-1)
      immune.pvalue <- c(immune.pvalue, p.fit1)
      nonimmune.pvalue <- c(nonimmune.pvalue, p.fit2)
    }
    a <- length(which(nonimmune.pvalue<unique(immune.pvalue)))
    # Plot the median non-immunotherapy p value survival plots.
    SampleRecord <- SampleRecord[[which(nonimmune.pvalue==sort(nonimmune.pvalue)[500])]]
    surv_object.new <- Surv(time = SampleRecord$OtherSurvival, event = SampleRecord$status)
    fit2 <- survfit(surv_object.new ~ shape, data = SampleRecord)
    p.survival <- ggsurvplot(fit2, data = SampleRecord, pval = TRUE)
    pdf(paste0('/Share/home/JiFansen/TCGA-', cluster_number, 
               '-', top.genes,'-non-immuno-survival-oneTime.pdf'))
    print(p.survival$plot)
    dev.off()
    # Plot the immunotherapy p value survival plot.
    fit1 <- survfit(surv_object ~ shape, data = survival.table)
    p.survival <- ggsurvplot(fit1, data = survival.table, pval = TRUE)
    pdf(paste0('/Share/home/JiFansen/TCGA-', cluster_number, 
               '-', top.genes,'-survival.pdf'))
    print(p.survival$plot)
    dev.off()
    cat('immunotherapy p.value: ', unique(immune.pvalue), '\n')
    cat('how many times non-immunotherapy is significant than immunotherapy: ', a, '\n')
    ############################################ Point out which group survive longer. ######################################################
    tcga.status <- c()
    for(h in 1:cluster_number){
      survival.subset <- survival.table[which(survival.table$shape==h), ]
      tcga.status <- c(tcga.status, mean(survival.subset$survival))
    }
    centriod.status <- rep('between', cluster_number)
    centriod.status[which(tcga.status==max(tcga.status))] <- 'high'
    centriod.status[which(tcga.status==min(tcga.status))] <- 'low'
    ##########################################################################################################################################
    groups <- as.character(groups)
    immuno_label <- as.character(dataset$group)
    immuno_label[which(immuno_label=='Clinical')] <- 'immunotherapy'
    immuno_label[which(immuno_label=='No-clinical')] <- 'other'
    annotation_col = data.frame(cluster = factor(groups), cancertype = cancertype, treatment = as.factor(immuno_label))
    rownames(annotation_col) <- colnames(expression.candidate)
    treatment = c("magenta2", "grey92")
    names(treatment) = c("immunotherapy", "other" )
    ann_colors = list(treatment = treatment)
    pdf(file = paste0('/Share/home/JiFansen/TCGA-', cluster_number, 
                      '-', top.genes,'.pdf'), width = 40, height = 20)
    #hv <- pheatmap(as.matrix(expression.candidate), scale = 'row', color = redblue(100)[sort(seq(1:100), decreasing = T)], 
    #               show_colnames = FALSE, treeheight_row = 0, annotation_col = annotation_col, cluster_cols = hc,cellheight = 20, fontsize = 20)
    pheatmap(as.matrix(expression.candidate), scale = 'row', color = redblue(100)[sort(seq(1:100), decreasing = T)], 
                   show_colnames = FALSE, treeheight_row = 0, annotation_col = annotation_col, annotation_colors = ann_colors, cluster_cols = hc,cellheight = 20, fontsize = 20)
    dev.off()
    genes.index <- hv$tree_row$order
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
    colnames(survival_42)[4] <- 'status'
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
    ############################################# Plot gse data. ####################################################################
    dataset.gse <- dataset.gse[match(candidate.genes.top20, rownames(dataset.gse)), ]
    dataset.gse <- dataset.gse[, -cutoffIndex(data = survival.gse, colTime = 4, colStatus = 5, cutoff = 365)]
    survival.gse <- survivalCutoff(data = survival.gse, colTime = 4, colStatus = 5, cutoff = 365)
    #dataset.gse <- dataset.gse[genes.index, ]
    hc<-hclust(dist(t(dataset.gse)))
    groups <- cutree(hc,k=cluster_number)
    groups <- as.character(groups)
    annotation_col = data.frame(cluster = factor(groups))
    rownames(annotation_col) <- colnames(dataset.gse)
    pdf(file = paste0('/Share/home/JiFansen/gse-', cluster_number, 
                      '-', top.genes,'.pdf'), width = 40, height = 20)
    hv <- pheatmap(as.matrix(dataset.gse), scale = 'row', color = redblue(100)[sort(seq(1:100), decreasing = T)], 
                   show_colnames = FALSE, treeheight_row = 0, annotation_col = annotation_col, cluster_cols = hc, 
                   cluster_row = FALSE,cellheight = 20, fontsize = 20)
    dev.off()
    label <- apply(dataset.gse, 2, FUN = function(x){
      Distance <- apply(expression.candidate,2, FUN = function(y){
        return(sqrt(sum((x-y)^2)))
      })
      if(length(as.numeric(which(Distance==min(Distance))))>1){
        return(as.numeric(which(Distance==min(Distance)))[1])
      }else(
        return(as.numeric(which(Distance==min(Distance))))
      )
    })
    predicted <- centriod.status[as.numeric(Groups)[label]]
    #################################################### gse data survival analysis. ###############################################################
    #survival.table <- data.frame(survival.gse, shape = as.factor(predicted))
    #surv_object <- Surv(time = survival.table$drug.survival.times, event = survival.table$status)
    #fit1 <- survdiff(surv_object ~ shape, data = survival.table)
    #p.fit1 <- 1-pchisq(fit1$chisq, length(fit1$n)-1)
    #fit1 <- pairwise_survdiff(surv_object ~ shape, data = survival.table)
    #p.fit1 <- as.numeric(fit1$p.value)
    #fit1 <- survfit(surv_object ~ shape, data = survival.table)
    #p.survival <- ggsurvplot(fit1, data = survival.table, pval = TRUE)
    #pdf(paste0('/Share/home/JiFansen/JiFansen/Plots/10-07/gse-', cluster_number, 
    #           '-', top.genes,'-survival.pdf'))
    #print(p.survival$plot)
    #dev.off()
    ############################################################# Calculate accuracy. ############################################################
    class <- rep('between', dim(survival.gse)[1])
    #quantile.number <- quantile(survival.gse$drug.survival.times, 
    #                            probs = seq(0, 1, by=1/cluster_number))
    class[which(survival.gse$status==0)] <- 'high'
    class[which(survival.gse$status==1)] <- 'low'
    if(cluster_number==2){
      #confusion.result <- confusionMatrix(predicted, class)
      print(paste0('cluster:', cluster_number, '   genes:', top.genes, '    gse'))
      #print(confusion.result$table)
      print(table(predicted, class))
    }else{
      #class[which(survival.gse$drug.survival.times<quantile.number[2])] <- 'low'
      #class[which(survival.gse$drug.survival.times>=quantile.number[3])] <- 'high'
      between.index <- which(predicted=='between')
      predicted <- predicted[-between.index]
      class <- class[-between.index]
      #confusion.result <- confusionMatrix(predicted, class)
      print(paste0('cluster:', cluster_number, '   genes:', top.genes, '    gse'))
      #print(confusion.result$table)
      print(table(predicted, class))
    }
    ############################################# Plot bms data. ####################################################################
    dataset.bms <- dataset.bms[match(candidate.genes.top20, rownames(dataset.bms)), ]
    dataset.bms <- dataset.bms[, -cutoffIndex(data = survival.bms, colTime = 4, colStatus = 5, cutoff = 365)]
    survival.bms <- survivalCutoff(data = survival.bms, colTime = 4, colStatus = 5, cutoff = 365)
    #dataset.bms <- dataset.bms[genes.index, ]
    hc<-hclust(dist(t(dataset.bms)))
    groups <- cutree(hc,k=cluster_number)
    groups <- as.character(groups)
    annotation_col = data.frame(cluster = factor(groups))
    rownames(annotation_col) <- colnames(dataset.bms)
    pdf(file = paste0('/Share/home/JiFansen/bms-', cluster_number, 
                      '-', top.genes,'.pdf'), width = 40, height = 20)
    hv <- pheatmap(as.matrix(dataset.bms), scale = 'row', color = redblue(100)[sort(seq(1:100), decreasing = T)], 
                   show_colnames = FALSE, treeheight_row = 0, annotation_col = annotation_col, cluster_cols = hc, 
                   cluster_row = FALSE,cellheight = 20, fontsize = 20)
    dev.off()
    label <- apply(dataset.bms, 2, FUN = function(x){
      Distance <- apply(expression.candidate,2, FUN = function(y){
        return(sqrt(sum((x-y)^2)))
      })
      if(length(as.numeric(which(Distance==min(Distance))))>1){
        return(as.numeric(which(Distance==min(Distance)))[1])
      }else(
        return(as.numeric(which(Distance==min(Distance))))
      )
    })
    predicted <- centriod.status[as.numeric(Groups)[label]]
    ################################################### bms data survival analysis. #################################################################
    #survival.table <- data.frame(survival.bms, shape = as.factor(predicted))
    #surv_object <- Surv(time = survival.table$drug.survival.times, event = survival.table$status)
    #fit1 <- survdiff(surv_object ~ shape, data = survival.table)
    #p.fit1 <- 1-pchisq(fit1$chisq, length(fit1$n)-1)
    #fit1 <- pairwise_survdiff(surv_object ~ shape, data = survival.table)
    #p.fit1 <- as.numeric(fit1$p.value)
    #fit1 <- survfit(surv_object ~ shape, data = survival.table)
    #p.survival <- ggsurvplot(fit1, data = survival.table, pval = TRUE)
    #pdf(paste0('/Share/home/JiFansen/JiFansen/Plots/10-07/bms-', cluster_number, 
    #           '-', top.genes,'-survival.pdf'))
    #print(p.survival$plot)
    #dev.off()
    class <- rep('between', dim(survival.bms)[1])
    class[which(survival.bms$status==0)] <- 'high'
    class[which(survival.bms$status==1)] <- 'low'
    if(cluster_number==2){
      #confusion.result <- confusionMatrix(predicted, class)
      print(paste0('cluster:', cluster_number, '   genes:', top.genes, '    gse'))
      #print(confusion.result$table)
      print(table(predicted, class))
    }else{
      #class[which(survival.gse$drug.survival.times<quantile.number[2])] <- 'low'
      #class[which(survival.gse$drug.survival.times>=quantile.number[3])] <- 'high'
      between.index <- which(predicted=='between')
      predicted <- predicted[-between.index]
      class <- class[-between.index]
      #confusion.result <- confusionMatrix(predicted, class)
      print(paste0('cluster:', cluster_number, '   genes:', top.genes, '    gse'))
      #print(confusion.result$table)
      print(table(predicted, class))
    }
    ############################################# Plot sample42 data. ####################################################################
    dataset_42 <- dataset_42[match(candidate.genes.top20, rownames(dataset_42)), ]
    if(length(cutoffIndex(data = survival_42, colTime = 2, colStatus = 4, cutoff = 365))>0){
      dataset_42 <- dataset_42[, -cutoffIndex(data = survival_42, colTime = 2, colStatus = 4, cutoff = 365)]
      survival_42 <- survivalCutoff(data = survival_42, colTime = 2, colStatus = 4, cutoff = 365)
    }
    #dataset_42 <- dataset_42[genes.index, ]
    hc<-hclust(dist(t(dataset_42)))
    groups <- cutree(hc,k=cluster_number)
    groups <- as.character(groups)
    annotation_col = data.frame(cluster = factor(groups))
    rownames(annotation_col) <- colnames(dataset_42)
    pdf(file = paste0('/Share/home/JiFansen/sample42-', cluster_number, 
                      '-', top.genes,'.pdf'), width = 40, height = 20)
    hv <- pheatmap(as.matrix(dataset_42), scale = 'row', color = redblue(100)[sort(seq(1:100), decreasing = T)], 
                   show_colnames = FALSE, treeheight_row = 0, annotation_col = annotation_col, cluster_cols = hc, 
                   cluster_row = FALSE, cellheight = 20, fontsize = 20)
    dev.off()
    label <- apply(dataset_42, 2, FUN = function(x){
      Distance <- apply(expression.candidate,2, FUN = function(y){
        return(sqrt(sum((x-y)^2)))
      })
      if(length(as.numeric(which(Distance==min(Distance))))>1){
        return(as.numeric(which(Distance==min(Distance)))[1])
      }else(
        return(as.numeric(which(Distance==min(Distance))))
      )
    })
    predicted <- centriod.status[as.numeric(Groups)[label]]
    ############################################################ sample42 survival analysis. ##################################################
    #survival.table <- data.frame(survival_42, shape = as.factor(predicted))
    #surv_object <- Surv(time = survival.table$overall_survival, event = survival.table$Status)
    #fit1 <- survdiff(surv_object ~ shape, data = survival.table)
    #p.fit1 <- 1-pchisq(fit1$chisq, length(fit1$n)-1)
    #fit1 <- pairwise_survdiff(surv_object ~ shape, data = survival.table)
    #p.fit1 <- as.numeric(fit1$p.value)
    #fit1 <- survfit(surv_object ~ shape, data = survival.table)
    #p.survival <- ggsurvplot(fit1, data = survival.table, pval = TRUE)
    #pdf(paste0('/Share/home/JiFansen/JiFansen/Plots/10-07/sample42-', cluster_number, 
    #           '-', top.genes,'-survival.pdf'))
    #print(p.survival$plot)
    #dev.off()
    class <- rep('between', dim(survival_42)[1])
    class[which(survival_42$status==0)] <- 'high'
    class[which(survival_42$status==1)] <- 'low'
    if(cluster_number==2){
      #confusion.result <- confusionMatrix(predicted, class)
      print(paste0('cluster:', cluster_number, '   genes:', top.genes, '    gse'))
      #print(confusion.result$table)
      print(table(predicted, class))
    }else{
      #class[which(survival.gse$drug.survival.times<quantile.number[2])] <- 'low'
      #class[which(survival.gse$drug.survival.times>=quantile.number[3])] <- 'high'
      between.index <- which(predicted=='between')
      predicted <- predicted[-between.index]
      class <- class[-between.index]
      #confusion.result <- confusionMatrix(predicted, class)
      print(paste0('cluster:', cluster_number, '   genes:', top.genes, '    gse'))
      #print(confusion.result$table)
      print(table(predicted, class))
    }
    library(clusterProfiler)
    candidate.genes.top20 <- candidate.genes[1:top.genes]
    expression.candidate <- expression.table[match(candidate.genes.top20, rownames(expression.table)), ]
    hc<-hclust(dist(t(expression.candidate)))
    groups <- cutree(hc,k=cluster_number)
    Groups <- groups
    library(org.Hs.eg.db)
    immuno.genes <- bitr(rownames(expression.table), fromType = 'SYMBOL', toType ="ENTREZID", OrgDb="org.Hs.eg.db")
    immuno.expression <- expression.table[immuno.genes$SYMBOL,]
    if(cluster_number==2){
      immuno.pvalues <- apply(immuno.expression,1,FUN = function(x){
        return(wilcox.test(as.numeric(x)~as.factor(Groups))$p.value)
      })
      fold.change <- apply(immuno.expression, 1, FUN = function(x){
        zero.samples <- median(x[which(Groups==0)])
        one.samples <- median(x[which(Groups==1)])
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
      write.table(kk2[,1:9], file = '/Share/home/JiFansen/JiFansen/TCGA-Project/differentialAnalysis/GSEA-no-correction-cluster-2-label.txt', row.names = F, sep = '\t', quote = F)
    }else{
      clinical <- read.csv(file = '/fshare2/JiFansen/Clinical/three_year_data.csv', header = T, sep = ',')
      clinical <- clinical[,c(1,2,3,10,12,81)]
      clinical$X.1 <- gsub('\\.','-',as.character(clinical$X.1))
      match.result <- match(clinical$X.1, names(Groups))
      big.index <- which(is.na(match.result)==FALSE)
      small.index <- match.result[big.index]
      clinical <- clinical[big.index, ]
      clinical$group <- Groups[small.index]
      surv_object <- Surv(time = clinical$times, event =clinical$vital_status)
      fit1 <- survfit(surv_object ~ group, data = clinical)
      p.survival <- ggsurvplot(fit1, data = clinical, pval = TRUE)
      pdf('~/detect-high-low-survival.pdf')
      p.survival$plot
      dev.off()
      library(dplyr)
      planes <- group_by(clinical, group)
      summarise(planes, sum = median(times))
      rm.index <- which(Groups == 2)
      immuno.pvalues <- apply(immuno.expression[, -rm.index],1,FUN = function(x){
        return(wilcox.test(as.numeric(x)~as.factor(Groups)[-rm.index])$p.value)
      })
      geneList <- data.frame(gene = immuno.genes$SYMBOL, entrez = immuno.genes$ENTREZID, p.value = immuno.pvalues)
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
      write.table(kk2[,1:9], file = '/Share/home/JiFansen/JiFansen/TCGA-Project/differentialAnalysis/GSEA-no-correction-cluster-3-label.txt', row.names = F, sep = '\t', quote = F)
    }
  }
}

