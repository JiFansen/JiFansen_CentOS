w <- 1
for(file in list.files('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/New-196-Mutation-Signature-Calculation/')){
  tmp <- read.table(file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/New-196-Mutation-Signature-Calculation/', file), header = T, sep = '\t', row.names = 1)
  colnames(tmp) <- unlist(strsplit(file,'\\.'))[1]
  if(w==1){
    background <- tmp
    w <- w+1
  }else{
    background <- cbind.data.frame(background, tmp)
  }
}
sampleNumber <- c()
w <- 1
for(file in list.files('/Share/home/JiFansen/JiFansen/Mutation-Correction/MAF-summary/', pattern = '.barcode.txt')){
  cancertype <- unlist(strsplit(file, '-'))[1]
  print(cancertype)
  tmp <- read.table(file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/MAF-summary/', file), header = T, sep = '\t')
  tmp$cancertype <- cancertype
  sample.barcode <- substr(as.character(tmp$Tumor_Sample_Barcode), 1, 16)
  sampleNumber <- c(sampleNumber, length(unique(sample.barcode)))
  if(w==1){
    data.initial <- tmp
    w <- w+1
  }else{
    data.initial <- rbind.data.frame(data.initial, tmp)
    w <- w+1
  }
}
total.background <- as.matrix(apply(background, 1, FUN = function(x){return(sum(x*sampleNumber))}))
write.table(total.background, file = '/Share/home/JiFansen/JiFansen/Mutation-Correction/Correction/total-cancertype-background.txt', sep = '\t', quote = FALSE, col.names = FALSE)
data.initial <- data.initial[which(data.initial$Variant_Classification=='Missense_Mutation'),]
match.result <- match(as.character(data.initial$change), rownames(total.background))
data.initial$true.value <- as.numeric(total.background[,1])[match.result]
data.initial$count <- 1
data.summary <- aggregate(x = data.initial[,13:14], by = list(data.initial$gene, data.initial$Protein_position), FUN = sum)
data.summary$absolute <- data.summary$true.value/data.summary$count
write.table(data.summary, file = '/Share/home/JiFansen/JiFansen/Mutation-Correction/Correction/Gene-Count-TrueValue.txt', sep = '\t', quote = FALSE, row.names = FALSE)
uniq.gene <- unique(as.character(data.summary$Group.1))
w <- 1
for(genename in uniq.gene){
  tmp <- data.summary[which(data.summary$Group.1==genename),]
  tmp$scale <- (tmp$absolute-min(tmp$absolute))/(max(tmp$absolute)-min(tmp$absolute))
  if(w == 1){
    score.correction <- tmp
    w <- w+1
  }else{
    score.correction <- rbind.data.frame(score.correction, tmp)
  }
}
sig.genes <- read.table(file = '/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/New-Overlap.txt', header = F, sep = '\t')
sig.genes <- as.character(sig.genes[,1])
match.result <- match(as.character(score.correction$Group.1), sig.genes)
sig.index <- which(is.na(match.result)==FALSE)
score.correction.sig <- score.correction[sig.index, ]
write.table(score.correction.sig, file = '/Share/home/JiFansen/JiFansen/Mutation-Correction/Correction/Total-Codon-Score.txt', sep = '\t', quote = FALSE, row.names = FALSE)

gene.codon <- paste0(as.character(score.correction.sig$Group.1), '-', as.character(score.correction.sig$Group.2))

for(file in list.files('/Share/home/JiFansen/JiFansen/Mutation-Correction/MAF-summary/', pattern = '.barcode.txt')){
  cancertype <- unlist(strsplit(file, '-'))[1]
  tmp <- read.table(file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/MAF-summary/', file), header = T, sep = '\t')
  tmp <- tmp[which(tmp$Variant_Classification=='Missense_Mutation'), ]
  tmp$cancertype <- cancertype
  sample.barcode <- substr(as.character(tmp$Tumor_Sample_Barcode), 1, 16)
  tmp$Tumor_Sample_Barcode <- sample.barcode
  tmp.gene.codon <- paste0(as.character(tmp$gene), '-',as.character(tmp$Protein_position))
  match.result <- match(tmp.gene.codon, gene.codon)
  big.index <- which(is.na(match.result)==FALSE)
  small.index <- match.result[big.index]
  tmp.sig <- tmp[big.index,]
  tmp.sig$score <- score.correction.sig$scale[small.index]
  tmp.sig <- tmp.sig[,c(1,10,13)]
  frameshift <- read.table(file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Frameshift/', cancertype, '-shift.txt-New'), header = T, sep = '\t', stringsAsFactors = FALSE)
  frameshift <- t(frameshift)
  colnames(frameshift) <- as.character(frameshift[1,])
  frameshift <- frameshift[-c(1,2),]
  frameshift <- frameshift[-nrow(frameshift),]
  frameshift <- as.data.frame(frameshift)
  for(m in 1:dim(frameshift)[1]){
    index <- which(as.numeric(as.matrix(frameshift[m,]))>0)
    frame.tmp <- data.frame(gene = rownames(frameshift)[m], Tumor_Sample_Barcode = colnames(frameshift)[index], score = 0)
    tmp.sig <- rbind.data.frame(tmp.sig, frame.tmp)
  }
  Sample.uniq <- unique(as.character(tmp.sig$Tumor_Sample_Barcode))
  p <- 1
  for(i in Sample.uniq){
    gene.score <- c()
    test.1 <- tmp.sig[which(tmp.sig$Tumor_Sample_Barcode==i),]
    Gene.uniq <- unique(as.character(test.1$gene))
    for(j in Gene.uniq){
      test.2 <- test.1[which(test.1$gene==j),]
      gene.score <- c(gene.score,prod(as.numeric(test.2$score)))
    }
    result <- data.frame(sample = i, gene = Gene.uniq, score = gene.score)
    if(p==1){
      result.cancertype <- result
      p <- p+1
    }else{
      result.cancertype <- rbind.data.frame(result.cancertype, result)
    }
  }
  write.table(result.cancertype, file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Correction/', cancertype, '-Score-Correction.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
}
