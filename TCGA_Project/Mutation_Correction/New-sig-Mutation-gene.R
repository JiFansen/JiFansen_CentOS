w=1
cancerlable <- c()
for(File in list.files('/Share/home/JiFansen/JiFansen/Mutation-Correction/MAF-summary/', pattern = '*-MAF-summary.txt')){
  total.summary <- read.table(file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/MAF-summary/', File), header = T, sep = '\t')
  non.sys.summary <- total.summary
  uniq.plus.genes <- unique(non.sys.summary$gene)
  final.results1 <- c();final.results2 <- c();final.results3 <- c()
  for (i in 1:length(uniq.plus.genes)){
    gene.name <- as.character(uniq.plus.genes[i])
    test.data <- non.sys.summary[which(non.sys.summary$gene==gene.name),]
    non.count <- length(which(test.data$Variant_Classification=="Missense_Mutation"))
    sys.count <- length(which(test.data$Variant_Classification=="Silent"))
    final.results1 <- c(final.results1,gene.name)
    final.results2 <- c(final.results2, non.count)
    final.results3 <- c(final.results3, sys.count)
  }
  final.count <- data.frame(geneName = final.results1, non.count = final.results2, sys.count = final.results3)
  genename.initial <- as.character(final.count$geneName)
  genename.initial[which(genename.initial=="7-Mar")] <- "LDOC1";genename.initial[which(genename.initial=="2-Mar")] <- "PEG10";
  genename.initial[which(genename.initial=="8-Mar")] <- "RTL8C";genename.initial[which(genename.initial=="6-Mar")] <- "RTL6";
  genename.initial[which(genename.initial=="10-Mar")] <- "MARCH10";genename.initial[which(genename.initial=="1-Mar")] <- "Mar1";
  genename.initial[which(genename.initial=="4-Mar")] <- "RTL4";genename.initial[which(genename.initial=="11-Mar")] <- "MARCH11";
  genename.initial[which(genename.initial=="8-Sep")] <- "SEPT8";genename.initial[which(genename.initial=="3-Sep")] <- "SEPT3";
  genename.initial[which(genename.initial=="6-Sep")] <- "SEPT6";genename.initial[which(genename.initial=="2-Sep")] <- "SEPT2";
  genename.initial[which(genename.initial=="7-Sep")] <- "SEPT7";genename.initial[which(genename.initial=="4-Sep")] <- "SEPT4";
  genename.initial[which(genename.initial=="5-Sep")] <- "SEPT5";genename.initial[which(genename.initial=="1-Sep")] <- "SEPT1";
  genename.initial[which(genename.initial=="10-Sep")] <- "SEPT10";genename.initial[which(genename.initial=="11-Sep")] <- "SEPT11";
  genename.initial[which(genename.initial=="12-Sep")] <- "SEPT12";genename.initial[which(genename.initial=="15-Sep")] <- "SEP15";
  genename.initial[which(genename.initial=="14-Sep")] <- "SEPT14";genename.initial[which(genename.initial=="1-Dec")] <- "DEC1";
  final.count$gene <- genename.initial
  cancertype <- unlist(strsplit(File,'-'))[1]
  cancertype.background <- read.table(file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/New-Background-Method/', 
                                                    cancertype, '-background-tmp.txt'), header = T, sep = '\t')
  match.result <- match(as.character(final.count$gene), as.character(cancertype.background$gene))
  na.index <- which(is.na(match.result)==TRUE)
  final.results <- cbind.data.frame(final.count[-na.index,], cancertype.background[match.result[-na.index],])
  final.results <- final.results[,c(4,2,3,7,8)]
  write.table(final.results, file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/GeneCount-Background/', 
                                           cancertype, '-final.result.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
}

w=1
cancerlable <- c()
for(File in list.files('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/GeneCount-Background/')){
  cancertype <- unlist(strsplit(File,'-'))[1]
  final.results <- read.table(file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/GeneCount-Background/', File), header = T, 
                             sep = '\t')
  final.results$n.dn <- final.results$non.count/final.results$non.background
  final.results$s.ds <- final.results$sys.count/final.results$sys.background
  final.results$correct.n.dn <- final.results$n.dn+1
  final.results$correct.s.ds <- final.results$s.ds+1
  final.results$non.ratio <- final.results$non.background/(final.results$non.background+final.results$sys.background)
  final.results$total.count <- final.results$non.count+final.results$sys.count
  final.results$dN.dS <- final.results$correct.n.dn/final.results$correct.s.ds
  tmp <- c()
  for(i in 1:dim(final.results)[1]){
    tmp[i] <- binom.test(as.numeric(final.results$non.count)[i],as.numeric(final.results$total.count)[i],as.numeric(final.results$non.ratio)[i],alternative = "greater")$p.value
  }
  final.results$p.value <- tmp
  write.table(final.results, file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/dN-dS/', cancertype,'-dN.dS.txt'), sep = '\t',
              quote = FALSE, row.names = FALSE)
  result.1 <- final.results[which(final.results$dN.dS>1),]
  print(paste0(cancertype, ':', length(which(result.1$p.value<=0.05))))
  result <- result.1[which(result.1$p.value<=0.05),]
  cancerlable <- c(cancerlable, rep(cancertype,length(which(result.1$p.value<=0.05))))
  if(w==1){
    final.sig <- result
    w <- w+1
  }else{
    final.sig <- rbind.data.frame(final.sig, result)
    w <- w+1
  }
}
final.sig$cancertype <- cancerlable
write.table(final.sig, file = '/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/Final-sig-gene.txt', 
            sep = '\t', quote = FALSE, row.names = FALSE)
save.image(file = '/Share/home/JiFansen/JiFansen/working-image/New-sig-Mutation-gene.RData')
