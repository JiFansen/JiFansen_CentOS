Args <- commandArgs()
p <- as.numeric(Args[6])
cds <- read.table(file = '/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/CCDS/sequences', header = F, sep = '\t')
header <- read.table(file = '/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/CCDS/newheader', header = F, sep = '\t')
maptable <- read.table(file = '/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/CCDS/useful.CCDS.information', header = F, sep = '\t')
colnames(maptable) <- c('chromosome', 'gene', 'geneID', 'CCDSID')
map.gene <- as.character(maptable[,2])
map.ccdsid <- apply(maptable,1,FUN = function(x){unlist(strsplit(as.character(x[4]),"\\."))[1]})
gene <- as.matrix(apply(header,1,FUN = function(x){
  a <- unlist(strsplit(as.character(x),"\\|"))[4]
  index <- which(map.ccdsid==a)
  return(as.character(map.gene[index]))
}))
header$gene <- gene
########################################################### Read the codon table. #################################################################
codon <- read.table(file="/Share/home/JiFansen/JiFansen/Mutation-Correction/MAF-summary/codon1.txt",sep="\t",header=F)
codon <- codon[,3:4]
colnames(codon) <- c("Amino Acid","Codon")
codon.table <- matrix(ncol = 2)
for (i in 1:dim(codon)[1]){
  test <- unlist(strsplit(as.character(codon$Codon[i]),split=","))
  matrix <- as.matrix(data.frame(x=rep(codon$`Amino Acid`[i],length(unlist(strsplit(
    as.character(codon$Codon[i]),split=",")))),y=unlist(strsplit(
      as.character(codon$Codon[i]),split=","))))
  codon.table <- rbind(codon.table,matrix)
}
codon.table <- as.data.frame(codon.table[-1,])
colnames(codon.table) <- c("Amino.Acid","Codon")
codon.table <- codon.table[-62,]
for (i in 1:dim(codon.table)[1]){
  char <- unlist(strsplit(as.character(codon.table$Codon[i]),""))
  for (j in 1:length(char)){
    if(char[j]=="U")
    {
      char[j] = "T"
    }
  }
  codon.table[i,3] <- paste(char[1],char[2],char[3],sep="")
}
colnames(codon.table)[3] <- "Coding.strand"
all.change <- c()
bases <- c('A','T','C','G')
for(i in bases){
  for(j in bases){
    all.change <- c(all.change,paste0(i,j))
  }
}
all.change <- all.change[-c(1,6,11,16)]
all.change <- unlist(strsplit(all.change,""))
mutation.condition <- function(x){
  sequence <- as.character(x[1])
  codon.number <- nchar(sequence)/3
  non.total.background <- c()
  sys.total.background <- c()
  for(w in 2:(codon.number-1)){
    string = substr(sequence, 3*w-2, 3*w)
    amino.acid <- codon.table[which(as.character(codon.table$Coding.strand)==string),1]
    string <- unlist(strsplit(string,""))
    if(string[1]=='A'){
      char1 <- data.frame(a=paste0(all.change[c(2,4,6)],string[2],string[3]),b=paste0('A',paste0(all.change[c(2,4,6)])))
    }
    if(string[1]=='T'){
      char1 <- data.frame(a=paste0(all.change[c(8,10,12)],string[2],string[3]),b=paste0('T',paste0(all.change[c(8,10,12)])))
    }
    if(string[1]=='C'){
      char1 <- data.frame(a=paste0(all.change[c(14,16,18)],string[2],string[3]),b=paste0('C',paste0(all.change[c(14,16,18)])))
    }
    if(string[1]=='G'){
      char1 <- data.frame(a=paste0(all.change[c(20,22,24)],string[2],string[3]),b=paste0('G',paste0(all.change[c(20,22,24)])))
    }
    if(string[2]=='A'){
      char2 <- data.frame(a=paste0(string[1], all.change[c(2,4,6)],string[3]),b=paste0('A',paste0(all.change[c(2,4,6)])))
    }
    if(string[2]=='T'){
      char2 <- data.frame(a=paste0(string[1],all.change[c(8,10,12)],string[3]),b=paste0('T',paste0(all.change[c(8,10,12)])))
    }
    if(string[2]=='C'){
      char2 <- data.frame(a=paste0(string[1],all.change[c(14,16,18)],string[3]),b=paste0('C',paste0(all.change[c(14,16,18)])))
    }
    if(string[2]=='G'){
      char2 <- data.frame(a=paste0(string[1],all.change[c(20,22,24)],string[3]),b=paste0('G',paste0(all.change[c(20,22,24)])))
    }
    if(string[3]=='A'){
      char3 <- data.frame(a=paste0(string[1],string[2],all.change[c(2,4,6)]),b=paste0('A',paste0(all.change[c(2,4,6)])))
    }
    if(string[3]=='T'){
      char3 <- data.frame(a=paste0(string[1],string[2],all.change[c(8,10,12)]),b=paste0('T',paste0(all.change[c(8,10,12)])))
    }
    if(string[3]=='C'){
      char3 <- data.frame(a=paste0(string[1],string[2],all.change[c(14,16,18)]),b=paste0('C',paste0(all.change[c(14,16,18)])))
    }
    if(string[3]=='G'){
      char3 <- data.frame(a=paste0(string[1],string[2],all.change[c(20,22,24)]),b=paste0('G',paste0(all.change[c(20,22,24)])))
    }
    char <- rbind.data.frame(char1,char2,char3)
    char4 <- c()
    for(j in 1:9){
      if(codon.table[which(as.character(codon.table[,3])==as.character(char[j,1])),1]==amino.acid){
        char4 <- c(char4,"sys")
      }
      else{
        char4 <- c(char4,"non")
      }
    }
    char$c <- char4
    char5 <- c()
    between <- unlist(strsplit(as.character(x[1]),""))
    for (z in 1:3){
      char5 <- c(char5, paste("(",between[3*w-3],between[3*w-1],")",sep=""))
    }
    for(z in 4:6){
      char5 <- c(char5, paste("(",between[3*w-2],between[3*w],")",sep=""))
    }
    for (z in 7:9){
      char5 <- c(char5, paste("(",between[3*w-1],between[3*w+1],")",sep=""))
    }
    char$d <- char5
    combination1 <- as.character(final.background$change)
    combination2 <- paste0(char$b,char$d)
    match.result <- match(combination2,combination1)
    char$e <- as.numeric(background[,3])[match.result]
    non.total.background <- c(non.total.background, sum(char[which(char[,3]=="non"),5]))
    sys.total.background <- c(sys.total.background, sum(char[which(char[,3]=="sys"),5]))
  }
  return(c(sum(non.total.background), sum(sys.total.background)))
}
for(p in 1:length(cancername)){
  background <- read.table(file = '/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/final_background.txt', header = T, sep = '\t')
  cancername <- colnames(background)[-c(1,2)]
  cancertype <- cancername[p]
  final.background <- read.table(file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/New-196-Mutation-Signature-Calculation/', cancertype,'.txt'),
                                 header = T, sep = '\t')
  final.result <- t(apply(as.matrix(cds),1,mutation.condition))
  colnames(final.result) <- c('non.background','sys.background')
  final.result <- cbind(as.matrix(header), final.result)
  write.table(final.result, file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/New-Background-Method/', cancertype, '-background-tmp.txt'), sep = '\t',
              quote = FALSE, row.names = FALSE)
  #final.result <- as.data.frame(final.result)
  #final.result$gene <- header$gene
  #write.table(final.result, file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Cancertype/', cancertype, '-background.txt'), sep = '\t',
  #            quote = FALSE, row.names = FALSE)
  print(cancertype)
  #save.image(file = '/Share/home/JiFansen/JiFansen/working-image/background-calculation.RData')
}
  
