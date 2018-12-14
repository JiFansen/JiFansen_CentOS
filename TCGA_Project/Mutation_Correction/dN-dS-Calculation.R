################################################Old Background Mutation Frequency Calculation. ###############################################
a <- read.csv(file = "/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/motiffreqtypes.csv", header = T)
col1<-as.vector(a[,1])
motif1<-strsplit(col1," ")
snp<-c()
conten<-c()
for(i in 1:96){
  a1<-motif1[[i]][1]
  a2<-motif1[[i]][2]
  snp<-c(snp,a1)
  conten<-c(conten,a2)
}
snp[which(snp=="CA")]="CAGT"
snp[which(snp=="CG")]="CGGC"
snp[which(snp=="CT")]="CTGA"
snp[which(snp=="TA")]="TAAT"
snp[which(snp=="TC")]="TCAG"
snp[which(snp=="TG")]="TGAC"
motif<-data.frame(snp,conten)
CHOH<-read.csv("TCGA.SKCM.varscan.txt",header = T, sep = '\t')
s_n<-length(unique(CHOH$Tumor_Sample_Barcode))
indata<-CHOH[,c(1,5,6,7,9,10,11,12,13,16,55,56,57,112)]
indata<-indata[which(indata$Variant_Type=="SNP"),]
tri_motif<-substr(indata$CONTEXT,5,7)
substr(tri_motif,2,2)<-"."
subsit<-paste0(indata$Reference_Allele,indata$Tumor_Seq_Allele2)
num<-c()
for (i in 1:96){
  ind1<-which(subsit==substr(motif[i,1],1,2))
  ind2<-which(subsit==substr(motif[i,1],3,4))
  subsit_match<-c(ind1,ind2)
  conten_match<-which(tri_motif==as.character(motif[i,2]))
  match_index<-intersect(subsit_match,conten_match)
  num<-c(num,length(match_index))
}
num<-num/sum(num)
out<-data.frame(snv,conten,num)
skcm <- c(out$num,out$num)
back <- read.table(file = "background_cancer.txt", sep = '\t', header = T)
back$SKCM <- skcm
write.table(back,file="final_background.txt",sep="\t")
################################################## New Methods Background Calculation. #####################################################
erik.background <- read.csv(file = '/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/Erik.csv', header = T, row.names = 1)
snv <- c()
content <- c()
for(i in 1:ncol(erik.background)){
  snv <- c(snv, paste0(unlist(strsplit(colnames(erik.background),""))[7*i-5], unlist(strsplit(colnames(erik.background),""))[7*i-1]))
  content <- c(content, paste0("(",unlist(strsplit(colnames(erik.background),""))[7*i-6], unlist(strsplit(colnames(erik.background),""))[7*i-4],")"))
}
snv.reverse <- rep('a',96)
snv.reverse[which(snv=='CA')] <- 'GT'
snv.reverse[which(snv=='CG')] <- 'GC'
snv.reverse[which(snv=='CT')] <- 'GA'
snv.reverse[which(snv=='TA')] <- 'AT'
snv.reverse[which(snv=='TC')] <- 'AG'
snv.reverse[which(snv=='TG')] <- 'AC'
content.reverse <- rep('a',96)
content.reverse[which(content=='(AA)')] <- '(TT)'
content.reverse[which(content=='(AC)')] <- '(GT)'
content.reverse[which(content=='(AG)')] <- '(CT)'
content.reverse[which(content=='(AT)')] <- '(AT)'
content.reverse[which(content=='(CA)')] <- '(TG)'
content.reverse[which(content=='(CC)')] <- '(GG)'
content.reverse[which(content=='(CG)')] <- '(CG)'
content.reverse[which(content=='(CT)')] <- '(AG)'
content.reverse[which(content=='(GA)')] <- '(TC)'
content.reverse[which(content=='(GC)')] <- '(GC)'
content.reverse[which(content=='(GG)')] <- '(CC)'
content.reverse[which(content=='(GT)')] <- '(AC)'
content.reverse[which(content=='(TA)')] <- '(TA)'
content.reverse[which(content=='(TC)')] <- '(GA)'
content.reverse[which(content=='(TG)')] <- '(CA)'
content.reverse[which(content=='(TT)')] <- '(AA)'
snv <- c(snv,snv.reverse)
content <- c(content, content.reverse)
background.table <- paste0(snv, content)
for(File in list.files("/Share/home/JiFansen/JiFansen/Mutation-Correction/MAF-summary/", pattern = '*-MAF-summary.txt')){
  CHOH <- read.csv(file = paste0("/Share/home/JiFansen/JiFansen/Mutation-Correction/MAF-summary/", File), sep = '\t',header = T)
  cancertype <- unlist(strsplit(File, '-'))[1]
  print(cancertype)
  Frequency <- c()
  CHOH.change <- as.character(CHOH$change)       
  for(i in 1:length(background.table)){
    Frequency <- c(Frequency, length(which(CHOH.change==background.table[i]))) 
  }
  Frequency <- Frequency/sum(Frequency)
  final.background <- data.frame(change = background.table, frequency = Frequency)
  write.table(final.background, file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/New-196-Mutation-Signature-Calculation/',cancertype,'.txt'),
              sep = '\t', quote = FALSE, row.names = FALSE)
}

################################################## Calculate abnormal mutations for each sample.############################################
for(File in list.files('/Share/home/lanxun5/Data/TCGA/MAF_TCGA/csv_format/MAF-File/')){
  cancertype <- unlist(strsplit(File,'\\.'))[2]
  if(cancertype=='SKCM'){
    CHOH <- read.csv(file = paste0('/Share/home/lanxun5/Data/TCGA/MAF_TCGA/csv_format/MAF-File/',File), header = T, sep = '\t', fileEncoding = "GBK")
  }else{
    CHOH <- read.csv(file = paste0('/Share/home/lanxun5/Data/TCGA/MAF_TCGA/csv_format/MAF-File/',File), header = T, fileEncoding = "GBK")
  }
  useful.name <- c("Hugo_Symbol","Variant_Classification","Variant_Type",
                   "Reference_Allele","Tumor_Seq_Allele2","Consequence",
                   "Amino_acids","Codons","VARIANT_CLASS","CONTEXT","Protein_position","Tumor_Sample_Barcode")
  data <- CHOH[,match(useful.name,colnames(CHOH))]
  data <- data[,c(1,4,5,7,8,10,11,2,3,6,9,12)]
  data$tumor.sample <- substr(as.character(data$Tumor_Sample_Barcode),1,16)
  data$tumor.sample <- substr(as.character(data$Tumor_Sample_Barcode),1,16)
  data.1 <- data[which(data$Variant_Classification=="Frame_Shift_Del"),]
  data.2 <- data[which(data$Variant_Classification=="Splice_Site"),]
  data.3 <- data[which(data$Variant_Classification=="Nonsense_Mutation"),]
  data.4 <- data[which(data$Variant_Classification=="Frame_Shift_Ins"),]
  data.5 <- data[which(data$Variant_Classification=="Splice_Region"),]
  data.8 <- data[which(data$Variant_Classification=="Nonstop_Mutation"),]
  data.9 <- data[which(data$Variant_Classification=="Translation_Start_Site"),]
  data.total <- rbind.data.frame(data.1,data.2,data.3,data.4,data.5,data.8,data.9)
  rm(data.1,data.2,data.3,data.4,data.5,data.8,data.9)
  test <- as.matrix(aggregate(as.character(data.total[,1]),by=list(data.total$tumor.sample,data.total$Variant_Classification), FUN=table))
  test <- as.data.frame(test)
  test$sum <- apply(test,1,FUN = function(x){sum(as.numeric(x[-c(1,2)]))})
  colnames(test)[1:2] <- c('Barcode','Mutation.Type')
  colnames(test)[-c(1,2,ncol(test))] <- unlist(lapply(strsplit(colnames(test)[-c(1,2,ncol(test))],'\\.'), FUN = function(x){x[2]}))
  write.table(test, file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Frameshift/',cancertype,'-shift.txt'), 
              quote = FALSE, row.names = FALSE)
}
################################################## MAF files summary preparation.#######################################################
strand.reverse <- function(x){
  result.2 <- as.character(x[5])
  result.3 <- paste0(as.character(x[2]),as.character(x[3]),'(',substr(as.character(x[6]),5,5),
                     substr(as.character(x[6]),7,7),")")
  results.4 <- as.character(x[6])
  forward <- unlist(strsplit(as.character(x[6]),""))
  reverse <- c()
  for (z in 1:length(forward)){
    if(forward[z]=="A"){
      reverse[z]="T"
    }
    if(forward[z]=="T"){
      reverse[z]="A"
    }
    if(forward[z]=="C"){
      reverse[z]="G"
    }
    if(forward[z]=="G"){
      reverse[z]="C"
    }
  }
  reverse <- paste(reverse[11],reverse[10],reverse[9],reverse[8],reverse[7],reverse[6],
                   reverse[5],reverse[4],reverse[3],reverse[2],reverse[1],sep="")
  mutation <- unlist(strsplit(as.character(x[5]),""))[which(unlist(
    strsplit(toupper(as.character(x[5])),""))==unlist(strsplit(as.character(x[5]),"")))]
  position <- which(unlist(strsplit(toupper(as.character(x[5])),""))==unlist(strsplit(as.character(x[5]),"")))[1]
  codon.mutation <- paste(mutation[1],mutation[3],sep="")
  return(c(result.2,result.3,results.4,reverse,codon.mutation,position,x[8],x[7],x[12]))
}
for(File in list.files('/Share/home/lanxun5/Data/TCGA/MAF_TCGA/csv_format/MAF-File/')){
  cancertype <- unlist(strsplit(File,'\\.'))[2]
  if(cancertype=='SKCM'){
    CHOH <- read.csv(file = paste0('/Share/home/lanxun5/Data/TCGA/MAF_TCGA/csv_format/MAF-File/',File), header = T, sep = '\t', fileEncoding = "GBK")
  }else{
    CHOH <- read.csv(file = paste0('/Share/home/lanxun5/Data/TCGA/MAF_TCGA/csv_format/MAF-File/',File), header = T, fileEncoding = "GBK")
  }
  useful.name <- c("Hugo_Symbol","Variant_Classification","Variant_Type",
                   "Reference_Allele","Tumor_Seq_Allele2","Consequence",
                   "Amino_acids","Codons","VARIANT_CLASS","CONTEXT","Protein_position","Tumor_Sample_Barcode")
  data <- CHOH[,match(useful.name,colnames(CHOH))]
  data <- data[which(data$VARIANT_CLASS=="SNV"),]
  data <- data[,c(1,4,5,7,8,10,11,2,3,6,9,12)]
  data.1 <- data[which(data$Variant_Classification=="Missense_Mutation"),]
  data.2 <- data[which(data$Variant_Classification=="Silent"),]
  missense.data <- rbind.data.frame(data.1,data.2)
  uniq.missence.genes <- as.character(unique(missense.data$Hugo_Symbol))
  final <- data.frame()
  for(genename in uniq.missence.genes){
    test <- missense.data[which(missense.data$Hugo_Symbol==genename),]
    result <- t(apply(test,1,strand.reverse))
    result <- data.frame(gene = genename, result)
    final <- rbind.data.frame(final,result)
  }
  colnames(final)[c(2,3,4,5,6,7)] <- c('codon','change','forward','reverse','codon.change','codon.position')
  results.10 <- character()
  for (i in 1:dim(final)[1]){
    reference.change <- substr(as.character(final[i,3]),1,2)
    codon.change <- as.character(final[i,6])
    if(reference.change==codon.change){
      results.10[i] <- "+"
    }
    else{
      results.10[i] <- "-"
    }
  }
  final$strand <- results.10
  write.table(final,file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/MAF-summary/', 
                                  cancertype,'-MAF-summary.barcode.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
}
################################################## Calculate significantly mutation genes. #####################################
codon <- read.table(file="/Share/home/JiFansen/JiFansen/Mutation-Correction/dN-dS/codon1.txt",sep="\t",header=F)
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
  string <- substr(as.character(x[11]),2,4)
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
  between <- unlist(strsplit(as.character(x[11]),""))
  for (z in 1:3){
    char5 <- c(char5, paste("(",between[1],between[3],")",sep=""))
  }
  for(z in 4:6){
    char5 <- c(char5, paste("(",between[2],between[4],")",sep=""))
  }
  for (z in 7:9){
    char5 <- c(char5, paste("(",between[3],between[5],")",sep=""))
  }
  char$d <- char5
  combination1 <- paste0(background$snv,background$context)
  combination2 <- paste0(char$b,char$d)
  match.result <- match(combination2,combination1)
  char$e <- as.numeric(background[,3])[match.result]
  return(c(sum(char[which(char[,3]=="non"),5]),sum(char[which(char[,3]=="sys"),5])))
}

#for(File in c('KIRP-MAF-summary.txt','LAML-MAF-summary.txt','LGG-MAF-summary.txt','LIHC-MAF-summary.txt','LUAD-MAF-summary.txt','LUSC-MAF-summary.txt')){
for(File in list.files('/Share/home/JiFansen/JiFansen/Mutation-Correction/dN-dS', pattern = '*MAF-summary.txt')){
  cancertype <- unlist(strsplit(File,'-'))[1]
  final <- read.table(file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/dN-dS/', File), header = T, sep = '\t')
  strand.plus <- final[which(final$strand=="+"),]
  strand.subtract <- final[which(final$strand=="-"),]
  index <- 2*(3-as.numeric(strand.plus$codon.position))+as.numeric(strand.plus$codon.position)
  m <- substr(as.character(strand.plus$forward),index,index+4)
  strand.plus$useful.context <- m
  index <- 2*(3-as.numeric(strand.subtract$codon.position))+as.numeric(strand.subtract$codon.position)
  m <- substr(as.character(strand.subtract$reverse),index,index+4)
  strand.subtract$useful.context <- m
  strand.plus <- rbind.data.frame(strand.plus,strand.subtract)
  background <- read.table(file = '/Share/home/JiFansen/JiFansen/Mutation-Correction/Background/final_background.txt', header = T)
  number <- which(colnames(background)==cancertype)
  background <- background[,c(1,2,number)]
  a <- unlist(strsplit(as.character(background[,2]),split="\\."))
  for (i in 1:192){
    b <- paste("(",a[2*i-1],a[2*i],")",sep = "")
    background[i,4] <- b
  }
  colnames(background)[4] <- "context"
  non.sys <- t(apply(strand.plus,1,mutation.condition))
  colnames(non.sys) <- c("non.background","sys.background")
  non.sys <- as.data.frame(non.sys)
  strand.plus <- cbind.data.frame(strand.plus, non.sys)
  combination3 <- paste0(background$snv,background$context)
  match.result <- match(as.character(strand.plus$change),combination3)
  value <- as.numeric(background[,3])[match.result]
  strand.plus$true.value <- value
  write.table(strand.plus, file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/OldMethods-non-sys/', 
                                         cancertype,'-non-sys.txt'), sep = '\t', quote = FALSE, row.names = FALSE)
  uniq.plus.genes <- unique(strand.plus$gene)
  final.results1 <- c();final.results2 <- c();final.results3 <- c();final.results4 <- c();final.results5 <- c()
  for (i in 1:length(uniq.plus.genes)){
    gene.name <- as.character(uniq.plus.genes[i])
    test.data <- strand.plus[which(strand.plus$gene==gene.name),]
    non.count <- length(which(test.data$Variant_Classification=="Missense_Mutation"))
    sys.count <- length(which(test.data$Variant_Classification=="Silent"))
    uniq.protein.position <- unique(test.data$Protein_position)
    non.background <- 0
    sys.background <- 0
    for (j in 1:length(uniq.protein.position)){
      protein.position <- uniq.protein.position[j]
      non.background <- non.background+test.data[which(test.data$Protein_position==protein.position)[1],12]
      sys.background <- sys.background+test.data[which(test.data$Protein_position==protein.position)[1],13]
    }
    final.results1 <- c(final.results1,gene.name)
    final.results2 <- c(final.results2, non.count)
    final.results3 <- c(final.results3, sys.count)
    final.results4 <- c(final.results4, non.background)
    final.results5 <- c(final.results5, sys.background)
  }
  final.results <- data.frame(final.results1,final.results2,final.results3,final.results4,final.results5)
  colnames(final.results) <- c("gene","non.count","sys.count","non.background","sys.background")
  final.results[,6] <- final.results[,4]/(final.results[,4]+final.results[,5])
  final.results[,7] <- final.results[,5]/(final.results[,4]+final.results[,5])
  final.results[,2] <- final.results[,2]+1
  final.results[,3] <- final.results[,3]+1
  final.results[,8] <- (final.results[,2]/final.results[,4])/(final.results[,3]/final.results[,5])
  colnames(final.results)[c(6,7,8)] <- c("non.ratio","sys.ratio","dN/ds")
  final.results[,9] <- final.results[,2]+final.results[,3]
  for(i in 1:dim(final.results)[1]){
    final.results[i,10] <- binom.test((final.results[i,2]-1),(final.results[i,9]-2),final.results[i,6],alternative = "greater")$p.value
  }
  
  colnames(final.results)[c(9,10)] <- c("total.count","p.value")
  write.table(final.results, file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/OldMethods-non-sys/', cancertype,'-final.results.txt'),
              sep = '\t', quote = FALSE, row.names = FALSE)
  result.1 <- final.results[which(final.results$`dN/ds`>1),]
  length(which(result.1$p.value<=0.05))
  result <- result.1[which(result.1$p.value<=0.05),]
  write.table(result, file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/OldMethods-non-sys/', cancertype,'-significant.results.txt'))
}