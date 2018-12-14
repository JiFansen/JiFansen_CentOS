Args <- commandArgs()
File <- as.character(Args[6])
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
