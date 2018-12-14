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
####################################### Read the mutation count. #####################################################################
i <- 1
for(mutation.sample in list.files('/Share/home/JiFansen/JiFansen/Mutation-Correction/Correction/', 
                                  pattern = '*-Score-Correction.txt')){
  mutation.count <- read.table(file = paste0('/Share/home/JiFansen/JiFansen/Mutation-Correction/Correction/', mutation.sample),
                               header = T, sep = '\t')
  if(i == 1){
    final.file <- mutation.count
    i <- i+1
  }else{
    final.file <- rbind.data.frame(final.file, mutation.count)
    i <- i+1
  }
  print(i)
}
final.file <- final.file[,-3]
stat.summary <- table(final.file$sample)
match.result <- match(substr(colnames(expression.table), 1,16), names(stat.summary))
big.index <- which(is.na(match.result)==FALSE)
small.index <- match.result[big.index]
expression.mutation <- expression.table[, big.index]
stat.summary <- stat.summary[small.index]
###################################################################################################################################
clinical <- read.csv(file = '/fshare2/JiFansen/Clinical/three_year_data.csv', header = T)
clinical <- clinical[,c(1,2,3,10,12,13)]
clinical$X.1 <- gsub('\\.','-',as.character(clinical$X.1))
index <- match(clinical$X.1, colnames(expression.mutation))
immuno.clinical <- clinical[which(is.na(index)==FALSE), ]
index <- index[which(is.na(index)==FALSE)]
immuno.expression <- expression.mutation[, index]
immuno.mutation <- stat.summary[index]
wilcox.test(as.numeric(immuno.mutation)~as.factor(immuno.clinical$vital_status), exact = FALSE)$p.value
data.summary <- data.frame(mutation.count = as.numeric(immuno.mutation), group = as.factor(immuno.clinical$vital_status))
library(ggplot2)
library(ggbeeswarm)
p <- ggplot(data.summary, aes(x=group, y=mutation.count)) + 
  geom_violin(aes(fill=group), trim=FALSE)+ theme_classic()
p + geom_quasirandom()+ theme_classic()

aggregate(x=data.summary[,1], by = list(data.summary$group), FUN = median)

p+ylim(c(0,100))
data.summary.new <- data.summary[which(data.summary$mutation.count<=100), ]
#################################################################################################################################################
nonimmuno.clinical <- read.csv(file = '/fshare2/JiFansen/Clinical/40candidatesgpredict_nonimmune_lable.csv', header = T)
nonimmuno.clinical$X <- gsub('\\.', '-', as.character(nonimmuno.clinical$X))
index <- match(nonimmuno.clinical$X, colnames(expression.mutation))
nonimmuno.clinical <- nonimmuno.clinical[which(is.na(index)==FALSE), ]
index <- index[which(is.na(index)==FALSE)]
nonimmuno.expression <- expression.mutation[, index]
nonimmuno.mutation <- stat.summary[index]
wilcox.test(as.numeric(nonimmuno.mutation)~as.factor(nonimmuno.clinical$lab), exact = FALSE)$p.value
data.summary <- data.frame(mutation.count = as.numeric(nonimmuno.mutation), group = as.factor(nonimmuno.clinical$lab))
aggregate(x=data.summary[,1], by = list(data.summary$group), FUN = median)
p <- ggplot(data.summary, aes(x=group, y=mutation.count)) + 
  geom_violin(aes(fill=group), trim=FALSE)+ theme_classic()
p + geom_quasirandom(alpha = 0.1)+ theme_classic()+ylim(c(0,150))
###################################################################################################################################################
allen_mutation <- read.csv(file = '/Share/home/JiFansen/JiFansen/TCGA-Project/mutation_load/TestData/Allen_mutation.csv', header = T)
mis.index <- which(allen_mutation$Variant_Classification%in%c('Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del', 'In_Frame_Ins', 'Missense_Mutation',
                                                               'Nonsense_Mutation', 'Nonstop_Mutation', 'Splice_Site', 'Start_Codon_Ins', 'Stop_Codon_Del'))
allen_mutation <- allen_mutation[mis.index,]
candidate.genes <- read.csv(file = '/Share/home/JiFansen/pathways/TBCP_candidgenes.csv', header = T)
candidate.genes <- as.character(candidate.genes[,2])
candidate.index <- which(allen_mutation$Hugo_Symbol%in%candidate.genes)
allen_mutation <- as.data.frame(table(allen_mutation$patient))
colnames(allen_mutation)[1] <- 'sample'
allen_clinical <- read.csv(file = '/Share/home/JiFansen/JiFansen/TCGA-Project/mutation_load/TestData/Allen_clinical.csv', header = T)
allen_clinical <- allen_clinical[,c(1,5)]
colnames(allen_clinical)[c(1,2)] <- c('sample', 'os')
colnames(allen_mutation)[2] <- 'non.sys'
#############################################################################################################################################
hugo_mutation <- read.csv(file = '/Share/home/JiFansen/JiFansen/TCGA-Project/mutation_load/TestData/Hugo_clinical_new.csv', header = T)
hugo_clinical <- read.csv(file = '/Share/home/JiFansen/JiFansen/TCGA-Project/mutation_load/TestData/Hugo_Clinical.csv', header = T)
hugo_mutation <- data.frame(sample = hugo_mutation$Patient.ID, non.sys = hugo_mutation$TotalNonSyn)
hugo_clinical <- data.frame(sample = hugo_clinical$Patient.ID, OS = hugo_clinical$Overall.Survival)
hugo_clinical <- hugo_clinical[-28, ]
colnames(hugo_clinical)[2] <- 'os'
hugo_clinical <- hugo_clinical[-8,]
hugo_mutation <- hugo_mutation[-8,]
#############################################################################################################################################
naiyer_mutation <- read.csv(file = '/Share/home/JiFansen/JiFansen/TCGA-Project/mutation_load/TestData/Naiyer_mutation.csv', header = T)
naiyer_clinical <- read.csv(file = '/Share/home/JiFansen/JiFansen/TCGA-Project/mutation_load/TestData/Naiyer_clinical.csv', header = T)
naiyer_mutation <- data.frame(sample = naiyer_clinical$Study.ID, non.sys = naiyer_clinical$Nonsyn.)
naiyer_clinical <- data.frame(sample = naiyer_clinical$Study.ID, OS = naiyer_clinical$PFS..mos.)
naiyer_clinical$OS <- 30*naiyer_clinical$OS
colnames(naiyer_clinical)[2] <- 'os'
##############################################################################################################################################
Synder_mutation <- read.csv(file = '/Share/home/JiFansen/JiFansen/TCGA-Project/mutation_load/TestData/Synder_mutation.csv', header = T)
Synder_clinical <- read.csv(file = '/Share/home/JiFansen/JiFansen/TCGA-Project/mutation_load/TestData/Synder_clinical.csv', header = T)
Synder_mutation <- as.data.frame(table(Synder_mutation$Sample))
colnames(Synder_mutation)[1] <- 'sample'
Synder_clinical <- Synder_clinical[which(Synder_clinical$Biopsy=='r'), ]
Synder_clinical$OS <- 365*Synder_clinical$OS
match.result <- match(Synder_mutation$sample, Synder_clinical$StudyID)
big.index <- which(is.na(match.result)==FALSE)
small.index <- match.result[big.index]
Synder_mutation <- Synder_mutation[big.index, ]
Synder_clinical <- Synder_clinical[small.index, ]
Synder_clinical <- Synder_clinical[,-3]
colnames(Synder_clinical)[c(1,2)] <- c('sample', 'os')
colnames(Synder_mutation)[2] <- 'non.sys'
################################################################################################################################################
allen_clinical$type <- 'allen'
hugo_clinical$type <- 'hugo'
naiyer_clinical$type <- 'naiyer'
Synder_clinical$type <- 'synder'
os.summary <- rbind.data.frame(allen_clinical, hugo_clinical, naiyer_clinical, Synder_clinical)
ggplot(os.summary, aes(x=os, col = type, fill = type)) +geom_histogram(binwidth=10)+theme_classic()
os.summary$status <- 0
os.summary$status[which(os.summary$os<365)] <- 1
os.summary$status <- as.factor(os.summary$status)
mutation.summary <- rbind.data.frame(allen_mutation, hugo_mutation, naiyer_mutation, Synder_mutation)
wilcox.test(as.numeric(mutation.summary$non.sys)~as.factor(os.summary$status))$p.value
data.summary <- cbind.data.frame(os.summary, non.sys = mutation.summary[,2])
p <- ggplot(data.summary, aes(x=status, y=non.sys)) + 
  geom_violin(aes(fill=status), trim=FALSE)+ theme_classic()
p + geom_quasirandom(alpha = 0.5)+ theme_classic()+ylim(c(0,2000))
