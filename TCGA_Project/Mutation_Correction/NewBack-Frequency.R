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
