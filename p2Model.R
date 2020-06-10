# phase 2 of model selection
library(tidyverse)
# read in df_File of model that we currently think is best

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0){
	stop("supply input file: chr##_ng_##_df.txt")
}

inFile <- args[1]
workDf <- read.table(inFile,header=T) %>% tail(,n=500)
alphSize <- mean(workDf$alphsize)
numGps <- length(grep("mixProp",colnames(workDf)))
# for each group in model
for (i in 0:(numGps-1)){ #0 indexed
	#isolate relevant entries for model
	consGp <- c('06','27')
	padGp <- str_pad(i,2,pad="0")
	padChar <- str_pad(0:(alphSize-1),2,pad="0")
	tempStr <- sprintf("alphaGp%s_",padGp)
	consHeadings <- paste0(tempStr,consGp)
	gpHeadings <- paste0(tempStr,padChar)

	consDf <- select(workDf,consHeadings)
	gpDf <- select(workDf,gpHeadings)

	plotDf <- data.frame(consProp=rowSums(consDf)/rowSums(gpDf),
		iter=1:nrow(workDf))
	# gpDf <- mutate(gpDf,iter=1:nrow(gpDf))
	# gpDf <- mutate(gpDf,consProp=sum(consHeadings)/sum(gpHeadings))

	ggplot(plotDf, aes(x=iter,y=consProp,group=1))+
		geom_line() + xlab('Iteration') + ylab('consProp') + 
		ggtitle(sprintf('consProp, Gp=%s',padGp)) + theme_minimal()

	ggsave(sprintf('consProp_Gp%s.pdf',padGp), plot = last_plot())
}
nameStr <- paste0(numGps,'_cons_graphs')
system(sprintf('mkdir %s', nameStr))
system(sprintf('mv consProp*.pdf ./%s',nameStr))