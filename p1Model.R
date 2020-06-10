# what this script should do:
# first phase of model selection. Making plots for AIC, BIC, DICV
library(tidyverse)
library("reshape2")
#find all the relevant files 
allFiles <- system('ls *_df.txt', intern = T)
#see how many there are
numFiles <- length(allFiles)
#create empty storage for output
df <- data.frame(numGp=integer(numFiles),AIC=double(numFiles),
	BIC=double(numFiles),DICV=double(numFiles))
#for each file
for (i in 1:numFiles){
	foo <- allFiles[i]
	#read in df last 500 iters
	workDf <- read.table(foo,header=T) %>% tail(,n=500)
	#compute mean lnLike
	lnL_bar <- mean(workDf$logLikelihood)
	#compute mean numCp
	numCp <- mean(workDf$numCP)
	#get alphsize
	alphSize <- mean(workDf$alphsize)
	#get numGps
	numGps <- length(grep("mixProp",colnames(workDf)))
	#get lengthSeq
	seqLen <- mean(workDf$seqLen)
	#calculate AIC sec 2.2 Oldmeadow2011
	AIC <- -2*lnL_bar+2*(numCp+numGps*(alphSize+1))
	#calculate BIC sec 2.3 Oldmeadow 2011
	BIC <- -2*lnL_bar+(numCp+numGps*(alphSize+1))*log(seqLen)
	#calculate DICV sec 2.4 Oldmeadow 2011
	DICV <- 0.5*var(-2*workDf$logLikelihood)+(-2*lnL_bar)
	#store into a row of df
	df[i,] <- c(numGps,AIC,BIC,DICV)
}

#arrange df by numGps
df <- arrange(df,numGp)
outDf <- select(df,numGp,DICV) 
outDf <- arrange(outDf,DICV)

write.table(outDf,file='DICVs_ordered.txt')



#plot them?
long_df <- melt(df,id='numGp')

pdf(file='infoCrit.pdf')
ggplot(data=df,
	aes(x=numGp,y=AIC)) +
	geom_line() + xlab('numGps') + ylab('Information Criterion') + 
	ggtitle('AIC') 

ggplot(data=df,
       aes(x=numGp,y=BIC)) +
  geom_line() + xlab('numGps') + ylab('Information Criterion') + 
  ggtitle('BIC') 

ggplot(data=df,
       aes(x=numGp,y=DICV)) +
  geom_line() + xlab('numGps') + ylab('Information Criterion') + 
  ggtitle('DICV') 

dev.off()