##### I am updating this Jan 2020
# This script is for determining whether burnin has completed
rm(list = ls())
library(ggplot2)
library(reshape2)
library(dplyr)
# Input the file you care about
system('mkdir logPlots')
# Read in the header
allFiles <- system('ls *.log', intern = TRUE)
for (foo in allFiles){
  fileName <- foo
  #create subdirectory to save plots to 
  subDir <- paste0('subdir',foo)
  system(sprintf('mkdir %s', subDir))
  command1 <- sprintf("head -n 1 %s > tempFile", fileName)
  system(command1)
  cPointCommand <- read.table(file = 'tempFile', header = FALSE)
  cPointCommand <- t(cPointCommand)
  # catch strings
  flags <- c("-n", "-ng")
  # Process header to retrieve numSims and numGp:
  for (i in 1:length(cPointCommand)){
    if (cPointCommand[i] == "-n"){
      numSims <- cPointCommand[i+1]
      numSims <- as.integer(numSims)
    } else if (cPointCommand[i] == "-ng") {
      
      numGp <- cPointCommand[i+1]
      numGp <- as.integer(numGp)
    }
  }
  # find out the alphabet size
  command2 <- sprintf("grep 'alphsize' %s | cut -d'=' -f2", fileName)
  alphsize <- as.integer(system(command2, intern=T))
  seqLenStr <- sprintf("grep 'Sequence length' %s | cut -d'=' -f2", fileName)
  seqLen <- as.integer(system(seqLenStr, intern=T)) 

  command3 <- sprintf("grep -n 'Beginning MCMC' %s | grep -Eo '^[^:]+'", 
    fileName)

  startLine <- system(command3, intern = TRUE)
  df <- read.table(fileName, header = FALSE, nrows = numSims, sep = " ", 
    skip = startLine)
  
  # Process the data frame
  # Drop first colum (only contains indecies)
  df <- df[-1]
  # Drop the last column
  df <- df[-length(df)]
  # name the columns
  var1 <- "numCP"
  mixProps <- sprintf("mixProp%02d", 0:(numGp-1))
  alphas <- rep(sprintf("alphaGp%02d_",0:(numGp-1)), alphsize)
  alphas <- matrix(alphas, nrow = alphsize, byrow = TRUE)
  alphabetChar <- sprintf("%02d",0:(alphsize-1))
  
  for (i in 1:numGp){
    alphas[,i] <- paste0(alphas[,i],alphabetChar)
  }
  
  temp <- c(alphas)
  
  lnLike <- "logLikelihood"
  
  headers <- c(var1, mixProps, temp, lnLike)
  
  colnames(df) <- headers
  # create a Df to be read for model selection later
  save_df <- mutate(df, alphsize=alphsize,seqLen=seqLen)
  tableName <- paste0(strsplit(fileName,'\\.')[[1]][1],'_df.txt')
  write.table(save_df,file=tableName)

  # The mixture proportions will be given by the
  # next ng entries.
  # plots
  iterSeq <- seq(1:numSims)
  # number of change points
  ggplot(data = df, aes(x = seq(1:numSims), y = numCP, group = 1)) + 
    geom_line() + xlab('Iteration') + ylab('Number of change points') + 
    ggtitle(sprintf('NumCp, ng=%s',numGp)) + theme_minimal()

  ggsave('NumCp.pdf', plot = last_plot(), path = subDir)
  # df[,2] to df[,2+(numGp-1)] will be the mix Props
  mixPropdf <- df[grep('mixProp', colnames(df))]
  mixPropMeans <- colMeans(tail(mixPropdf,n=500))
  write.table(mixPropMeans, file = paste(fileName,'MixPropMeans.txt',sep = ''))
  mixPropdf$iter <- iterSeq
  mixPropdf <- melt(mixPropdf, id.vars = 'iter')
  
  ggplot(mixPropdf, aes(x = iter, y = value, colour = variable)) + 
    geom_line() + xlab('Iteration') + ylab('Mixture Proportions') +
    ggtitle(sprintf('Mixture proportions Time Series, ng=%s', numGp)) + 
    theme_minimal()

  ggsave('MixProps.pdf', plot = last_plot(), path = subDir)
  
  # trace plot for alpha vector per group
  gpIndex = sprintf('%02d_',0:(numGp-1))
  for (i in gpIndex){
    print(i)
    gpDf <- df[grep(sprintf('alphaGp%s',i), colnames(df))]

    gpDf$iter <- iterSeq

    gpDf <- melt(gpDf, id.vars = 'iter')
    
    ggplot(gpDf, aes(x = iter, y = value, colour = variable)) +
      geom_line() + xlab('Iteration') + ylab('Param Value') + 
      ggtitle(sprintf('Pram Time Series, ng=%s', numGp)) + theme_minimal() +
      theme(legend.position = "none")

    ggsave(sprintf('%sGp.pdf',i), plot = last_plot(), path = subDir)
  }
  # trace plot for logLikelihood 
  ggplot(data = df, aes(x = iterSeq, y = logLikelihood)) +
    geom_line() + xlab('Iteration') + ylab('Log Likelihood') + 
    ggtitle(sprintf('Ln Likelihood Time Series, ng=%s',numGp)) + theme_minimal()
  
  logPlotName <- sprintf("LogLikelihood_%02d.pdf", as.numeric(numGp))
  ggsave(logPlotName, plot = last_plot(), path = subDir)
  moveLogPlot <- sprintf('cp ./%s/%s ./logPlots',subDir,logPlotName)
  system(moveLogPlot)
}

system('mkdir logPlotDirs')
system('mv subdir* ./logPlotDirs')