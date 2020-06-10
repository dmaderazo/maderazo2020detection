args <- commandArgs(trailingOnly = TRUE)

if (length(args)==0) {
		stop('Must provide input classified file')
	} 

library(jagsUI)
library(coda)
modelString = "
model{
    for(i in 1:N){ # N is the number of individuals
      for(j in 1:K){ # K is the number of classifiers
        
        # I
        class[j,i] ~ dbern(params[j,diseasepl[i]])
      }
      # Same as disease but with offset 1, instead of 0
      diseasepl[i] <- disease[i]+1
      # Actual disease status ie true classification
      disease[i] ~ dbern(pd)
    }
    # Distributions of classifier parameters
    for(j in 1:K){
      temp[j,1] ~ dunif(0,1)
      temp[j,2] ~ dunif(0,1)

      # Parameter for the distribution of the probability of getting a false negative
      params[j,1] <- 1 - sqrt(temp[j,1]);
      # Parameter for the distribution of the probability of getting a true positive
      params[j,2] <- 1 - sqrt(temp[j,1]) + temp[j,2]*sqrt(temp[j,1]);
      
      # Sensitivity
      sens[j] <- params[j,2]
      # Specificity
      spec[j] <- 1-params[j,1]
    }
    
    # Distribution of prevalence
    pd ~ dunif(0,1)
}
# ... end of JAGS model specification
"

# Write the modelString to a file, using R commands:
writeLines(modelString,con="binaryClassModel.txt")

fn <- args[1]

mainInfo <- unlist(strsplit(fn,'_'))[1:2]
mainStr <- paste(mainInfo[1],mainInfo[2],sep='_')

estimateName <- paste(mainStr,'estimate.txt',sep='_')
outName <- paste(mainStr,'out.txt',sep='_')
mydf <- read.csv(fn)

dataMat <-t( data.matrix(mydf))
numRows <- nrow(mydf) #1000
numCols <- ncol(mydf) #3 

tfbsClassifications <- list(N=numRows,K=numCols,class=dataMat)


#data = list("swineFlu$N", "swineFlu$K", "swineFlu$class")
# inits

inits <- function(){list(
              pd=0.5,temp=structure(.Data = rep(0.5,2*numCols),
               .Dim=c(numCols,2)))}

parameters = c("params", "sens", "spec","pd")

#This runs the model
jagsModel <- autojags(model.file = "binaryClassModel.txt", 
	data=tfbsClassifications,
	inits=inits,n.chains=4, iter.increment=502500,parameters.to.save=parameters,
	n.burnin=2500,n.thin=2000,parallel =TRUE, max.iter = 1000000)

sink(estimateName)
jagsModel
sink()

save.image(file=paste(mainStr,'dat.RData',sep='_'))
clean_file = 'rm -rf binaryClassModel.txt'
system(clean_file)

# Now evaluate the classifiers

# Creaete an array to record the combinations that make up each classifier. 
# Will be a 2^N*N matrix of 1s and 0s, where a 0 in the {i,j}th position indicates to 
# intersect with classifier j in combination i.

comboarray <- NULL

# you better hope these have at least 1000 rows
# or this will shit itself
sens <- jagsModel$sims.list$sens
spec <- jagsModel$sims.list$spec
N <- ncol(sens)


for(i in 1:2^N-1){
	newrow <- NULL
	for(j in 1:N){
		newrow <- cbind(newrow,(1-(-1)^floor(i/(2^(j-1))))/2)
	}
	comboarray <- rbind(comboarray,newrow)
}

# Calculate the sensitivity and specificity of 'classifiers' 
# created by intersecting actual classifiers with each other
# Initialise sensitivity and specificity vectors
sensitivity <- array(1,dim=c(2^N,500))
specificity <- array(0,dim=c(2^N,500))

for(i in 1:2^N){
	for(j in 1:N){
		for(k in 501:1000){
			sensitivity[i,(k-500)] <- sensitivity[i,(k-500)]*((1-sens[k,j])^comboarray[i,j])*(sens[k,j]^(1-comboarray[i,j]))
		}
	}

	for(j in 1:N){
		for(k in 501:1000){
			specificity[i,(k-500)] <- specificity[i,(k-500)]+((1-spec[k,j])^comboarray[i,j])*(spec[k,j]^(1-comboarray[i,j]))-specificity[i,(k-500)]*((1-spec[k,j])^comboarray[i,j])*(spec[k,j]^(1-comboarray[i,j]))
		}
	}
}

# Calculate sensitivities and specificities of unions of the intersections of classifiers
# Create another comboarray, where a 0 in position [i,j] indicates not to include intersection j in union i.
comboarray2 <- NULL
for(i in 1:2^(2^N)-1){
	newrow <- NULL
	for(j in 1:2^N){
		newrow <- cbind(newrow,(1-(-1)^floor(i/(2^(j-1))))/2)
	}
	comboarray2 <- rbind(comboarray2,newrow)
}

# Initialise arrays for the sensitivity and specificity of these unions of intersections of classifiers
sensitivity2 <- array(0,dim=c(2^(2^N),500))
specificity2 <- array(1,dim=c(2^(2^N),500))

# Calculate the sensitivities and specificities
for(k in 1:500){
	for(i in 1:(2^(2^N))){
		for(j in 1:(2^N)){
			sensitivity2[i,k] <- sensitivity2[i,k]+sensitivity[j,k]*comboarray2[i,j]
		}

		for(j in 1:2^N){
			specificity2[i,k] <- specificity2[i,k]+(specificity[j,k]-1)*comboarray2[i,j]
		}
	}
}

# Pick out the best combination based on each of four different ranking methods

# First method: Product of Sensitivity and Specificity
bestcombo <- array(0,dim=c(1,500))
bestrank1 <- array(0,dim=c(1,500))

for(i in 1:(2^(2^N))){
	for(k in 1:500){
		if(sensitivity2[i,k]*specificity2[i,k] > bestcombo[k]){
			bestcombo[k] <- sensitivity2[i,k]*specificity2[i,k]
			bestrank1[k] <- i
		}
	}
}

# Second method: Sum of Squares of Sensitivity and Specificity
bestcombo <- array(0,dim=c(1,500))
bestrank2 <- array(0,dim=c(1,500))

for(i in 1:(2^(2^N))){
	for(k in 1:500){
		if(sensitivity2[i,k]^2+specificity2[i,k]^2 > bestcombo[k]){
			bestcombo[k] <- (sensitivity2[i,k]^2+specificity2[i,k]^2)
			bestrank2[k] <- i
		}
	}
}

# Third method: Sum of Absolute Values of Sensitivity and Specificity
bestcombo <- array(0,dim=c(1,500))
bestrank3 <- array(0,dim=c(1,500))

for(i in 1:(2^(2^N))){
	for(k in 1:500){
		if(abs(sensitivity2[i,k])+abs(specificity2[i,k]) > bestcombo[k]){
			bestcombo[k] <- (abs(sensitivity2[i,k])+abs(specificity2[i,k]))
			bestrank3[k] <- i
		}
	}
}

# Fourth method: Minimum of Sensitivity and Specificity
bestcombo <- array(0,dim=c(1,500))
bestrank4 <- array(0,dim=c(1,500))

for(i in 1:(2^(2^N))){
	for(k in 1:500){
		if(min(sensitivity2[i,k],specificity2[i,k]) > bestcombo[k]){
			bestcombo[k] <- min(sensitivity2[i,k],specificity2[i,k])
			bestrank4[k] <- i
		}
	}
}

# Record the combinations (this is how you calculate mode in R)
best1 <- as.numeric(names(which.max(table(bestrank1))))
best2 <- as.numeric(names(which.max(table(bestrank2))))
best3 <- as.numeric(names(which.max(table(bestrank3))))
best4 <- as.numeric(names(which.max(table(bestrank4))))


comboname1 <- ""
for(i in 1:ncol(comboarray2)){
	comboname1 <- paste(comboname1,comboarray2[best1,i],sep="")
}

comboname2 <- ""
for(i in 1:ncol(comboarray2)){
	comboname2 <- paste(comboname2,comboarray2[best2,i],sep="")
}

comboname3 <- ""
for(i in 1:ncol(comboarray2)){
	comboname3 <- paste(comboname3,comboarray2[best3,i],sep="")
}

comboname4 <- ""
for(i in 1:ncol(comboarray2)){
	comboname4 <- paste(comboname4,comboarray2[best4,i],sep="")
}


# Determine the probability that the mode is correct for each method
chancecorrect1 <-  as.numeric(table(bestrank1)[names(which.max(table(best1)))])/500
chancecorrect2 <-  as.numeric(table(bestrank2)[names(which.max(table(best2)))])/500
chancecorrect3 <-  as.numeric(table(bestrank3)[names(which.max(table(best3)))])/500
chancecorrect4 <-  as.numeric(table(bestrank4)[names(which.max(table(best4)))])/500

# Print the results

sink(outName)
print(paste("According to Method 1 (Product of Sensitivity and Specificity) there is a probability of ",chancecorrect1, " that combination ", best1, " (", comboname1, ")", " is the best combination, with median sensitivity ", median(sensitivity2[best1,]), " and median specificity ", median(specificity2[best1,]), ".", sep=""))
print(paste("mean sensitivity ", mean(sensitivity2[best1,]), " and mean specificity ", mean(specificity2[best1,]), ".", sep=""))
print("This combination is a union of the following intersections (given as a binary code.)")
for(i in 1:2^N){
	if(comboarray2[best1,i] == 1){
		print(comboarray[i,])
	}
}

print(paste("According to Method 2 (Sum of Squares of Sensitivity and Specificity) there is a probability of ",chancecorrect2, " that combination ", best2, " (", comboname2, ")", " is the best combination, with median sensitivity ", median(sensitivity2[best2,]), " and median specificity ", median(specificity2[best2,]), ".", sep=""))
print(paste("mean sensitivity ", mean(sensitivity2[best2,]), " and mean specificity ", mean(specificity2[best2,]), ".", sep=""))
print("This combination is a union of the following intersections (given as a binary code.)")
for(i in 1:2^N){
	if(comboarray2[best2,i] == 1){
		print(comboarray[i,])
	}
}

print(paste("According to Method 3 (Sum of Absolute Values of Sensitivity and Specificity) there is a probability of ",chancecorrect3, " that combination ", best3, " (", comboname3, ")", " is the best combination, with median sensitivity ", median(sensitivity2[best3,]), " and median specificity ", median(specificity2[best3,]), ".", sep=""))
print(paste("mean sensitivity ", mean(sensitivity2[best3,]), " and mean specificity ", mean(specificity2[best3,]), ".", sep=""))
print("This combination is a union of the following intersections (given as a binary code.)")
for(i in 1:2^N){
	if(comboarray2[best3,i] == 1){
		print(comboarray[i,])
	}
}

print(paste("According to Method 4 (Minimum of Sensitivity and Specificity) there is a probability of ",chancecorrect4, " that combination ", best4, " (", comboname4, ")", " is the best combination, with median sensitivity ", median(sensitivity2[best4,]), " and median specificity ", median(specificity2[best4,]), ".", sep=""))
print(paste("mean sensitivity ", mean(sensitivity2[best4,]), " and mean specificity ", mean(specificity2[best4,]), ".", sep=""))
print("This combination is a union of the following intersections (given as a binary code.)")
for(i in 1:2^N){
	if(comboarray2[best4,i] == 1){
		print(comboarray[i,])
	}
}
sensOut <- paste(mainStr,'sens2.txt',sep='_')
specOut <- paste(mainStr,'spec2.txt',sep='_')
write.table(sensitivity2,file = sensOut,sep = ",")
write.table(specificity2,file = specOut,sep = ",")
