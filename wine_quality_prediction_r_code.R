graphics.off() 
rm(list=ls()) 
library(ggplot2)
library(ggpubr)
library(ks)
library(rjags)
library(runjags)
library(readr) 
setwd("C:/Users/Admin/Downloads/Master Of Analytics/Second Sem 2024/Applied Bayesian Statistics/Assignment3")
source("DBDA2E-utilities.R") 

#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================

smryMCMC_HD = function(  codaSamples , compVal = NULL,  saveName=NULL) {
  summaryInfo = NULL
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  paramName = colnames(mcmcMat)
  for ( pName in paramName ) {
    if (pName %in% colnames(compVal)){
      if (!is.na(compVal[pName])) {
        summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] , 
                                                          compVal = as.numeric(compVal[pName]) ))
      }
      else {
        summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] ) )
      }
    } else {
      summaryInfo = rbind( summaryInfo , summarizePost( paramSampleVec = mcmcMat[,pName] ) )
    }
  }
  rownames(summaryInfo) = paramName
  
  # summaryInfo = rbind( summaryInfo , 
  #                      "tau" = summarizePost( mcmcMat[,"tau"] ) )
  if ( !is.null(saveName) ) {
    write.csv( summaryInfo , file=paste(saveName,"SummaryInfo.csv",sep="") )
  }
  return( summaryInfo )
}

#===============================================================================

plotMCMC_HD = function( codaSamples , data , xName="x" , yName="y", preds = FALSE ,
                        showCurve=FALSE ,  pairsPlot=FALSE , compVal = NULL, 
                        saveName=NULL , saveType="jpg" ) {
  # showCurve is TRUE or FALSE and indicates whether the posterior should
  #   be displayed as a histogram (by default) or by an approximate curve.
  # pairsPlot is TRUE or FALSE and indicates whether scatterplots of pairs
  #   of parameters should be displayed.
  #-----------------------------------------------------------------------------
  y = data[,yName]
  x = as.matrix(data[,xName])
  mcmcMat = as.matrix(codaSamples,chains=TRUE)
  chainLength = NROW( mcmcMat )
  zbeta0 = mcmcMat[,"zbeta0"]
  zbeta  = mcmcMat[,grep("^zbeta$|^zbeta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { zbeta = matrix( zbeta , ncol=1 ) }
  beta0 = mcmcMat[,"beta0"]
  beta  = mcmcMat[,grep("^beta$|^beta\\[",colnames(mcmcMat))]
  if ( ncol(x)==1 ) { beta = matrix( beta , ncol=1 ) }
  if (preds){
    pred = mcmcMat[,grep("^pred$|^pred\\[",colnames(mcmcMat))]
  } # Added by Demirhan
  guess = mcmcMat[,"guess"]
  #-----------------------------------------------------------------------------
  # Compute R^2 for credible parameters:
  YcorX = cor( y , x ) # correlation of y with each x predictor
  Rsq = zbeta %*% matrix(YcorX,ncol=1)
  Rsq <- Rsq[,1]
  #-----------------------------------------------------------------------------
  if ( pairsPlot ) {
    # Plot the parameters pairwise, to see correlations:
    openGraph()
    nPtToPlot = 1000
    plotIdx = floor(seq(1,chainLength,by=chainLength/nPtToPlot))
    panel.cor = function(x, y, digits=2, prefix="", cex.cor, ...) {
      usr = par("usr"); on.exit(par(usr))
      par(usr = c(0, 1, 0, 1))
      r = (cor(x, y))
      txt = format(c(r, 0.123456789), digits=digits)[1]
      txt = paste(prefix, txt, sep="")
      if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
      text(0.5, 0.5, txt, cex=1.25 ) # was cex=cex.cor*r
    }
    pairs( cbind( beta0 , beta , tau )[plotIdx,] ,
           labels=c( "beta[0]" , 
                     paste0("beta[",1:ncol(beta),"]\n",xName) , 
                     expression(tau) ) , 
           lower.panel=panel.cor , col="skyblue" )
    if ( !is.null(saveName) ) {
      saveGraph( file=paste(saveName,"PostPairs",sep=""), type=saveType)
    }
  }
  #-----------------------------------------------------------------------------
  # Marginal histograms:
  
  decideOpenGraph = function( panelCount , saveName , finished=FALSE , 
                              nRow=2 , nCol=3 ) {
    # If finishing a set:
    if ( finished==TRUE ) {
      if ( !is.null(saveName) ) {
        saveGraph( file=paste0(saveName,ceiling((panelCount-1)/(nRow*nCol))), 
                   type=saveType)
      }
      panelCount = 1 # re-set panelCount
      return(panelCount)
    } else {
      # If this is first panel of a graph:
      if ( ( panelCount %% (nRow*nCol) ) == 1 ) {
        # If previous graph was open, save previous one:
        if ( panelCount>1 & !is.null(saveName) ) {
          saveGraph( file=paste0(saveName,(panelCount%/%(nRow*nCol))), 
                     type=saveType)
        }
        # Open new graph
        openGraph(width=nCol*7.0/3,height=nRow*2.0)
        layout( matrix( 1:(nRow*nCol) , nrow=nRow, byrow=TRUE ) )
        par( mar=c(4,4,2.5,0.5) , mgp=c(2.5,0.7,0) )
      }
      # Increment and return panel count:
      panelCount = panelCount+1
      return(panelCount)
    }
  }
  
  # Original scale:
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( beta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(beta[0]) , main="Intercept", compVal = as.numeric(compVal["beta0"] ))
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
    if (!is.na(compVal[paste0("beta[",bIdx,"]")])){
      histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx],
                           compVal = as.numeric(compVal[paste0("beta[",bIdx,"]")]))
    } else{
      histInfo = plotPost( beta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(beta[.(bIdx)]) , main=xName[bIdx])
    }
  }
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(R^2) , main=paste("Prop Var Accntd") , finished=FALSE )
  
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
  histInfo = plotPost( guess , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab="Guess parameter" , main=paste("Prop Var Accntd") , finished=TRUE )
  
  panelCount = 1
  if ( pred){
    
    for ( pIdx in 1:ncol(pred) ) {
      panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMarg") )
      histInfo = plotPost( pred[,pIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                           xlab=bquote(pred[.(pIdx)]) , main=paste0("Prediction ",pIdx) ) 
    }
  }# Added by Demirhan
  # Standardized scale:
  panelCount = 1
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( zbeta0 , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(z*beta[0]) , main="Intercept" )
  for ( bIdx in 1:ncol(beta) ) {
    panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
    histInfo = plotPost( zbeta[,bIdx] , cex.lab = 1.75 , showCurve=showCurve ,
                         xlab=bquote(z*beta[.(bIdx)]) , main=xName[bIdx] )
  }
  panelCount = decideOpenGraph( panelCount , saveName=paste0(saveName,"PostMargZ") )
  histInfo = plotPost( Rsq , cex.lab = 1.75 , showCurve=showCurve ,
                       xlab=bquote(R^2) , main=paste("Prop Var Accntd") )
  panelCount = decideOpenGraph( panelCount , finished=TRUE , saveName=paste0(saveName,"PostMargZ") )
  
  #-----------------------------------------------------------------------------
}

#===============PRELIMINARY FUNCTIONS FOR POSTERIOR INFERENCES====================
whitewine <- read_csv("whitewinetrain.csv")

whitewine$Column1 <- as.numeric(as.factor(whitewine$Column1)) - 1

summary(whitewine)
str(whitewine)
sum(is.na(whitewine))


logistic_model <- glm(Column1 ~ `fixed acidity` + `volatile acidity` + `citric acid` +
                        `residual sugar` + chlorides + `free sulfur dioxide` +
                        `total sulfur dioxide` + density + pH + sulphates + alcohol, 
                      data = whitewine, family = binomial)
logistic_model

# THE DATA.
y = as.numeric(whitewine$Column1)
x = as.matrix(whitewine[,2:12])
cat("\nCORRELATION MATRIX OF PREDICTORS:\n ")
show( round(cor(x),3) )
cat("\n")

PredData = read_csv("whitewinetest.csv")
head(PredData)
sum(is.na(PredData))
summary(PredData)

xPred = as.matrix(PredData[-1])

Nx = ncol(x)

# Specify the data in a list, for later shipment to JAGS:
dataList <- list(
  x = x ,
  y = y ,
  xPred = xPred ,
  Ntotal = length(y),
  Nx = Nx, 
  Npred = nrow(xPred)
)


# THE MODEL
modelString = "
data {
  for ( j in 1:Nx ) {
    xm[j]  <- mean(x[,j])
    xsd[j] <-   sd(x[,j])
    for ( i in 1:Ntotal ) {
      zx[i,j] <- ( x[i,j] - xm[j] ) / xsd[j]
    }
  }
}

model {
  # Informative run
  for ( i in 1:Ntotal ) {
    # In JAGS, ilogit is logistic:
    y[i] ~ dbern( mu[i] )
      mu[i] <- ( guess*(1/2) + (1.0-guess)*ilogit(zbeta0+sum(zbeta[1:Nx]*zx[i,1:Nx])) )
  }
   zbeta0 ~ dnorm( 0, 1/2^2 )
zbeta[1] ~ dnorm(0.0968, 1/(1/xsd[1]^2))     # fixed acidity 
zbeta[2] ~ dnorm(-2.3552, 1/(1/xsd[2]^2))    # volatile acidity 
zbeta[3] ~ dnorm(-0.8745, 1/(1/xsd[3]^2))    # citric acid 
zbeta[4] ~ dnorm(0.3952, 1/(1/xsd[4]^2))     # residual sugar 
zbeta[5] ~ dnorm(-2.3361, 1/(1/xsd[5]^2))    # chlorides 
zbeta[6] ~ dnorm(0.0084, 1/(1/xsd[6]^2))     # free sulfur dioxide 
zbeta[7] ~ dnorm(-0.001314, 1/(1/xsd[7]^2))  # total sulfur dioxide 
zbeta[8] ~ dnorm(-478.1, 1/(1/xsd[8]^2))     # density 
zbeta[9] ~ dnorm(2.1594, 1/(1/xsd[9]^2))     # pH 
zbeta[10] ~ dnorm(5.37, 1/(1/xsd[10]^2))     # sulphates 
zbeta[11] ~ dnorm(0.2856, 1/(1/xsd[11]^2))   # alcohol
   guess ~ dbeta(1.2,2.8)
   
  #Transform to original scale:
  beta[1:Nx] <- zbeta[1:Nx] / xsd[1:Nx]
  beta0 <- zbeta0 - sum( zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx] )

  # Compute predictions at every step of the MCMC
  for ( k in 1:Npred){
    pred[k] <- ilogit(beta0 + sum(beta[1:Nx] * xPred[k,1:Nx]))
  }
}
" # close quote for modelString
# Write out modelString to a text file
writeLines( modelString , con="TEMPmodel.txt" )


parameters = c( "zbeta0" , "beta0")
# parameters = c( "zbeta0")
# parameters = c( "beta0")
for ( i in 1:Nx){
  parameters = c(parameters, paste0("zbeta[",i,"]"), paste0("beta[",i,"]"))
  #parameters = c(parameters, paste0("zbeta[",i,"]"))
  # parameters = c(parameters, paste0("beta[",i,"]"))
}
for ( i in 1:nrow(xPred)){
  parameters = c(parameters, paste0("pred[",i,"]")) 
}

parameters = c(parameters, "guess")
adaptSteps = 5000  # Number of steps to "tune" the samplers
burnInSteps = 10000
nChains = 2
thinSteps = 50
numSavedSteps = 3000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

startTime = proc.time()
runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel.txt" ,
                        monitor=parameters  ,
                        data=dataList ,
                        #inits=initsList ,
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps ,
                        sample=numSavedSteps ,
                        thin=thinSteps , summarise=FALSE , plots=FALSE )
stopTime = proc.time()
duration = stopTime - startTime
show(duration)
codaSamples = as.mcmc.list( runJagsOut )


diagMCMC( codaSamples , parName="zbeta0" )
for ( i in 1:Nx){
  diagMCMC( codaSamples , parName=paste0("zbeta[",i,"]") )
}

diagMCMC( codaSamples , parName="beta0" )
for ( i in 1:Nx){
  diagMCMC( codaSamples , parName=paste0("beta[",i,"]") )
}
for ( i in 1:nrow(xPred)){
  diagMCMC( codaSamples , parName=paste0("pred[",i,"]") )
}
gelman.plot(as.mcmc.list(lapply(codaSamples, function(chain) chain[, "zbeta0", drop = FALSE])), 
            main = "zbeta0")
for (i in 1:Nx) {
  gelman.plot(as.mcmc.list(lapply(codaSamples, function(chain) chain[, paste0("zbeta[", i, "]"), drop = FALSE])), 
              main = paste0("zbeta[", i, "]"))
}
compVal <- data.frame("beta0" = 0, "beta[1]" = 0, "beta[2]" = 0, "beta[3]" = 0, "beta[4]" =  0,  "beta[5]" =  0, 
                      "beta[6]" =  0, "beta[7]" =  0, "beta[8]" =  0, "beta[9]" =  0, "beta[10]" =  0,"beta[11]" =  0,check.names=FALSE)

summaryInfo <- smryMCMC_HD( codaSamples = codaSamples , compVal = compVal )
print(summaryInfo)

plotMCMC_HD( codaSamples = codaSamples , data = whitewine, xName=c("fixed acidity",
                                                                   "volatile acidity",
                                                                   "citric acid",
                                                                   "residual sugar",
                                                                   "chlorides",
                                                                   "free sulfur dioxide",
                                                                   "total sulfur dioxide",
                                                                   "density",
                                                                   "pH",
                                                                   "sulphates",
                                                                   "alcohol") ,
             yName="Column1", compVal = compVal, preds = FALSE)

preds <- data.frame(ID = PredData[,1], PredProb = summaryInfo[26:1495,3], Quality = PredData[,1])

threshold <- 0.5 # summaryInfo[427,3]
preds[which(preds[,2]<threshold),3] <- 0
preds[which(preds[,2]>threshold),3] <- 1

predsSorted <- preds[order(preds$ID),]

actual <- PredData[,"Column1"]

# ============ Predictive check ============

confusionMatrix <- function(resp, pred){
  classRes <- data.frame(response = resp , predicted = pred)
  conf = xtabs(~ predicted + response, data = classRes)
  
  accuracy = sum(diag(conf))/sum(conf)
  accuracy
  precision = conf[1,1]/(conf[1,1]+conf[1,2])
  precision
  recall = conf[1,1]/(conf[1,1]+conf[2,1])
  recall
  Fscore = 2*((precision*recall)/(precision+recall))
  Fscore
  return(list(accuracy = accuracy, precision = precision, recall = recall, Fscore = Fscore, conf = conf))
}

# Confusion Matrix: approximately 80% of the predictions are accurate:
confusionMatrix(resp = actual$Column1, pred = predsSorted[,3])
Conf_Mat <- data.frame(confusionMatrix(resp = actual$Column1, pred = predsSorted[,3]))

# Accurate prediction %
(Conf_Mat[1,7]+Conf_Mat[4,7])/length(actual$Column1)

# Large MCMC run
parameters = c(parameters, "guess")
adaptSteps = 5000  # Number of steps to "tune" the samplers
burnInSteps = 10000
nChains = 2
thinSteps = 50
numSavedSteps = 5000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )

# Iteration 5 settings (smaller MCMC run)
parameters = c(parameters, "guess")
adaptSteps = 5000  # Number of steps to "tune" the samplers
burnInSteps = 10000
nChains = 2
thinSteps = 50
numSavedSteps = 3000
nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )


# JAGS model string for non-informative run
modelString = "
data {
  for ( j in 1:Nx ) {
    xm[j]  <- mean(x[,j])
    xsd[j] <-   sd(x[,j])
    for ( i in 1:Ntotal ) {
      zx[i,j] <- ( x[i,j] - xm[j] ) / xsd[j]
    }
  }
}

model {
  # Non-informative run with standardisation
  for ( i in 1:Ntotal ) {
    # In JAGS, ilogit is logistic:
    y[i] ~ dbern( mu[i] )
      mu[i] <- ( guess*(1/2) + (1.0-guess)*ilogit(zbeta0+sum(zbeta[1:Nx]*zx[i,1:Nx])) )
  }
  # Priors vague on standardized scale:
  zbeta0 ~ dnorm( 0 , 1/2^2 )
  # non-informative run
  for ( j in 1:Nx ) {
    zbeta[j] ~ dnorm( 0 , 1/2^2 )
  }
  #Transform to original scale:
  beta[1:Nx] <- zbeta[1:Nx] / xsd[1:Nx]
  beta0 <- zbeta0 - sum( zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx] )
  guess ~ dbeta(1,9)

  # Compute predictions at every step of the MCMC
  for ( k in 1:Npred){
    pred[k] <- ilogit(beta0 + sum(beta[1:Nx] * xPred[k,1:Nx]))
  }
}
" # close quote for modelString
# Write out modelString to a text file
writeLines( modelString , con="TEMPmodel.txt" )

# parameters = c( "zbeta0" , "beta0")
# parameters = c( "zbeta0")
parameters = c( "beta0")
for ( i in 1:Nx){
  # parameters = c(parameters, paste0("zbeta[",i,"]"), paste0("beta[",i,"]"))
  # parameters = c(parameters, paste0("zbeta[",i,"]"))
  parameters = c(parameters, paste0("beta[",i,"]"))
}
for ( i in 1:nrow(xPred)){
  parameters = c(parameters, paste0("pred[",i,"]")) 
}

# Further Thinning ####
#furtherThin <- 3
#thiningSequence <- seq(1,nrow(codaSamples[[1]]), furtherThin)
#newCodaSamples <- mcmc.list()
#for ( i in 1:nChains){
#  newCodaSamples[[i]] <- as.mcmc(codaSamples[[i]][thiningSequence,])
#}
# Run JAGS ####
startTime = proc.time()
runJagsOut <- run.jags( method="parallel" ,
                        model="TEMPmodel.txt" ,
                        monitor=parameters  ,
                        data=dataList ,
                        # inits=initsList ,
                        n.chains=nChains ,
                        adapt=adaptSteps ,
                        burnin=burnInSteps ,
                        sample=numSavedSteps ,
                        thin=thinSteps , summarise=FALSE , plots=FALSE )
stopTime = proc.time()
duration = stopTime - startTime
show(duration)
codaSamples = as.mcmc.list( runJagsOut )


# Diagnostics checks ####
diagMCMC( codaSamples , parName="zbeta0" )
for ( i in 1:Nx){
  diagMCMC( codaSamples , parName=paste0("zbeta[",i,"]") )
}
diagMCMC( codaSamples , parName="beta0" )
for ( i in 1:Nx){
  diagMCMC( codaSamples , parName=paste0("beta[",i,"]") )
}
#for ( i in 1:nrow(xPred)){
diagMCMC( codaSamples , parName=paste0("pred[",i,"]") )
#}


graphics.off()

# Posterior Distributions ####

compVal <- data.frame("beta0" = 0, "beta[1]" = 0, "beta[2]" = 0, "beta[3]" = 0, "beta[4]" =  0,  "beta[5]" =  0, 
                      "beta[6]" =  0, "beta[7]" =  0, "beta[8]" =  0, "beta[9]" =  0, "beta[10]" =  0,"beta[11]" =  0,check.names=FALSE)

summaryInfo <- smryMCMC_HD( codaSamples = codaSamples , compVal = compVal )
print(summaryInfo)

plotMCMC_HD( codaSamples = codaSamples , data = whitewine, xName=c("fixed acidity",
                                                                   "volatile acidity",
                                                                   "citric acid",
                                                                   "residual sugar",
                                                                   "chlorides",
                                                                   "free sulfur dioxide",
                                                                   "total sulfur dioxide",
                                                                   "density",
                                                                   "pH",
                                                                   "sulphates",
                                                                   "alcohol") ,
             yName="Column1", compVal = compVal, preds = FALSE)

preds <- data.frame(ID = PredData[,1], PredProb = summaryInfo[14:1483,3], Quality = PredData[,1])

threshold <- 0.5 # summaryInfo[427,3]
preds[which(preds[,2]<threshold),3] <- 0
preds[which(preds[,2]>threshold),3] <- 1

predsSorted <- preds[order(preds$ID),]

# ============ Predictive check ============

confusionMatrix <- function(resp, pred){
  classRes <- data.frame(response = resp , predicted = pred)
  conf = xtabs(~ predicted + response, data = classRes)
  
  accuracy = sum(diag(conf))/sum(conf)
  accuracy
  precision = conf[1,1]/(conf[1,1]+conf[1,2])
  precision
  recall = conf[1,1]/(conf[1,1]+conf[2,1])
  recall
  Fscore = 2*((precision*recall)/(precision+recall))
  Fscore
  return(list(accuracy = accuracy, precision = precision, recall = recall, Fscore = Fscore, conf = conf))
}

# Confusion Matrix: approximately 80% of the predictions are accurate:
confusionMatrix(resp = actual$Column1, pred = predsSorted[,3])
Conf_Mat <- data.frame(confusionMatrix(resp = actual$Column1, pred = predsSorted[,3]))

# Model for sensitivity analysis for non-informative (this model is run with the same settings as the initial non-informative run)
modelString = "
data {
  for ( j in 1:Nx ) {
    xm[j]  <- mean(x[,j])
    xsd[j] <-   sd(x[,j])
    for ( i in 1:Ntotal ) {
      zx[i,j] <- ( x[i,j] - xm[j] ) / xsd[j]
    }
  }
}

model {
  # Non-informative run with standardisation
  for ( i in 1:Ntotal ) {
    # In JAGS, ilogit is logistic:
    y[i] ~ dbern( mu[i] )
      mu[i] <- ( guess*(1/2) + (1.0-guess)*ilogit(zbeta0+sum(zbeta[1:Nx]*zx[i,1:Nx])) )
  }
  # Priors vague on standardized scale:
  zbeta0 ~ dnorm( 0 , 1/12^2 )
  # non-informative run
  for ( j in 1:Nx ) {
    zbeta[j] ~ dnorm( 0 , 1/12^2 )
  }
  #Transform to original scale:
  beta[1:Nx] <- zbeta[1:Nx] / xsd[1:Nx]
  beta0 <- zbeta0 - sum( zbeta[1:Nx] * xm[1:Nx] / xsd[1:Nx] )
  guess ~ dbeta(1,9)

  # Compute predictions at every step of the MCMC
  for ( k in 1:Npred){
    pred[k] <- ilogit(beta0 + sum(beta[1:Nx] * xPred[k,1:Nx]))
  }
}
" # close quote for modelString
# Write out modelString to a text file
writeLines( modelString , con="TEMPmodel.txt" )
