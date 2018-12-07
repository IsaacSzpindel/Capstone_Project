
##################
# Pre-Processing #
##################

#Set the directory.
setwd('C:/Users/Zachary/Documents/Capstone/Analysis/Code')

#call on the following libraries.
library(caret)
library(pROC)

#Set the number of randomizations to run.
nrun = 100

#Read in the data.
treedata <- read.csv("Tree_Data_Final3.csv")

#Create empty data frames for the data.
RanData <- data.frame(Teak=integer(),
                      Eucalyptus=integer(),
                      Oak=integer(),
                      Pine=integer(),
                      Other_Conifers=integer(),
                      Mangrove_Tree=integer(),
                      Mangrove_Shrub=integer(),
                      Carbon_Absorb=integer(),
                      stringsAsFactors=FALSE)


CalcData <- data.frame(Teak=integer(),
                       Eucalyptus=integer(),
                       Oak=integer(),
                       Pine=integer(),
                       Other_Conifers=integer(),
                       Viable=integer(),
                       stringsAsFactors=FALSE)


#Calculate the amount of mangroves that are in salt water.
MTArea <- treedata$SWater[1]/2
MSArea <- treedata$SWater[1]/2

#Read how much land is left for forest due to water.
Left <- (1-treedata$SWater[1]-treedata$FWater[1])

for (i in 1:nrun){
  #Calculate the amount of land that is occupied by each species of tree.
  TArea <- runif(1,treedata$Teak[2],treedata$Teak[3])
  EArea <- runif(1,treedata$Eucalyptus[2],treedata$Eucalyptus[3])
  OArea <- runif(1,treedata$Oak[2],treedata$Oak[3])
  PArea <- runif(1,treedata$Pine[2],treedata$Pine[3])
  OCArea <- runif(1,treedata$Other_Conifers[2],treedata$Other_Conifers[3])
  
  five_ran <- TArea+EArea+OArea+PArea+OCArea
  
  TArea <- (TArea/five_ran)*Left
  EArea <- (EArea/five_ran)*Left
  OArea <- (OArea/five_ran)*Left
  PArea <- (PArea/five_ran)*Left
  OCArea <- (OCArea/five_ran)*Left
  
      #Calculate if the species is out of the ideal range.
      if (TArea>treedata$Teak[5]){
        TAreacalc = 0
      } else if (TArea<treedata$Teak[4]){
        TAreacalc = 0
      } else {
        TAreacalc = 1
      }
      
      if (EArea>treedata$Eucalyptus[5]){
        EAreacalc = 0
      } else if (EArea<treedata$Eucalyptus[4]){
        EAreacalc = 0
      } else {
        EAreacalc = 1
      }
      
      if (OArea>treedata$Oak[5]){
        OAreacalc = 0
      } else if (TArea<treedata$Oak[4]){
        OAreacalc = 0
      } else {
        OAreacalc = 1
      }
      
      if (PArea>treedata$Pine[5]){
        PAreacalc = 0
      } else if (TArea<treedata$Pine[4]){
        PAreacalc = 0
      } else {
        PAreacalc = 1
      }
      
      if (OCArea>treedata$Other_Conifers[5]){
        OCAreacalc = 0
      } else if (TArea<treedata$Other_Conifers[4]){
        OCAreacalc = 0
      } else {
        OCAreacalc = 1
      }
  
      #Calculate if the forest will be viable or not.
      if (TAreacalc+EAreacalc+OAreacalc+PAreacalc+OCAreacalc < 3){
        Viable = 0
      } else{
        Viable = 1
      }
  
  #Calculate the amount of carbon dioxide that is absorbed.
  Carbon_Absorb <- (TArea*treedata$Teak[1]+EArea*treedata$Eucalyptus[1]+OArea*treedata$Oak[1]+PArea*treedata$Pine[1]+OCArea*treedata$Other_Conifers[1]+MTArea*treedata$Mangrove_Tree[1]+MSArea*treedata$Mangrove_Shrub[1])*treedata$Years[1]
  
  #Add the data into its respective dataframe.
  TempGroup1 = c(TArea,EArea,OArea,PArea,OCArea,MTArea,MSArea,Carbon_Absorb)
  TempGroup2 = c(TAreacalc,EAreacalc,OAreacalc,PAreacalc,OCAreacalc,Viable)
  
  RanData <- rbind(RanData,TempGroup1)
  CalcData <- rbind(CalcData,TempGroup2)
}

#Set the names for the dataframes.
RanNames <- c("Teak", "Eucalyptus", "Oak", "Pine", "Other_Conifers", "Mangrove_Tree", "Mangrove_Shrub", "Carbon_Absorb")
colnames(RanData) <- RanNames

CalcNames <- c("Teak", "Eucalyptus", "Oak", "Pine", "Other_Conifers","Viable")
colnames(CalcData) <- CalcNames

####################
# Analyze the Data #
####################

set.seed(123)

#Test-train split
train.rows<- createDataPartition(y= CalcData$Viable, p=0.7, list = FALSE)
trainData<- CalcData[train.rows,]
testData<- CalcData[-train.rows,]

source("Jags-Ydich-XmetMulti-Mlogistic.R")

varList = colnames(CalcData)[1:5]

graphFileType = "pdf"

hitRates = list()

bayesLog <- function(trainData, testDataX, testDataY, parameter, Titles){
  print(parameter)
  parameterMatrix = as.matrix(parameter)
  testDataMatrix = as.matrix(testDataX)
  yHat = testDataMatrix%*%parameterMatrix[-1] + parameterMatrix[1]
  prob = 1/(1+exp(-yHat))
  roc_data <- data.frame(prob, testDataY)
  dev.new()
  plot(roc(roc_data$prob, roc_data$testDataY), main = Titles)
  binaryPred = ifelse(prob>0.5,1,0)
  hitRate = mean(binaryPred==testDataY)
  print(table(binaryPred,testDataY))
  print(cbind(hitRate))
  return(hitRate)
}

for(i in varList){
  
  yName = "Viable" ; xName = c(i)
  fileNameRoot = "output"
  numSavedSteps=15000 ; thinSteps=2
  
  mcmcCoda = genMCMC( data=trainData , xName=xName , yName=yName , 
                      numSavedSteps=numSavedSteps , thinSteps=thinSteps , 
                      saveName=fileNameRoot )
  #------------------------------------------------------------------------------- 
  # Display diagnostics of chain, for specified parameters:
  parameterNames = varnames(mcmcCoda) # get all parameter names
  for ( parName in parameterNames ) {
    diagMCMC( codaObject=mcmcCoda , parName=parName , 
              saveName=fileNameRoot , saveType=graphFileType )
  }
  #------------------------------------------------------------------------------- 
  # Get summary statistics of chain:
  summaryInfo = smryMCMC( mcmcCoda , 
                          saveName=fileNameRoot )
  show(summaryInfo)
  # Display posterior information:
  plotMCMC( mcmcCoda , data=trainData , xName=xName , yName=yName , 
            pairsPlot=TRUE , showCurve=FALSE ,
            saveName=fileNameRoot , saveType=graphFileType )
  #------------------------------------------------------------------------------- 
  
  parameter_file = read.csv(paste0(fileNameRoot,"SummaryInfo.csv"),row.names = 1)
  
  parameter = parameter_file[1:2,"Mode"]
  
  hitRates[i] = bayesLog(trainData, testData[[xName]], testData[[yName]], parameter, xName)
  
}

print(hitRates)

#Export out the area calculations for each species of tree, as well as the
#viability for each species of tree and the forest as a whole, in two separate
#csv files. The area calculations file also has the total carbon dioxide that
#would be absorbed by each forest.
write.csv(RanData, file = "AreaCalculations_Final3.csv")
write.csv(CalcData, file = "ViabilityCalculations_Final3.csv")
