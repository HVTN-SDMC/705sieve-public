# --------------------------------------------------------------------------- #
# Perform pre-screening for all binary features using minimum variability filter
# The criteria is with ≥ 4 cases at each feature level in their single analyzed AA sequence
# Input:   705 primary endpoints,  sequence features specified in the sieve SAP that are to be adjusted in multiple comparison
# Output:  Tables with marks pass the minimum variability filter


# refresh the workspace
rm(list=ls(all=TRUE))


# Setting directory paths -------------------------------------------------
here::i_am("README.md")
repoDir <- here::here()
dataDir <- file.path(repoDir, "data")
codeDir <- file.path(repoDir, "code/0_preScreen")
figureDir <- file.path(repoDir, "figures")
tableDir <- file.path(repoDir, "tables")

# Initialization
library(tidyverse)
library(plyr)
# load supplement
source(file.path(codeDir, "screenUtils.R"))
source(file.path(repoDir, "code/common.R"))

# load input data
sieveData <- read.csv(file.path(dataDir, datFile))

######################################################################################################
# Screen residues at AA locations in tier 1 
posIsAA <- colnames(sieveData)[grepl(".is.",colnames(sieveData))&grepl("hxb2",colnames(sieveData)) &grepl("tier1",colnames(sieveData))& (!grepl("is.sequon.tier1",colnames(sieveData)))]
posIsAAlist <- strsplit(posIsAA, split = ".", fixed = TRUE)
posIsPosLetter <- plyr::ldply(posIsAAlist , function(list){c(list[2],list[4])})
uniquePos <- unique(posIsPosLetter$V1)
posFreq <- plyr::ldply(uniquePos, function(x){sum(posIsPosLetter$V1==x)})
posFreq <- data.frame(HXB2.Position = uniquePos, freq = posFreq$V1)
subInd <- rep(TRUE, length(posIsAA))
for(i in 1:length(posIsAA)){
  pos = posIsPosLetter$V1[i]
  residue = posIsPosLetter$V2[i]
  #remove single dominating residues and if there are two residues at a specific position, remove one of the residues
  posFreqsub <- filter(posFreq, HXB2.Position==pos)
  residues <- filter(posIsPosLetter, V1 == pos)
  if(posFreqsub$freq==1){
    subInd[i] = FALSE
  }
  if(posFreqsub$freq==2){
    if(residue == residues$V2[2]){
      subInd[i] = FALSE
    }
    
  }
  
}
subPosIsAA <- posIsAA[subInd]
subPosIsAAlist <- strsplit(subPosIsAA, split = ".", fixed = TRUE)
posIsAANewName <- plyr::laply(subPosIsAAlist , function(list){paste0(list[2], ".",list[4])})
posIsAAmatrix <- data.frame(plyr::laply(subPosIsAAlist , function(list){c(list[2],list[4])}))
colnames(posIsAAmatrix) <- c("position","residue")
featureScreen.is.aa.print(sieveData, subPosIsAA , posIsAANewName, "posIsAAtier1",tableDir, "")



# Screen residues for AA locations in tier 2 
posIsAA <- colnames(sieveData)[grepl(".is.",colnames(sieveData))&grepl("hxb2",colnames(sieveData)) & (!grepl("tier1",colnames(sieveData)))& (!grepl("is.sequon.tier1",colnames(sieveData)))]
posIsAAlist <- strsplit(posIsAA, split = ".", fixed = TRUE)
#remove residues that are vaccine insert in vaccine insert region
posIsPosLetter <- plyr::ldply(posIsAAlist , function(list){c(list[2],list[4])})

uniquePos <- unique(posIsPosLetter$V1)
posFreq <- plyr::ldply(uniquePos, function(x){sum(posIsPosLetter$V1==x)})
posFreq <- data.frame(HXB2.Position = uniquePos, freq = posFreq$V1)
subInd <- rep(TRUE, length(posIsAA))
for(i in 1:length(posIsAA)){
  pos = posIsPosLetter$V1[i]
  residue = posIsPosLetter$V2[i]
  #remove exactly complimentary residues and single dominating residues
  posFreqsub <- filter(posFreq, HXB2.Position==pos)
  residues <- filter(posIsPosLetter, V1 == pos)
  if(posFreqsub$freq==1){
    subInd[i] = FALSE
  }
  if(posFreqsub$freq==2){
    if(residue == residues$V2[2]){
      subInd[i] = FALSE
    }
    
  }
  
}
subPosIsAA <- posIsAA[subInd]
subPosIsAAlist <- strsplit(subPosIsAA, split = ".", fixed = TRUE)
posIsAANewName <- plyr::laply(subPosIsAAlist , function(list){paste0(list[2], ".",list[4])})
posIsAAmatrix <- data.frame(plyr::laply(subPosIsAAlist , function(list){c(list[2],list[4])}))
colnames(posIsAAmatrix) <- c("position","residue")
featureScreen.is.aa.print(sieveData, subPosIsAA , posIsAANewName, "posIsAAtier2", tableDir, "")

######################################################################################################
# Screen sequons
seqonAbsence <- colnames(sieveData)[grepl("is.sequon.tier1",colnames(sieveData))]
seqonAbsenceList <- strsplit(seqonAbsence, split = ".", fixed = TRUE)
seqonAbsenceNewName <- plyr::laply(seqonAbsenceList , function(list){paste0(list[2])})
featureScreen.presenceVSabsence(sieveData, seqonAbsence, seqonAbsenceNewName, tableDir, figureDir, "")
