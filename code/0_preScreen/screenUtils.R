#produce a character string with screened sites 
cellPaste <- function(screenedIn, celllength = 65){
  
  positions <- screenedIn$position
  digits <- nchar(as.character(positions))
  totalDigits <- sum(digits)
  cumSumDigits <- cumsum(digits)
  
  digitSeq <- seq(1, totalDigits, celllength)
  if(max(digitSeq)< totalDigits){digitSeq <- c(digitSeq, totalDigits)}
  seq <- 1
  for(i in 2:length(digitSeq)){
    seq <- c(seq, seq(1, length(positions),1)[cumSumDigits == max(cumSumDigits[cumSumDigits<=digitSeq[i]])])
  }
  outputString <- paste(positions[1:(seq[2])],collapse=", ")
  if(length(seq) > 2){
    for(i in 2:(length(seq)-1)){
      outputString <- paste0(outputString,  paste(positions[(seq[i]+1):(seq[i+1])],collapse=", "), sep= ", ")
    }
  }
  if(max(seq) < length(positions)){
    outputString <- paste0(outputString,  paste(positions[(seq[i]+1):(length(positions))],collapse=", "), sep = ", ")
    
  }
  return(outputString)
}



#produce tables for pre-screening presence/absence of residues
featureScreen.is.aa.print <- function(sieveData, featureVars, newFeatureName, insert, tableDir, fileTag){
  
  featureData <- subset(dplyr::select(sieveData, all_of(c("subjid","armdesc", "hiv1event", featureVars))), hiv1event == 1)
  colnames(featureData ) <- c("subjid","armDesc", "event", newFeatureName)
  
  featureLongFormat <- tidyr::pivot_longer(featureData, cols = all_of(newFeatureName), 
                                           names_to = "position", 
                                           values_to = "present",
                                           values_drop_na = TRUE)
  
  freqsTable <- plyr::ddply(featureLongFormat, .(armDesc, position, present), 
                            function(df){
                              length(unique(df$subjid))
                            })
  
  table1 <- tibble(armDesc = freqsTable$armDesc,
                   position = freqsTable$position,
                   present = freqsTable$present,
                   "nCases" = freqsTable$V1)
  
  table2 <- tibble("position" = character(), "present" = numeric(), "nTxPoolCases" = numeric())
  
  d_ply(table1, .(position), function(df){
    for(pre.k in c(0,1)){
      nTxPoolCases <- ifelse(is.null(subset(df, present == pre.k)), 0, sum(df$nCases[df$present == pre.k]))
      table2 <<- add_row(.data = table2, "position" = df$position[1], "present" = pre.k,
                         "nTxPoolCases" = nTxPoolCases)
    }
  })
  
  cutoff = 4
  table2$ind <- as.factor(ifelse(table2$nTxPoolCases>=cutoff, paste0(">= ", cutoff), paste0("< ", cutoff)))
  table2$ind <- factor(table2$ind, levels = c(paste0(">= ", cutoff), paste0("< ", cutoff)))
  table2$plotY <- factor(table2$position, levels = rev(newFeatureName))
  
  table2 <- mutate(table2,matchString = case_when(
    present == 1 ~ "Present",
    present == 0 ~ "Absent",
  ))

  table2$matchString <- factor(table2$matchString, levels = c("Present", "Absent"))
  
  #Screened in comparison table
  posScreenedIn <- subset(ddply(table2,.(position), function(df){sum(df$nTxPoolCases>=cutoff)==2}), V1 == TRUE)
  posScreenedIn$position <- sort(factor(posScreenedIn$position, levels = newFeatureName))
  table3 <- tibble("Insert" = character(), "m" = numeric(), "Positions Screened in" = character())
  
  table3 <- add_row(.data = table3, "Insert" = insert, "m" = length(posScreenedIn$position), 
                    "Positions Screened in" = cellPaste(posScreenedIn,celllength = 80))
  write.csv(table3, file.path(tableDir, paste0(insert, "posScreenedIn", fileTag, ".csv")), row.names = FALSE)
  
  x = data.frame(vars = featureVars, newName = newFeatureName)
  write.csv(x$vars[x$newName %in% posScreenedIn$position], file.path(tableDir, paste0(insert, "posVarScreenedIn", fileTag, ".csv")), row.names = FALSE)
  
}

#produce tables for pre-screening presence/absence of sequons
featureScreen.presenceVSabsence <- function(sieveData, featureVars, newFeatureName, tableDir, figureDir, fileTag){
  #browser()
  featureData <- subset(dplyr::select(sieveData, all_of(c("subjid","armdesc", "hiv1event",featureVars))), hiv1event == 1)
  colnames(featureData ) <- c("subjid","armDesc", "event",newFeatureName)
  
  featureLongFormat <- tidyr::pivot_longer(featureData, cols = all_of(newFeatureName), 
                                           names_to = "position", 
                                           values_to = "presence",
                                           values_drop_na = TRUE)
  
  freqsTable <- plyr::ddply(featureLongFormat, .(armDesc, position, presence), 
                            function(df){
                              length(unique(df$subjid))
                            })
  
  table1 <- tibble(armDesc = freqsTable$armDesc,
                   position = freqsTable$position,
                   presence = freqsTable$presence,
                   "nCases" = freqsTable$V1)
  
  table2 <- tibble("position" = character(), "presence" = numeric(), "nTxPoolCases" = numeric())
  
  d_ply(table1, .(position), function(df){
    for(pre.k in c(0,1)){
      nTxPoolCases <- ifelse(is.null(subset(df, presence == pre.k)), 0, sum(df$nCases[df$presence == pre.k]))
      table2 <<- add_row(.data = table2, "position" = df$position[1], "presence" = pre.k,
                         "nTxPoolCases" = nTxPoolCases)
    }
  })
  
  cutoff = 4
  table2$ind <- as.factor(ifelse(table2$nTxPoolCases>=cutoff, paste0(">= ", cutoff), paste0("< ", cutoff)))
  table2$ind <- factor(table2$ind, levels = c(paste0(">= ", cutoff), paste0("< ", cutoff)))
  table2$plotY <- factor(table2$position, levels = rev(newFeatureName))
  
  table2 <- mutate(table2,presenceString = case_when(
    presence == 1 ~ "Presence",
    presence == 0 ~ "Absence",
  ))
  
  table2$presenceString <- factor(table2$presenceString, levels = c("Presence", "Absence"))
  
  #Screened in comparison table
  posScreenedIn <- subset(ddply(table2,.(position), function(df){sum(df$nTxPoolCases>=cutoff)==2}), V1 == TRUE)
  posScreenedIn$position <- sort(factor(posScreenedIn$position, levels = newFeatureName))
  table3 <- tibble("m" = numeric(),"Positions Screened in" = character())
  table3 <- add_row(.data = table3, "m" = length(posScreenedIn$position), "Positions Screened in" = paste0(posScreenedIn$position, collapse=","))
  
  
  write.csv(table3, file.path(tableDir, paste0("sequonPosScreenedIn", fileTag, ".csv")), row.names = FALSE)
  x = data.frame(vars = featureVars, newName = newFeatureName)
  write.csv(x$vars[x$newName %in% posScreenedIn$position], file.path(tableDir, paste0("sequonPosVarScreenedIn",fileTag, ".csv")), row.names = FALSE)
}
