# Calculate incidence rate
data.summary <- function(markData, var){

  fit.data.cases <- subset(markData,eventInd==1)
  table.cases <- table(fit.data.cases[,var], fit.data.cases$tx)
  #the order is 0-1 for mark and 0-1 for armLabel
  #number of cases per person-years
  vaccinePersonYears <- sum(markData$eventTime[markData$tx==1])/365.25
  placeboPersonYears <- sum(markData$eventTime[markData$tx==0])/365.25
  
  vaccineIncRate <- round(table.cases[, "1"]/vaccinePersonYears*100,1)
  placeboIncRate <- round(table.cases[, "0"]/placeboPersonYears*100,1)
  
  vaccineCases <- table.cases[,"1"]
  placeboCases <- table.cases[,"0"]
  
  return(list(vaccineCases = vaccineCases, placeboCases = placeboCases, 
              vaccineIncRate = vaccineIncRate, 
              placeboIncRate = placeboIncRate))
}

# Identify the location in the string that matches the specified pattern
loc <- function(string, pattern){
  seq = seq(1:length(string))
  output = NULL
  for(i in 1:length(pattern)){
    output[i] = seq[string == pattern[i]]
  }
  return(output)
}

# Format VE estimate or the CI of VE to show percentages with one digit after zero
format.VE <- function(VE){
  if(as.numeric(VE) < (-10)){return("-Inf")}
  else{
    return(as.character(format(round(as.numeric(VE)*100,1),nsmall=1)))
  }}

#' Convert a p-value to a string and if it is less than 0.001, return '<0.001'
#' @param p An p-value
#' @example paste.p(0.1)
format.p <- function(p){if(as.numeric(p) < 0.001){return("< 0.001")}else{return(as.character(format(as.numeric(p),digits=2)))}}

#' locate the position of a pattern in a vector of string
loc <- function(stringVec, pattern){
  seq = 1:length(stringVec)
  return(seq[stringVec %in% pattern])
}



# Produce VE summary table for residue present/absent
table.VE.residue <- function(markResultList, marks, WestfallYoungAdjPvalues, forestplot = FALSE){
  #browser()
  if(forestplot){
    VE <- tibble(mark = character(), position = character(),residue = character(), cases = character(), VE = character(), mean = numeric(),
                 ll = numeric(), ul = numeric(), P = character(), 
                 diffP = character(), FWER = character(), Qvalue = character() )
    for(mark in marks){
      VEtable <- markResultList[[mark]]$VEtable
      y = strsplit(mark,split="\\.")[[1]]
      pos <- y[2]
      residue <- y[4]
      unadjp <- markResultList[[mark]]$diffP
      fdr = subset(WestfallYoungAdjPvalues, X==mark)$"p.FDR"
      holm = subset(WestfallYoungAdjPvalues, X==mark)$"p.FWER"
      Qvalue = format.p(fdr)
      if(fdr <= 0.2 & unadjp <= 0.05){Qvalue = paste0(Qvalue,"*")}
      
      VE <- add_row(VE, 
                    mark = mark,
                    position = pos,
                    residue = "",
                    cases = "",
                    VE = "",
                    mean = NA, ll = NA, ul = NA, P = "",
                    diffP = format.p( markResultList[[mark]]$diffP), 
                    FWER = format.p(holm), Qvalue = Qvalue)
      
      VE <- add_row(VE, 
                    mark = mark,
                    position = "",
                    residue = residue,
                    cases = VEtable$inc[1],
                    VE = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                    mean = VEtable$VE[1]*100, ll = VEtable$LB[1]*100, ul = VEtable$UB[1]*100,P = format.p(VEtable$p[1]), 
                    diffP = "", 
                    FWER = "",Qvalue = "")
      VE <- add_row(VE, 
                    mark = mark,
                    position = "",
                    residue = paste0 ("Not ", residue),
                    cases = VEtable$inc[2],
                    VE = paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")") ,
                    mean = VEtable$VE[2]*100, ll = VEtable$LB[2]*100, ul = VEtable$UB[2]*100,P = format.p(VEtable$p[2]), 
                    diffP = "", 
                    FWER = "", Qvalue = "")
  }}else{
    VE <- tibble(position = character(), residue = character(), cases1 = character(), VE1 = character(),P1 = character(), 
                 ws1 = character(),
                 cases0 = character(), VE0 = character(),P0 = character(), 
                 ws0 = character(),
                 diffP = character(), FWER = character(), Qvalue = character() )
    for(mark in marks){
      VEtable <- markResultList[[mark]]$VEtable
      y = strsplit(mark,split="\\.")[[1]]
      pos <- y[2]
      residue <- y[4]
      
      unadjp <- markResultList[[mark]]$diffP
      sig.unadj <- ifelse( unadjp<=0.05, 1, 0) 
      annotatedPosition = ifelse(sig.unadj, paste(pos,"$^{\\P}$",sep=""), pos)
      
      fdr = subset(WestfallYoungAdjPvalues, X==mark)$"p.FDR"
      holm = subset(WestfallYoungAdjPvalues, X==mark)$"p.FWER"
      sig.fdr <- ifelse(fdr<=0.2 & unadjp <=0.05, 1, 0) 
      sig.holm <- ifelse(holm<=0.05, 1, 0) 
      annotatedPosition = ifelse(sig.fdr, paste(annotatedPosition,"$^{\\S}$",sep=""), annotatedPosition)
      annotatedPosition = ifelse(sig.holm, paste(annotatedPosition,"$^{\\dag}$",sep=""), annotatedPosition)
      
      VE <- add_row(VE, 
                    position = annotatedPosition,
                    residue = residue,
                    cases1 = VEtable$inc[1],
                    VE1 = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                    P1 = format.p(VEtable$p[1]), 
                    ws1 = "",
                    cases0 = VEtable$inc[2], 
                    VE0 =  paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")"),
                    P0 = format.p(VEtable$p[2]), 
                    ws0 = "",
                    diffP = format.p( markResultList[[mark]]$diffP), 
                    FWER = format.p(holm), Qvalue = format.p(fdr)
      )
    }
  }
  
  return(VE) 
}

# Produce VE summary table for sequon present/absent
table.VE.sequon <- function(markResultList, marks, WestfallYoungAdjPvalues, forestplot = FALSE){
  #browser()
  if(forestplot){
    VE <- tibble(mark = character(), position = character(), feature = character(), cases = character(), VE = character(),mean = numeric(), ll = numeric(), ul = numeric(),
                 P = character(), diffP = character(), FWER = character(), Qvalue = character() )
    #browser()
    for(mark in marks){
      VEtable <- markResultList[[mark]]$VEtable
      pos <- strsplit(mark,split="\\.")[[1]][2]
      unadjp <- markResultList[[mark]]$diffP
      fdr = subset(WestfallYoungAdjPvalues, X==mark)$"p.FDR"
      holm = subset(WestfallYoungAdjPvalues, X==mark)$"p.FWER"
      Qvalue = format.p(fdr)
      if(fdr <= 0.2 & unadjp <= 0.05){Qvalue = paste0(Qvalue,"*")}
      
      
      VE <- add_row(VE, mark = mark, position = pos, feature = "", cases = "", VE = "" , mean = NA, ll = NA, ul = NA, P = "",
                    diffP = format.p(unadjp), FWER = format.p(holm), Qvalue = Qvalue)
      
      VE <- add_row(VE, mark = mark, position = "", feature  = paste0("PNGS"), cases = VEtable$inc[1],
                    VE = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                    
                    mean = VEtable$VE[1]*100, ll = VEtable$LB[1]*100, ul = VEtable$UB[1]*100,
                    P = format.p(VEtable$p[1]), diffP = "", FWER = "", Qvalue = "")
      VE <- add_row(VE,mark = mark,  position = "", feature = paste0("No PNGS"), cases = VEtable$inc[2],
                    VE = paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")") ,
                    mean = VEtable$VE[2]*100, ll = VEtable$LB[2]*100, ul = VEtable$UB[2]*100,
                    P = format.p(VEtable$p[2]), diffP = "", FWER = "", Qvalue = "")
    }
  }else{
    VE <- tibble(position = character(), cases1 = character(), VE1 = character(),P1 = character(), 
                 ws1 = character(),
                 cases0 = character(), VE0 = character(),P0 = character(), 
                 ws0 = character(),
                 diffP = character(), FWER = character(), Qvalue = character() )
    #browser()
    for(mark in marks){
      VEtable <- markResultList[[mark]]$VEtable
      pos <- strsplit(mark,split="\\.")[[1]][2]
      unadjp <- markResultList[[mark]]$diffP
      sig.unadj <- ifelse( unadjp<=0.05, 1, 0) 
      annotatedPosition = ifelse(sig.unadj, paste(pos,"$^{\\P}$",sep=""), pos)
      
      fdr = subset(WestfallYoungAdjPvalues, X==mark)$"p.FDR"
      holm = subset(WestfallYoungAdjPvalues, X==mark)$"p.FWER"
      sig.fdr <- ifelse(fdr<=0.2 & unadjp <=0.05, 1, 0) 
      sig.holm <- ifelse(holm<=0.05, 1, 0) 
      annotatedPosition = ifelse(sig.fdr, paste(annotatedPosition,"$^{\\S}$",sep=""), annotatedPosition)
      annotatedPosition = ifelse(sig.holm, paste(annotatedPosition,"$^{\\dag}$",sep=""), annotatedPosition)
      
      VE <- add_row(VE, 
                    position = annotatedPosition,
                    cases1 = VEtable$inc[1],
                    VE1 = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                    P1 = format.p(VEtable$p[1]), 
                    ws1 = "",
                    cases0 = VEtable$inc[2], 
                    VE0 =  paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")"),
                    P0 = format.p(VEtable$p[2]), 
                    ws0 = "",
                    diffP = format.p( markResultList[[mark]]$diffP), 
                    FWER = format.p(holm), Qvalue = format.p(fdr)
      )
    }
  }
  
  return(VE) 
}

# Produce VE summary table for sequon present/absent without adjusted p-values
table.VE.sequon.wo.adj <- function(markResultList, marks, forestplot = FALSE){
  #browser()
  if(forestplot){
    VE <- tibble(mark = character(), position = character(), feature = character(), cases = character(), VE = character(),mean = numeric(), ll = numeric(), ul = numeric(),
                 P = character(), diffP = character())
    #browser()
    for(mark in marks){
      VEtable <- markResultList[[mark]]$VEtable
      pos <- strsplit(mark,split="\\.")[[1]][2]
      unadjp <- markResultList[[mark]]$diffP
  
      VE <- add_row(VE, mark = mark, position = pos, feature = "", cases = "", VE = "" , mean = NA, ll = NA, ul = NA, P = "",
                    diffP = format.p(unadjp))
      
      VE <- add_row(VE, mark = mark, position = "", feature  = paste0("PNG Present"), cases = VEtable$inc[1],
                    VE = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                    
                    mean = VEtable$VE[1]*100, ll = VEtable$LB[1]*100, ul = VEtable$UB[1]*100,
                    P = format.p(VEtable$p[1]), diffP = "")
      VE <- add_row(VE,mark = mark,  position = "", feature = paste0("PNG Absent"), cases = VEtable$inc[2],
                    VE = paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")") ,
                    mean = VEtable$VE[2]*100, ll = VEtable$LB[2]*100, ul = VEtable$UB[2]*100,
                    P = format.p(VEtable$p[2]), diffP = "")
    }
  }else{
    VE <- tibble(position = character(), cases1 = character(), VE1 = character(),P1 = character(), 
                 ws1 = character(),
                 cases0 = character(), VE0 = character(),P0 = character(), 
                 ws0 = character(),
                 diffP = character())
    #browser()
    for(mark in marks){
      VEtable <- markResultList[[mark]]$VEtable
      pos <- strsplit(mark,split="\\.")[[1]][2]
      unadjp <- markResultList[[mark]]$diffP
      sig.unadj <- ifelse( unadjp<=0.05, 1, 0) 
      annotatedPosition = ifelse(sig.unadj, paste(pos,"$^{\\P}$",sep=""), pos)
      
     
      VE <- add_row(VE, 
                    position = annotatedPosition,
                    cases1 = VEtable$inc[1],
                    VE1 = paste0(format.VE(VEtable$VE[1]), " (", format.VE(VEtable$LB[1]),", ",format.VE(VEtable$UB[1]),")") ,
                    P1 = format.p(VEtable$p[1]), 
                    ws1 = "",
                    cases0 = VEtable$inc[2], 
                    VE0 =  paste0(format.VE(VEtable$VE[2]), " (", format.VE(VEtable$LB[2]),", ",format.VE(VEtable$UB[2]),")"),
                    P0 = format.p(VEtable$p[2]), 
                    ws0 = "",
                    diffP = format.p( markResultList[[mark]]$diffP)
      )
    }
  }
  
  return(VE) 
}

