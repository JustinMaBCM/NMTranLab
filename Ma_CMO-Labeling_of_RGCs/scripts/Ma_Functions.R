CR0.8_Assign <- function(DF, CMOnum = 6){
  DFtemp <- DF[,c(2:(CMOnum+1),(CMOnum+3):(CMOnum+6), CMOnum+2)]
  DFtemp$CR0.8 <- 0
  for (i in 1:nrow(DFtemp)){
    if (DFtemp$Assignment_Probability[i] >= 0.8) {DFtemp$CR0.8[i] <- colnames(DFtemp)[which(match(DFtemp[i,1:(CMOnum+2)],DFtemp[i,(CMOnum+4)]) == 1)]
    } else {DFtemp$CR0.8[i] <- "Unassigned"}
  }
  return(DFtemp$CR0.8)
}


