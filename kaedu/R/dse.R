rm(list=ls())
  

DSE <- function(n,P,TMT){
  message("Diesner Simple Exchange")
  
  N <- length(P)
  periods <- length(TMT)/n
  
  period_A <- P[1:(N/periods)]
  period_B <- P[((N/periods)+1):N]
  differenceA <- NULL
  differenceB <- NULL

  treatmentA <- NULL
  treatmentB <- NULL
  
  UE_total <- NULL
  dif = 0
  #calculating treatments A and B
  for(i in 1: (N/periods) ){
      if(TMT[i] == "a"){
        treatmentA <- c(treatmentA, P[i])
        treatmentB <- c(treatmentB, P[i+n])
      }
      
      if(TMT[i] == "b"){
        treatmentB <- c(treatmentB, P[i])
        treatmentA <- c(treatmentA, P[i+n])
      }
    UE_total <- c(UE_total,P[i]+P[i+n])
  }
  
  #calculating vectors of diffeces
  for(i in 1: (N/periods) ){
    if(TMT[i] == "a"){
      dif = P[i]-P[n+i]
      differenceA <- c(differenceA, dif)
    }else{
      dif = P[i]-P[n+i]
      differenceB <- c(differenceB, dif)
    }
  }
  
  cat("treatment A:", treatmentA, "\n")
  cat("treatment B:", treatmentB, "\n")
  
  cat("diff A:", differenceA, "\n")
  cat("diff B:", differenceB, "\n")
  
  cat("sum P:", sum(P), "\n")
  cat("period A:", period_A, "\n")
  cat("period B:", period_B, "\n")
  cat("ue. total:", UE_total)
  treat = length(summary(factor(TMT)))
  #cat("sum: ", sum(period_A)," ",sum(period_B)," ",sum(P)," ",sum(treatmentA)," ",sum(treatmentB)," ",sum(differenceA)," ",sum(differenceB))
  GL_i = n-1
  GL_p = periods-1
  GL_tmt = treat-1
  GL_t =  N-1
  GL_e = (N-1)-(GL_i+GL_p+GL_tmt)
  #sum of squares
  SC_p = N
  SC_t = (sum(P^2) - ((sum(P)^2)/N))
  SC_i = ((sum(UE_total^2)/periods) - ((sum(P)^2)/N))
  SC_tmt =((  (sum(sum(treatmentA)^2+sum(treatmentB)^2))/n ) - ((sum(P)^2)/N)) 
  SC_e = (SC_t - SC_i-SC_tmt-SC_p)
  #square means
  CM_i = SC_i / GL_i
  CM_p = SC_p / GL_p
  CM_tmt = SC_tmt / GL_tmt
  CM_e = SC_e / GL_e
  # F 
  F_i = CM_i / CM_e
  F_p = CM_p / CM_e
  F_tmt = CM_tmt / CM_e
  
  #####################################################################
  message("Design Simple Exchange (DSE)  \n By -Yeison Eduardo")
  cat("--------------------------------ANOVA--------------------------------","\n")
  cat("   F.V.      d.f.        SC      CM              Fc    ","\n")
  cat(" I.......   ",GL_i,"\t",round(SC_i,4),"\t",round(CM_i,4),"\t",F_i,"\n")
  cat(" Periods    ",GL_p,"\t\t",round(SC_p,4),"\t",round(CM_p,4),"\t\t",F_p,"\n")
  cat(" Treatment  ",GL_tmt,"\t\t",round(SC_tmt,4),"\t",round(CM_tmt,4),"\t\t",F_tmt,"\n")
  cat(" Error      ",GL_e,"\t",round(SC_e,4),"\t",round(CM_e,4),"\t","--------","\n")
  cat(" Total      ",GL_t,"\t",round(SC_t,4),"\n")
  cat("---------------------------------------------------------------------","\n")
}

p <- c(21,13,24,13,14,11,17,23,24,28,8,18,11,33,28,20,10,26,+
       22,24,18,14,8,17,26,31,20,18,16,12,28,27,24,30,19,24)
tmt <- c("a","b","b","b","a","a","b","b","b","a","b","a","b","a","a","b","a","a","b","a","a","a","b","b","a","a","a","b","a","b","a","b","b","a","b","b")           
DSE(18,p,tmt)



