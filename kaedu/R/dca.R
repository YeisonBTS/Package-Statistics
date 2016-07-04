rm(list=ls())

DCAsmD <- function(Y,TMT,SM){
  message("DCA with Sub samping ...")
  cat("suma de y:", sum(Y),"\n")
  cat("suma de tratamiento: ", sum(TMT),"\n")
  cat("suma de tratamiento: ",Y[8],"\n")

  nij <- c(SM)
  ni  <- c(rep(0,max(TMT)))
  Yij <- c(rep(0,length(TMT)))
  Yij.<- c(rep(0,length(TMT)))
  Yi..<- c(rep(0,max(TMT)))
  Yi..m <- c(rep(0,max(TMT)))

  #create the matrix with -1 all elements
  Matx  <- matrix(-1, length(TMT), max(SM))

  k = 1
  #conver the vectors to matrix
  for(i in 1: length(TMT)){
    for (j in 1: max(SM)) {
      if(j <= SM[i]){
        Matx[[i,j]] <- Y[k]
        k = k+1;
      }
    }
  }

 #calculating ni
  for(j in 1: max(TMT)){#j 3
    for (i in 1: length(TMT)) {#i 18
          if(TMT[i] == j){
            ni[j] = ni[j]+nij[i]
          }

    }
  }

 for (i in 1: max(TMT)) {
   cat("ni ", ni[i], "\n")
 }

  #calculating Yij and Yij.
  for(i in 1: length(TMT)){#18
    for (j in 1: max(SM)) {#4
      if(Matx[[i,j]] != -1){
        Yij[i] = Yij[i] + Matx[[i,j]]
      }
    }
    Yij.[i] = Yij[i]/nij[i]
  }

  for (i in 1: length(TMT)) {
    cat("Yij ", Yij[i]," Yij. ",Yij.[i],"\n")
  }

  #calculating Yi.. and Yi..m
  for (j in 1: max(TMT)) {#3
    for (i in 1: length(TMT)) {#18
      if(TMT[i] == j ){
        Yi..[j] = Yi..[j] + Yij[i]
      }
    }
    Yi..m[j] = Yi..[j]/ni[j]
  }

  for (i in 1: max(TMT)) {
    cat("Yi..", Yi..[i]," Yi..m ", Yi..m[i],"\n")
  }

  n.. = sum(nij)
  Y... = sum(Y)
  Y...M = Y.../n..
  cat("n.. ",n..,"\n")
  cat("Y... ",Y...,"\n")

  ##################ANOVA#############################
  #calculating degrees of freedom
  t = max(TMT)
  GLtrat = t-1
  GLerror = length(TMT)-t
  GLmuestreo = sum(nij)-length(TMT)
  GLtotal = sum(nij)-1

  cat("GLtrat ",GLtrat,"\n")
  cat("GLerror ",GLerror,"\n")
  cat("GLmuestreo ",GLmuestreo,"\n")
  cat("GLtotal ",GLtotal,"\n")

  #sum of squares
  sc1 = 0
  sc2 = 0
  for (i in 1: t) {
    sc1 = sc1 + ((Yi..[i]^2)/ni[i])
  }

  for (i in 1: length(SM)) {
    sc2 = sc2 + ((Yij[i]^2)/nij[i])
  }

  SCtrat = sc1 - ((Y...^2) / n..)
  SCerror = sc2 - sc1
  SCmuestreo = sum(Y^2) - sc2
  SCtotal = sum(Y^2) - ((Y...^2) / n..)

  cat("SCtrat ",SCtrat,"\n")
  cat("SCerror ",SCerror,"\n")
  cat("SCmuestreo ",SCmuestreo,"\n")
  cat("SCtotal ",SCtotal,"\n")

  #calculating square means
  CMtrat = SCtrat/GLtrat
  CMerror = SCerror/GLerror
  CMmuestreo = SCmuestreo/GLmuestreo


  cat("CMtrat ",CMtrat,"\n")
  cat("CMerror ",CMerror,"\n")
  cat("CMmuestreo ",CMmuestreo,"\n")


  #if the number to samples fro U.E are different
  #proceed to calculate the coefficients
  ni. = c(ni)

  #calculating A
  A = 0
  a1 = NULL
  for (j in 1: max(TMT)) {
    a1 = 0
    for (i in 1: length(SM)) {
      if(TMT[i] == j){
        a1 = a1 +nij[i]^2
      }
    }
    A = A+(a1/ni.[j])
  }

  #calculating B
  B = sum(nij^2)

  #calculating D
  D = sum(ni.^2)


  N = sum(nij)

  cat("A ", A,"\n")
  cat("B ", B,"\n")
  cat("D ", D,"\n")

  c1 = (1/(t-1))*(A-(B/N))
  c2 = (1/(t-1))*(N-(D/N))
  c3 = (1/(length(SM)-t))*(N-A)

  cat("c1 ", c1,"\n")
  cat("c2 ", c2,"\n")
  cat("c3 ", c3,"\n")
  ##############################################if sub sample is .....
  subsample_different = TRUE
    M = 0
  if(subsample_different){
    a1 = c1/c3
    a2 = 1-a1
    M = (a1*CMerror) + (a2*CMmuestreo)
    cat("M ", M,"\n")
  }else{
    c1 = 36####################################
    c3 = 36
  }

  #calculating degrees of freedom for M
  v = (M^2) /( (((a1*CMerror)^2)/GLerror) + (((a1*CMmuestreo)^2)/GLmuestreo))
  v = round(v+0.1)
  cat("V ", v,"\n")

  #F_aprox
  F_prox = CMtrat/M;
  #####################################################################
  message("DCA with Sub samping unbalanced \n By -Yeison Eduardo \n    -Sergio ")
  cat("--------------------------------ANOVA--------------------------------","\n")
  cat("   F.V.      d.f.          SC              CM              Fc    ","\n")
  cat(" Treatment   ",GLtrat,    "\t",round(SCtrat,4),"\t",round(CMtrat,4),"\t",F_prox,"\n")
  cat(" Error Exp.  ",GLerror,   "\t",round(SCerror,4),"\t",round(CMerror,4),"\t",CMerror/CMmuestreo,"\n")
  cat(" Sampling    ",GLmuestreo,"\t",round(SCmuestreo,4),"\t",round(CMmuestreo,4),"\t","--------","\n")
  cat(" Total       ",GLtotal,   "\t",round(SCtotal,4),"\n")
  cat("---------------------------------------------------------------------","\n")

}
