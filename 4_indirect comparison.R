
## IV.Indirect comparison
####Author's Statement: The code for this study has been modified based on the work of Hatswell et al., 
####and specific technical details and annotations can be referred to at: Link to the research paper 
####(https://onlinelibrary.wiley.com/doi/10.1002/jrsm.1511),and Link to the GitHub repository
####(https://github.com/remiroazocar/population_adjustment_simstudy).

rm(list=ls())
load(file="survival_settings.RData")
source('3functions_chars.R')


####1.  settings
set.seed(10010)
load(file="survival_settings.RData")
scenarios<-nrow(pc)
if(!require(dplyr)) {install.packages("dplyr"); library(dplyr)}
if(!require(tidyr)) {install.packages("tidyr"); library(tidyr)}
if(!require(parallel)) {install.packages("parallel"); library(parallel)}
if(!require(doSNOW)) {install.packages("doSNOW"); library(doSNOW)}
if(!require(survival)) {install.packages("survival"); library(survival)}

#Import IPD data for all scenarios
IPD.A.all<- vector(mode="list", scenarios)
IPD.B.all<- vector(mode="list", scenarios)
IPD.PA.all<- vector(mode="list", scenarios)
IPD.PB.all<- vector(mode="list", scenarios)
AgD.B.all<- vector(mode="list", scenarios)

#Name all scenarios
for (i in 1:scenarios) {
  file.id <- paste0("N_A", pc$N_A[i], "N_B", pc$N_B[i],"N_P", pc$N_P[i],"b_X", round(pc$b_X[i], digits=2), 
                    "probX_A", pc$probX_A[i], "corX", pc$corX[i]) 
  load(paste0("Data/IPD_A_", file.id, ".RData"))
  load(paste0("Data/IPD_B_", file.id, ".RData"))
  load(paste0("Data/IPD_PA_", file.id, ".RData"))
  load(paste0("Data/IPD_PB_", file.id, ".RData"))
  load(paste0("Data/AgD_B_", file.id, ".RData"))
  IPD.A.all[[i]] <- IPD.A
  IPD.B.all[[i]] <- IPD.B
  IPD.PA.all[[i]] <- IPD.PA
  IPD.PB.all[[i]] <- IPD.PB
  AgD.B.all[[i]] <- AgD.B
}
gc()


####2.  wrapper function
#naive indirect comparison
naive.wrapper<-function(data.A, data.B){
  data.AB<- rbind(data.A,data.B)
  data.AB$trt<- factor(data.AB$trt, levels=c("B","A"))
  d.AB.naive<-summary(coxph(Surv(time, status)~trt, data=data.AB))$coef[1]
  var.d.AB.naive<-vcov(coxph(Surv(time, status)~trt, data=data.AB))[1]
  list(d.AB.naive, var.d.AB.naive)
}

#compute the truth function
true.wrapper<-function(data.A, data.B){
  data.AB<- rbind(data.A,data.B)
  data.AB$trt<- factor(data.AB$trt, levels=c("PB","PA"))
  d.AB.true<-summary(coxph(Surv(time, status)~trt, data=data.AB))$coef[1]
  var.d.AB.true<-vcov(coxph(Surv(time, status)~trt, data=data.AB))[1]
  list(d.AB.true, var.d.AB.true)
}

#matching-adjusted indirect comparison
maic.wrapper <- function(data.A, data.Bsum, data.Bipd, vars) { # ems indexes the position (columns) of the effect modifiers
  A.vars <- data.A[,1+vars] # column 1 of IPD is treatment indicator 
  weights <- maic(A.X=A.vars, B.summary=data.Bsum[vars]) # maic weights through method of moments
  maic.aess <- approx.ess(weights) # approximate effective sample size
  # fit weighted Cox proportional hazards model using robust=TRUE
  data.Aall<-cbind(data.A,weights)
  data.Ball<-mutate(data.Bipd,weights=c(1))
  data.AB<-rbind(data.Aall,data.Ball)
  data.AB$trt<-factor(data.AB$trt,levels=c("B","A"))
  outcome.fit.maic <- coxph(Surv(time, status)~trt, robust=TRUE, weights=data.AB$weights, data=data.AB)
  d.AB.maic <- summary(outcome.fit.maic)$coef[1]
  var.d.AB.maic <- vcov(outcome.fit.maic)[[1]] 
  list(d.AB.maic, var.d.AB.maic, maic.aess,data.AB$weights)
}


####3. parallel computing
num.cores <- detectCores()
cluster <- makeCluster(num.cores, type="SOCK", outfile="")
registerDoSNOW(cluster)
pb <- txtProgressBar(max=replicates, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)
comb <- function(x, ...) {
  lapply(seq_along(x),
         function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}


####4. calculation
for(i in 1:scenarios) {
  IPD.A <- IPD.A.all[[i]]
  IPD.B <- IPD.B.all[[i]]
  IPD.PA <- IPD.PA.all[[i]]
  IPD.PB <- IPD.PB.all[[i]]
  AgD.B <- AgD.B.all[[i]]
  vars_i <- vars
  chars_i <- chars
  file.id <- paste0("N_A", pc$N_A[i],"N_B", pc$N_B[i],"N_P", pc$N_P[i], "b_X", round(pc$b_X[i], digits=2), 
                    "probX_A", pc$probX_A[i], "corX", pc$corX[i]) 
  # NIC
  naive.results <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                           .init=list(list(), list()), .options.snow=opts,
                           .packages=c("dplyr","tidyr", "survival")) %dopar% {
                             results <- naive.wrapper(IPD.A[[j]], IPD.B[[j]])
                             return(results)
                           }
  close(pb)
  means <- unlist(naive.results[[1]])
  variances <- unlist(naive.results[[2]])
  save(means, file=paste0("Results/Naive/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/Naive/variances_", file.id, ".RData"))  
  
  # marginal relative treatment effects
  true.results <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                           .init=list(list(), list()), .options.snow=opts,
                           .packages=c("dplyr","tidyr", "survival")) %dopar% {
                             results <- true.wrapper(IPD.PA[[j]], IPD.PB[[j]])
                             return(results)
                           }
  close(pb)
  means <- unlist(true.results[[1]])
  save(means, file=paste0("Results/True/means_", file.id, ".RData"))

  # MAIC
  maic.results <- foreach(j=1:replicates, .combine='comb', .multicombine=TRUE,
                          .init=list(list(), list(), list(),list()), .options.snow=opts,
                          .packages=c("dplyr","tidyr", "survival")) %dopar% {
                            results <- maic.wrapper(IPD.A[[j]], AgD.B[[j]],IPD.B[[j]],vars=vars_i)
                            return(results)
                          }
  close(pb)
  means <- unlist(maic.results[[1]])
  variances <- unlist(maic.results[[2]])
  approx.ess.maic <- unlist(maic.results[[3]])
  weights.all<- unlist(maic.results[[4]])
  save(means, file=paste0("Results/MAIC/means_", file.id, ".RData"))
  save(variances, file=paste0("Results/MAIC/variances_", file.id, ".RData"))  
  save(approx.ess.maic, file=paste0("Results/MAIC/aess_", file.id, ".RData"))
  save(weights.all, file=paste0("Results/MAIC/weights_", file.id, ".RData"))
}  

stopCluster(cluster)
