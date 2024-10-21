
## II.Data generation
####Author's Statement: The code for this study has been modified based on the work of Antonio Remiro-Azócar et al., 
####and specific technical details and annotations can be referred to at: Link to the research paper 
####(https://onlinelibrary.wiley.com/doi/10.1002/jrsm.1511),and Link to the GitHub repository
####(https://github.com/remiroazocar/population_adjustment_simstudy).


####1.  Load the software package and data
rm(list = ls())
if(!require(simstudy)){install.packages("simstudy");library(simstudy)} 
if(!require(dplyr)){install.packages("dplyr");library(dplyr)}
if(!require(survival)){install.packages("survival");library(survival)}
load(file="survival_settings.RData")


####2.  Generate data matrix
set.seed(204204204)##(This seed number is derived from the office room number I worked in during my master's degree)
gen.data<- function(no.vars, no.chars, N_A, N_B, b_trt, 
                    b_X, probX_A, probX_B, sdX, weib_scale, weib_shape, 
                    cens_rate, corX){
  #Set the correlation matrix
  R<- matrix(corX, nrow=no.vars, ncol=no.vars) 
  diag(R)<- rep(1, no.vars)
  states.data<- .Random.seed
  #Record the sample size
  N_A<- floor(N_A)
  N_B<- floor(N_B)
  N_P<-10000 #the large sample used to calculate the marginal relative treatment effect
  #The covariate matrix
  X_A<- as.matrix(genCorGen(n=N_A, nvar=no.vars, 
                            params1= as.numeric(rep(probX_A, no.vars)),
                            dist= "binary", corstr= "cs", corMatrix= R, wide=TRUE)[,-1]) #去掉第一列
  
  X_B<- as.matrix(genCorGen(n=N_B, nvar=no.vars, 
                            params1= as.numeric(rep(probX_B, no.vars)),
                            dist= "binary", corstr= "cs", corMatrix= R, wide=TRUE)[,-1])
  
  X_P<- as.matrix(genCorGen(n=N_P/2, nvar=no.vars, 
                            params1= as.numeric(rep(probX_B, no.vars)),
                            dist= "binary", corstr= "cs", corMatrix= R, wide=TRUE)[,-1])
  #set HR values for patients
  betaX_AB_A<- rep(b_trt, N_A)
  betaX_AB_B<- rep(0, N_B)
  betaX_AB_PA<- rep(b_trt, N_P/2)
  betaX_AB_PB<- rep(0, N_P/2)
  #Naming column names
  col.names<- NULL
  for (k in 1:no.vars) {
    col.names<-c(col.names, paste0("X", k))
    betaX_AB_A<- betaX_AB_A+ b_X * X_A[,k] #Patients in the drug A trial received drug A
    betaX_AB_B<- betaX_AB_B+ b_X * X_B[,k] #Patients in the drug B trial received drug B
    betaX_AB_PA<- betaX_AB_PA+ b_X * X_P[,k] #patients in the drug B trial were given drug A (simulated counterfactual results in a large population P)
    betaX_AB_PB<- betaX_AB_PB+ b_X * X_P[,k] #Patients in drug B trials were given drug B (simulated normal results in large population P)
  }
  
  #IPD framework for trial A was generated
  X_A= as.data.frame(X_A)
  colnames(X_A)<- col.names
  U= runif(n=N_A)
  Tlat =(-log(U)/(weib_scale*exp(betaX_AB_A)))^(1/weib_shape)
  C= rexp(n=N_A, rate=cen_rate)
  time= pmin(Tlat, C)#Survival times following a Weibull-Cox distribution were generated
  status= as.numeric(Tlat<=C)
  trt<- c(rep("A", N_A))
  IPD.A<- as.data.frame(cbind(trt, X_A, time, status))

  #IPD framework for trial B was generated
  X_B= as.data.frame(X_B)
  colnames(X_B)<- col.names
  U= runif(n=N_B)
  Tlat =(-log(U)/(weib_scale*exp(betaX_AB_B)))^(1/weib_shape)
  C= rexp(n=N_B, rate=cen_rate)
  time= pmin(Tlat, C)
  status= as.numeric(Tlat<=C)
  trt<- c(rep("B", N_B))
  IPD.B<- as.data.frame(cbind(trt, X_B, time, status))
  
  #IPD framework for a large sample population was generated
  X_P= as.data.frame(X_P)
  colnames(X_P)<- col.names
  U= runif(n=N_P/2)
  Tlat =(-log(U)/(weib_scale*exp(betaX_AB_PA)))^(1/weib_shape)
  C= rexp(n=N_P/2, rate=cen_rate)
  time= pmin(Tlat, C)
  status= as.numeric(Tlat<=C)
  trt<- c(rep("PA", N_P/2))
  IPD.PA<- as.data.frame(cbind(trt, X_P, time, status))
  U= runif(n=N_P/2)
  Tlat =(-log(U)/(weib_scale*exp(betaX_AB_PB)))^(1/weib_shape)
  C= rexp(n=N_P/2, rate=cen_rate)
  time= pmin(Tlat, C)
  status= as.numeric(Tlat<=C)
  trt<- c(rep("PB", N_P/2))
  IPD.PB<- as.data.frame(cbind(trt, X_P, time, status))
  #IPD.P<- rbind(trt, X_P, time, status)
  
  #The large sample IPD was summarized into AgD
  AgD.B<- as.data.frame(cbind(
    summarise(IPD.B, mean.X1=mean(X1), mean.X2=mean(X2), mean.X3=mean(X3), mean.X4=mean(X4))))
  list(IPD.A, IPD.B, AgD.B, IPD.PA, IPD.PB, states.data)
}


####3.  IPDs were generated for trials A and B
for (i in 1:scenarios) {
  print(i)
  IPD.A<- vector(mode="list", replicates)
  IPD.B<- vector(mode="list", replicates)
  AgD.B<- vector(mode="list", replicates)
  IPD.PA<- vector(mode="list", replicates)
  IPD.PB<- vector(mode="list", replicates)
  states.data<- vector(mode="list", replicates+1)
    for (j in 1:replicates) {
    gen.datasets<- gen.data(no.vars=no.vars, no.chars=no.chars, N_A=pc$N_A[i], N_B=pc$N_B[i],  b_trt=b_trt, 
                            b_X=pc$b_X[i], probX_A=pc$probX_A[i], probX_B=probX_B,
                            weib_scale=weib_scale, weib_shape=weib_shape, cens_rate=cens_rate, corX=pc$corX[i])
    IPD.A[[j]]<-gen.datasets[[1]]
    IPD.B[[j]]<-gen.datasets[[2]]
    AgD.B[[j]]<-gen.datasets[[3]]
    IPD.PA[[j]]<-gen.datasets[[4]]
    IPD.PB[[j]]<-gen.datasets[[5]]
    states.data[[j]]<- gen.datasets[[6]]
  }
  states.data[[replicates+1]]<- .Random.seed
  file.id<- paste0("N_A", pc$N_A[i], "N_B", pc$N_B[i],"N_P", pc$N_P[i],"b_X", round(pc$b_X[i],digits=2), 
                   "probX_A", pc$probX_A[i], "corX", pc$corX[i])
  save(IPD.A, file=paste0("Data/IPD_A_", file.id, ".RData"))
  save(IPD.B, file=paste0("Data/IPD_B_", file.id, ".RData"))
  save(IPD.PA, file=paste0("Data/IPD_PA_", file.id, ".RData"))
  save(IPD.PB, file=paste0("Data/IPD_PB_", file.id, ".RData"))
  save(AgD.B, file=paste0("Data/AgD_B_", file.id, ".RData"))
  save(states.data, file=paste0("Data/States_", file.id, ".RData"))
}
