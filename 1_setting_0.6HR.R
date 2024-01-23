
## I.Survival parameter setting
####Author's Statement: The code for this study has been modified based on the work of Hatswell et al., 
####and specific technical details and annotations can be referred to at: Link to the research paper 
####(https://onlinelibrary.wiley.com/doi/10.1002/jrsm.1511),and Link to the GitHub repository
####(https://github.com/remiroazocar/population_adjustment_simstudy).


#####1.  Clear the history of all variables and set the working directory
rm(list = ls())
setwd("C:/Users/Administrator/Desktop/the last chars HR=1")


#####2. Set variable values and calculate all cross scenarios
replicates<- 1000  
no.vars<- 4
no.chars<- 4
vars<- 1:4
chars<- 1:4
N_A<- c(50,150,300)
N_B<- c(50,150,300)
b_trt<- c(log(0.6))#log0.3\log0.6\1
b_X<- c(log(1.1),log(1.5),log(2))
cen_prob<- c(0.3)
probX_A<- c(0.4,0.65,0.8)
probX_B<- c(0.9)
corX<- c(0,0.2,0.4)
param.combinations<- expand.grid(N_A=N_A, N_B=N_B, b_trt=b_trt, b_X=b_X, cen_prob=cen_prob, probX_A=probX_A, corX=corX)
pc<- param.combinations
scenarios<- nrow(param.combinations)


#####3.  Set the value of the survival function distribution parameter
weib_scale<- 8.5
weib_shape<- 1.3
getwd()


#####4.  Set the rate function of patient deletion
optim.function<- function(param, scale, shape, b_trt, cen_prob, N){
  cens_rate<- param
  U<- runif(N)
  Tlat<- -log(U)/(scale*exp(b_trt)^(1/shape))
  C<- rexp(n=N, rate=cens_rate)
  prop_cens<- sum(Tlat>C)/N
  fit<- sum((prop_cens-cen_prob)^2)
  return(fit)
}
cen_rate<-optim(par=1, fn=optim.function, scale=weib_scale, shape=weib_shape, 
                b_trt=b_trt, cen_prob=cen_prob, N=1000000, method="Brent",
                lower=0, upper=10)$par
save.image(file="survival_settings.RData")