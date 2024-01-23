
## VI.Plotting
####Author's Statement: The code for this study has been modified based on the work of Hatswell et al., 
####and specific technical details and annotations can be referred to at: Link to the research paper 
####(https://onlinelibrary.wiley.com/doi/10.1002/jrsm.1511),and Link to the GitHub repository
####(https://github.com/remiroazocar/population_adjustment_simstudy).

rm(list=ls())
setwd('C:/Users/Administrator/Desktop/alldata')
source('3functions_chars.R')
scenarios.df <- read.csv(file="C:/Users/Administrator/Desktop/alldata/alldata.csv", header=TRUE, sep=",")#need to put in all the data from the combined 729 scenarios


####1. Nestedloop
# bias
nldata <- nestedloop(scenarios.df,
                     varnames=c("probX_A","b_X","HR","N_A","N_B","corX"),
                     varlabels=c("Overlap of covariates",
                                 "Strength of covariates",
                                 "Strength of relative treatment effect",
                                 "Sample size of individual patient data",
                                 "Sample size of aggregate data",
                                 "Covariate correlation"),
                     sign=c(1, 1, 1, 1, 1, 1))
pd.bias <- nldata
pd.bias$probX_A <- factor(pd.bias$probX_A, levels=c(0.4,0.65,0.8), 
                          labels=c("6%~20%","36%~57%","72%~85%"))
pd.bias$b_X <- factor(pd.bias$b_X, levels=c(0.09531018,0.405465108,0.693147181), 
                      labels=c("ln(1.1)","ln(1.5)","ln(2)"))

# MSE
nldata <- nestedloop(scenarios.df,
                     varnames=c("probX_A","b_X","HR","N_A","corX","N_B"),
                     varlabels=c("Overlap of covariates",
                                 "Strength of covariates",
                                 "Strength of relative treatment effect",
                                 "Sample size of individual patient data",
                                 "Covariate correlation",
                                 "Sample size of aggregate data"),
                     sign=c(1, 1, 1, 1, 1, 1))
pd.mse <- nldata
pd.mse$probX_A <- factor(pd.mse$probX_A,levels=c(0.4,0.65,0.8),
                         labels=c("6%~20%","36%~57%","72%~85%"))
pd.mse$b_X <- factor(pd.mse$b_X, levels=c(0.09531018,0.405465108,0.693147181), 
                     labels=c("ln(1.1)","ln(1.5)","ln(2)"))


# ESE
nldata <- nestedloop(scenarios.df,
                     varnames=c("probX_A","N_A","corX","b_X","HR","N_B"),
                     varlabels=c("Overlap of covariates",
                                 "Sample size of individual patient data",
                                 "Covariate correlation",
                                 "Strength of covariates",
                                 "Strength of relative treatment effect",
                                 "Sample size of aggregate data"),
                     sign=c(1, 1, 1, 1, 1, 1))
pd.ese <- nldata
pd.ese$probX_A <- factor(pd.ese$probX_A,levels=c(0.4,0.65,0.8),
                         labels=c("6%~20%","36%~57%","72%~85%"))
pd.ese$b_X <- factor(pd.ese$b_X, levels=c(0.09531018,0.405465108,0.693147181), 
                     labels=c("ln(1.1)","ln(1.5)","ln(2)"))

# coverage
nldata <- nestedloop(scenarios.df,
                     varnames=c("b_X","probX_A","HR","N_A","corX","N_B"),
                     varlabels=c("Strength of covariates",
                                 "Overlap of covariates",
                                 "Strength of relative treatment effect",
                                 "Sample size of individual patient data",
                                 "Covariate correlation",
                                 "Sample size of aggregate data" ),
                     sign=c(1, 1, 1, 1, 1, 1))
pd.cover <- nldata
pd.cover$probX_A <- factor(pd.cover$probX_A,levels=c(0.4,0.65,0.8),
                           labels=c("6%~20%","36%~57%","72%~85%"))
pd.cover$b_X <- factor(pd.cover$b_X, levels=c(0.09531018,0.405465108,0.693147181), 
                       labels=c("ln(1.1)","ln(1.5)","ln(2)"))


####2. TIFF
# bias
tiff("Analysis/bias.tiff", res=800, width=18, height=15, units='in')
par(pty="m")
plot(pd.bias$maic.bias, type="n", ylim=c(-1.3,1), bty="n", xlab="Scenario", ylab="Bias",
     las=1, xaxt="n", cex.axis=1.5, cex.lab=1.5)
abline(h=0, col="grey") # no bias
lines(pd.bias, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE")) # add vertical lines
lines(pd.bias, which="r",ymin.refline=0.35, ymax.refline=1, cex.ref=1.2) # add reference lines
lines(pd.bias$maic.bias, col="#FF8C00", lty=5,type="s", lwd=2) # performance measures 
lines(pd.bias$naive.bias, col="#4B0082", lty=2, type="s", lwd=1.8)
legend("bottom",lwd=c(2,1.8),col=c("#FF8C00", "#4B0082"),
       lty=c(5,2),cex=1.2,bty="n", c("MAIC", "NIC")) # legend
dev.off()

# MSE
tiff("Analysis/mse.tiff", res=800, width=18, height=15, units='in')
par(pty="m")
plot(pd.mse$maic.mse,
     type="n",
     ylim=c(0,2), bty="n",
     xlab="Scenario",
     ylab="Mean square error (MSE)",
     las=1, xaxt="n",
     cex.axis=1.5, cex.lab=1.5) 
lines(pd.mse, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.mse, which="r", ymin.refline=1.6, ymax.refline=2, cex.ref=1.2)
lines(pd.mse$maic.mse, col="#FF8C00", lty=5,type="s", lwd=2)
lines(pd.mse$naive.mse, col="#4B0082", lty=2, type="s", lwd=1.8)
legend("right",lwd=c(2,1.8),col=c("#FF8C00","#4B0082"), lty=c(5,2),
       cex=1.2,bty="n",c("MAIC", "NIC"))
dev.off()

# ESE
tiff("Analysis/EmpSE.tiff", res=800, width=18, height=15, units='in')
par(pty="m")
plot(pd.ese$maic.empse, type="n",ylim=c(0, 1.4), bty="n", xlab="Scenario",
     ylab="Empirical standard error (ESE)",las=1, xaxt="n",
     cex.axis=1.5, cex.lab=1.5)
lines(pd.ese, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.ese, which="r", ymin.refline=1.0, ymax.refline=1.4, cex.ref=1.2)
lines(pd.ese$maic.empse, col="#FF8C00", lty=5,type="s",lwd=2)    
lines(pd.ese$naive.empse, col="#4B0082", lty=2, type="s", lwd=1.8) 
legend("right", lwd=c(2,1.8), col=c("#FF8C00","#4B0082"), lty=c(5,2), 
       cex=1.2, bty="n", c("MAIC", "NIC"))
dev.off()

# coverage
tiff("Analysis/coverage.tiff", res=800, width=18, height=15, units='in')
par(pty="m")
plot(pd.cover$maic.cov*100, type="n", ylim=c(0, 140), bty="n", xlab="Scenario",
     ylab="Coverage of 95% confidence intervals (%)", las=1, xaxt="n",
     cex.axis=1.5, cex.lab=1.5)
abline(h=95, col="grey") # nominal 
lines(pd.cover, col=c("#BBBBBB", "#CCCCCC", "#DDDDDD", "#EEEEEE"))
lines(pd.cover, which="r", ymin.refline=100, ymax.refline=140, cex.ref=1.2)
lines(pd.cover$maic.cov*100, col="#FF8C00", lty=5,type="s", lwd=2)    
lines(pd.cover$naive.cov*100, col="#4B0082", lty=2, type="s", lwd=1.8) 
legend("bottomright", lwd=c(2,1.8), col=c("#FF8C00", "#4B0082"), 
       lty=c(5,2), cex=1.2, bty="n",c("MAIC", "NIC"))
dev.off()
