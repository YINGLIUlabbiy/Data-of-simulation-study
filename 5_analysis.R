
## V.Data analysis
####Author's Statement: The code for this study has been modified based on the work of Hatswell et al., 
####and specific technical details and annotations can be referred to at: Link to the research paper 
####(https://onlinelibrary.wiley.com/doi/10.1002/jrsm.1511),and Link to the GitHub repository
####(https://github.com/remiroazocar/population_adjustment_simstudy).

setwd("C:/Users/Administrator/Desktop/the last chars HR=1")
getwd()
rm(list=ls())
source('3functions_chars.R')
load(file="survival_settings.RData")


####1.  Performance calculation function
#means
maic.means.list <- vector(mode="list", scenarios)
maic.variances.list <- vector(mode="list", scenarios)
naive.means.list <- vector(mode="list", scenarios)
naive.variances.list <- vector(mode="list", scenarios)
true.means.list <- vector(mode="list", scenarios)
#ate
maic.ate <- rep(NA, scenarios) 
naive.ate <- rep(NA, scenarios)
maic.ate.mcse <- rep(NA, scenarios) 
naive.ate.mcse <- rep(NA, scenarios)
#CI
maic.lci <- vector(mode="list", scenarios)
naive.lci <- vector(mode="list", scenarios)
maic.lci.mean <- rep(NA, scenarios)
naive.lci.mean <- rep(NA, scenarios)
maic.lci.mcse <- rep(NA, scenarios)
naive.lci.mcse <- rep(NA, scenarios)
maic.uci <- vector(mode="list", scenarios)
naive.uci <- vector(mode="list", scenarios)
maic.uci.mean <- rep(NA, scenarios)
naive.uci.mean <- rep(NA, scenarios)
maic.uci.mcse <- rep(NA, scenarios)
naive.uci.mcse <- rep(NA, scenarios)
maic.ciwidth <- vector(mode="list", scenarios)  
naive.ciwidth <- vector(mode="list", scenarios)
#bias
maic.bias <- rep(NA, scenarios) 
naive.bias <- rep(NA, scenarios)
maic.bias.mcse <- rep(NA, scenarios) 
naive.bias.mcse <- rep(NA, scenarios)
#mae
maic.abs.err <- vector(mode="list", scenarios)
naive.abs.err <- vector(mode="list", scenarios)
maic.mae <- rep(NA, scenarios)
naive.mae <- rep(NA, scenarios)
maic.mae.mcse <- rep(NA, scenarios)
naive.mae.mcse <- rep(NA, scenarios)
#mse
maic.mse <- rep(NA, scenarios)
naive.mse <- rep(NA, scenarios)
maic.mse.mcse <- rep(NA, scenarios)
naive.mse.mcse <- rep(NA, scenarios)
#vr
maic.vr <- rep(NA, scenarios)
naive.vr <- rep(NA, scenarios)
maic.vr.mcse <- rep(NA, scenarios)
naive.vr.mcse <- rep(NA, scenarios)
#ese
maic.empse <- rep(NA, scenarios)
naive.empse <- rep(NA, scenarios)
maic.empse.mcse <- rep(NA, scenarios)
naive.empse.mcse <- rep(NA, scenarios)
#coverage
maic.cov <- rep(NA, scenarios)
naive.cov <- rep(NA, scenarios)
maic.cov.mcse <- rep(NA, scenarios)
naive.cov.mcse <- rep(NA, scenarios)
#proportion of worse
maic.error.worse.than <- rep(NA,scenarios)
maic.error.worse.than.mcse <- rep(NA,scenarios)
#std.bias
maic.std.bias <- rep(NA,scenarios) 
naive.std.bias <- rep(NA,scenarios)
#ess
maic.aess.mean <- rep(NA,scenarios)
maic.aess.mcse <- rep(NA,scenarios)


####2.  Head of table
scenarios.df <- data.frame()
for (i in 1:scenarios) {
  file.id <- paste0("N_A", pc$N_A[i], "N_B", pc$N_B[i],"N_P", pc$N_P[i],"b_X", round(pc$b_X[i], digits=2), 
                    "probX_A", pc$probX_A[i], "corX", pc$corX[i])   
  
  
  # calculate the true marginal relative treatment effects
  load(paste0("Results/True/means_", file.id, ".RData"))
  true.means.list[[i]] <- means
  
  # MAIC calculate
  load(paste0("Results/MAIC/means_", file.id, ".RData"))
  load(paste0("Results/MAIC/variances_", file.id, ".RData"))
  #number of replicates for which MAIC weighted regression has separability issues 
  maic.separability.list[i] <- sum(abs(means)>max.effect)
  # discard these replicates 
  maic.no.sep <- abs(means)<max.effect
  means <- means[maic.no.sep]
  variances <- variances[maic.no.sep]
  # MAIC results
  maic.means.list[[i]] <- means
  maic.variances.list[[i]] <- variances
  maic.bias[i] <- bias(maic.means.list[[i]], true.means.list[[i]])  
  maic.bias.mcse[i] <- bias.mcse(maic.means.list[[i]])
  maic.mae[i] <- mae(maic.means.list[[i]], true.means.list[[i]])
  maic.abs.err[[i]] <- maic.means.list[[i]] - (true.means.list[[i]])
  maic.mae.mcse[i] <- mcse.estimate(maic.abs.err[[i]])
  maic.mse[i] <- mse(maic.means.list[[i]], true.means.list[[i]]) 
  maic.mse.mcse[i] <- mse.mcse(maic.means.list[[i]], true.means.list[[i]]) 
  maic.ate[i] <- mean(maic.means.list[[i]])
  maic.ate.mcse[i] <- mcse.estimate(maic.means.list[[i]])
  # construct confidence interval using normal distribution
  maic.lci[[i]] <- maic.means.list[[i]] + qnorm(0.025)*sqrt(maic.variances.list[[i]])
  maic.uci[[i]] <- maic.means.list[[i]] + qnorm(0.975)*sqrt(maic.variances.list[[i]])
  maic.ciwidth[[i]] <- maic.uci[[i]] - maic.lci[[i]]
  maic.lci.mean[i] <- mean(maic.lci[[i]])
  maic.lci.mcse[i] <- mcse.estimate(maic.lci[[i]])
  maic.uci.mean[i] <- mean(maic.uci[[i]])
  maic.uci.mcse[i] <- mcse.estimate(maic.uci[[i]])
  maic.cov[i] <- coverage(maic.lci[[i]], maic.uci[[i]], true.means.list[[i]])
  maic.cov.mcse[i] <- coverage.mcse(maic.cov[i], length(maic.lci[[i]]))
  maic.empse[i] <- empse(maic.means.list[[i]])
  maic.empse.mcse[i] <- empse.mcse(maic.empse[i], length(maic.means.list[[i]]))
  maic.vr[i] <- var.ratio(maic.means.list[[i]], sqrt(maic.variances.list[[i]]))
  maic.vr.mcse[i] <- var.ratio.mcse(avg.se=mean(sqrt(maic.variances.list[[i]])), 
                                    emp.se=maic.empse[i],
                                    var.avg.se=mcse.estimate(sqrt(maic.variances.list[[i]]))^2,
                                    var.emp.se=maic.empse.mcse[i]^2)
  maic.std.bias[i] <- (maic.bias[i]*100)/maic.empse[i]
  load(paste0("Results/MAIC/aess_", file.id, ".RData"))
  maic.aess.mean[i] <- mean(approx.ess.maic)
  maic.aess.mcse[i] <- mcse.estimate(approx.ess.maic) 
 
  
  ### NIC calculate
  load(paste0("Results/Naive/means_", file.id, ".RData"))
  load(paste0("Results/Naive/variances_", file.id, ".RData")) 
  naive.means.list[[i]] <- means
  naive.variances.list[[i]] <- variances
  naive.bias[i] <- bias(naive.means.list[[i]], true.means.list[[i]])  
  naive.bias.mcse[i] <- bias.mcse(naive.means.list[[i]])
  naive.mae[i] <- mae(naive.means.list[[i]], true.means.list[[i]])
  naive.abs.err[[i]] <- naive.means.list[[i]] - (true.means.list[[i]])
  naive.mae.mcse[i] <- mcse.estimate(naive.abs.err[[i]])
  naive.mse[i] <- mse(naive.means.list[[i]], true.means.list[[i]]) 
  naive.mse.mcse[i] <- mse.mcse(naive.means.list[[i]], true.means.list[[i]]) 
  naive.ate[i] <- mean(naive.means.list[[i]]) 
  naive.ate.mcse[i] <- mcse.estimate(naive.means.list[[i]])
  naive.lci[[i]] <- naive.means.list[[i]] + qnorm(0.025)*sqrt(naive.variances.list[[i]])
  naive.uci[[i]] <- naive.means.list[[i]] + qnorm(0.975)*sqrt(naive.variances.list[[i]])
  naive.ciwidth[[i]] <- naive.uci[[i]] - naive.lci[[i]]
  naive.lci.mean[i] <- mean(naive.lci[[i]])
  naive.lci.mcse[i] <- mcse.estimate(naive.lci[[i]])
  naive.uci.mean[i] <- mean(naive.uci[[i]])
  naive.uci.mcse[i] <- mcse.estimate(naive.uci[[i]])
  naive.cov[i] <- coverage(naive.lci[[i]], naive.uci[[i]], true.means.list[[i]])
  naive.cov.mcse[i] <- coverage.mcse(naive.cov[i], length(naive.lci[[i]]))
  naive.empse[i] <- empse(naive.means.list[[i]])
  naive.empse.mcse[i] <- empse.mcse(naive.empse[i], replicates)
  naive.vr[i] <- var.ratio(naive.means.list[[i]], sqrt(naive.variances.list[[i]]))
  naive.vr.mcse[i] <- var.ratio.mcse(avg.se=mean(sqrt(naive.variances.list[[i]])), 
                                     emp.se=naive.empse[i],
                                     var.avg.se=mcse.estimate(sqrt(naive.variances.list[[i]]))^2,
                                     var.emp.se=naive.empse.mcse[i]^2) 
  naive.std.bias[i] <- (naive.bias[i]*100)/naive.empse[i]
  maic.error.worse.than[i] <- sum(abs(naive.means.list[[i]][maic.no.sep]-true.means.list[[i]])<abs(maic.means.list[[i]]-true.means.list[[i]]))/sum(maic.no.sep)
  maic.error.worse.than.mcse[i] <- coverage.mcse(maic.error.worse.than[i], sum(maic.no.sep))
  
  ## Combining data
  tmp.scenarios <- cbind(i, pc$N_A[i],pc$N_B[i],pc$N_P[i],pc$b_X[i], pc$probX_A[i], pc$corX[i])
  maic.tmp.metrics <- cbind(maic.ate[i], maic.ate.mcse[i], maic.lci.mean[i],
                            maic.lci.mcse[i], maic.uci.mean[i], maic.uci.mcse[i],
                            maic.bias[i], maic.bias.mcse[i], maic.mse[i], maic.mse.mcse[i],
                            maic.mae[i], maic.mae.mcse[i], maic.cov[i], maic.cov.mcse[i],
                            maic.empse[i], maic.empse.mcse[i], maic.vr[i], maic.vr.mcse[i],
                            maic.error.worse.than[i], maic.error.worse.than.mcse[i], 
                            maic.std.bias[i], maic.aess.mean[i], maic.aess.mcse[i])
  naive.tmp.metrics <- cbind(naive.ate[i], naive.ate.mcse[i], naive.lci.mean[i],
                             naive.lci.mcse[i], naive.uci.mean[i], naive.uci.mcse[i],
                             naive.bias[i], naive.bias.mcse[i], naive.mse[i], naive.mse.mcse[i],
                             naive.mae[i], naive.mae.mcse[i], naive.cov[i], naive.cov.mcse[i],
                             naive.empse[i], naive.empse.mcse[i], naive.vr[i], naive.vr.mcse[i], 
                             naive.std.bias[i]) 
  tmp.scenarios<- cbind(tmp.scenarios, maic.tmp.metrics, naive.tmp.metrics)
  scenarios.df <- rbind(scenarios.df, tmp.scenarios)
}


####3.  data table
colnames(scenarios.df) <- c("Scenario", "N_A","N_B","b_X","probX_A", "corX",
                            "maic.ate", "maic.ate.mcse", "maic.lci.mean",
                            "maic.lci.mcse", "maic.uci.mean", "maic.uci.mcse",
                            "maic.bias", "maic.bias.mcse", "maic.mse", "maic.mse.mcse",
                            "maic.mae", "maic.mae.mcse", "maic.cov", "maic.cov.mcse",
                            "maic.empse", "maic.empse.mcse", "maic.vr", "maic.vr.mcse",
                            "maic.error.worse.than", "maic.error.worse.than.mcse", 
                            "maic.std.bias", "maic.aess.mean", "maic.aess.mcse",
                            
                            "naive.ate", "naive.ate.mcse", "naive.lci.mean",
                            "naive.lci.mcse", "naive.uci.mean", "naive.uci.mcse",
                            "naive.bias", "naive.bias.mcse", "naive.mse", "naive.mse.mcse",
                            "naive.mae", "naive.mae.mcse", "naive.cov", "naive.cov.mcse",
                            "naive.empse", "naive.empse.mcse", "naive.vr", 
                            "naive.vr.mcse", "naive.std.bias")
write.csv(scenarios.df, "Analysis/scenarios.csv", row.names = FALSE)
