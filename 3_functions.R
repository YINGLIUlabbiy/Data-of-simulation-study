
## III.Functions
####Author's Statement: The code for this study has been modified based on the work of Hatswell et al., 
####and specific technical details and annotations can be referred to at: Link to the research paper 
####(https://onlinelibrary.wiley.com/doi/10.1002/jrsm.1511),and Link to the GitHub repository
####(https://github.com/remiroazocar/population_adjustment_simstudy).

rm(list = ls())


#####1.  MAIC correlation function
maic<- function(A.X, B.summary){
  M<-length(B.summary)
  for (i in 1:M) {
    A.X[,i]<- A.X[,i]-B.summary[,i]
  }
  A.X<-as.matrix(A.X)
     objective.function <- function(beta.pars,X){
    return(sum(exp(X %*% beta.pars)))
  }
   beta.start<- rep(1,M)
  out.optim <- optim(fn=objective.function, X=A.X, par=beta.start, method = "CG")
  beta.pars <- out.optim$par
  temp.w<-rep(0,nrow(A.X))
  for (i in 1:M){
    temp.w<-temp.w+beta.pars[i]*A.X[,i]
  }
  w<-exp(temp.w)
  return(w)
}
#ESS
approx.ess <- function(w) {
  ess <- sum(w)^2/sum(w^2)
  return(ess)
}

#####2.  Performance evaluation correlation functions

###1.bias
bias <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  est <- sum(theta.hat)/nsim - mean(theta)
  return(est)
}
# Monte Carlo SE of bias estimate
bias.mcse <- function(theta.hat) {
  nsim <- length(theta.hat)
  tmp <- sum((theta.hat - mean(theta.hat))^2)
  mcse <- sqrt(1/(nsim*(nsim-1))*tmp)
  return(mcse)
}

###2.coverage
coverage <- function(theta.hat.low, theta.hat.upp, theta) {
  nsim <- length(theta.hat.low)
  est <- sum(ifelse(theta>=theta.hat.low & theta<=theta.hat.upp,1,0))/nsim
  return(est)
}
# Monte Carlo SE of coverage estimate
coverage.mcse <- function(coverage, nsim) {
  mcse <- sqrt((coverage*(1-coverage))/nsim)
  return(mcse)
}

###3.MSE
mse <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  est <- sum((theta.hat-theta)^2)/nsim
  return(est)
}
# Monte Carlo SE of MSE estimate
mse.mcse <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  tmp <- (theta.hat-theta)^2
  mse.est <- sum(tmp)/nsim
  mcse <- sqrt(sum((tmp - mse.est)^2)/(nsim*(nsim-1)))
  return(mcse)
}

###4.MAE
mae <- function(theta.hat, theta) {
  nsim <- length(theta.hat)
  est <- sum(abs(theta.hat-theta))/nsim
  return(est)
}
# Monte Carlo SE of any continuous estimate
mcse.estimate <- function(perf.measure) {
  nsim <- length(perf.measure)
  perf.measure.mean <- sum(perf.measure)/nsim
  mcse <- sqrt(sum((perf.measure-perf.measure.mean)^2)/(nsim*(nsim-1)))
  return(mcse)
}

###5.ESE
empse <- function(theta.hat) {
  nsim <- length(theta.hat)
  tmp <- sum((theta.hat - mean(theta.hat))^2)
  est <- sqrt(tmp/(nsim-1))
  return(est)
}
# EmpSE MCSE
empse.mcse <- function(empse, nsim) {
  mcse <- empse/(sqrt(2*(nsim-1)))
  return(mcse)
} 

###6.Variability ratio
var.ratio <- function(theta.hat, std.err) {
  nsim <- length(theta.hat)
  num <- sum(std.err)/nsim
  denom <- sqrt(sum((theta.hat-mean(theta.hat))^2)/(nsim-1))
  est <- num/denom
  return(est)    
}
# Variability ratio MCSE
var.ratio.mcse <- function(avg.se, emp.se, var.avg.se, var.emp.se) {
  # approximation of ratio variance based on independence of avg. se and emp.se
  # see Wolter, K., 2007. Introduction to variance estimation. 
  mcse <- sqrt((1/emp.se^2)*var.avg.se + (((avg.se^2)/(emp.se^4))*var.emp.se))
  return(mcse)
}


#####3.  Graphing
# nestedloop graph
nestedloop <- function(x,
                       varnames, sign=rep(1, length(varnames)),
                       varlabels=NULL){
  ##
  if (!inherits(x, "data.frame"))
    stop("Argument 'x' must be a data.frame.")
  ##
  mo <- matrix(sign,
               nrow=dim(x)[[1]], ncol=length(varnames),
               byrow=TRUE)
  xo <- x[,varnames]
  ##
  ## Re-ordering:
  res <- x[do.call(order, mo*xo),]
  ##
  attr(res, "varnames")  <- varnames
  attr(res, "varlabels") <- varlabels
  attr(res, "sign") <- sign
  ##
  class(res) <- c("nestedloop", class(res))
  ##
  res
}

# Reference lines
lines.nestedloop <- function(x,
                             varnames=attributes(x)$varnames,
                             varlabels=attributes(x)$varlabels,
                             which="v",
                             col=if (which=="r") "#999999" else "black",
                             ymin.refline, ymax.refline,
                             cex.ref=0.9,
                             log=TRUE,
                             ...){
  ##
  nvar <- length(varnames)
  ##
  if (length(col)==1)
    col <- rep(col, nvar)
  ##
  if (which=="v"){
    ##
    nlen <- rep(NA, nvar)
    ##
    for (i in 1:nvar)
      nlen[i] <- length(unique(x[,varnames[i]]))
    ##
    cnlen <- cumprod(nlen)
    ##
    for (i in (nvar-1):1)
      abline(v=cnlen[nvar]*(0:cnlen[i])/cnlen[i]+1,
             col=col[i])
  }
  else if (which=="r"){
    ##
    if (is.null(varlabels))
      varlabels <- varnames
    ##
    labels.varnames <- rep("", nvar)
    ##
    for (i in 1:length(varnames)){
      if (is.factor(x[,varnames[i]]))
        varvals <- unique(x[,varnames[i]])
      else{
        varvals <- format(unique(x[,varnames[i]]))
        varvals <- sub("^[[:space:]]*(.*?)[[:space:]]*$",
                       "\\1",
                       varvals,
                       perl=TRUE)
      }
      ##
      labels.varnames[i] <- paste(varlabels[i],
                                  " (",
                                  paste(varvals,
                                        collapse=", "),
                                  ")", sep="")
    }
    ##
    if (log){
      ymax <- log(ymax.refline)
      ymin <- log(ymin.refline)
    }
    else{
      ymax <- ymax.refline
      ymin <- ymin.refline
    }
    ##
    distance <- (ymax-ymin)/nvar
    ##
    ypos <- ymax-0.2*distance-(1/nvar)*(0:(nvar-1))*(ymax-ymin)
    ypos.ref.max <- ypos-0.20*distance
    ypos.ref.min <- ypos-0.75*distance
    ##
    if (log){
      ypos <- exp(ypos)
      ypos.ref.max <- exp(ypos.ref.max)
      ypos.ref.min <- exp(ypos.ref.min)
    }
    ##
    for (i in 1:nvar){
      text(1, ypos[i], labels.varnames[i], adj=0, cex=cex.ref)
      ##
      xvar <- x[,varnames[i]]
      if (is.factor(xvar))
        xvar <- as.numeric(xvar)
      xvar <- ypos.ref.min[i] +
        (xvar-min(xvar))/(max(xvar)-min(xvar))*
        (ypos.ref.max[i]-ypos.ref.min[i])
      lines(xvar, col=col[i], type="s", lwd=1)
    }
  }
  ##
  invisible(NULL)
}
