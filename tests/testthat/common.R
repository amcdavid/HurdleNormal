Indep <- HurdleStructure(G=-13*diag(3), H=diag(3)*5, K=diag(3), gibbs=FALSE)
Gdep <- Indep
Gdep$true$G[1,2] <- Gdep$true$G[2,1] <- 1.5
Kdep1 <- Indep
Kdep1$true$K[1,2] <- Kdep1$true$K[2,1] <- -.7
Kdep2 <- Kdep1
Kdep2$true$H[1,2] <- Kdep2$true$H[2,1] <- -2.25
Hupdep <- Indep
Hupdep$true$H[1,2] <- 1
Hupdep$true$G[2,2] <- Hupdep$true$G[1,1] <- -16
Hlodep <- Indep
Hlodep$true$H[2,1] <- 1
Hlodep$true$G[2,2] <- Hlodep$true$G[1,1] <- -17
thin <- .01
diag(Indep$true$G) <- -12.5

Indep <- getGibbs(Indep, thin=thin)
Gdep <- getGibbs(Gdep, thin=thin)
Kdep2 <- getGibbs(Kdep2, thin=thin)
Hupdep <- getGibbs(Hupdep, thin=thin)
Hlodep <- getGibbs(Hlodep, thin=thin)


library(plyr)
library(nloptr)
getLimit <- function(hs, engine='R'){
    burnin=1000
    series <- c(100, 400, 2000)
    fitseries <- matrix(NA, nrow=length(series), ncol=3)
    colnames(fitseries) <- c("G", "H", "K")
    for(i in seq_along(series)){
        sub <- series[i]
        n <- sub/thin+burnin
        hs <- suppressMessages(getGibbs(hs, Nt=n, burnin=burnin, thin=thin))
        #hs$gibbs <- addPseudoCounts(hs$gibbs)
        fit <- getConditionalMLE(hs, testGrad=TRUE, engine=engine) 
        gj <- getJoint(fit)
        err <- sapply(c('G', 'H', 'K'), function(i) sum((hs$true[[i]]-gj[[i]])^2))
        fitseries[i,] <- err
    }
    fitseries
}

library(numDeriv)
checkGrad <- function(hs, j=1, theta, lambda=0, engine='R'){
    sample <- hs$gibbs
    par <- parmap(ncol(sample))
    if(missing(theta))     theta <- setNames(rep(2, length(par)), par)
    glR <- generatelogLik(sample[,j], sample[,-j, drop=FALSE], lambda=lambda)
    grR <- generatelogLik(sample[,j], sample[,-j, drop=FALSE], lambda=lambda, returnGrad=TRUE)
    if(engine=='R'){
        ## numeric
        gnum <- grad(glR, theta)
        ganalytic <- grR(theta)
} else if(engine=='cpp'){
    hl <- wrapHurdleLikelihood(sample[,j], sample[,-j, drop=FALSE], theta=theta, lambda=lambda)
    gnum <- grad(wrapLLall, x=theta, hl=hl,penalize=FALSE)
    ganalytic <- wrapGradAll(hl, theta, penalize=FALSE)
    expect_equal(wrapLLall(hl, theta, penalize=FALSE), glR(theta), tolerance=1e-6, check.names=FALSE)
} else{
    stop('unreachable')
}
    
    expect_equal(gnum, ganalytic,tolerance=1e-4, check.names=FALSE)
}
