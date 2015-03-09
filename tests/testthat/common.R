Indep <- HurdleStructure(G=-13*diag(3), H=diag(3)*5, K=diag(3))
Gdep <- Indep
Gdep$true$G[1,2] <- Gdep$true$G[2,1] <- 1.5
Kdep1 <- Indep
Kdep1$true$K[1,2] <- Kdep1$true$K[2,1] <- -.7
Kdep2 <- Kdep1
Kdep2$true$H[1,2] <- Kdep2$true$H[2,1] <- -2.25
Hupdep <- Indep
Hupdep$true$H[1,2] <- 1.5
Hupdep$true$G[2,2] <- Hupdep$true$G[1,1] <- -18.5
Hlodep <- Indep
Hlodep$true$H[2,1] <- 1.5
Hlodep$true$G[2,2] <- Hlodep$true$G[1,1] <- -18.5

Indep <- getGibbs(Indep)
Gdep <- getGibbs(Gdep)
Kdep2 <- getGibbs(Kdep2)
Hupdep <- getGibbs(Hupdep)
Hlodep <- getGibbs(Hlodep)


library(plyr)
library(nloptr)
getLimit <- function(hs, engine='R'){
    thin=.1
    burnin=1000
    series <- c(50, 250, 1250)
    fitseries <- matrix(NA, nrow=length(series), ncol=3)
    colnames(fitseries) <- c("G", "H", "K")
    for(i in seq_along(series)){
        sub <- series[i]
        n <- sub/thin+burnin
        hs <- suppressMessages(getGibbs(hs, Nt=n, burnin=burnin, thin=thin))
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
    hl <- HurdleLikelihood(sample[,j], sample[,-j, drop=FALSE], theta=theta, lambda=lambda)
    gnum <- grad(hl$LLall, theta)
    ganalytic <- hl$gradAll(theta)
    expect_equal(hl$LLall(theta), glR(theta), tolerance=1e-6, check.names=FALSE)
    expect_equal(hl$gradAll(theta), grR(theta), tolerance=1e-6, check.names=FALSE)
} else{
    stop('unreachable')
}
    
    expect_equal(gnum, ganalytic,tolerance=1e-4, check.names=FALSE)
}
