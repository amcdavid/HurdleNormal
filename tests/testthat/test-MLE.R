context('Simulate from trivial distr')
G <- diag(c(-14, -11, -15))
H <- diag(c(5, 5, 8))
K <- diag(c(1, 1, 2))
## G <- diag(c(-13, -22))
## H <- diag(c(5, 10))
## K <- diag(c(1, 2))

#no thinning necessary because the gibbs sampler always moves??
library(Biobase)
library(plyr)
library(nloptr)
getLimit <- function(rgh, G, H, K){
    colnames(rgh) <- LETTERS[seq_len(nrow(G))]
    eset <- ExpressionSet(t(rgh))
    sc <- as(eset, 'SingleCellAssay')
    series <- lapply(c(3e2, 15e2, 75e2), seq, from=1)
    gj <- fit <- NULL
    fitseries <- laply(series, function(sub) {
        #browser()
        fit <<- fitZifNetwork(sc[sub,], response='cg.mle', additive.effects='0', gene.predictors='zero.inflated', precenter=FALSE, debug=FALSE, min.freq=0)
        gj <<- getJoint210(fit)
        #browser(expr=length(sub)>2e3)
        c(Gerror=sum((gj$G-G)^2), Herror=sum((gj$H-H)^2), Kerror=sum((gj$K-K))^2)
    })
    structure(fitseries, joint=gj, fit=fit)
}

library(numDeriv)
checkGrad <- function(sample, j=1, theta, lambda=0){
    gl <- generatelogLik(sample[,j], sample[,-j, drop=FALSE], lambda=lambda)
    gr <- generatelogLik(sample[,j], sample[,-j, drop=FALSE], lambda=lambda, returnGrad=TRUE)
    par <- get('par', environment(gl))
    if(missing(theta))     theta <- setNames(rep(2, length(par)), par)
    expect_equal(grad(gl, theta), gr(theta),tolerance=1e-4)
}

checkConverge <- function(sample, limit){
    fit <- attr(limit, 'fit')
    for(j in 1:nrow(fit)){
        checkGrad(sample, j, fit[[j,1]]$coef)
    }
}


suff.stat <- function(sample){
    s1 <- sample>0
    colnames(s1) <- paste0('1', colnames(sample))
    s1x1y <- crossprod(s1)
    sxy <- crossprod(sample)
    s1yx <- crossprod(s1, sample)
    list( s1x1y=s1x1y, sxy=sxy, s1yx=s1yx)
}

context('MV hurdle MLE')
## test_that('Finite difference gradients approximate analytic', {
##     ee <- exprs(vbetaT)    
##     gl <- generatelogLik(ee[,1], ee[,-1, drop=FALSE], lambda=1)
##     par <- get('par', environment(gl))
##     th <- setNames(rep(2, length(par)), par)

## })

rgh.adjust <- function(rgh){
    rgh2 <- addPseudoCounts(rgh)
    colnames(rgh2) <- LETTERS[seq_len(ncol(rgh2))]
    eset <- ExpressionSet(t(rgh2))
    sc <- as(eset, 'SingleCellAssay')
    fit <-   fitZifNetwork(sc, response='cg.mle', additive.effects='0', gene.predictors='zero.inflated', precenter=FALSE, debug=FALSE, min.freq=0)
    gj <- getJoint210(fit)
    rgh3 <- rGibbsHurdle(gj$G, gj$H, gj$K, Nt=nrow(rgh)*4+1000, burnin=1000)
    rg3 <- rgh3[seq(from=1, to=nrow(rgh), by=4),]
    list(rgh=rgh3, parm=gj)
}

## set.seed(1234)
## test_that('Converge under independence', {
##     rgh <- rGibbsHurdle(G, H, K, 1e4, 1e3)
##     checkGrad(rgh)
##     err <- getLimit(rgh, G, H, K)
##     expect_false(any(diff(err)>0))
## })


## test_that('Converge under K dependence', {
##     K[1,2] <- K[2,1] <- -.2
##     G[1,1] <- G[1,1]-7
##     G[2,2] <- G[2,2]-7
##     rgh <- rGibbsHurdle(G, H, K, 9e3, 1e3)
##     browser()
##     refit <- rgh.adjust(rgh)
##     err <- getLimit(refit$rgh, refit$parm$G, refit$parm$H, refit$parm$K)
##     expect_false(any(diff(err)>0))
## })

## test_that('Converge under G dependence', {
##     G[1,2] <- G[2,1] <- 3
##     G[1,1] <- G[1,1] - 4
##     G[2,2] <- G[2,2] -4
##     rgh <- rGibbsHurdle(G, H, K, 1e4, 1e3)
##     refit <- rgh.adjust(rgh)
##     err <- getLimit(refit$rgh, refit$parm$G, refit$parm$H, refit$parm$K)
##     expect_false(any(diff(err)>0))
## })

test_that('Converge under Hupper dependence', {
    H[1,2] <- -2
    G[2,2] <- G[1,1] <- G[1,1]+4
    rgh <- rGibbsHurdle(G, H, K, 1e4, 1e3)
    refit <- rgh.adjust(rgh)
    err <- getLimit(refit$rgh, refit$parm$G, refit$parm$H, refit$parm$K)
    expect_false(any(diff(err)>0))
})

test_that('Converge under Hlower dependence', {
    H[2,1] <- -2
    G[2,2] <- G[2,2]+3
    G[1,1] <- G[1,1]+3
    rgh <- rGibbsHurdle(G, H, K, 1e4, 1e3)
    refit <- rgh.adjust(rgh)
    err <- getLimit(refit$rgh, refit$parm$G, refit$parm$H, refit$parm$K)
    expect_false(any(diff(err)>0))
})

context("Penalized likelihood")
test_that('Penalized Gradient matches', {
    rgh <- rGibbsHurdle(G, H, K, 1e4, 1e3)
    checkGrad(rgh, j=1, lambda=0)
    checkGrad(rgh, j=1, lambda=1)
    checkGrad(rgh, j=3, lambda=3)
    checkGrad(rgh, j=3, lambda=.5)    
})



library(numDeriv)
context('MV hurdle MLE')
test_that('Finite difference gradients approximate analytic', {
    ee <- exprs(vbetaT)    
    gl <- generatelogLik(ee[,1], ee[,-1, drop=FALSE])
    dr <- generatelogLik(ee[,1], ee[,-1, drop=FALSE], returnGrad=TRUE)
    par <- get('par', environment(gl))
    th <- setNames(rep(2, length(par)), par)
    expect_equal(grad(gl, th), dr(th), tolerance=1e-4)

    gl <- generatelogLik(ee[,1], ee[,-1, drop=FALSE], lambda=1)
    dr <- generatelogLik(ee[,1], ee[,-1, drop=FALSE], returnGrad=TRUE, lambda=1)
    par <- get('par', environment(gl))
    th <- setNames(rep(2, length(par)), par)
    expect_equal(grad(gl, th), dr(th),tolerance=1e-4)
})
