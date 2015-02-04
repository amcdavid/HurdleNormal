context('Test convergence of gradient and MLE')
source('common.R', local=TRUE)
## for scenario indep, Gdep, Hdep, Kdep
## check gradient
## check convergence
## check penalized gradient
library(plyr)
library(nloptr)
getLimit <- function(hs){
    thin=.1
    burnin=1000
    fitseries <- laply(c(50, 250, 1250), function(sub) {
        n <- sub/thin+burnin
        hs <- suppressMessages(getGibbs(hs, Nt=n, burnin=burnin, thin=thin))
        fit <- getConditionalMLE(hs, testGrad=TRUE) 
        gj <- getJoint(fit)
        err <- sapply(c('G', 'H', 'K'), function(i) sum((hs$true[[i]]-gj[[i]])^2))
        err
    })
    fitseries
}

library(numDeriv)
checkGrad <- function(hs, j=1, theta, lambda=0){
    sample <- hs$gibbs
    gl <- generatelogLik(hs$gibbs[,j], sample[,-j, drop=FALSE], lambda=lambda)
    gr <- generatelogLik(sample[,j], sample[,-j, drop=FALSE], lambda=lambda, returnGrad=TRUE)
    par <- get('par', environment(gl))
    if(missing(theta))     theta <- setNames(rep(2, length(par)), par)
    expect_equal(grad(gl, theta), gr(theta),tolerance=1e-4)
}


## suff.stat <- function(sample){
##     s1 <- sample>0
##     colnames(s1) <- paste0('1', colnames(sample))
##     s1x1y <- crossprod(s1)
##     sxy <- crossprod(sample)
##     s1yx <- crossprod(s1, sample)
##     list( s1x1y=s1x1y, sxy=sxy, s1yx=s1yx)
## }

set.seed(1234)
test_that('Converge under independence', {
    checkGrad(Indep)
    err <- getLimit(Indep)
    expect_false(any(diff(err)>0))
})


test_that('Converge under K dependence', {
    checkGrad(Kdep2)
    err <- getLimit(Kdep2)
    expect_false(any(diff(err)>0))
})

test_that('Converge under G dependence', {
    checkGrad(Gdep)
    err <- getLimit(Gdep)
    expect_false(any(diff(err)>0))

})

test_that('Converge under Hupper dependence', {
    checkGrad(Hupdep)
    err <- getLimit(Hupdep)
    expect_false(any(diff(err)>0))
})

test_that('Converge under Hlower dependence', {
            checkGrad(Hlodep)
    err <- getLimit(Hlodep)
    expect_false(any(diff(err)>0))

})

context("Penalized likelihood")
test_that('Penalized Gradient matches', {
    checkGrad(Indep, j=1, lambda=0)
    checkGrad(Indep, j=1, lambda=1)
    checkGrad(Indep, j=3, lambda=3)
    checkGrad(Indep, j=3, lambda=.5)    
})
