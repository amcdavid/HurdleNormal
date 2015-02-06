context('Test convergence of gradient and MLE')
source('common.R', local=TRUE)
## for scenario indep, Gdep, Hdep, Kdep
## check gradient
## check convergence
## check penalized gradient


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

