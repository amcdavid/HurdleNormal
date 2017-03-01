context("Cpp")
source('common.R', local=TRUE)
set.seed(1234)

test_that('unpack is inverse of pack', {
    theta <- 1:11
    ptheta <- unpackTheta(theta)
    pptheta <- unpackTheta(ptheta, 'pack')
    expect_equal(theta, pptheta)

})

test_that('Cpp: Converge under independence', {
    checkGrad(Indep, engine='cpp')
    err <- getLimit(Indep, engine='cpp')
    ## error is monotone (within tolerance defined by final error)
    expect_false(any(t(diff(err))>err[3,]))
})

## test_that('Unpenalized likelihood matches', {
##     sample <- Indep$gibbs
##     par <- parmap(ncol(sample))
##     theta <- setNames(rep(2, length(par)), par)
##     hl <- wrapHurdleLikelihood(sample[,1], sample[,-1, drop=FALSE], theta=theta, lambda=1)
##     hl0 <- wrapHurdleLikelihood(sample[,1], sample[,-1, drop=FALSE], theta=theta, lambda=0)
##     expect_equal(hl$LLall(theta, penalize=FALSE), hl0$LLall(theta, penalize=TRUE))   

## })


test_that('Cpp: Converge under K dependence', {
    checkGrad(Kdep2, engine='cpp')
    err <- getLimit(Kdep2, engine='cpp')
    expect_false(any(t(diff(err))>err[3,]))
})

test_that('Cpp: Converge under G dependence', {
    checkGrad(Gdep, engine='cpp')
    err <- getLimit(Gdep, engine='cpp')
    expect_false(any(t(diff(err))>err[3,]))

})

test_that('Cpp: Converge under Hupper dependence', {
    checkGrad(Hupdep, engine='cpp')
    err <- getLimit(Hupdep, engine='cpp')
    expect_false(any(t(diff(err))>err[3,]))
})

test_that('Cpp: Converge under Hlower dependence', {
    checkGrad(Hlodep, engine='cpp')
    err <- getLimit(Hlodep, engine='cpp')
    #browser()
    expect_false(any(t(diff(err))>err[3,]))

})

## test_that('Cpp: Penalized Gradient matches', {
##     checkGrad(Indep, j=1, lambda=0, engine='cpp')
##     checkGrad(Indep, j=1, lambda=1, engine='cpp')
##     checkGrad(Indep, j=3, lambda=3, engine='cpp')
##     checkGrad(Indep, j=3, lambda=.5, engine='cpp')
## })
