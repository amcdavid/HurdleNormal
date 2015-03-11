context("Cpp")
source('common.R', local=TRUE)
set.seed(1234)
test_that('Cpp: Converge under independence', {
    checkGrad(Indep, engine='cpp')
    err <- getLimit(Indep, engine='cpp')
    expect_false(any(diff(err)>0))
})

test_that('Unpenalized likelihood matches', {
    sample <- Indep$gibbs
    par <- parmap(ncol(sample))
    theta <- setNames(rep(2, length(par)), par)
    hl <- HurdleLikelihood(sample[,1], sample[,-1, drop=FALSE], theta=theta, lambda=1)
    hl0 <- HurdleLikelihood(sample[,1], sample[,-1, drop=FALSE], theta=theta, lambda=0)
    expect_equal(hl$LLall(theta, penalize=FALSE), hl0$LLall(theta, penalize=TRUE))   

})


test_that('Cpp: Converge under K dependence', {
    checkGrad(Kdep2, engine='cpp')
    err <- getLimit(Kdep2, engine='cpp')
    expect_false(any(diff(err)>0))
})

test_that('Cpp: Converge under G dependence', {
    checkGrad(Gdep, engine='cpp')
    err <- getLimit(Gdep, engine='cpp')
    expect_false(any(diff(err)>0))

})

test_that('Cpp: Converge under Hupper dependence', {
    checkGrad(Hupdep, engine='cpp')
    err <- getLimit(Hupdep, engine='cpp')
    expect_false(any(diff(err)>0))
})

test_that('Cpp: Converge under Hlower dependence', {
    checkGrad(Hlodep, engine='cpp')
    err <- getLimit(Hlodep, engine='cpp')
    expect_false(any(diff(err)>0))

})

test_that('Cpp: Penalized Gradient matches', {
    checkGrad(Indep, j=1, lambda=0, engine='cpp')
    checkGrad(Indep, j=1, lambda=1, engine='cpp')
    checkGrad(Indep, j=3, lambda=3, engine='cpp')
    checkGrad(Indep, j=3, lambda=.5, engine='cpp')
})
