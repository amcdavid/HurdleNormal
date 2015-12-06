G <- diag(c(-14, -14, -15))
H <- diag(c(5, 5, 8))
K <- diag(c(1, 1, 2))
G[1,3] <- G[3,1] <- -1


rgh <- rGibbsHurdle(G, H, K, 8e3, 1e3, thin=.5)
context("Blocks")
this.model <- makeModel(rgh[,-1, drop=FALSE])
test_that('Can Generate Generically', {
    blk <- Block(this.model)
    expect_equal(unique(blk$map$block), 1:3)
    expect_equal(blk$map[block==1,lambda], c(0, 0, 0))
})


context('Fit one penalized')
test_that('Can fit', {
    y.zif <- rgh[,1]
    #this.model <- rgh[,-1,drop=FALSE]
    theta <- c(gbb=-10, 0, 0, hbb=10, 0, 0, 0, 0, 0, 0, 1)+.02
    names(theta) <- parmap(3)
    lambda <- c(2.8, 1, .8, .7, .6, .5, .3, .1, .1, .05, .05, .01, 0)

    paths1 <- cgpaths(y.zif, this.model, Block(this.model), lambda=lambda, control=list(tol=1e-3, maxrounds=400, debug=1, method='proximal', stepsize=1, updatehess=200, stepexpand=.05), penaltyFactor='full')
    #paths2 <- cgpaths(y.zif, makeModel(rgh[,-1, drop=FALSE]), lambda=lambda, control=list(tol=1e-5, maxrounds=2000, debug=0, method='block'), standardize=FALSE)
    hm <- HurdleLikelihood(y.zif, this.model, theta=paths1$path[13,], lambda=0)
    thetareg <- optim(paths1$path[13,], hm$LLall, hm$gradAll, method='L-BFGS-B', hessian=TRUE, control=list(pgtol=1-8, maxit=1e4))
    distreg <- (paths1$path[13,]-thetareg$par)
    manoblis <- crossprod(distreg, thetareg$hess) %*% distreg
    expect_less_than(manoblis, sqrt(1e-6))
})
