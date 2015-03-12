G <- diag(c(-14, -14, -15))
H <- diag(c(5, 5, 8))
K <- diag(c(1, 1, 2))
G[1,3] <- G[3,1] <- -1

context('Fit one penalized')
rgh <- rGibbsHurdle(G, H, K, 8e3, 1e3, thin=.5)
test_that('Can fit', {
    ## rgh <- apply(rgh, 2, function(col){
    ##     xI <- abs(col)>1e-7
    ##     col[xI] <- scale(col[xI])
    ##     col
    ## })
    y.zif <- rgh[,1]
    this.model <- rgh[,-1,drop=FALSE]
    theta <- c(gbb=-10, 0, 0, hbb=10, 0, 0, 0, 0, 0, 0, 1)+.02
    names(theta) <- parmap(3)
    lambda <- c(2.8, 1, .8, .7, .6, .5, .3, .1, .05, .01, 0)
    paths1 <- cgpaths(y.zif, this.model, lambda=lambda, control=list(tol=1e-9, maxrounds=2000, debug=0, method='proximal', stepsize=3, updatehess=200), standardize=FALSE)
    paths2 <- cgpaths(y.zif, this.model, lambda=lambda, control=list(tol=1e-5, maxrounds=2000, debug=0, method='block'), standardize=FALSE)
    hm <- HurdleLikelihood(y.zif, this.model, theta=paths1$path[13,], lambda=0)
    thetareg <- optim(theta, hm$LLall, hm$gradAll, method='BFGS', hessian=TRUE)
    distreg <- (paths1$path[13,]-thetareg$par)
    manoblis <- crossprod(distreg, thetareg$hess) %*% distreg
    distreg2 <- paths2$path[13,]-thetareg$par
    manoblis2 <- crossprod(distreg2, thetareg$hess) %*% distreg2
    distreg3 <- t(paths2$path-paths1$path)
    manoblis3 <- diag(crossprod(distreg3, thetareg$hess) %*% distreg3)
    expect_equivalent(paths1$path, paths2$path)
})
