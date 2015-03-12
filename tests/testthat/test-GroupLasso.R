G <- diag(c(-16, -16, -15))
H <- diag(c(5, 5, 8))
K <- diag(c(1, 1, 2))
G[1,3] <- G[3,1] <- 1

context('Fit one penalized')
rgh <- rGibbsHurdle(G, H, K, 8e3, 1e3)
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
    lambda <- c(2.8, 1, .8, .8, .6, .6, .3, .3, .1, .1, .05, .01, 0)
    browser()
    paths1 <- cgpaths(y.zif, this.model, lambda=lambda, control=list(tol=1e-5, maxrounds=300, debug=2, method='proximal', stepsize=3, precondition=TRUE), standardize=FALSE)
    paths2 <- cgpaths(y.zif, this.model, lambda=lambda, control=list(tol=1e-5, maxrounds=300, debug=2, method='proximal', stepsize=1, precondition=FALSE), standardize=FALSE)
    paths2 <- cgpaths(y.zif, this.model, lambda=lambda, control=list(tol=1e-5, maxrounds=300, debug=1, method='block'), standardize=FALSE)
    hm <- HurdleLikelihood(y.zif, this.model, theta=theta, lambda=0)
    thetareg <- optim(theta, hm$LLall, hm$gradAll, method='BFGS')
    expect_equivalent(paths1$path[9,], thetareg$par)

})
