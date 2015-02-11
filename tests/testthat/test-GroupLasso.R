G <- diag(c(-15, -25, -15))
H <- diag(c(5, 5, 8))
K <- diag(c(1, 1, 2))
K[1,2] <- 0
H[1,2] <- 2

context('Fit one penalized')
rgh <- rGibbsHurdle(G, H, K, 5e3, 1e3)
test_that('Can fit', {
    theta <- c(gbb=-10, 0, 0, hbb=10, 0, 0, 0, 0, 0, 0, 1)+.02
    names(theta) <- parmap(3)
    y.zif <- rgh[,1]
    this.model <- rgh[,-1,drop=FALSE]
    cgpaths(theta, y.zif, this.model, lambda=c(8, 1, .5, .2, .05, 1e-6), control=list(tol=1e-8, maxrounds=200, debug=3))

    hm <- HurdleLikelihood(y.zif, this.model, theta=theta, lambda=0)
    thetareg <- optim(theta, hm$LLall, hm$gradAll, method='BFGS')
    expect_equivalent(cgpaths[5,], thetareg$par)

})
