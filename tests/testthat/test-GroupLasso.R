G <- diag(c(-15, -25, -15))
H <- diag(c(5, 5, 8))
K <- diag(c(1, 1, 2))
K[1,2] <- 0
H[1,2] <- 2

context('Fit one penalized')
rgh <- rGibbsHurdle(G, H, K, 5e3, 1e3)
test_that('Can fit', {
    theta <- c(gbb=-10, 0, 0, hbb=10, 0, 0, 0, 0, 0, 0, 1)
    names(theta) <- parmap(3)
    y.zif <- rgh[,1]
    this.model <- rgh[,-1,drop=FALSE]
    thetapen <- solvePen(theta, 0, y.zif, this.model, control=list(tol=1e-6, maxrounds=100, maxit=200, debug=TRUE))
    ll <- generatelogLik(y.zif, this.model, 0)
    grad <- generatelogLik(y.zif, this.model, returnGrad=TRUE)
    thetareg <- optim(theta, ll, grad, method='BFGS')

    thetapen <- solvePen(theta, .1, y.zif, this.model, control=list(tol=1e-6, maxrounds=100, maxit=200, debug=3))
    cgpaths(theta, y.zif, this.model, lambda=c(.3, .2, .1, .05, 0), control=list(tol=1e-8, maxrounds=100, maxit=50, debug=3))
})
