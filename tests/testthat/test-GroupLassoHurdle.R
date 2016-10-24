G <- diag(c(-14, -14, -15))
H <- diag(c(5, 5, 8))
K <- diag(c(1, 1, 2))
G[1,3] <- G[3,1] <- -1


rgh <- rGibbsHurdle(G, H, K, 8e3, 1e3, thin=.5)
colnames(rgh) <- LETTERS[1:ncol(rgh)]
context("Blocks")
this.model1 <- makeModel(rgh[,-1, drop=FALSE], conditionalCenter=TRUE,scale=FALSE)
test_that('Can Generate generically', {
    m <- ncol(this.model1)
    blk1 <- Block(this.model1)
    expect_equal(unique(blk1$map$block), 1:3)
    ## block 1 is intercepts (non penalized)
    expect_equal(blk1$map[block==1,lambda], c(0, 0, 0))
    ## block 1 parameters should occur at position 1, m + 1, 2*m +1 (precision kbb here)
    expect_equal(blk1$map[block==1,paridx], c(1, m+1, 2*m+1))
    ## block 2 parameters should occur at position 2, 3, m+2, m+3
    expect_equal(blk1$map[block==2,paridx], c(2, 3, m+2, m+3))
    expect_equal(blk1$map[,nodeId], c('(fixed)', 'B', 'B', 'C', 'C', '(fixed)', 'B', 'B', 'C', 'C', '(fixed)'))
})

test_that('Can Generate with fixed columns', {
    this.modelF <- makeModel(rgh[,-1, drop=FALSE], fixed=cbind(1, Z=rnorm(nrow(rgh))))
    m <- ncol(this.modelF)
    blk1 <- Block(this.modelF)    
    expect_equal(unique(blk1$map$block), 1:3)
    ## block 1 is intercepts (non penalized)
    expect_equal(blk1$map[block==1,lambda], c(0, 0, 0, 0, 0))
    ## block 1 parameters should occur at position 1, 2, m + 1, m+2, 2*m +1 (precision kbb)
    expect_equal(blk1$map[block==1,paridx], c(1, 2, 1+m, 2+m, 2*m+1))
    expect_equal(blk1$map[,nodeId], c('(fixed)', '(fixed)', 'B', 'B', 'C', 'C', '(fixed)', '(fixed)', 'B', 'B', 'C', 'C', '(fixed)'))
})


context('Fit one penalized')
test_that('Can fit', {
    y.zif <- rgh[,1]
    #this.model <- rgh[,-1,drop=FALSE]
    theta <- c(gbb=-10, 0, 0, hbb=10, 0, 0, 0, 0, 0, 0, 1)+.02
    names(theta) <- parmap(3)
    lambda <- c(2.8, 1, .8, .7, .6, .5, .3, .1, .1, .05, .05, .01, 0)
    paths1 <- cgpaths(y.zif, this.model1, Block(this.model1), lambda=lambda, control=list(tol=5e-4, maxrounds=1000, debug=1, stepsize=1, stepexpand=.01, newton0=TRUE), penaltyFactor='full')
    expect_true(inherits(paths1, 'SolPath'))
    #paths2 <- cgpaths(y.zif, makeModel(rgh[,-1, drop=FALSE]), lambda=lambda, control=list(tol=1e-5, maxrounds=2000, debug=0, method='block'), standardize=FALSE)
    hm <- HurdleLikelihood(y.zif, this.model1, theta=paths1$path[13,], lambda=0)
    thetareg <- optim(paths1$path[13,], hm$LLall, hm$gradAll, method='L-BFGS-B', hessian=TRUE, control=list(pgtol=1-8, maxit=1e4, factr=1e5))
    distreg <- (paths1$path[13,]-thetareg$par)
    manoblis <- crossprod(distreg, thetareg$hess) %*% distreg
    expect_lt(manoblis, sqrt(1e-6))
})
