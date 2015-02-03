context('Log Density')
library(plyr)
library(reshape)


Gdep <- HurdleStructure(G=matrix(c(-2, 1, 1, -2), nrow=2),
                         H =matrix(c(2, 0, 0, 2), nrow=2),
                         K =diag(2))
Hlowdep <- HurdleStructure(G=diag(2)*-4, H=matrix(c(2, 1, 0, 2), nrow=2),K=diag(2))
Hupdep <- HurdleStructure(G=diag(2)*-4, H=matrix(c(2, 0, 1, 2), nrow=2),K=diag(2))



EPS <- .1
fun <- getDensity(Gdep, tol=EPS/2)
grid <- .evalGridHurdleMeasure(c(-1, -1), c(4,4), EPS, fun=fun, tol=EPS/2)$feval

test_that('Symmetric density is symmetric', {
    expect_equivalent(grid, t(grid))
})

context('Conditional agrees with joint density')
## get a sample from the conditional
## compare the odds r
test_that('Conditional sampler for G dependence', {
    set.seed(123)
    EA <- 1e3
    x <- matrix(rep(c(0, EPS), each=EA), nrow=2*EA)
    TOL <- EPS/2
    y <- with(Gdep$true, rCondHurdle210(x, 2, G, H, K, tol=TOL))
    yz <- mean(abs(y[seq_len(EA)])>TOL)
    ynz <- mean(abs(y[EA+seq_len(EA)])>TOL)
    logOR <- -log(yz/ynz*(1-ynz)/(1-yz))
    eOR <- with(Gdep$true, 2*G[1,2]+((H[1,1]+H[2,1])^2 - H[1,1]^2)/(2*K[1,1]))
    expect_equal(logOR, eOR, tolerance=10/sqrt(EA))

    p10 <- 1-grid['0','0']/sum(grid['0',])
    p11 <- 1-grid['0', '0.1']/sum(grid['0.1',])
    expect_equal(p10, yz, tolerance=10/sqrt(EA))
    expect_equal(p11, ynz, tolerance=10/sqrt(EA))
})

test_that('Higher order conditional sampler', {
    G <- matrix(c(2, 1, 1,
                  1,-2, 1,
                  1, 1, -2), nrow=3)
    H <- matrix(c(5, -1, -1,
                  -1, 2, -1,
                  -1, -1, 2), nrow=3)
    K <- diag(3)
    x <- matrix(rep(c(0, EPS), each=200), nrow=200)
    y <- rCondHurdle210(x, 3, G, H, K, tol=EPS/2)
})

context('Gibbs approximates joint')
test_that('Marginals for Gdep', {
    set.seed(1234)
    Gdep <- getGibbs(Gdep)
    samp <- Gdep$gibbs
    df <- data.frame(samp, nz=abs(samp)>.05)
    ## compare 2d contour plots to 2d density estimate
    ## ggplot(df, aes(x=X1, y=X2))+geom_density2d()+  geom_contour(aes(x=x1, y=x2, z=ap), data.frame(x, ap))
    cs <- colSums(grid)
    Fn <- cumsum(cs)
    Z <- Fn[length(Fn)]
    F <- data.frame(F=Fn/Z, q=as.numeric(names(Fn)))
    f <- data.frame(f=cs/Z, q=as.numeric(names(Fn)))
    f$counts <- f$f*nrow(df)
    ## graphical comparison
    ##browser()
    ##ggplot(df, aes(x=X1))+stat_ecdf() + geom_step(aes(x=q, y=F), data=F)

    
    suppressWarnings(gibbs1d <- chisq.test(table(cut(df$X1, breaks=c(f$q, max(f$q)+EPS), right=FALSE)), p=f$f))
    expect_more_than(gibbs1d$p.value, .05)
})

adjustedCondDistr <- function(hs){
    offt1 <- with(hs$true, -.5*log(K[1,1]/(2*pi)) + (H[1,1]-K[1,2]*hs$gibbs[,2]+H[2,1]*(hs$gibbs[,2]>0))^2/(2*K[1,1]))
    offt2 <- with(hs$true, -.5*log(K[2,2]/(2*pi)) + (H[2,2]-K[2,1]*hs$gibbs[,1]+H[1,2]*(hs$gibbs[,1]>0))^2/(2*K[2,2]))
    gibbs <- as.data.frame(hs$gibbs)
    g12 <- glm(V1 >0 ~ V2+I(V2>0)+offset(offt1), data=gibbs, family='binomial')
    g21 <- glm(V2 >0 ~ V1+I(V1>0)+offset(offt2), data=gibbs, family='binomial')
    b12 <- lm(V1 ~ V2 + I(V2>0), data=gibbs, subset=V1>0)
    b21 <- lm(V2 ~ V1 + I(V1>0), data=gibbs, subset=V2>0)
    ldply(list(g12=g12, g21=g21, b12=b12, b21=b21), function(m){
      coef(m)/sqrt(diag(vcov(m)))  
    })
}

test_that('Conditionals for Hlowdep', {
    Hlowdep <- getGibbs(Hlowdep)
    adjustedCondDistr(Hlowdep)
    Hupdep <- getGibbs(Hupdep)
    adjustedCondDistr(Hupdep)    
})
