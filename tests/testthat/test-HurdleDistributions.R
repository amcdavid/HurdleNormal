context('Log Density')
library(plyr)
library(reshape)


G <- matrix(c(-2, 1, 1, -2), nrow=2)
H <- matrix(c(2, -1, -1, 2), nrow=2)
K <- diag(2)
EPS <- .01
xidx <- sort(c(seq(-2, 5, by=.2), EPS))
dx <- diff(xidx)
dx[dx==.01] <- 1
dx <- c(dx, .2)
x <- as.matrix(expand.grid(x1=xidx, x2=xidx))
ap <- aaply(x, 1, dHurdle210, G=G, H=H, K=K, .inform=TRUE, tol=EPS/2)
grid <- as.matrix(cast(data.frame(x, ap), x1 ~ x2, value='ap'))


test_that('Symmetric density is symmetric', {
    expect_equivalent(grid, t(grid))
})

context('Conditional agrees with joint density')
## get a sample from the conditional
## compare the odds r
test_that('Conditional sampler', {
    set.seed(123)
    EA <- 1e3
    x <- matrix(rep(c(0, EPS), each=EA), nrow=2*EA)
    TOL <- EPS/2
    y <- rCondHurdle210(x, 2, G, H, K, tol=TOL)
    yz <- mean(abs(y[seq_len(EA)])>TOL)
    ynz <- mean(abs(y[EA+seq_len(EA)])>TOL)
    logOR <- -log(yz/ynz*(1-ynz)/(1-yz))
    eOR <- 2*G[1,2]+((H[1,1]+H[2,1])^2 - H[1,1]^2)/(2*K[1,1])
    expect_equal(logOR, eOR, tolerance=10/sqrt(EA))

    p10 <- 1-grid['0','0']/sum(grid['0',]*dx)
    p11 <- 1-grid['0', '0.01']/sum(grid['0.01',]*dx)
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
    x <- matrix(rep(c(0, EPS), each=2*EA), nrow=2*EA)
    y <- rCondHurdle210(x, 3, G, H, K, tol=EPS/2)
})

context('Gibbs approximates joint')
test_that('Can sample', {
    set.seed(1234)
    samp <- rGibbsHurdle(G, H, K, 1e4, thin=.5, tol=EPS/2)
    df <- data.frame(samp, it=rep(seq_len(nrow(samp)/nchain), each=nchain), chain=rep(seq_len(nchain), times=nrow(samp)/nchain), nz=abs(samp)>.05)
    ## compare 2d contour plots to 2d density estimate
    ## ggplot(df, aes(x=X1, y=X2))+geom_density2d()+  geom_contour(aes(x=x1, y=x2, z=ap), data.frame(x, ap))
    cs <- colSums(grid*dx)
    Fn <- cumsum(cs*dx)
    Z <- Fn[length(Fn)]
    F <- data.frame(F=Fn/Z, q=as.numeric(names(Fn)))
    f <- data.frame(f=cs*dx/Z, q=as.numeric(names(Fn)))
    f$counts <- f$f*nrow(df)
    ## graphical comparison
     ggplot(df, aes(x=X1))+stat_ecdf() + geom_step(aes(x=q, y=F), data=F)

    
    suppressWarnings(gibbs1d <- chisq.test(table(cut(df$X1, breaks=c(f$q, max(f$q)+max(dx)), right=FALSE)), p=f$f))
    expect_more_than(ct$p.value, .05)
})
