expit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

.checkArgs <- function(x, G, H, K){
    stopifnot(is.matrix(H))
    stopifnot(ncol(H)==nrow(H))
    stopifnot(dim(G)==dim(H))
    stopifnot(dim(K)==dim(H))
    if(!missing(x)) stopifnot(ncol(x)==ncol(H))
}

## joint density function
## x is row vector!!!
## But should be a column vector...
dHurdle210 <- function(x, G, H, K, tol=5e-2){
    if(length(dim(x))<2) x <- t(as.matrix(x))
    .checkArgs(x, G, H, K)
    xI <- (abs(x)>tol)*1
    logDens <- xI %*% G %*% t(xI) + (xI*x)%*%H%*%t(xI) - .5*(xI*x)%*%K%*%t(xI*x)
    exp(logDens)
}

## simulate from the conditional distribution of a hurdle model
## where i gives the indices being provided in x
## rCondHurdle210 <- function(x, i, G, H, K, tol=5e-2){
##     stopifnot(length(i) == ncol(x), length(i) == ncol(G)-1)
##     .checkArgs(cbind(1, x), G, H, K)
##     xI <- (abs(x)>tol)*1
##     midx <- seq_len(ncol(G))
##     noti <- setdiff(midx, i)
    
##     Gba <- G[noti,noti,drop=TRUE]+2*xI%*%G[i,noti,drop=FALSE] + (xI*x)%*%t(H[noti, i,drop=FALSE])
##     Hba <- H[noti, noti]+xI%*%H[i,noti,drop=FALSE] - (xI*x)%*%K[i,noti,drop=FALSE]
##     Kba <- K[noti, noti,drop=TRUE]

##     logitP <- Gba-.5*log(Kba/(2*pi))+Hba^2/(2*Kba)
##     mu <- Hba/Kba
##     yI <- runif(nrow(x))<expit(logitP)
##     y <- rnorm(nrow(x), sd=1/sqrt(Kba))+mu
##     y*yI
## }


##' Sample from the full conditional distribution
##'
##' Return a sample from f(x_i | x_{-i}), that is, a sample of the conditional distribution of the \code{i}th coordinates, given \code{x} aside from \code{x_i}.
##' 
##' \code{x} is assumed in the same order as the parameter vectors, aside from its missing ith element.
##' This function is vectorized over x.
##' @param x vector, or matrix, to be conditioned upon.
##' @param i coordinate to be returned
##' @param G binary dependence matrix
##' @param H location dependence matrix
##' @param K precision matrix
##' @param tol tolerance for calling a point to be numerically equal to zero.
##' @return numeric sample from x_i | x_{-i}
##' @export
rCondHurdle210 <- function(x, i, G, H, K, tol=5e-4){
    stopifnot(length(j) == 1 && ncol(x)  == ncol(G)-1)
    .checkArgs(G=G, H=H, K=K)
    if(is.matrix(x)) aaply(x, 1, .rCondHurdle, j=j, G=G, H=H, K=K, tol=tol) else .rCondHurdle(x, j, G, H, K, tol)
}


rv <- function(x){
    dim(x) <- c(1, length(x))
    x
}

## rGibbsHurdle <- function(G, H, K, Nt, burnin=floor(Nt/2), tol=5e-2){
##     p <- ncol(G)
##     .checkArgs(matrix(ncol=p), G, H, K)
##     ## samp is internally arrayed backwards from what rCondHurdle uses,
##     ## but it is much easier to have it in column-major order for the purposes of gibbs sampling
##     samp <- rep(NA, p*Nt)
##     samp[seq_len(p)] <-0
##     ## pointer to previous p-1 positions in vectorized sample
##     notp <- seq_len(p-1)+1
##     for(i in seq(p+1, Nt*p)){
##             ## and coordinates of these samples
##             notpidx <- ((notp-1) %% p)+1
##             samp[i] <- rCondHurdle210(rv(samp[notp]), i=notpidx, G=G, H=H, K=K, tol=tol)
##             if(is.na(samp[i])) stop("The gibbs sampler has gotten sad.")
##             notp <- notp+1
##         }
##     ## we are in column-major order, but will want to be in row-major for consistency elsewhere
##     t(matrix(samp, nrow=p))[-seq_len(burnin),]
## }

rGibbsHurdle <- function(G, H, K, Nt, burnin=floor(Nt/2), thin=1, tol=5e-4){
    p <- ncol(G)
    .checkArgs(matrix(ncol=p), G, H, K)
    ## samp is internally arrayed backwards from what rCondHurdle uses,
    ## but it is much easier to have it in column-major order for the purposes of gibbs sampling
    tmat <- .rGibbsHurdle(G, H, K, Nt, tol)
    mat <- t(tmat)[-seq_len(burnin),]
    keep <- seq(from=1, to=nrow(mat), by=round(1/thin))
    mat[keep,]
}


addPseudoCounts <- function(rgh){
    ## (1) At least 3 totally positive counts (so mean and covariance are estimible)
    ## (2) and every pair of coordinates has full 2x2 contingency table
    ## (so that binary dependence is estimible)
    ## Adding another p observations (singly positive) fills the other three cells in the contigency tables
    colmeans <- apply(rgh, 2, function(x) mean(x[abs(x)>0]))
    colmeans[is.na(colmeans)] <- 5
    pos <- matrix(rnorm((3+1)*ncol(rgh)), nrow=ncol(rgh))
    pos <- pos+colmeans
    allpos <- pos[,1:3]
    contigentPos <- pos[,4]
    contigentPos <- diag(contigentPos)
    rbind(rgh, t(allpos), t(contigentPos))
}
