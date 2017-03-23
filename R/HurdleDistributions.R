## Quell warnings about missing imports from base packages
##' @importFrom stats approx coef cor.test median na.omit optim optimize rnorm runif setNames var glm.fit
##' @importFrom utils file_test
##' @importFrom methods new
NULL

expit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

.checkArgs <- function(x, G, H, K){
    stopifnot(is.matrix(H))
    stopifnot(ncol(H)==nrow(H))
    stopifnot(dim(G)==dim(H))
    stopifnot(dim(K)==dim(H))
    stopifnot(all.equal(t(G), G, check.attributes=FALSE))
    stopifnot(all.equal(t(K), K, check.attributes=FALSE))
    eK <- eigen(K, only.values=TRUE)$values
    stopifnot(all(eK>0) &&  all(Im(eK)==0))
    if(!missing(x)) stopifnot(length(x) == ncol(H) || nrow(x)==ncol(H))
}

## joint density function
## if x is a matrix, it should be cbind'd column vectors
dHurdle210 <- function(x, G, H, K, tol=5e-2){
    if(length(dim(x))<2) x <- as.matrix(x)
    .checkArgs(x, G, H, K)
    xI <- (abs(x)>tol)*1
    x[xI==0] <- 0
    logDens <- t(xI) %*% G %*% xI + t(xI)%*%H%*%x - .5*t(x)%*%K%*%(x)
    exp(logDens)
}

##' Evaluate a Riemann sum on a grid of points that includes zero
##'
##' @param lower vector of lower limits
##' @param upper upper limits
##' @param dx scalar gridsize
##' @param fun function to evaluate.  First argument should be the vector point at which it's being eval'd
##' @param ... passed to fun
##' @param tol passed to fun
##' @return list of: reimann-weighted function on the grid/lattice (feval) and
##' breaks (lattice points)
##' @rdname eval_grid_hurdle_measure
.evalGridHurdleMeasure <- function(lower, upper, dx, fun, ..., tol=5e-2){
    dxli <- limitli <- vector('list', length=length(lower))
    for(i in seq_along(lower)){
        limitli[[i]] <- unique(sort(c(0, seq(from=lower[i], to=upper[i], by=dx))))
        dxli[[i]] <- diff(limitli[[i]])
        ## cut off last interval since it's open
        limitli[[i]] <- limitli[[i]][-length(limitli[[i]])]
        dxli[[i]][abs(limitli[[i]])<tol] <- 1
    }
    x <- as.matrix(do.call(expand.grid, limitli))
    dx <- do.call(expand.grid, dxli)
    prodDx <- apply(dx, 1, prod)
    feval <- plyr::aaply(x, 1, fun, ...)*prodDx
    dim(feval) <- c(sapply(limitli, length))
    dimnames(feval) <- limitli
    list(feval=feval, breaks=limitli)
}

.additiveVstate <- function(p){
    rv <- cv <- states <- rep(NA, length=p*(p+1)/2)
    i <- 1
    for(r in 0:(p-1)){
        for(c in r:(p-1)){
            states[i] <- if(r==c) 2^r else 2^r + 2^c
            print(paste0('state = ', states[i], ' r=', r, ' c=', c))
            rv[i] <- r+1
            cv[i] <- c+1
            i <- i+1
        }
    }
    list(states=states, indices=cbind(rv, cv))
}

.calcNormalizing <- function(statevec, subsample=Inf, v, G, H, K, returnCovariance=FALSE, log=TRUE){
    .checkArgs(G=G, H=H, K=K)
    p <- ncol(H)
    
    if(!missing(v)){
        statevec <- sum(2^(v-1))
    }
    
    if(missing(statevec)){
        nstate <- 2^p-1
        statevec <- if(subsample<Inf) stop('Not implemented') else  seq(1L, nstate)
    } else{
        nstate <- length(statevec)
    }
    if(nstate>=1e6) stop('Too hazardous')
    logP <- N <- rep(NA, nstate)
    #N[1] <- 1
    state <- matrix(FALSE, nrow=p, ncol=nstate)
    mu <- matrix(0, nrow=p, ncol=nstate)
    if(returnCovariance) Sigma <- array(0, dim=c(p, p, nstate)) else Sigma <- NULL
    for(i in seq_len(nstate)){
        v <- statevec[i]
        bits <- as.logical(intToBits(v)[1:p])
        state[,i] <- bits
        subK <- K[bits,bits,drop=FALSE]
        cov <- solve(subK)
        if(returnCovariance) Sigma[bits,bits,i] <- cov
        subH <- t(crossprod(bits, H))[bits,,drop=F]      #H order is wrong throughout, since 1y is supposed to be a selector, eg eq 6?  Can we fix by inverting 6?
        submu <- cov %*% subH
        mu[bits,i] <- submu
        N[i] <- .5*(t(subH) %*% submu -  determinant(subK/(2*pi), logarithm=TRUE)$modulus)
        logP[i] <- sum(G[bits, bits]) + N[i]
        if(!log) N[i] <- exp(N[i])

    }
    list(N=N, state=state, mu=mu, logP=logP,Sigma=Sigma)
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
##' @param j coordinate to be returned
##' @param G binary dependence matrix
##' @param H location dependence matrix
##' @param K precision matrix
##' @param tol tolerance for calling a point to be numerically equal to zero.
##' @return numeric sample from x_i | x_{-i}
##' @export
rCondHurdle210 <- function(x, j, G, H, K, tol=5e-4){
    stopifnot(length(j) == 1 && ncol(x)  == ncol(G)-1)
    .checkArgs(G=G, H=H, K=K)
    if(is.matrix(x)) plyr::aaply(x, 1, .rCondHurdle, j=j-1, G=G, H=H, K=K, tol=tol) else .rCondHurdle(x, j-1, G, H, K, tol)
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


#' Sample from a multivariate hurdle model
#'
#' @param G symmetric discrete interaction matrix
#' @param H unstructured location matrix
#' @param K Symmetric positive definite conditional precision matrix
#' @param Nt Number of unthinned samples, including burnin
#' @param burnin how many samples to discard from burn in
#' @param thin how many samples to thin
#' @param tol Numeric tolerance for zero
#' @param Nkeep (optional) number of samples, post-burnin and thinning
#' @return matrix of (Nt-burnin)*thin samples
#' @export
#' @examples
#' G = matrix(c(-15, 1, 0,
#' 1, -15, 1.5,
#' 0, 1.5, -15), nrow=3)
#' H = diag(5, nrow=3)
#' K = diag(1, nrow=3)
#' y = rGibbsHurdle(G, H, K, 2000, thin = .2, burnin=1000)
rGibbsHurdle <- function(G, H, K, Nt, burnin=500, thin=.1, tol=5e-4, Nkeep=500){
    p <- ncol(G)
    .checkArgs(matrix(ncol=p), G, H, K)
    ## coords X samples
    if(!missing(Nkeep)){
        if(!missing(Nt)) stop("Provide only one of `Nt` and `Nkeep`")
        Nt <- (Nkeep)*round(1/thin) + burnin
    }
    if(burnin>Nt) stop("More burnin samples than total samples!")
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
