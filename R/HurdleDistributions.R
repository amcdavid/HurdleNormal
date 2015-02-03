expit <- function(x) exp(x)/(1+exp(x))
logit <- function(x) log(x/(1-x))

.checkArgs <- function(x, G, H, K){
    stopifnot(is.matrix(H))
    stopifnot(ncol(H)==nrow(H))
    stopifnot(dim(G)==dim(H))
    stopifnot(dim(K)==dim(H))
    stopifnot(all.equal(t(G), G))
    stopifnot(all.equal(t(K), K))
    eK <- eigen(K, only.values=TRUE)$values
    stopifnot(all(eK>0) &&  all(Im(eK)==0))
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
    feval <- aaply(x, 1, fun, ...)*prodDx
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

.calcNormalizing <- function(statevec, H, K){
    .checkArgs(G=diag(1, nrow(H)), H=H, K=K)
    p <- ncol(H)
    if(missing(statevec)){
        nstate <- 2^p
        statevec <- seq(from=1, to=nstate-1)
    } else{
        nstate <- length(statevec)
    }
    if(nstate>=1024) stop('Too hazardous')
    N <- rep(NA, nstate)
    #N[1] <- 1
    state <- matrix(FALSE, nrow=nstate, ncol=ncol(H))
    for(i in seq_len(nstate)){
        v <- statevec[i]
        bits <- as.logical(intToBits(v)[1:p])
        state[i,] <- bits
        subK <- K[bits,bits,drop=FALSE]
        subH <- t(crossprod(bits, H))[bits,,drop=F]      #H order is wrong throughout, since 1y is supposed to be a selector, eg eq 6?  Can we fix by inverting 6?
        N[i] <- det(subK/(2*pi))^(-.5)*exp(.5*t(subH) %*% solve(subK) %*% subH)
    }
    list(N=N, state=state)
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
    if(is.matrix(x)) aaply(x, 1, .rCondHurdle, j=j-1, G=G, H=H, K=K, tol=tol) else .rCondHurdle(x, j-1, G, H, K, tol)
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

loopMLE <- function(samp){
    p <- ncol(samp)
    separm <- parm <- matrix(0, nrow=p, ncol=length(parmap(p)))
    colnames(separm) <- colnames(parm) <- parmap(p)
    parm[,'kbb'] <- 1
    for(j in 1:p){
        y <- samp[,j]
        x <- samp[,-j]
        ll <- generatelogLik(y, x)
        grad <- generatelogLik(y, x, returnGrad=TRUE)
        O <- optim(parm[j,], ll, grad, method='BFGS', hessian=TRUE)
        parm[j,] <- O$par
        separm[j,] <- try(sqrt(diag(solve(O$hessian))))
    }
    list(parm, separm)
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
