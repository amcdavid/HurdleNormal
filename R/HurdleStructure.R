HurdleStructure <- function(G=NULL, H=NULL, K=NULL, sample=NULL, gibbs=TRUE, estimate=NULL){
    if(!is.null(G)){
        .checkArgs(G=G, H=H, K=K)
    } else if(!is.null(sample) || !is.null(gibbs)){
        p <- if(!is.null(sample)) ncol(sample) else ncol(gibbs)
        G <- H <- K <- matrix(NA, nrow=ncol(sample), ncol=ncol(sample))
    } else if(!is.null(estimate)){
        G <- estimate$G
        K <- estimate$K
        H <- estimate$H
    } else{
        stop('Empty constructor')
    }
    true <- list(G=G, H=H, K=K)
    if(isTRUE(gibbs)) gibbs <- getGibbs(list(true=true))$gibbs

    trueFactor <- sign(K)+2^sign(H)+4^sign(G)
    Kt <- signToChar(K, 'K')
    Ht <- signToChar(H, 'H')
    Gt <- signToChar(G, 'G')
    trueChar <- paste0(Kt, Ht, Gt)
    trueFactor <- matrix(trueFactor, nrow=nrow(G))
    trueChar <- matrix(trueChar, nrow=nrow(G))

    li <- list(true=true, sample=sample, gibbs=gibbs, estimate=estimate, norm=NULL, trueFactor=trueFactor, trueChar=trueChar)
    class(li) <- c('list', 'HurdleStructure')
    li
}

signToChar <- function(v, char) {
    char <- c(paste0('-', char),
              '',
              paste0('+', char))
    char[sign(v)+2]
}

getNormalizing <- function(hs, subsample=Inf){
    ##sample(seq(from=1L, to=nstate-1L), subsample) else
    
    norm <- .calcNormalizing(G=hs$true$G, H=hs$true$H, K=hs$true$K, log=TRUE)
    maxP <- max(norm$logP)
    state0 <- 0-maxP
    logP0 <- -(log(sum(exp(norm$logP-maxP))+
                    exp(state0)          #state 0
                  )+maxP)
    P <- exp(norm$logP +logP0)
    margins <- colSums(t(norm$state)*P)
    list(logP0=logP0, P=c(exp(logP0), P), margins=margins, lMargin=colSums(t(norm$stat)*norm$logP), mu=norm$mu)
}

getGibbs <- function(hs, Nt=2000, ...){
    if(!is.null(hs$gibbs)) message("Replacing gibbs sample")
    gibbs <- with(hs$true, rGibbsHurdle(G, H, K, Nt=Nt, ...))
    ct <- cor.test(gibbs[-1,1], gibbs[-nrow(gibbs),
                                  1], conf.level=.99)$conf.int[1]
    if(!is.na(ct) && ct >.1) warning('Significant auto correlation found')
    hs$gibbs <- gibbs
    
    hs
}

getDensity <- function(hs, ...){
    function(x) with(hs$true, dHurdle210(x, G, H, K, ...))
}



##' Convert from old style to new style parameter vector or vice versa
##'
##' This just used to bridge between the R and C++ likelihood calcs, hence is legacy code now.
##' @param theta parameter vector
##' @param method \code{character} one of 'unpack', 'pack', or 'gradient'
##' @return ordered parameter
unpackTheta <- function(theta, method='unpack'){
    p <- (length(theta)-3)/4+1
    component <- c(rep(0, p), #gbb, gba
                   rep(2, p), #hbb, hba
                   rep(1, p-1),#hab
                   rep(3, p-1),#kba
                   4) #kbb
    scaling <- c(1,   #gbb
                 rep(2, p-1),           #gba
                 rep(1, 2*p-1),             #hbb, hba, hab,
                 rep(-1, p-1),              #kba
                 1)                  #kbb
    grporder <- order(component, coordmap(p))
    method <- match.arg(method, c('unpack', 'pack', 'gradient'))
    if(method=='unpack'){
        return(theta[grporder]*scaling[grporder])
    } else if(method=='pack') {
        return(theta[order(grporder)]/scaling[order(grporder)])
    } else if(method=='gradient'){
        return(theta[order(grporder)]*scaling[order(grporder)]) #apply jacobian
    }
    
}



wrapHurdleLikelihood <- function(y, x, theta, lambda){
    x <- cbind(1, abs(x)>0, x)
    if(!missing(theta)){
        theta <- unpackTheta(theta)
        HurdleLikelihood(y, x, ,theta, lambda)
    } else{
           HurdleLikelihood(y, x, lambda=lambda)
    }
}

wrapLLall <- function(hl, theta, penalize){
    theta <- unpackTheta(theta)
    hl$LLall(theta, penalize)
}

wrapGradAll <- function(hl, theta, penalize){
    unpackTheta(hl$gradAll(unpackTheta(theta), penalize), 'gradient')
}


##' Loop over indices and solve condiiton MLE
##'
##' @param hs HurdleStructure
##' @param using 'gibbs' or 'sample'
##' @param testGrad Should we check the gradient for convergence.
##' @param engine \code{character}.  If `R` then the R likelihood/gradient will be used. Else C++.
##' @param ... passed to optim
##' @return list of length two giving parameter values and standard errors. j is the index of the response, i is the index of the coefficient.
##' @export
getConditionalMLE <- function(hs,using='gibbs', testGrad=FALSE, engine='R', ...){
    samp <-  if(using=='gibbs') hs$gibbs else hs$sample
    p <- ncol(samp)
    coefIndex <- separm <- parm <- matrix(0, nrow=p, ncol=length(parmap(p)), dimnames=list(1:p, parmap(p)))
    parm[,'kbb'] <- 1
    for(j in 1:p){
        y <- samp[,j]
        x <- samp[,-j,drop=FALSE]
        hl <- wrapHurdleLikelihood(y, x, theta=parm[j,], lambda=0)
        if(engine=='R'){
            ll <- generatelogLik(y, x, lambda=0)
            grad <- generatelogLik(y, x, lambda=0, returnGrad=TRUE)
        } else{
            ll <- function(x) wrapLLall(hl, theta=x, penalize=FALSE)
            grad <- function(x) wrapGradAll(hl, theta=x, penalize=FALSE)
        }
        O <- try(hushWarning(optim(parm[j,], ll, method='BFGS', hessian=TRUE, ...), fixed('NaNs produced')))
        if(inherits(O, 'try-error'))
        if(testGrad && !all(abs(grad(O$par))<.1)){
            stop('Gradient not zero at putative solution.  Max |grad| =  ', max(abs(grad(O$par))))
        }

        ## O2 <- optim(parm[j,], ll, gr=grad, method='BFGS', hessian=TRUE)
        ## hl$LLall(O1$par)
        ## hl$grad(O2$par[c(1, 4, 11)], grp=-1)
        ## hl$gradAll(O1$par)
        parm[j,] <- O$par
        try(separm[j,] <- sqrt(diag(solve(O$hessian*nrow(samp)))), silent=TRUE) #gradient is scaled by 1/N
        ## what index in data do these components refer to?
        ## remap coordmap by current permutation of indices
        coefIndex[j,] <- c(j, setdiff(1:p, j))[coordmap(p)]
    }
    Parm <- melt(parm)
    Separm <- melt(separm)
    names(Parm) <- names(Separm) <-  c('j', 'par', 'value')
    Parm <- cbind(Parm, i=as.numeric(coefIndex))
    Separm <- cbind(Separm, i=as.numeric(coefIndex))
    list(Parm=Parm, Separm=Separm)
}

momentChar <- function(hs, v){
    norm <- .calcNormalizing(v=v, H=hs$true$H, K=hs$true$K, log=TRUE)
    logP <- sum(hs$true$G[v,v]) +  norm$N
    cov <- solve(norm$subK)
    mu <- norm$mu
    list(logP=logP, cov=cov, mu=mu)
    
}

getJoint <- function(fit){
    C <- lapply(reshape::cast(fit$Parm, i~j | par), function(x) as.matrix(x[,-1]))
    G <- C$gba
    diag(G) <- diag(C$gbb)
    H <- C$hba + t(C$hab)
    diag(H) <- diag(C$hbb)
    K <- C$kba
    diag(K) <- diag(C$kbb)
    G <- (G+t(G))/2
    K <- (K+t(K))/2
    list(G=G, H=H, K=K)
}


## Used to simulate Erdos-Renyi networks
simulateHurdle210 <- function(N, p, dependence='G', structure='independence', structureArgs=list(sparsity=.1, groupwise=FALSE), intensityArgs=list(G=1, Hupper=1, Hlower=1, K=-.3, gamma=1), Gdiag=-14, Hdiag=5, Kdiag=1, tweak=NULL, tweakArgs=list(lH0=4, lG0 =2)){
    dependence <- match.arg(dependence, c('G', 'Hupper', 'Hlower', 'K'), several.ok=TRUE)
    structure <- match.arg(structure, c('independence', 'sparse', 'chain'))
    if(!is.null(tweak)) tweak <- match.arg(tweak, c("G", "GH"))
    Hlower <- Hupper <- diag(Hdiag, p)/2
    K <- diag(Kdiag, p)/2
    G <- matrix(0, p, p)

    ## loop over types of dependence
    for(d in dependence){
        ## empty matrix to be filled with dependences and added to diagonal matrix
        depMat <- matrix(0, p, p)
        ## get a offdiagonal location and fill with a nonzero entry (using intensityArgs)
        offdiag <- p*(p-1)/2
        nnz <- ceiling(structureArgs$sparsity*offdiag)
        if(structure=='sparse'){
            depidx <- if(nnz>1) sample(offdiag, nnz) else 1 #let's use the 1,2 entry unless we need more than one dependence
            depMat[upper.tri(depMat)][depidx] <- intensityArgs[[d]]*sign(intensityArgs[['gamma']]-runif(length(depidx)))
        } else if(structure=='chain'){
            i <- rep(seq_len(p), length.out=nnz)
            depidx <- cbind(i=i, j=(i+ceiling(seq_len(nnz)/p)))
            depidx[depidx[,'j']>p, 'j'] <- 1
            depMat[depidx] <- intensityArgs[[d]]
        }
        ## else independence, and we just add an empty matrix.
        ourMat <- get(d)+depMat
        assign(d, ourMat)
    }
    K <- K+t(K)
    eK <- min(eigen(K, only.values=TRUE)$values-.1, 0)
    diag(K) <- diag(K)-diag(K)*eK #Let's keep it PD folks
    
    H <- Hupper + t(Hlower)
    G <- G+t(G)
    diag(G) <- Gdiag
    tweakmat <- G
    if(!is.null(tweak)){
        if(tweak=='GH') tweakmat <- K
        lH0 <- tweakArgs[['lH0']]
        lG0 <- tweakArgs[['lG0']]
    } else{
        lG0 <- lH0 <- 0
    }
    knz <- abs(tweakmat)>0
    knz <- knz & t(knz)
    diag(knz) <- FALSE
    anynz <- abs(K) |abs(G)>0 | abs(H)>0
    diag(anynz) <- FALSE

    tweakOffDiag <- function(HS, lH, lG){
        HS$true$H[knz] <- lH*tweakmat[knz]
        if(!is.null(tweak) && tweak=='GH'){
            HS$true$G[knz] <- lG*tweakmat[knz]
        }
        HS <- suppressWarnings(getGibbs(HS))
        margins <- round(colMeans(abs(HS$gibbs)>0), 2)
        message(paste(margins, collapse=','))
        G <- HS$true$G
        message('Before: ', paste(round(diag(G), 2), collapse=','))
        diag(G) <- diag(G)-log((margins+.01)/(1.01-margins))
        message('After: ', paste(round(diag(G),2), collapse=','))
        HS$true$G <- G
        HS
    }
    hs <- HurdleStructure(G, H, K, gibbs=FALSE)
    if(!is.null(tweak)){
        warning('Beware of crunchy code')
        hs <- tweakOffDiag(hs, lH0, lG0)
        jamPcor <- function(lH=lH0, lG=lG0){
            message('**lH= ', round(lH, 2), 'lG = ', round(lG, 2), '**')
            hs <- tweakOffDiag(hs, lH, lG)
            rho <- as.matrix(pcor.shrink(hs$gibbs, lambda=.15))
            diag(rho) <- 0
            obj <- sum( (2*anynz-1)*abs(rho))
            message('--obj=', obj, '--')
            obj
        }
        O2 <- optimize(jamPcor, c(-2, lH0+3), tol=.15, lG=lG0)
        lH0 <- O2$min
        if(tweak=='GH'){
            O <- optimize(jamPcor, c(-lG0-1, lG0), tol=.15, lH=lG0)
            lG0 <- O$min
        }
        message('lH=', lH0, 'lG= ', lG0)
        hs <- tweakOffDiag(hs, lH0, lG0)
    }
    hs <- getGibbs(hs, N/.1+2000, burnin=2000, thin=.1)
    margins <- colMeans(abs(hs$gibbs)>0)
    message(paste(margins, collapse=','))
    hs
}



##Kdep <- simulateHurdle210(200, 3, 'K', 'sparse', Gdiag=c(-19, -19, -13.41))
## Gdep <- simulateHurdle210(200, 3, 'G', 'sparse')
## Hupperdep <- simulateHurdle210(2000, 3, 'Hupper', 'sparse', intensityArgs=list(Hupper=1, Hlower=1))
## 
## 
##  #discrete significant
## 

## Hlowerdep <- simulateHurdle210(2000, 3, 'Hlower', 'sparse', intensityArgs=list(Hupper=1, Hlower=1))
## g12 <- glm(samp[,1] >0 ~ samp[,2]+I(samp[,2]>0), family='binomial') 
## b12 <- lm(samp[,1] ~ samp[,2] + I(samp[,2]>0), subset=samp[,1]>0) #discrete significant
## g21 <- glm(samp[,2] >0 ~ samp[,1]+I(samp[,1]>0), family='binomial') #continuous significant
## b21 <- lm(samp[,2] ~ samp[,1] + I(samp[,1]>0), subset=samp[,2]>0) 

blockHS <- function(mat, times){
    m2 <- mat
    diag(m2) <- 0
    htile <- do.call(cbind, rep(list(m2), times=times))
    vtile <- do.call(rbind, rep(list(htile), times=times))
    diag(vtile) <- rep(diag(mat), times)
    vtile
}

hcenter <- function(samp){
   apply(samp, 2, function(x){
             xI <- abs(x)>0
             x[xI] <- scale(x[xI], scale=FALSE)
             x
         })
}
