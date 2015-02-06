HurdleStructure <- function(G=NULL, H=NULL, K=NULL, sample=NULL, gibbs=NULL, estimate=NULL){
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

    li <- list(true=list(G=G, H=H, K=K), sample=sample, gibbs=gibbs, estimate=estimate)
    class(li) <- c('list', 'HurdleStructure')
    li
}

getGibbs <- function(hs, Nt=2000, ...){
    if(!is.null(hs$gibbs)) message("Replacing gibbs sample")
    gibbs <- with(hs$true, rGibbsHurdle(G, H, K, Nt=Nt, ...))
    hs$gibbs <- gibbs
    hs
}

getDensity <- function(hs, ...){
    function(x) with(hs$true, dHurdle210(x, G, H, K, ...))
}

##' Loop over indices and solve condiiton MLE
##'
##' @param hs HurdleStructure
##' @param using 'gibbs' or 'sample'
##' @param ... passed to optim
##' @return list of length two giving parameter values and standard errors. j is the index of the response, i is the index of the coefficient.
##' @import reshape
##' @export
getConditionalMLE <- function(hs,using='gibbs', testGrad=FALSE, engine='R', ...){
    samp <-  if(using=='gibbs') hs$gibbs else hs$sample
    p <- ncol(samp)
    coefIndex <- separm <- parm <- matrix(0, nrow=p, ncol=length(parmap(p)), dimnames=list(1:p, parmap(p)))
    parm[,'kbb'] <- 1
    for(j in 1:p){
        y <- samp[,j]
        x <- samp[,-j,drop=FALSE]
        if(engine=='R'){
            ll <- generatelogLik(y, x, lambda=0)
            grad <- generatelogLik(y, x, lambda=0, returnGrad=TRUE)
        } else{
            hl <- HurdleLikelihood(y, x, theta=parm[j,], lambda=0)
            ll <- hl$LLall
            grad <- hl$gradAll
        }
        O <- hushWarning(optim(parm[j,], ll, method='BFGS', hessian=TRUE, ...), fixed('NaNs produced'))
        if(testGrad && !all(abs(grad(O$par))<.1)){
            stop('Gradient not zero at putative solution.  Max |grad| =  ', max(abs(grad(O$par))))
        }

        ## O2 <- optim(parm[j,], ll, gr=grad, method='BFGS', hessian=TRUE)
        ## hl$LLall(O1$par)
        ## hl$grad(O2$par[c(1, 4, 11)], grp=-1)
        ## hl$gradAll(O1$par)
        ## parm[j,] <- O$par
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

getJoint <- function(fit){
    C <- lapply(cast(fit$Parm, i~j | par), function(x) as.matrix(x[,-1]))
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


simulateHurdle210 <- function(N, p, dependence='G', structure='independence', structureArgs=list(sparsity=.1, groupwise=FALSE), intensityArgs=list(G=1, Hupper=1, Hlower=1, K=-.3), Hdiag=5, Kdiag=1, G0odds=0){
    dependence <- match.arg(dependence, c('G', 'Hupper', 'Hlower', 'K'), several.ok=TRUE)
    structure <- match.arg(structure, c('independence', 'sparse', 'chain'))
    Hlower <- Hupper <- diag(Hdiag, p)/2
    K <- diag(Kdiag, p)/2
    G <- matrix(0, p, p)

    ## loop over types of dependence
    for(d in dependence){
        ## empty matrix to be filled with dependences and added to diagonal matrix
        depMat <- matrix(0, p, p)
        ## get a offdiagonal location and fill with a nonzero entry (using intensityArgs)
        if(structure=='sparse'){
            offdiag <- p*(p-1)/2
            nnz <- ceiling(structureArgs$sparsity*offdiag)
            depidx <- if(nnz>1) sample(offdiag, nnz) else 1 #let's use the 1,2 entry unless we need more than one dependence
            depMat[upper.tri(depMat)][depidx] <- intensityArgs[[d]]
        } else if(structure=='chain'){
            depidx <- cbind(i=seq_len(p-1), j=2:p)
            depMat[depidx] <- intensityArgs[[d]]
        }
        ## else independence, and we just add an empty matrix.
        ourMat <- get(d)+depMat
        assign(d, ourMat)
    }
    
    K <- K+t(K)
    H <- Hupper + t(Hlower)
    G <- G+t(G)
    addIndices <- .additiveVstate(ncol(H))
    norm <- .calcNormalizing(addIndices$states, H=H, K=K)
    browser()
    G[addIndices$indices] <- -log(norm$N)

    Gprime <- matrix(c(-13.4, (-37.5 + 26.8)/2, 0,
                       (-37.5 + 26.8)/2, -13.4, 0,
                       0, 0, -13.4), nrow=3, ncol=3)
    samp <- rGibbsHurdle(Gprime, H, K, Nt=20*N, thin=.1)
    pairs(samp)

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
