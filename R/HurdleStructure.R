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
