HurdleStructure <- function(G=NULL, H=NULL, K=NULL, sample=NULL, gibbs=NULL, estimate=NULL){
    if(!is.null(G)){
        .checkArgs(G, H, K)
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
    class(li, 'HurdleStructure')
}

getGibbs <- function(hs, Nt=2000, ...){
    if(!is.null(hs$gibbs)) message("Replacing gibbs sample")
    gibbs <- with(hs$true, rGibbsHurdle(G, H, K, Nt=Nt, ...))
    hs$gibbs <- gibbs
}

simulateHurdle210 <- function(N, p, dependence='G', structure='independence', structureArgs=list(sparsity=.1, groupwise=FALSE), intensityArgs=list(G=1, Hupper=1, Hlower=1, K=-.3), Hdiag=5, Kdiag=1, GmarginalProb=.5, Gdiag){
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
    ## now try to guess at the marginal distribution by replacing xI and x
    ## in the conditional distributions with their expectations given independence
    ex <- diag(H)/diag(K)
    exI <- rep(.75, p)
    K0 <- K
    Hlower0 <-Hlower
    Hupper0<- Hupper
    diag(K0) <- diag(Hlower0) <- diag(Hupper0) <- 0
    
    guessMarginals <- function(Gdiag){
        G0 <- G
        diag(G0) <- 0
        for(i in seq_len(10)){
            gbaNoDiag <- crossprod(G0, exI)+ crossprod(Hlower0, ex*exI)
            hba <- diag(H) + crossprod(Hupper0, exI) + #only H's diagonal
                -crossprod(K0, ex*exI)/2 # not K's diagonal
            logitP <- gbaNoDiag-.5*log(diag(K)/(2*pi))+hba^2/(2*diag(K))
            ex <- hba/diag(K)
            Gdiag <- pmax(pmin(logit(GmarginalProb)-logitP, 1000), -1000)
            exI <- expit(logitP+Gdiag)
            print(paste0('ex=', ex))
            print(paste0('exI=',t(exI)))
            print(paste0('Gdiag=', t(Gdiag)))
        }
        if(any(abs(exI-GmarginalProb) > .05)) warning("Wasn't able to match marginal probabilities")
        Gdiag
    }
    if(missing(Gdiag)){
        diag(G) <- guessMarginals(diag(G))
    } else{
        diag(G) <- Gdiag
    }
    browser()
    samp <- rGibbsHurdle(G, H, K, Nt=20*N, thin=.1)
    pairs(samp)

}


Kdep <- simulateHurdle210(200, 3, 'K', 'sparse', Gdiag=c(-19, -19, -13.41))
Gdep <- simulateHurdle210(200, 3, 'G', 'sparse')
Hupperdep <- simulateHurdle210(2000, 3, 'Hupper', 'sparse', intensityArgs=list(Hupper=1, Hlower=1))
g12 <- glm(samp[,1] >0 ~ samp[,2]+I(samp[,2]>0), family='binomial')
b12 <- lm(samp[,1] ~ samp[,2] + I(samp[,2]>0), subset=samp[,1]>0)
g21 <- glm(samp[,2] >0 ~ samp[,1]+I(samp[,1]>0), family='binomial') #discrete significant
b21 <- lm(samp[,2] ~ samp[,1] + I(samp[,1]>0), subset=samp[,2]>0) 

Hlowerdep <- simulateHurdle210(2000, 3, 'Hlower', 'sparse', intensityArgs=list(Hupper=1, Hlower=1))
g12 <- glm(samp[,1] >0 ~ samp[,2]+I(samp[,2]>0), family='binomial') 
b12 <- lm(samp[,1] ~ samp[,2] + I(samp[,2]>0), subset=samp[,1]>0) #discrete significant
g21 <- glm(samp[,2] >0 ~ samp[,1]+I(samp[,1]>0), family='binomial') #continuous significant
b21 <- lm(samp[,2] ~ samp[,1] + I(samp[,1]>0), subset=samp[,2]>0) 
