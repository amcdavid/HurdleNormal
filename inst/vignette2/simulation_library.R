library(HurdleNormal)
library(plyr)
library(data.table)
library(parallel)
library(huge)
library(reshape2)


##' \code{rbind} a list of data.tables adding column `L1` with names
##'
##' @param li list(possibly named)
##' @param nm optional vector of names
##' @param nm.col what should the column of names be called?
##' @param fill should columns be filled with NAs if column names don't match
##' @return data.table
namedrbindlist <- function(li, nm=NULL, nm.col='L1', fill=TRUE){
    if(!is.null(nm)) names(li) <- nm
    li <- lapply(names(li), function(linm){
        x <- as.data.table(li[[linm]])
        x[,(nm.col):=linm]
        x
               })
    rbindlist(li, fill=fill)

}

##' Plot a heatmap of the true adjacency matrix
##'
##' @param HS `HurdleStructure`
##' @return plots a heatmap of signed adjacency.  Sign is the sum of the weigh matrices (??)
##' @author Andrew
plotTrue <- function(HS){
    S <- with(HS$true, {
    sign(G) + sign(H) + sign(K)
    })
              heat2(S)
    }

##' Calculate true positives, false negatives and false positives
##'
##' @param slice adjacency matrix (possibly not symmetrized)
##' @param slicefun function to apply to adjacency matrix (ie, to symmetrize) 
##' @param trueAdj true adjacency matrix
##' @return named numeric vector with tp, fp, fn
getRoc <- function(slice, slicefun=function(sl) (sl & t(sl)), trueAdj){
    sl <- as.matrix(slicefun(slice))
    diag(sl) <- FALSE
    sl1 <- sl*1
    true1 <- trueAdj*1
    bitflag <- true1 -sl1 + 2*(true1* sl1)
#    heat2(bitflag)
    sl <- onlyTri(sl)
    trueAdj <- onlyTri(trueAdj)
    tp <- sum(sl & trueAdj)
    fp <- sum(sl & !trueAdj)
    fn <- sum(!sl & trueAdj)
    tn <- sum(!sl & !trueAdj)
    c(tp=tp, fp=fp, fn=fn)
}

## Contaminate positive values with t distributed noise scaled to match the standard deviation of each coordinate
contaminate <- function(x, df=8, tol=0){
    for(i in seq_len(ncol(x))){
        xi = x[,i]
        xpos <- abs(xi)>tol
        sx <- sd(xi[xpos])
        x[xpos,i] <- xi[xpos] + rt(sum(xpos), df=df)*sx
    }
    x
}

crunch <- function(x, v){
    pmin(abs(x), abs(v)) * sign(v)
}

ecoli_genenetweaver <- function(N_vertex=Inf, vertex_seq, Hdiag=0, Gdiag=0, Kdiag=1, Klim=.4, Glim=1){
    signed_nel <- fread(system.file('vignette2', 'Ecoli_goldstandard_signed.tsv', package='HurdleNormal'), header=FALSE)
    signed_nel[, weight:=ifelse(V3=='+', 1, -1)]
    g  <- graph_from_edgelist(el=as.matrix(signed_nel[,.(V1, V2)]), directed=FALSE)
    E(g)$weights <- signed_nel$weight
    if(is.finite(N_vertex) || !missing(vertex_seq)){
        if(missing(vertex_seq)) vertex_seq <- sample(V(g), N_vertex)
        ss <- setdiff(V(g), vertex_seq)
        g <- delete_vertices(g, ss)
    }
    g <- simplify(g, remove.multiple = FALSE, remove.loops = TRUE)
    el <- as.data.table(as_edgelist(g, names=FALSE))
    el$weights <- E(g)$weights
    ## Set edge weights
    el$G <- crunch(rnorm(nrow(el))+Glim/3, el$weights*Glim)
    el$K <- crunch(rnorm(nrow(el))/2+Klim/3, el$weights*Klim)
    el[runif(nrow(el))<.5, K:=0]
    el[runif(nrow(el))<.25 & abs(K)>0, G:=0]
    
    ## ensure PD precision and reweight normalization constant
    el2 <- copy(el)
    el <- rbind(el, el2[,':='(V1=V2, V2=V1)])
    Kdiag <- el[,.(K=sum(abs(K)), G=sum(G)-1),keyby=V1]
    Kdiag[,K:=pmax(1, K)]
    Kdiag[, V2:=V1]
    el$Hoff <- 0#2*el$K
    el <- rbind.fill(el, Kdiag)
    dn = list(names(V(g)), names(V(g)))
    Gmat <- sparseMatrix(i=el$V1, j=el$V2, x=el$G, dims=c(length(V(g)), length(V(g))), dimnames=dn)
    Kmat <- sparseMatrix(i=el$V1, j=el$V2, x=el$K, dims=c(length(V(g)), length(V(g))), dimnames=dn)
    Hmat <- sparseMatrix(i=el$V1, j=el$V2, x=el$Hoff, dims=c(length(V(g)), length(V(g))), dimnames=dn)
    Hmat <- Hmat + t(Hmat)
    #diag(Gmat) <- Gdiag
    diag(Kmat) <- pmax(1, diag(Kmat))
    diag(Hmat) <- Hdiag
    
    hs = HurdleStructure(G=as.matrix(Gmat), K=as.matrix(Kmat), H=as.matrix(Hmat), gibbs=FALSE)
    hs = getGibbs(hs, Nt=4e4)
    
    ## Shift marginal data
    x = hs$gibbs
    for(i in seq_len(ncol(x))){
       xi = x[,i]
       xpos <- abs(xi)>0
       x[xpos,i] <- xi[xpos] + rnorm(1, mean=4, sd=1)
       }
    hs$gibbs <- addPseudoCounts(x)
    colnames(hs$gibbs) = dn[[1]]
    hist(colMeans(hs$gibbs))
    hist(colMeans(abs(hs$gibbs)>0))
    
    ## Fit to get estimates after shifting data
    hh = fitHurdle(hs$gibbs, parallel=TRUE, nlambda=25, lambda.min=.2, makeModelArgs = list(conditionalCenter=FALSE, center=TRUE, scale=TRUE), control=list(tol=1e-2, newton0=TRUE, refit=FALSE), keepNodePaths = TRUE)
    
    result <- attr(hh, 'nodePaths')
    lall <- list(G=neighborhoodToArray(result, summaryFun=summaryG, nobs=nrow(hs$gibbs), self_edges=TRUE),
                 H=neighborhoodToArray(result, summaryFun=summaryHij,nobs=nrow(hs$gibbs), self_edges=TRUE))
    for(par in names(lall)){ 
        diag(hs$true[[par]]) <- diag(lall[[par]][[1]][[1]])
    }
    hs = getGibbs(hs, Nt=1e2)
    hs
}


##' Simulate new data from a model and fit via various algorithms
##'
##' 1. Sample a data set from model
##' 2. Fit various hurdle and other methods
##' 3. Calculate true positives, false positives, false negatives as a function of the number of edges
##' 4. Calculate the maximum sensitivity (TP/Positives) at 10% FDR
##' 5. Find the range of FP and interpolate TP over this set for each method
##' 6. Glue together as data.table.  `L1` gives method.
##' @param model `HurdleStructure`
##' @param n number of data points to simulate
##' @param makeModelArgs arguments passed  to `makeModel`
##' @param maxFDR FDR threshold at which to report sensitivity
##' @return data.table
fitAndSummarize <- function(model, n=1000, makeModelArgs, maxFDR=.1, parallel=FALSE){
    model <- getGibbs(model, n*20, thin=.1)
    P <- ncol(model$gibbs)
    true <- (abs(model$true$K)+abs(model$true$G)+abs(model$true$H))>0
    diag(true) <- FALSE
    mgibbs <- apply(model$gibbs, 2, function(col){
        col[abs(col)>.01] <-  scale(col[abs(col)>.01], scale=FALSE)
        col
    })
    
    ## Common parameters
    tol <- 5e-3
    newton0 <- TRUE
    nlambda <- 50
    lambda.min.ratio <- .4
    hfitReg <- fitHurdle(model$gibbs, lambda.min.ratio=lambda.min.ratio, nlambda=nlambda, makeModelArgs=list(center=TRUE, scale=FALSE, conditionalCenter=FALSE, newton0=newton0, tol=tol), parallel=parallel, penaltyFactor='identity')
    hfitCreg <- fitHurdle(mgibbs, lambda.min.ratio=lambda.min.ratio, nlambda=nlambda, makeModelArgs=list(center=TRUE, scale=FALSE, conditionalCenter=FALSE, newton0=newton0, tol=tol), parallel=parallel, penaltyFactor='identity')
    #hfitCCreg <- fitHurdle(model$gibbs, lambda.min.ratio=lambda.min.ratio, nlambda=nlambda, makeModelArgs=list(center=TRUE, scale=FALSE, conditionalCenter=TRUE, newton0=newton0, tol=tol), parallel=parallel, penaltyFactor='identity')
    hfitCfull <- fitHurdle(mgibbs, lambda.min.ratio=lambda.min.ratio, nlambda=nlambda, makeModelArgs=list(center=TRUE, scale=FALSE , conditionalCenter=FALSE, newton0=newton0, tol=tol), parallel=parallel, penaltyFactor='full')
    #hfitCCfull <- fitHurdle(model$gibbs, lambda.min.ratio=lambda.min.ratio, nlambda=nlambda, makeModelArgs=list(center=TRUE, scale=FALSE, conditionalCenter=TRUE, newton0=newton0, tol=tol), parallel=parallel, penaltyFactor='full')
    hfitFull <- fitHurdle(model$gibbs, lambda.min.ratio=lambda.min.ratio, nlambda=nlambda, makeModelArgs=list(center=TRUE, scale=FALSE , conditionalCenter=FALSE, newton0=newton0, tol=tol), parallel=parallel, penaltyFactor='full')
    glasso = huge(model$gibbs,method = "mb", nlambda=100,lambda.min.ratio = 0.01)
    S.npn = huge.npn(model$gibbs, npn.func="truncation")
    npn <- huge(S.npn,method = "mb", nlambda=100,lambda.min.ratio = 0.01)
    logiArr <- autoLogistic(model$gibbs, lambda.min.ratio=.02, nlambda=100)
    adjGlasso <- laply(glasso$path, getRoc, trueAdj=true)
    adjNpn <- laply(npn$path, getRoc, trueAdj=true)
    adjLogi <- laply(logiArr$adjMat, getRoc, trueAdj=true)
    #hurdlelist <- list(reg=hfitReg,Creg=hfitCreg, Cfull=hfitCfull, full=hfitFull, CCfull=hfitCCfull, CCreg=hfitCCreg)
    hurdlelist <- list(reg=hfitReg,Cfull=hfitCfull,  Creg=hfitCreg, full=hfitFull)
    hurdles <- lapply(hurdlelist, function(hfit){
        rocOr <- laply(hfit$adjMat, getRoc, slicefun=function(sl) sl | Matrix::t(sl), trueAdj=true)
        rocAnd <- laply(hfit$adjMat, getRoc, trueAdj=true)
        cbind(rocOr, timing=attr(hfit, 'timing')['elapsed'])
    })

    Mmetric <- namedrbindlist(c(list('Gaussian'=adjGlasso, 'NPN'=adjNpn, 'Logistic'=adjLogi), hurdles), fill=TRUE)
    Mmetric[, type:= c('Hurdle "Or"'='Hurdle',  'Hurdle "And"'='Hurdle', 'Gaussian'='Gaussian', 'NPN'='Gaussian', 'Logistic'='Logistic')[L1]]

    Mmetric[,FDR:=pmax(fp/(tp+fp), 0, na.rm=T), keyby=L1]
    setkey(Mmetric, L1, FDR)
    supFDR <- Mmetric[FDR<maxFDR, .(sensitivity=tp[1]/(tp[1]+fn[1])), key=L1]
    
    fprange <- range(Mmetric$fp)
    fpgrid <- fprange[1]:fprange[2]
    fpInterpolate <- Mmetric[,{
                ap <- approx(fp, tp, fpgrid, rule=2)
                list(fpI=ap$x,tpI=ap$y, Method=L1[1])
            }, keyby=L1]
    rbind(supFDR, fpInterpolate, fill=TRUE)
}
