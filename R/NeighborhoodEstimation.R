##' Define parameter groups
##'
##' Define parameter groups and relative penalization.
##' If \code{this.model} is specified, then, neither \code{blist} nor \code{nlist} should be provided.
##' If \code{this.model} is not specified then both of the preceeding must be provided.
##' When \code{this.model} is provided, it is assumed to have been created by with \code{makeModel}.
##'  Column indices indicated in the attribute `fixedCols` are unpenalized intercepts, then followed by pairs of interaction terms.
##' The grouping corresponds to the parameter set for a node.
##' Otherwise, blist, mlist and nlist should be set.  See details.
##' @details
##' There are four components that all need to be mapped between each other.
##' In increasing abstraction, with variable prefixes in parenthesis:
##' 1. (p)arameter vector. The map is given in these terms. 
##' 2. (mm) model matrix--columns from the covariate matrix.
##' 3. (b)locks -- penalty groups
##' 4. (n)odes. -- the graph-theoretic structure.
##' 5. lambda -- penalties as a function of blocks.
##' All of these components are provided in the `map`
##' @param this.model (optional) a matrix created by \code{\link{makeModel}}
##' @param blist a list of parameter indices, one per block.  By default the first block is assumed to be unpenalized.
##' @param mlist a list of parameter indices, one per column of the model.matrix. If omitted, assumed to equal to the identity.
##' @param nlist a named list of block indices, one per node
##' @param lambda currently ignored
##' @param group \code{character}: one of components, or none.
##' @param penalty.scale optional list containing elements `scale` and `group`.
##' `group` should be one of 'block' or 'none'.  `scale` should be \code{numeric} of length `blist` or the sum of the `blist` lengths.
##' @return a list containing a data.table `map` giving the mapping between parameters, groups and penalty scales and some other components
##' @export
Block <- function(this.model, blist, mlist, nlist, group='components', lambda, penalty.scale=NULL){
    ## only one of this.model, bidx and bset should be non-missing
    ## group={components, none}
    ## penalty scale should be provided as a list of indices and group={groups, components,none}


    group <- match.arg(group, c('components', 'none'))
    if(!is.null(penalty.scale)){
        if(!is.list(penalty.scale)) stop("`penalty.scale` should be list with components `scale` and `group`")
        penalty.scale.group <- match.arg(penalty.scale$group, c('components', 'none', 'groups'))
        if(!is.numeric(penalty.scale$scale)) stop("`penalty.scale$scale` should give integer indices into groups or parameters")
        penalty.scale.scale <- penalty.scale$scale
    }
    ## this.model provided
    if(!missing(this.model)){
        nfixed <- length(attr(this.model, 'fixedCols'))
        nodeId <- c('(fixed)', attr(this.model, 'nodeId'))

        ## used to assign nodeIds and blocks when group=='components'
        p <- ncol(this.model)
        nc <- (p-nfixed)/2 #penalized blocks
        ## block index into model matrix
        bidxMM <- c(rep(1, nfixed), (rep(1:nc, each=2)
                        +1)) #offset from non-pen block
        ## node index into parameter vector
        nidxPar <- c(bidxMM, bidxMM, 1) #variance

        ## group lasso
        if(group=='components'){
            ## block index into parameter
            bidxPar <- nidxPar
            ## node -> block map
        } else{
            nc <- (p-nfixed)               #regular lasso
            bidxMM <- c(rep(1, nfixed), ((1:nc)+1))
            bidxMM2 <- c(rep(1, nfixed), ((1:nc)+nc+1))
            bidxPar <- c(bidxMM, bidxMM2, 1)
        }
        nodeMap <- data.table(nodeId, nidx=seq_along(nodeId))
        map <- data.table(paridx=seq_len(2*p+1), block=bidxPar, mmidx=c(1:length(bidxMM), 1:length(bidxMM), NA), nidx=nidxPar)
        map <- merge(map, nodeMap, by='nidx')
    }
    ## end this.model provided
    else{
        bvec <- data.table(paridx=unlist(blist), block=rep(seq_along(blist), times=sapply(blist, length)))
        if(any(duplicated(bvec$paridx))) stop('blist had duplicated parameters!')
        if(is.list(nlist)){
            nvec <- setNames(melt(nlist), c('block', 'nodeId'))
        } else{
            nvec <- data.table(block=seq_along(nlist), nodeId=nlist)
        }
        if(any(duplicated(nvec$block))) stop('nlist had duplicated blocks!')
        if(missing(mlist)) mlist <- as.list(bvec$paridx)
        mvec <- data.table(paridx=unlist(mlist), mmidx=rep(seq_along(mlist), times=sapply(mlist, length)))
        if(any(duplicated(mvec$mmidx))) stop('mlist had duplicated mmidx!')
        map <- merge(bvec, nvec, by='block')
        map <- merge(map, mvec, by='paridx')        
    }
    
    if(is.null(penalty.scale)){
        blocks <- sort(unique(map$block))
        penalty.map <- data.table(block=blocks, lambda=c(0, rep(1, length(blocks)-1)))
        penalty.scale.group <- 'groups'
    }    
    if(penalty.scale.group=='groups'){
        map <- merge(map, penalty.map, by='block')
    } else if(penalty.scale.group=='components'){
        stop('Not implemented')
    } else{
        penalty.map <- data.table(paridx=seq_along(penalty.scale.scale), lambda=penalty.scale.lambda)
        map <- merge(map, penalty.map, by='paridx')
    }

    setkey(map, paridx)
    out <- list(map=map, nonpenpar=map[lambda<=0, paridx], nonpenMM=unique(map[lambda<=0 & !is.na(mmidx), mmidx]), nonpengrp=unique(map[lambda<=0, block]), lambdablock=map[,list(lambda=lambda[1]), keyby=block]$lambda)
    class(out) <-  'Block'
    out
}

globalVariables(c('penalty.scale.lambda', 'paridx', 'mmidx', 'block'))

##' @export
##' @import data.table
##' @import Matrix
##' @describeIn fitHurdle Fit an auto-model (Ising or Gaussian) to \code{samp} using glmnet.  checkpointDir is currently ignored.
##' @param nlambda number of lambda values on grid (default 200).
##' @param lambda.min.ratio minimum lambda ratio (as a function of lambda0, where the first predictor enters; default .1)
##' @param family in the case of \code{autoGLM} one of "gaussian" or "logistic"
##' @aliases autoLogistic
autoGLM <- function(samp, fixed=NULL, parallel=FALSE, keepNodePaths=FALSE, checkpointDir = 'ignored', nlambda=200, lambda.min.ratio=.1, family='binomial'){
    if(is.null(colnames(samp))) colnames(samp) <- seq_len(ncol(samp))
    colnames(samp) = make_unique(colnames(samp))
    samp0 <- if(family=='binomial') (abs(samp)>0)*1 else samp
    applyfun <- if(parallel) function(X, FUN) parallel::mclapply(X, FUN, mc.preschedule=TRUE) else lapply
    if(is.null(fixed)) fixed <- matrix(1, nrow=nrow(samp0))
    if(any(fixed[,1] != 1)) stop('Column 1 of `fixed` covariates must be intercept!')
    nid <- colnames(samp)
    blist <- singletonMap(ncol(samp)+ncol(fixed), ncol(fixed))
    
    if(any(constant <- apply(samp, 2, var)<.Machine$double.eps)) warning(sum(constant), ' constant gene(s) found.  Dangerous ground; this may trigger an error in future versions.')
    
    timing <- system.time(result <- applyfun(seq_len(ncol(samp)), function(i){
        model <- cbind(fixed, samp0[,-i])
        thisId <- nid[i]
        blk <- Block(blist=blist, nlist=c('(fixed)', setdiff(nid, thisId)))
        posobs <- sum(samp0[,i]>0)
        
        if( posobs > 2 && (family=='gaussian' | (nrow(samp0)-posobs)>2)){
            net <- glmnet::glmnet(model[,-1], #glmnet has a bug in which it always include the intercept `column` in its internal design, hence supplying our own
                                  ## and setting penalty.factor accordingly fails.
                                  samp0[,i], family=family,lambda.min.ratio=lambda.min.ratio, nlambda=nlambda, penalty.factor=blk$map$lambda[-1], standardize=FALSE)            
            path <- Matrix::t(coef(net))
        } else{
            net <- list(df=c(1), lambda=c(0))
            path <- Matrix::Matrix( c(1, rep(0, ncol(model)-1)), nrow=1, sparse=TRUE)
        }
        rownames(path) <- net$lambda
        refit <- refitGLMVector(path, model, samp0[,i], blk, family)
        res <- list(path=path, blocks=blk, lambda=net$lambda, df=net$df, nodeId=thisId, path_np=refit$path_np, loglik_np=refit$loglik_np, nobs=nrow(samp))
        class(res) <- 'SolPath'
        res
    }))
    arr <- neighborhoodToArray(result, vnames=colnames(samp), nobs=nrow(samp))
    if(keepNodePaths){
        return(structure(arr, timing=timing, nodePaths=result))
    }
    return(structure(arr, timing=timing))
}

##Legacy synonym
autoLogistic <- autoGLM

## Refit the model with an unpenalized regression on the appropriate subsets
refitGLMVector <- function(path, this.model, y.zif, blocks, family, fuzz=0){
    path_np = path
    loglik_np = rep(NA_real_, nrow(path))
    for(i in seq_len(nrow(path))){
        activetheta <- which(abs(path[i,])>fuzz)
        activemm <- unique(na.omit(blocks$map[list(paridx=activetheta),mmidx, on='paridx']))
        refit <- glm.fit(this.model[,activemm,drop=FALSE], y.zif, family=do.call(family, list()))
        path_np[i,activetheta] = refit$coef

        p <- refit$rank
        if (family %in% c("gaussian", "Gamma", "inverse.gaussian")) 
            p <- p + 1
        loglik_np[i] <- p - refit$aic/2
    }
    ## This should only arise when we fit a constant gene somehow
    ## Just zero out the log-likelihood
    if(any(!is.finite(loglik_np))) loglik_np[] = 0
    list(path_np=path_np, loglik_np=loglik_np)
}

accessDirOrDie <- function(dir){
    if(is.null(dir) || file_test("-d", dir)) return(TRUE)
    if(!dir.create(dir)) stop('Cannot create/access', dir)
    return(TRUE)
}

#' Center data so that it has conditional column mean zero
#'
#' @param samp matrix of data
#' @return centered matrix
#' @export
#' @examples
#' dat = data.frame(x=c(1:5, 0, 0, 0), y=c(6:10, 0, 0, 0))
#' conditionalCenter(dat)
#' #compare to
#' scale(dat, scale=FALSE)
conditionalCenter <- function(samp) {
    apply(samp, 2, function(x) {
        xI <- abs(x) > 0
        x[xI] <- scale(x[xI], scale = FALSE)
        x
    })
}

make_unique = function(x){
    uniq = make.unique(x)
    if( length(not_uniq <- which(x != uniq))>0) warning('Element(s) ', paste(x[not_uniq], collapse = ', '), ' have been renamed to be unique.  (Are there unintended duplicates?)')
    uniq
}

##' Fit the hurdle model coordinate-by-coordinate on a sample
##'
##' @param samp matrix of data, columns are variables
##' @param fixed data.frame of fixed covariates (to be conditioned upon)
##' @param parallel parallelize over variables using "mclapply"?
##' @param keepNodePaths return node-wise output (solution paths and diagnostics for each node) as attribute `nodePaths`
##' @param checkpointDir (optional) directory to save the fit of each node, useful for large problems.  If it exists, then completed nodes will be automatically loaded.
##' @param makeModelArgs (optional) arguments passed to the model matrix function
##' @param indices (optional) subset of indices to fit, useful for cluster parallelization.
##' @param ... passed to cgpaths
##' @return list of fits, one per coordinate and an attribute "timing"
##' @seealso neighborhoodToArray, autoGLM, interpolateEdges
##' @export
fitHurdle <- function(samp, fixed=NULL, parallel=TRUE, keepNodePaths=FALSE, checkpointDir=NULL, makeModelArgs=NULL,  indices, ...){
    applyfun <- if(parallel) function(X, FUN) parallel::mclapply(X, FUN, mc.preschedule=FALSE) else lapply
    allindices <- seq_len(ncol(samp))
    if(is.null(colnames(samp))) colnames(samp) <- seq_len(ncol(samp))
    colnames(samp) = make_unique(colnames(samp))
    indices <- if(missing(indices))  allindices else indices
    if(length(setdiff(indices, allindices))>0) stop('`indices` out of range')
    accessDirOrDie(checkpointDir)
    timing <- system.time(result <- applyfun(indices, function(i){
        message('node=', i, ' nodeId=', colnames(samp)[i])
        if(!is.null(checkpointDir) && file.exists(fname <- file.path(checkpointDir, paste0('gene', i, '.rds')))){
            res <- readRDS(fname)
            ##res <- NA
        } else{
            if(!is.null(fixed)){
                makeModelArgs <- c(list(zif=samp[,-i]), makeModelArgs, list(fixed=fixed))
            } else{
                makeModelArgs <- c(list(zif=samp[,-i]), makeModelArgs)
            }
            mm <- do.call(makeModel, makeModelArgs)
            blk <- Block(mm)
            res <- cgpaths(samp[,i], mm, Blocks=blk,  nodeId=colnames(samp)[i], ...)
            if(!is.null(checkpointDir)) saveRDS(res, fname)
    }
        res
    }))
    
    arr <- neighborhoodToArray(result, vnames=colnames(samp), nobs=nrow(samp))
    if(keepNodePaths){
        return(structure(arr, timing=timing, nodePaths=result))
    }
    return(structure(arr, timing=timing))
}


##' Set up the subsamples that will be taken for stability selection
##'
##' Stability selection is samples without replacement.  This initializes the sampling indices so that we can parallelize without worrying about the state of the random number generator (and recover from errors)
##' @param obs a matrix of observations from which rows will be sampled
##' @param strata an optional vector of stratifying covariates, over which the sampling will be balanced
##' @param B Number of resamples
##' @param seed random seed
##' @return list length \code{B} of indices
##' @export
setupStabilityIndex <- function(obs, strata=rep(1, nrow(obs)), B=50, seed=12345){
    if(!is.null(seed)) set.seed(12345)
    
    n <- nrow(obs)
    ns <- split(seq_len(n), strata)
    samplelist <- replicate(B, {
        x <- unlist(lapply(ns, function(x){
        s1 <- sample(seq_along(x), floor(length(x)/2))
        list(x[s1], x[setdiff(seq_along(x), s1)])
        }), recursive=FALSE)
        dim(x) <- c(length(ns), 2)
        y <- list(do.call(c, x[,1]), do.call(c, x[,2])) #glue together strata in each antithetical sample
        y
    })
    samplelist
}

##' Refit models for stability selection
##' 
##' The function \code{method} (which needs to follow the API of \code{autoGLM}) is called on subsampled data.
##' Exceptions are caught and output saved to disk since this can be quite computationally expensive.
##' @param obs a matrix of observations from which rows will be sampled
##' @param fixed a matrix of covariates
##' @param stabIndex output from \code{\link{setupStabilityIndex}}
##' @param step indices of components to run from \code{stabIndex}. Defaults to all.
##' @param method what method, eg, \code{fitHurdle} or \code{autoGLM}
##' @param stabilityCheckpointDir path to save output from each stability iteration
##' @param checkpointDir path to save intermediate output \emph{within} each stability iteration
##' @param ... arguments passed to \code{method}
##' @return list of output from \code{method}, eg, adjacency matrices.
##' @export
stability <- function(obs, fixed, stabIndex, step=seq_along(stabIndex), method, stabilityCheckpointDir=NULL, checkpointDir=NULL, ...){
    stabout <- list()
    accessDirOrDie(stabilityCheckpointDir)
    scheckpointDir <- sfixed <- NULL

    for(i in step){
        ss <- obs[stabIndex[[i]],,drop=FALSE]
        if(!is.null(fixed)){
            sfixed <- fixed[stabIndex[[i]],,drop=FALSE]
            ## test for rank-deficient fixed columns
            qrc <- qr(sfixed[,-1])
            nc <- qrc$rank
            lind <- 1:nc
            if(nc < (ncol(sfixed)-1)) sfixed <- cbind(1, sfixed[,qrc$pivot[lind]])
        } 
        if(!is.null(checkpointDir)){
            scheckpointDir <- file.path(dirname(checkpointDir), paste0(basename(checkpointDir), '_s', sprintf('%02d', i)))
        }
        stabi <- NULL
        tt <- try({
            stabi <- method(samp=ss, fixed=sfixed, checkpointDir=checkpointDir, ...)
        })
        if(inherits(tt, 'try-error')) warning(tt)
        if(!is.null(stabilityCheckpointDir) & !inherits(tt, 'try-error')){
            rfilename <- file.path(stabilityCheckpointDir, paste0('chk_s', sprintf('%02d', i)))
            saveRDS(stabi, file=paste0(rfilename, '.incomplete'))
            file.rename(paste0(rfilename, '.incomplete'), paste0(rfilename, '.rds'))
        }
        stabout[[i]] <- stabi
    }
    stabout
}

##' Process stability selected networks
##'
##' Stability selection attempts to estimate the probability
##' that an edge would have been included when we sample from the population.
##' Following Shah and Samsworth, we process the stability selection to give stability coefficients for each edge.
##' From Shah and Samsworth: p = #edges = m^2; q = max(iknot) and theta = q/p.
##' We bound the set of population transient edges, which have probabilities of selection below
##' theta (the same theta as above for convenience) when sampling from the population.
##' We call an edge \emph{empirically transient} when the stability coefficient is below tau.
##' Then for a sparsity theta < .1 and stability coefficient tau > .7, the ratio of empirical transient edges to population edges is less than 1%, and in fact typically more like .1%.
##' (Table 2 of Shah and Samworth)
##' @param stabout output of \code{\link{stability}}
##' @param theta sparsity of the solution in terms of the number of parameters. default .1.
##' @param stabilityCheckpointDir path to a directory to write checkpoint files
##' @return A list with components.  You probably want component 3, giving the stability coefficients for each edge.
### \code{stabCoefEdge}, a matrix of flattened, "vectorized" adjacency matrices; \code{estEdges} giving the number of edges estimated for each column \code{stabCoefEdge}, 
collectStability <- function(stabout, theta=.1, stabilityCheckpointDir=NULL){
    if(!is.null(stabilityCheckpointDir)){
        chkfiles <- list.files(stabilityCheckpointDir, pattern='chk_s.*.rds', full.names=TRUE)
        Get <- function(i) readRDS(chkfiles[i])
        ni <- length(chkfiles)
    } else{
        Get <- function(i) stabout[[i]]
        ni <- sum(!sapply(stabout, is.null))
    }
    x1 <- Get(1)
    m <- ncol(x1$adjMat[[1]])
    iknot <- floor(seq(0, 1.5*theta*m^2, length.out=50))
    e1 <- interpolateEdges(x1, knot=iknot)

    stabFlat <- Matrix(0, nrow=m^2, ncol=length(e1$trueEdges), sparse=TRUE)

    si <- 0
    for(i in seq_len(ni)){
        tt <- try({
        xx <- Get(i)
        ee <- interpolateEdges(xx, knot=iknot)
        hasedge <- (abs(sparseCbind(ee$edgeInterp))>0)*1
        #print(sum(hasedge))
        stabFlat <- stabFlat+ hasedge
        #print(table(as.vector(stabFlat)))
        si <- si+1
        })
        if(inherits(tt, 'try-error')) warning('Error reading trial ', i, ' in ', stabilityCheckpointDir)
    }
    stabFlat <- stabFlat/si
    fdri <- which.min(abs(iknot-floor(theta*m^2)))
    fdrStab <- stabFlat[,fdri]
    dim(fdrStab) <- c(m, m)
    list(stabCoefEdge=stabFlat, estEdges=iknot, highProbEdge=fdrStab)
}
