singletonMap <- function(nc, nf){
    blist <- as.list((nf+1):(nc+1))
    blist <- c(list(1:nf), blist)
    blist
}


##' @export
##' @import reshape
##' @import data.table
##' @import Matrix
##' @describeIn fitHurdle Fit an auto-model (Ising or Gaussian) to \code{samp} using glmnet
##' @param family in the case of \code{autoLogistic} one of "gaussian" or "logistic"
autoLogistic <- function(samp, fixed=NULL, nlambda=200, lambda.min.ratio=.1, parallel=FALSE, family='binomial', returnNodePaths=FALSE){
    samp0 <- if(family=='binomial') (abs(samp)>0)*1 else samp
    applyfun <- if(parallel) function(X, FUN) parallel::mclapply(X, FUN, mc.preschedule=TRUE) else lapply
    if(is.null(fixed)) fixed <- matrix(1, nrow=nrow(samp0))
    if(any(fixed[,1] != 1)) stop('Column 1 of `fixed` covariates must be intercept!')
    nid <- colnames(samp)
    blist <- singletonMap(ncol(samp)+ncol(fixed), ncol(fixed))
    
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
        list(path=path, blocks=blk, lambda=net$lambda, df=net$df, nodeId=thisId)
    }))
    arr <- neighborhoodToArray(result, vnames=colnames(samp))
    if(returnNodePaths){
        return(structure(arr, timing=timing, nodePaths=result))
    }
    return(structure(arr, timing=timing))
}

##' Fit the hurdle model coordinate-by-coordinate on a sample
##'
##' @param samp matrix of data, columns are variables
##' @param parallel parallelize over variables using "mclapply"?
##' @param checkpoint_dir (optional) directory to save the fit of each gene, useful for large problems.  If it exists, then completed genes will be automatically loaded.
##' @param makeModelArgs (optional) arguments passed to the model matrix function
##' @param returnNodePaths return node-wise output (solution paths and diagnostics for each node) as attribute `nodePaths`
##' @param indices (optional) subset of indices to fit, useful for cluster parallelization.
##' @param ... passed to cgpaths
##' @return list of fits, one per coordinate and an attribute "timing"
##' @export
fitHurdle <- function(samp, parallel=TRUE, checkpoint_dir=NULL, makeModelArgs=NULL, returnNodePaths=FALSE, indices, ...){
    applyfun <- if(parallel) function(X, FUN) parallel::mclapply(X, FUN, mc.preschedule=FALSE) else lapply
    allindices <- seq_len(ncol(samp))
    indices <- if(missing(indices))  allindices else indices
    if(length(setdiff(indices, allindices))>0) stop('`indices` out of range')
    timing <- system.time(result <- applyfun(indices, function(i){
        message('node=', i, ' nodeId=', colnames(samp)[i])
        if(!is.null(checkpoint_dir) && file.exists(fname <- file.path(checkpoint_dir, paste0('gene', i, '.rds')))){
            res <- readRDS(fname)
            ##res <- NA
        } else{
            mm <- do.call(makeModel, c(list(samp[,-i]), makeModelArgs))
            blk <- Block(mm)
            blk$id <- i
            res <- cgpaths(samp[,i], mm, Blocks=blk,  nodeId=colnames(samp)[i], ...)
            if(!is.null(checkpoint_dir)) saveRDS(res, fname)
    }
        res
    }))
    
    arr <- neighborhoodToArray(result, vnames=colnames(samp))
    if(returnNodePaths){
        return(structure(arr, timing=timing, nodePaths=result))
    }
    return(structure(arr, timing=timing))
}

##' Convert neighborhood estimates fit on differing lambda paths into adjacency matrices
##'
##' Join a series of paths fit node-wise into a list of adjacency matrices.
##' @param pathList a list of paths (coefficients in columns, paths in rows)
##' @param vnames a vector of names to be applied to the resulting adjacency matrix
##' @return list of (sparse) adjacency matrices, the number of non-zero elements for each matrix, and the lambda for each matrix.
neighborhoodToArray <- function(pathList, nknots, vnames=NULL){
    lambdaRange <- t(sapply(pathList, function(x){
        lambda <- if('lambda' %in% names(x))  x$lambda else 0
        c(range=range(lambda), n=length(lambda))
    }
                            ))
    nsol <- if(missing(nknots)) floor(max(log10(length(pathList)), 1)*median(lambdaRange[,'n'], na.rm=TRUE)) else nknots
    lambdaRange <- lambdaRange[lambdaRange[,'n']>3,]
    lmin <- min(lambdaRange[,'range1'], na.rm=TRUE)
    lmax <- max(lambdaRange[,'range2'], na.rm=TRUE)
    lpath <- exp(seq(log(lmin), log(lmax), length.out=nsol))
    if(is.null(vnames)){
        vnames <- sapply(pathList, '[[', 'nodeId')
        P <- length(vnames)
    } else{
        P <- length(vnames)
    }
    gridlist <- list()
    safeApprox <- getSafeApprox(lpath)
    
    for(i in seq_len(P)){
        gridlist[[i]] <- interpolateSummarizeCoefs(pathList[[i]], safeApprox)
    }
    
    allgrid <- rbindlist(gridlist)
    nodeId <- data.table(i=seq_along(vnames), nodeId1=vnames)
    setkey(allgrid, x, nodeId1,nodeId)
    allgrid <- merge(allgrid, nodeId, by='nodeId1')
    setnames(nodeId, c('i', 'nodeId1'), c('j', 'nodeId'))
    allgrid <- merge(allgrid, nodeId, by='nodeId')
    setkey(allgrid, x, i,j)
    adjMat <- list()
    nnz <- list()
    for(lidx in seq_along(lpath)){
        l <- lpath[lidx]
        ag <- allgrid[list(x=l),,nomatch=0]
        adjMat[[lidx]] <- Matrix::sparseMatrix(i=ag[,i], j=ag[,j], x=ag[,y], dims=c(P, P), dimnames=list(vnames, vnames))
        nnz[[lidx]] <- sum(abs(adjMat[[lidx]])>0)
    }
    list(adjMat=adjMat, trueEdges=unlist(nnz), lambda=lpath)
}


##' Return a function interpolating over a grid 
##'
##' In the interior of x, interpolation is used, otherwise the value of y at the endpoint is used.
##' Returns an empty list if length(x) = 0.  
##' @param lpath grid of points over which interpolation is required
##' @return function of (x,y) providing the interpolations for x=lpath
getSafeApprox <- function(lpath){
    fun <- function(x, y){
        lx <- length(x)
        if(lx>0){
            approx(x, y, lpath, rule=2, yright=0, method=ifelse(lx<2, 'constant', 'linear'))
        } else{
            list(x=numeric(0), y=numeric(0))
        }
    }
    return(fun)
}

##' Interpolate a solution path using approxFun and reduce the coefficients in a block to their signed L2 norm
##' @param sol solution (a list containing a sparse neighborhood matrix and the lambda over which it was evaluated)
##' @param blk a data table mapping between columns of sol$path and nodes
##' @param approxFun a function interpolating or otherwise approximating the solution between nodes
##' @return data.table with columns `y` giving normed, interpolated values, `x` giving lamba values and `block` giving the node in question
interpolateSummarizeCoefs <- function(sol, approxFun){
    ## mapping from parameters to blocks and nodes.  Only consider penalized blocks.
    blk <-  sol$blocks$map[lambda>0,.(paridx, nodeId)]
    ## path at node
    lambda <- sol$lambda
    trip <- toSparseTriples(sol$path)
    ldt <- data.table(i=seq_along(lambda), lambda=lambda)
    trip <- merge(trip, ldt, by='i')
    setnames(trip, c('j', 'x'), c('paridx', 'Coef'))
    trip <- merge(trip, blk, by='paridx')
    interpolate <- trip[,approxFun(x=lambda, y=Coef), keyby=list(paridx, nodeId)]
    summarized <- interpolate[,list(y=mean(y^2)*sign(y[which.max(abs(y))]) ), keyby=list(x,nodeId)]
    summarized[,nodeId1:=sol$nodeId]
    summarized[abs(y)>0,]
}

toSparseTriples <- function(mat){
    if(!inherits(mat, 'sparseMatrix')){
        mat <- Matrix::Matrix(mat, sparse=TRUE)
    }
    as.data.table(Matrix::summary(mat))
}

onlyTri <- function(mat, diag=FALSE, upper=TRUE)if(upper) mat[upper.tri(mat)] else mat[lower.tri(mat)]


##' Interpolate set of adjacency matrix across a set of edges
##'
##' @param Mpath list of (possibly signed, weighted) adjacency matrices
##' @param lambda penalty parameter that generated each matrix
##' @param knot optional vector of number of edges over which to conduct the interpolation.  If missing then we use something on the order of P*log(P), where P is the number of nodes.
##' @param nknot number knots over which to interpolate
##' @export
##' @return a list of `edgeInterp`: contains a list of length `nknot` of sparse adjacency matrices, `estEdges` the desired number of edges to be interpolated, `trueEdges` the actual number of edges
interpolateEdges <- function(Mpath, lambda, knot, nknot=100){
    lo <- order(lambda)
    Mpath <- Mpath[lo]
    lambda <- lambda[lo]
    Mpath <- lapply(Mpath, function(x) (x+Matrix::t(x))/2)
    nnz <-  sapply(Mpath, function(x) sum((abs(x)+abs(Matrix::t(x)))>0))
    P <- ncol(Mpath[[1]])
    if(missing(knot)){
        maxknot <- if(P<200) 2*P*sqrt(P) else P*log(P)/2 + 400*sqrt(200)
        knot <- seq(min(maxknot, max(nnz)), 2, length.out=nknot)
    }
    stopifnot(length(Mpath)==length(lambda))
    lknots <- approx(nnz, lambda, knot)$y
    edgeInterp <- list()
    for(i in seq_along(lknots)){
        lambdai <- lknots[i] #lambda corresponding to current edge count
        browser(expr=all(lambdai<lambda) || all(lambdai>lambda))
        lidx <- findInterval(lambdai, lambda, all.inside=TRUE) #greatest lambda smaller than lambdai
        gamma <- (lambda[lidx]-lambdai)/(lambda[lidx]-lambda[lidx+1]) #how much closer is lambdai to lambda[lidx] than lambda[lidx+1]
        edgeInterp[[i]] <- Mpath[[lidx]]*(1-gamma)+Mpath[[lidx+1]]*(gamma)
        if(gamma>.5){
            ## adjust zero entries
            nzi <- Matrix::which(abs(Mpath[[lidx]])>0)
            nziplus1 <- Matrix::which(abs(Mpath[[lidx+1]])>0)
            change <- setdiff(nzi, nziplus1)
            edgeInterp[[i]][change] <- 0
        }
    }
    truenz <- sapply(edgeInterp, function(x) sum(abs(x)+abs(Matrix::t(x))>0))
    list(edgeInterp=edgeInterp, estEdges=knot, trueEdges=truenz)
}
