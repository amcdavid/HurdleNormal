## Mapping from parameter indices to blocks
## Only used in simulation so far...redundant with Block()?
nativeMap <- function(p){
    covariates <- rep(seq(2, p), each=2)
    block <- c(1,                       #intercept
               covariates)
    total <- c(block, block, 1) #kbb
    total
}



singletonMap <- function(nc, nf){
    blist <- as.list((nf+1):(nc+1))
    blist <- c(list(1:nf), blist)
    blist
}

##' Recursively flatten and bind together sparse matrices
##'
##' @param x a matrix, or a list of matrices.
##' @param sparse should a sparse result be returned?
##' @return a Matrix, each column being a flattened matrix.
##' @export
sparseCbind <- function(x, sparse=TRUE){
    if(is.list(x)){
        if(inherits(x[[1]], 'sparseMatrix')){
            if(ncol(x[[1]])==1)             return(Matrix(do.call(cBind, x), sparse=sparse))
            return(sparseCbind(lapply(x, function(y){
                dim(y) <- c(nrow(y)^2, 1)
                y
            }), sparse=sparse))
        } else{
            return(lapply(x, sparseCbind, sparse=sparse))
        }
    }
    stop('Unreachable code')
}


##' Convert neighborhood estimates fit on differing lambda paths into adjacency matrices
##'
##' Join a series of paths fit node-wise into a list of adjacency matrices.
##' @param pathList a list of paths (coefficients in columns, paths in rows)
##' @param nknots number of lambda knots.  Be default, use the number of lambda in each neighborhood in pathList
##' @param vnames a vector of names to be applied to the resulting adjacency matrix
##' @param summaryFun  a function summarizing the node-level interaction parameters.  Should be a function of the vector of interaction parameters and the nodeId that returns a single numeric value.
##' @return list of (sparse) adjacency matrices, the number of non-zero elements for each matrix, and the lambda for each matrix.
##' @export
neighborhoodToArray <- function(pathList, nknots, vnames=NULL, summaryFun=summarySignedL1){
    lambdaRange <- t(sapply(pathList, function(x){
        lambda <- if('lambda' %in% names(x))  x$lambda else 0
        c(range=range(lambda), n=length(lambda))
    }
                            ))
    lambdaRange <- lambdaRange[lambdaRange[,'n']>3,]
    nsol <- if(missing(nknots)) floor(max(log10(length(pathList)), 1)*median(lambdaRange[,'n'], na.rm=TRUE)) else nknots
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
        if(inherits(pathList[[i]], 'SolPath')){
            gridlist[[i]] <- interpolateSummarizeCoefs(pathList[[i]], safeApprox, summaryFun)
        }
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
##' @param approxFun a function interpolating or otherwise approximating the solution between nodes
##' @param summaryFun a function summarizing the node-level interaction parameters.  Should be a function of the vector of interaction parameters and the nodeId that returns a single numeric value.
##' @param blk a data table mapping between columns of sol$path and nodes
##' @return data.table with columns `y` giving normed, interpolated values, `x` giving lamba values and `block` giving the node in question
interpolateSummarizeCoefs <- function(sol, approxFun, summaryFun){
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
    summarized <- interpolate[,list(y=summaryFun(y, nodeId)), keyby=list(x,nodeId)]
    summarized[,nodeId1:=sol$nodeId]
    summarized[abs(y)>0,]
}

summaryMaxL2 <- function(y, nodeId){
    sqrt(mean(y^2))*sign(y[which.max(abs(y))])
}

summaryL2 <- function(y, nodeId){
    sqrt(mean(y^2))
}

summarySignedL1 <- function(y, nodeId){
    pospart <- sum(y[y>0])
    negpart <- sum(-y[y<0])
    if(pospart >= negpart) pospart else -negpart
}

summaryHij <- function(y, nodeId){
    y[1]
}

summaryG <- function(y, nodeId){
    y[2]
}

summaryK <- function(y, nodeId){
    y[3]
}

summaryHji <- function(y, nodeId){
    y[4]
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
    lknots <- approx(nnz, lambda, knot, rule=2)$y
    edgeInterp <- list()
    for(i in seq_along(lknots)){
        lambdai <- lknots[i] #lambda corresponding to current edge count
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
