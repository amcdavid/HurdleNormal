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
##' @param summaryFun function to reduce a vector-valued parameter set at each node to a scalar (representing an edge weight).  Defaults to the "signed" L2 norm.
##' @param nobs number of observations the model was fit; used to calculate BIC
##' @param self_edges should self edges (loops) be returned in the adjacency matrix; this allows inference of intercept quantities.
##' @return list of (sparse) adjacency matrices, the number of non-zero elements (edges, before enforcing symmetry) for each matrix, the lambda for each matrix, the non-penalized, refitted pseudo log-likelihood, the number of parameters per edge, and the BIC
##' @seealso fitHurdle, autoGLM, interpolateEdges
##' @title neighborhoodToArray
##' @aliases print.HurdleNormalFit
##' @export
neighborhoodToArray <- function(pathList, nknots, vnames=NULL, summaryFun=summarySignedL1, nobs, self_edges=FALSE){
    lambdaRange <- t(sapply(pathList, function(x){
        lambda <- if('lambda' %in% names(x))  x$lambda else 0
        c(range=range(lambda), n=length(lambda))
    }
                            ))
    lambdaRange <- lambdaRange[lambdaRange[,'n']>3,]
    nsol <- if(missing(nknots)) floor(max(log10(length(pathList)), 1)*median(lambdaRange[,'n'], na.rm=TRUE)) else nknots
    if(is.na(nsol)) nsol <- 2
    lmin <- min(lambdaRange[,'range1'], 1e6, na.rm=TRUE)
    lmax <- max(lambdaRange[,'range2'], 1e-6, na.rm=TRUE)
    lpath <- exp(seq(log(lmin), log(lmax), length.out=nsol))
    if(is.null(vnames)){
        vnames <- sapply(pathList, '[[', 'nodeId')
        P <- length(vnames)
    } else{
        P <- length(vnames)
    }
    gridlist <- list()
    loglikmatrix <- matrix(0, nrow=P, ncol=length(lpath))
    safeApprox <- getSafeApprox(lpath)
    safeApproxPath <- function(x, y) safeApprox(x, y, yright=0)
    fail = 0
    for(i in seq_along(pathList)){
        if(inherits(pathList[[i]], 'SolPath')){
            gridlist[[i]] <- interpolateSummarizeCoefs(pathList[[i]], safeApproxPath, summaryFun, self_edges)
            llnp <- pathList[[i]]$loglik_np
            ## Take previous loglik on path as lower bound
            ## Needed especially when lambda is disjoint for some coordinates
            loglikmatrix[i,] <- safeApprox(pathList[[i]]$lambda, llnp, method='constant', f=1)$y
        } else{
            fail = fail + 1
        }
    }
    if(fail == length(pathList)) stop("No solution paths found in `pathList`")
    if(fail > 0) warning(fail, " failures in pathList, these nodes will be excluded from graph.")
    
    allgrid <- rbindlist(gridlist)
    nodeId <- data.table(i=seq_along(vnames), nodeId1=vnames)
    setkey(allgrid, x, nodeId1,nodeId)
    allgrid <- merge(allgrid, nodeId, by='nodeId1')
    setnames(nodeId, c('i', 'nodeId1'), c('j', 'nodeId'))
    allgrid <- merge(allgrid, nodeId, by='nodeId')
    setkey(allgrid, x, i,j)
    adjMat <- list()
    n_param_per_edge <- nnz <- list()
    for(lidx in seq_along(lpath)){
        l <- lpath[lidx]
        ag <- allgrid[list(x=l),,nomatch=0]
        adjMat[[lidx]] <- Matrix::sparseMatrix(i=ag[,i], j=ag[,j], x=ag[,y], dims=c(P, P), dimnames=list(vnames, vnames))
        nnz[[lidx]] <- sum(abs(adjMat[[lidx]])>0)
        n_param_per_edge[[lidx]] <- sum(ag$npar)/
            (if(nnz[[lidx]]>0) nnz[[lidx]] else 1) #don't divide by zero
    }
    out = list(adjMat=adjMat, trueEdges=unlist(nnz), lambda=lpath, pseudo_loglik_np=colSums(loglikmatrix, na.rm = TRUE), n_param_per_edge=unlist(n_param_per_edge))
    out$BIC = -2*out$pseudo_loglik_np + log(nobs)*out$trueEdges*out$n_param_per_edge
    BIC_etc <- data.table(trueEdges=out$trueEdges, lambda=out$lambda, pseudo_loglik_np=out$pseudo_loglik_np, BIC=out$BIC, adjMat_idx=seq_along(adjMat))
    out$BIC_etc <- BIC_etc
    class(out) <- 'HurdleNormalFit'
    out
}


##' @export
print.HurdleNormalFit <- function(x, ...){
    cat(sprintf('A fitted network on %d nodes containing a solution path over %d graphs.
Whose edge counts range from %d-%d edges.
Components are %s.
',
nrow(x$adjMat[[1]]), length(x$adjMat), min(x$trueEdges), max(x$trueEdges), paste(names(x), collapse=', ')))
}

##' Return a function interpolating over a grid 
##'
##' In the interior of x, interpolation is used, otherwise the value of y at the endpoint is used.
##' Returns an empty list if length(x) = 0.  
##' @param lpath grid of points over which interpolation is required
##' @return function of (x,y) providing the interpolations for x=lpath
getSafeApprox <- function(lpath){
    fun <- function(x, y, method, ...){
        lx <- length(unique(x))
        if(sum(!is.na(y))<2 && lx >= 2 ){           #approx errors out with less than two non-NA values ...
            return(list(x=lpath, y=rep(NA_real_, length(lpath))))
        }
        if(lx>0){
            if(missing(method)){
                method <- ifelse(lx<2, 'constant', 'linear')
            }
            return(approx(x, y, lpath, rule=2, method=method, ...))
        } else{
            return(list(x=numeric(0), y=numeric(0)))
        }
    }
    return(fun)
}

##' Interpolate a solution path using approxFun and reduce the coefficient vector in a block to a scalar
##' @param sol a SolPath solution
##' @param approxFun a function interpolating or otherwise approximating the solution between nodes
##' @inheritParams neighborhoodToArray
##' @return data.table with columns `y` giving normed, interpolated values, `x` giving lamba values and `block` giving the node in question
interpolateSummarizeCoefs <- function(sol, approxFun, summaryFun, self_edges){
    ## mapping from parameters to blocks and nodes.  Only consider penalized blocks.
    #blk <-  sol$blocks$map[lambda>0,.(paridx, nodeId)]
    blk <-  sol$blocks$map[,.(paridx, nodeId)]
    ## path at node
    lambda <- sol$lambda
    trip <- toSparseTriples(sol$path)
    ldt <- data.table(i=seq_along(lambda), lambda=lambda)
    trip <- merge(trip, ldt, by='i')
    setnames(trip, c('j', 'x'), c('paridx', 'Coef'))
    trip <- merge(trip, blk, by='paridx')
    interpolate <- trip[,approxFun(x=lambda, y=Coef), keyby=list(paridx, nodeId)]
    summarized <- interpolate[,list(y=summaryFun(y, nodeId), npar=.N), keyby=list(x,nodeId)]
    summarized[,nodeId1:=sol$nodeId]
    if(self_edges){
        summarized[nodeId=='(fixed)', nodeId:=sol$nodeId]
    } else{
        summarized <- summarized[nodeId!='(fixed)']
    }
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
    ## intercept of th_c (second block) or first value of th_d
    ## Fixme (maybe) in makeModel
    if(nodeId=='(fixed)') y[2] else y[1]
}

summaryG <- function(y, nodeId){
    ## intercept of th_d (first block) or second value of th_d
    ## Fixme (maybe) in makeModel
    if(nodeId=='(fixed)') y[1] else y[2]
}

summaryK <- function(y, nodeId){
    ## precision or first value of th_c
    ## Fixme (maybe) in makeModel
    y[3]
}

summaryHji <- function(y, nodeId){
    ## intercept of th_c 
    if(nodeId=='(fixed)') y[2] else y[4]
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
##' @param array \code{neighborhoodToArray} output
##' @param knot optional vector of number of edges over which to conduct the interpolation.  If missing then we use something on the order of P*log(P), where P is the number of nodes.
##' @param nknot number knots over which to interpolate
##' @export
##' @return a list of `edgeInterp`: contains a list of length `nknot` of sparse adjacency matrices, `estEdges` the desired number of edges to be interpolated, `trueEdges` the actual number of edges, `BIC`, the interpolated BIC.
interpolateEdges <- function(array, knot, nknot=100){
    lambda <- array$lambda
    Mpath <- array$adjMat
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

    safeApprox <- getSafeApprox(lknots)
    BIC <- safeApprox(lambda, array$BIC, f=1, method='constant')$y
    truenz <- sapply(edgeInterp, function(x) sum(abs(x)+abs(Matrix::t(x))>0))
    list(edgeInterp=edgeInterp, estEdges=knot, trueEdges=truenz, BIC=BIC)
}
