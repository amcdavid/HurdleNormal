##' @export
##' @import reshape
##' @import data.table
##' @import Matrix
##' @describeIn fitHurdle Fit an auto-logistic (Ising) model to \code{samp} using glmnet
##' @param fixed unpenalized covariates to include.  NOT IMPLEMENTED for autoLogistic yet.
##' @param nlambda number of knots in solution path
##' @param lambda.min.ratio minimum lambda as function lambda0, the smallest lambda such that all coordinates are zero
autoLogistic <- function(samp, fixed=NULL, nlambda=200, lambda.min.ratio=.1, parallel=FALSE){
    samp0 <- (abs(samp)>0)*1
     applyfun <- if(parallel) function(X, FUN) parallel::mclapply(X, FUN, mc.preschedule=TRUE) else lapply
    timing <- system.time(result <- applyfun(seq_len(ncol(samp)), function(i){
        model <- cbind(fixed, samp0[,-i])
        penalty.factor <- rep(1, ncol(model))
        if(!is.null(fixed)) penalty.factor[seq_len(fixed)] <- 0

        posobs <- sum(samp0[,i])
        if( posobs > 2 && (nrow(samp0)-posobs)>2){
            net <- glmnet::glmnet(model, samp0[,i], family='binomial',lambda.min.ratio=lambda.min.ratio, nlambda=nlambda, penalty.factor=penalty.factor)
            path <- Matrix::t(coef(net))
        } else{
            net <- list(df=c(1), lambda=c(0))
            path <- Matrix::Matrix( c(1, rep(0, ncol(model)-1)), nrow=1, sparse=TRUE)
        }
        rownames(path) <- net$lambda
        list(path=path, lambda=net$lambda, df=net$df)
    }))
    arr <- neighborhoodToArray(result, 1:ncol(samp), colnames(samp))
    structure(arr, timing=timing)
}

##' Fit the hurdle model coordinate-by-coordinate on a sample
##'
##' @param samp matrix of data, columns are variables
##' @param parallel parallelize over variables using "mclapply"?
##' @param checkpoint_dir (optional) directory to save the fit of each gene, useful for large problems.  If it exists, then completed genes will be automatically loaded.
##' @param makeModelArgs (optional) arguments passed to the model matrix function
##' @param returnNodePaths return node-wise output (solution paths and diagnostics for each node) as attribute `nodePaths`
##' @param indices (optional) subset of indices to fit, useful for cluster parallelization.  NOT IMPLEMENTED.
##' @param ... passed to cgpaths
##' @return list of fits, one per coordinate and an attribute "timing"
##' @export
fitHurdle <- function(samp, parallel=TRUE, checkpoint_dir=NULL, makeModelArgs=NULL, returnNodePaths=FALSE, indices, ...){
    applyfun <- if(parallel) function(X, FUN) parallel::mclapply(X, FUN, mc.preschedule=FALSE) else lapply
    timing <- system.time(result <- applyfun(seq_len(ncol(samp)), function(i){
        message('i=', i)
        if(!is.null(checkpoint_dir) && file.exists(fname <- file.path(checkpoint_dir, paste0('gene', i, '.rds')))){
            res <- readRDS(fname)
            ##res <- NA
        } else{
            mm <- do.call(makeModel, c(list(samp[,-i]), makeModelArgs))
            blk <- Block(mm)
            res <- cgpaths(samp[,i], mm, Blocks=blk,  ...)
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
##' @details
##' The \code{Coordmap} can be missing when \code{cgpaths} has been run, as the mapping is provided in the output.  Otherwise, \code{Coordmap} should be a vector, giving the node index for each coefficient in the model (to accodomdate models with several coefficients per node).  The node index `1` is assumed to correspond to intercept (nuisance) terms.
##' @param pathList a list of paths (coefficients in columns, paths in rows)
##' @param Coordmap A mapping (vector of length ncol(pathList[[1]])) between coefficient-indices and variables.  See details.
##' @param vnames a vector of names to be applied to the resulting adjacency matrix
##' @return list of (sparse) adjacency matrices, the number of non-zero elements for each matrix, and the lambda for each matrix.
neighborhoodToArray <- function(pathList, Coordmap, vnames=NULL){
    lambdaRange <- t(sapply(pathList, function(x){
        lambda <- if('lambda' %in% names(x))  x$lambda else 0
        c(range=range(lambda), n=length(lambda))
    }
        ))
    nsol <- floor(max(log10(length(pathList)), 1)*median(lambdaRange[,'n'], na.rm=TRUE))
    lambdaRange <- lambdaRange[lambdaRange[,'n']>3,]
    lmin <- min(lambdaRange[,'range1'], na.rm=TRUE)
    lmax <- max(lambdaRange[,'range2'], na.rm=TRUE)
    lpath <- exp(seq(log(lmin), log(lmax), length.out=nsol))
    if(missing(Coordmap)){
        blk <- pathList[[1]]$blocks$map
        setnames(blk, 'lambda', 'lambdaPenalty')
        setkey(blk, paridx)
        Coordmap <- blk[,block]
    } else{
        blk <- data.table(block=Coordmap, paridx=seq_along(Coordmap))
    }
    P <- max(blk$block)
    gridlist <- list()
    safeApprox <- function(x, y){
        lx <- length(x)
        if(lx>0){
            approx(x, y, lpath, rule=2, yright=0, method=ifelse(lx<2, 'constant', 'linear'))
        } else{
            list(x=numeric(0), y=numeric(0))
        }
    }
    for(i in seq_len(P)){
        #coefIndex <- c(i, setdiff(1:P, i))[Coordmap]
        lambda <- pathList[[i]]$lambda
        trip <- toSparseTriples(pathList[[i]]$path)
        ldt <- data.table(i=seq_along(lambda), lambda=lambda)
        trip <- merge(trip, ldt, by='i')
        setnames(trip, c('j', 'x'), c('paridx', 'Coef'))
        trip <- merge(trip, blk, by='paridx')
        trip <- trip[block!=1]
        #M <- data.table(lambda=pathList[[i]]$lambda, Coef=as.numeric(pathList[[i]]$path), group=rep(coefIndex, each=nrow(pathList[[i]]$path)))
        
        inModel <- trip[,list(inModel=mean(Coef^2)*sign(Coef[which.max(abs(Coef))]) ), keyby=list(lambda,block)]
        
        grid <- inModel[,safeApprox(x=lambda, y=inModel),keyby=list(block)]
        
        grid[,':='(i=i)]#, j=coefIndex)]
        gridlist[[i]] <- grid[abs(y)>0,]
    }
    allgrid <- rbindlist(gridlist)
    allgrid <- allgrid[,j:=ifelse(block<=i, block-1, block)]
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

toSparseTriples <- function(mat){
    if(!inherits(mat, 'sparseMatrix')){
        mat <- Matrix::Matrix(mat, sparse=TRUE)
    }
    as.data.table(Matrix::summary(mat))
}

onlyTri <- function(mat, diag=FALSE, upper=TRUE)if(upper) mat[upper.tri(mat)] else mat[lower.tri(mat)]


interpolateEdges <- function(Mpath, lambda, knot, nknot=100){
    if(all(sort(lambda) == rev(lambda))){
        lo <- order(lambda)
        Mpath <- Mpath[lo]
        lambda <- lambda[lo]
    }
    nnz <-  sapply(Mpath, function(x) sum(abs(x+Matrix::t(x))>0))
    P <- ncol(Mpath[[1]])
    if(missing(knot)){
        maxknot <- if(P<200) 2*P*sqrt(P) else P*log(P)/2 + 400*sqrt(200)
        knot <- seq(2, min(maxknot, max(nnz)), length.out=nknot)
    }
    stopifnot(length(Mpath)==length(lambda))
    lknots <- approx(nnz, lambda, knot)$y
    edgeInterp <- list()
    for(i in seq_along(lknots)){
        lambdai <- lknots[i] #lambda corresponding to current edge count
        browser(expr=all(lambdai<lambda) || all(lambdai>lambda))
        lidx <- findInterval(lambdai, lambda, all.inside=TRUE) #smallest lambda greater than lambdai
        gamma <- (lambda[lidx]-lambdai)/(lambda[lidx]-lambda[lidx+1]) #how much closer is lambdai to lambda[lidx+1] than lambda[lidx]
        edgeInterp[[i]] <- (Mpath[[lidx]]*(1-gamma)+Mpath[[lidx+1]]*gamma)
        nzi <- Matrix::which(abs(Mpath[[lidx]])>0)
        edgeInterp[[i]][nzi] <- round(edgeInterp[[i]][nzi]/Mpath[[lidx]][nzi], 0)*edgeInterp[[i]][nzi]
    }
    truenz <- sapply(edgeInterp, function(x) sum(abs(x+Matrix::t(x))>0))
    list(edgeInterp=edgeInterp, estEdges=knot, trueEdges=truenz)
}
