##' @export
##' @import reshape
##' @import glmnet
##' @import data.table
autoLogistic <- function(samp){
    samp0 <- (abs(samp)>0)*1
     lapply(seq_len(ncol(samp)), function(i){
        net <- glmnet(samp0[,-i], samp0[,i], family='binomial',lambda.min.ratio=.1, nlambda=200)
        path <- t(as.matrix(coef(net)))
        rownames(path) <- net$lambda
        list(path=path, lambda=net$lambda, df=net$df)
    })
}

fitHurdle <- function(samp, nlambda=100, control=list(tol=1e-3, maxrounds=3000, maxit=200, debug=0, method='proximal', updatehess=2500, stepsize=1), lambda.min.ratio=.1, makeModelArgs=NULL, parallel=TRUE, ...){
    applyfun <- if(parallel) mclapply else lapply
    applyfun(seq_len(ncol(samp)), function(i){
        message('i=', i)
        mm <- do.call(makeModel, c(list(samp[,-i]), makeModelArgs))
        cgpaths(samp[,i], mm, Blocks=Block(mm), nlambda=nlambda, control=control, lambda.min.ratio=.1, ...)
    })
}

## each element of pathlist should have its lambda as rownames
## coordmap gives mapping of elements in the pathlist to coordinates (includes intercept terms)
nsol <- 100
neighborhoodToArray <- function(pathList, Coordmap, vnames=NULL){
    require(data.table)
lambdaRange <- t(sapply(pathList, function(x) range(as.numeric(rownames(x$path)))))
lmin <- min(lambdaRange[,1])
lmax <- max(lambdaRange[,2])
lpath <- exp(seq(log(lmin), log(lmax), length.out=nsol))
P <- max(Coordmap) 
sol <- array(NA, dim=c(P, P, nsol))
gridlist <- list()
for(i in seq_len(P)){
    coefIndex <- c(i, setdiff(1:P, i))[Coordmap]
    M <- data.table(lambda=as.numeric(rownames(pathList[[i]]$path)), Coef=as.numeric(pathList[[i]]$path), group=rep(coefIndex, each=nrow(pathList[[i]]$path)))
    inModel <- M[,list(inModel=mean(Coef^2)*sign(Coef[which.max(abs(Coef))]) ), keyby=list(lambda,group)]
    grid <- inModel[,approx(x=lambda, y=inModel, xout=lpath, rule=2),keyby=list(group)]
    grid[,':='(i=i)]
    gridlist[[i]] <- grid
}
allgrid <- rbindlist(gridlist)
setkey(allgrid, x, i,group)
array(allgrid$y, dim=c(P, P, nsol), dimnames=list(vnames, vnames, lpath))
}

onlyTri <- function(mat, diag=FALSE, upper=TRUE)if(upper) mat[upper.tri(mat)] else mat[lower.tri(mat)]

