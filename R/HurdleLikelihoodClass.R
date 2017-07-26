validateTheta <- function(theta, p){
    stopifnot(length(theta)==2*p+1)
    stopifnot(is.numeric(theta))
}

HurdleLikelihood <- function(y, x, grp, theta, lambda=0, tol=1e-4){
    if(!is.numeric(x) || !is.matrix(x)) stop('`x` must be numeric matrix')
    if(!is.numeric(y) || length(y) != nrow(x)) stop('`y` must be numeric and length `nrow(x)`.')
    if(missing(grp)){
        p <- ncol(x)
        grp <- c(seq_len(p), seq_len(p))
    } else{
        stop("Grouping not implemented yet")
    }
    if(missing(theta)){
        theta <- c(rep(0, 2*length(grp)), kbb=1)
    }
    if(length(theta) != 2*ncol(x)+1) stop('2*ncol(x) +1  not equal theta length')
    if(!is.numeric(theta)) stop("`theta` must be numeric")
    if(theta[length(theta)]<=0) stop('kbb non-positive')
    if(length(lambda) ==1) lambda <- rep(lambda, ncol(x))
    if(length(lambda) != ncol(x)) stop('lambda must match theta length')
    if(!is.numeric(lambda)) stop("`lambda` must be numeric")
    if(any(floor(grp)!=grp)) stop("`grp` must be integer")
    if(!is.numeric(tol)) stop("`tol` must be numeric")
    new(HurdleLikelihoodMod, y, x, grp, theta, lambda, tol)
}
