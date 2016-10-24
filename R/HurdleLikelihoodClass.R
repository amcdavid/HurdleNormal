setClass("HurdleLikelihood", representation(pointer = "externalptr", p='integer'))

HurdleLikelihood_method <- function(name) {
    paste( "HurdleLikelihood", name, sep = "__" )
}

validateTheta <- function(theta, p){
    stopifnot(length(theta)==2*p+1)
    stopifnot(is.numeric(theta))
}

HurdleLikelihood_gradAll <- function(x, theta, penalize=TRUE){
    stopifnot(length(theta)==2*x@p+1)
    as.numeric(.Call(HurdleLikelihood_method('gradAll'), x@pointer, theta, penalize))
}

HurdleLikelihood_grad <- function(x, theta, grp, penalize=TRUE){
    stopifnot(grp>=1 && grp<=x@p)
    stopifnot((grp==1 && length(theta)==3) || (grp>1 && length(theta)==4))
    .Call(HurdleLikelihood_method('grad'), x@pointer, theta, grp-2, penalize)
}


HurdleLikelihood_LL <- function(x, theta, grp, penalize=TRUE){
    stopifnot(grp>=1 && grp<=x@p)
    stopifnot((grp==1 && length(theta)==3) || (grp>1 && length(theta)==4))
    .Call(HurdleLikelihood_method('LL'), x@pointer, theta, grp-2, penalize)
}

HurdleLikelihood_setLambda <- function(x, lambda){
    stopifnot(is.numeric(lambda))
    if(length(lambda)==1) lambda <- rep(lambda, x@p-1)
    stopifnot(length(lambda)==x@p-1)
    x@pointer <- .Call(HurdleLikelihood_method('setLambda'), x@pointer, lambda)
    NULL
}



## syntactic sugar to allow object$method( ... )
setMethod( "$", "HurdleLikelihood", function(x, name ) {
    if(name == 'gradAll'){
        function(theta, penalize=TRUE){
            if(missing(theta)){
                as.numeric(.Call(HurdleLikelihood_method('gradAllFixed'), x@pointer, penalize))
            } else{
                HurdleLikelihood_gradAll(x, theta, penalize)
            }
        }
    }else if(name == 'grad'){
        function(theta, grp, penalize=TRUE) HurdleLikelihood_grad(x, theta, grp, penalize)
    }else if(name == 'LL'){
        function(theta, grp, penalize=TRUE) HurdleLikelihood_LL(x, theta, grp, penalize)
    }  else if(name == 'setLambda'){
        function(lambda) HurdleLikelihood_setLambda(x,lambda)
    }else if(name=='LLall'){
        function(theta, penalize=TRUE){
            validateTheta(theta, x@p)
            .Call(HurdleLikelihood_method('LLall'), x@pointer, theta, penalize)
        }
    }else {
        function(...) .Call(HurdleLikelihood_method(name), x@pointer, ... )
    }
}
          )

HurdleLikelihood <- function(y, x, grp, theta, lambda=0, tol=1e-4){
    if(!is.numeric(x) || !is.matrix(x)) stop('`x` must be numeric matrix')
    if(!is.numeric(y) || length(y) != nrow(x)) stop('`y` must be numeric and length `nrow(x)`.')
    if(missing(grp)){
        p <- ncol(x)
        ## pminus1 <- (ncol(x)-1)/2
        ## if(abs(pminus1-floor(pminus1))>.1) stop('Expecting odd number of columns in x (that you included intercept).  If `x` is correct you will need to manually specify `grp`')
        ## ##grp <- c(0, seq_len(pminus1), seq_len(pminus1))
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
    new('HurdleLikelihood', pointer=.Call(HurdleLikelihood_method("new"), y, x, grp, theta, lambda, tol), p=as.integer(ncol(x)))
}
