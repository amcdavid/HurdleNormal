setClass("HurdleLikelihood", representation( pointer = "externalptr", p='integer') )

HurdleLikelihood_method <- function(name) {
paste( "HurdleLikelihood", name, sep = "__" )
 }

HurdleLikelihood_gradAll <- function(x, theta, penalize=TRUE){
    grp <- coordmap(x@p)-2
    stopifnot(length(theta)==length(grp))
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
    .Call(HurdleLikelihood_method('setLambda'), x@pointer, lambda)
}



## syntactic sugar to allow object$method( ... )
setMethod( "$", "HurdleLikelihood", function(x, name ) {
    if(name == 'gradAll'){
        function(theta, penalize=TRUE) HurdleLikelihood_gradAll(x, theta, penalize)
    }else if(name == 'grad'){
        function(theta, grp, penalize=TRUE) HurdleLikelihood_grad(x, theta, grp, penalize)
    }else if(name == 'LL'){
        function(theta, grp, penalize=TRUE) HurdleLikelihood_LL(x, theta, grp, penalize)
    }  else if(name == 'setLambda'){
        function(lambda) HurdleLikelihood_setLambda(x,lambda)
    }else if(name=='LLall'){
        function(theta, penalize) .Call(HurdleLikelihood_method('LLall'), x@pointer, theta, penalize=TRUE)
    }else {
        function(...) .Call(HurdleLikelihood_method(name), x@pointer, ... )
    }
}
          )

HurdleLikelihood <- function(y, x, grp, theta, lambda=0, tol=1e-4){
    if(!is.numeric(x) || !is.matrix(x)) stop('`x` must be numeric matrix')
    if(!is.numeric(y) || length(y) != nrow(x)) stop('`y` must be numeric and length `nrow(x)`.')
    
    if(missing(grp)) grp <- coordmap(ncol(x)+1)-2
    if(missing(theta)){
        rep(0, length(grp))
        theta[length(grp)] <- 1
    }
    if(length(theta) != length(grp)) stop('grp length not equal theta length')
    if(!is.numeric(theta)) stop("`theta` must be numeric")
    if(theta[length(grp)]<0) stop('kbb negative')
    if(length(lambda) ==1) lambda <- rep(lambda, ncol(x))
    if(length(lambda) != ncol(x)) stop('lambda must match theta length')
    if(!is.numeric(lambda)) stop("`lambda` must be numeric")
    if(any(floor(grp)!=grp)) stop("`grp` must be integer")
    if(!is.numeric(tol)) stop("`tol` must be numeric")
    new('HurdleLikelihood', pointer=.Call(HurdleLikelihood_method("new"), y, x, grp, theta, lambda, tol), p=as.integer(ncol(x)+1))
}
