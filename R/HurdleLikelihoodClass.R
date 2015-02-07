setClass("HurdleLikelihood", representation( pointer = "externalptr", p='integer') )

HurdleLikelihood_method <- function(name) {
paste( "HurdleLikelihood", name, sep = "__" )
 }

HurdleLikelihood_gradAll <- function(x, theta, penalize=TRUE){
    grp <- coordmap(x@p)-2
    stopifnot(length(theta)==length(grp))
    grad <- theta
    for(g in unique(grp)){
        grad[grp==g] <- .Call(HurdleLikelihood_method('grad'), x@pointer, theta[grp==g], g, penalize)
    }
    grad
}

HurdleLikelihood_grad <- function(x, theta, grp, penalize=TRUE){
    stopifnot(grp>=1 && grp<=x@p)
    stopifnot((grp==1 && length(theta)==3) || (grp>1 && length(theta)==4))
    .Call(HurdleLikelihood_method('grad'), x@pointer, theta, grp-2, penalize)
}


HurdleLikelihood_LL <- function(x, theta, grp){
    stopifnot(grp>=1 && grp<=x@p)
    stopifnot((grp==1 && length(theta)==3) || (grp>1 && length(theta)==4))
    .Call(HurdleLikelihood_method('LL'), x@pointer, theta, grp-2)
}



## syntactic sugar to allow object$method( ... )
setMethod( "$", "HurdleLikelihood", function(x, name ) {
    if(name == 'gradAll'){
        function(theta, penalize=TRUE) HurdleLikelihood_gradAll(x, theta, penalize)
    }else if(name == 'grad'){
        function(theta, grp, penalize=TRUE) HurdleLikelihood_grad(x, theta, grp, penalize)
    }else if(name == 'LL'){
        function(theta, grp) HurdleLikelihood_LL(x, theta, grp)
    }else {
        function(...) .Call(HurdleLikelihood_method(name), x@pointer, ... )
    }
}
          )

HurdleLikelihood <- function(y, x, grp, theta, lambda=0, tol=1e-4){
    if(missing(grp)) grp <- coordmap(ncol(x)+1)-2
    if(missing(theta)){
        rep(0, length(grp))
        theta[length(grp)] <- 1
    }
    if(length(theta) != length(grp)) stop('grp length not equal theta length')
    if(theta[length(grp)]<0) stop('kbb negative')
    if(length(lambda) ==1) lambda <- rep(lambda, ncol(x))
    if(length(lambda) != ncol(x)) stop('lambda must match theta length')
    new('HurdleLikelihood', pointer=.Call(HurdleLikelihood_method("new"), y, x, grp, theta, lambda, tol), p=as.integer(ncol(x)+1))
}
