setClass("HurdleLikelihood", representation( pointer = "externalptr" ) )

HurdleLikelihood_method <- function(name) {
paste( "HurdleLikelihood", name, sep = "__" )
 }

HurdleLikelihood_gradAll <- function(x, theta, penalize=TRUE){
    grp <- coordmap((length(theta)+1)/4)-2
    grad <- theta
    for(g in unique(grp)){
        grad[grp==g] <- .Call(HurdleLikelihood_method('grad'), x@pointer, theta[grp==g], g, penalize)
    }
    grad
}

## syntactic sugar to allow object$method( ... )
setMethod( "$", "HurdleLikelihood", function(x, name ) {
    if(name == 'gradAll'){
        function(theta) HurdleLikelihood_gradAll(x, theta)
    }else{
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
    new('HurdleLikelihood', pointer=.Call(HurdleLikelihood_method("new"), y, x, grp, theta, lambda, tol))
}
