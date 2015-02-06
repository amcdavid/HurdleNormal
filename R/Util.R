reraise <- function(err, convertToWarning=FALSE, silent=FALSE){
    if(silent) return(err)
    if(convertToWarning){
        warning(simpleWarning(message=err$message, call=err$call))
    } else{
        stop(simpleError(message=err$message, call=err$call))
    }
    return(err)
}

##' Selectively muffle warnings based on output
##'
##' @param expr an expression
##' @param regexp a regexp to be matched (with str_detect)
##' @return the result of expr
##' @import stringr
##' @export
hushWarning <- function(expr, regexp){
    withCallingHandlers(expr, warning=function(w){
        if(str_detect(conditionMessage(w), regexp)) invokeRestart("muffleWarning")
    })
}
