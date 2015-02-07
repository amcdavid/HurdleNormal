## for lambda = 0, ..., large
## gll, grad = generatelogLik(lambda, data)
## While not converged(theta):
## ## For blocks b=1,...,P
## ## ## If b != 1, test if 0 is minimizer.
## ## ## ## If so, next.
## ## ## ## If not, make sure we're at a differentiable point
## ## ## Let gb = Curry(gll, -b)
## ## ## theta(b) = Solve gb
## Done


##' Get a solution path for the CG model
##'
##' @param theta initial guess for parameter
##' @param y.zif (zero-inflated) response
##' @param this.model covariates
##' @param lambda penalty path desired
##' @param control optimization control parameters
##' @return matrix of parameters, one row per lambda
##' @useDynLib HurdleNormal
##' @importFrom Rcpp sourceCpp
cgpaths <- function(theta, y.zif, this.model, lambda=exp(seq(1, -5, length.out=10)), control=list(tol=1e-6, maxrounds=300, maxit=500, debug=1)){
    out <- matrix(NA, nrow=length(lambda), ncol=length(theta))
    p <- (length(theta)+1)/4
    colnames(out) <- parmap(p)
    rownames(out) <- lambda
    jerr <- rep(0, length(lambda))
    for(l in seq_along(lambda)){
        sp <-solvePen(theta, lambda[l], y.zif, this.model, control)
        out[l,] <- sp
        if(control$debug>0) message('Lambda = ', lambda[l], 'rounds = ', attr(sp, 'flag')['round'])
        ## update jerr?
    }
    out
}

solvePen <- function(theta, lambda, y.zif, this.model, control){
    hl <- HurdleLikelihood(y.zif, this.model, theta=theta, lambda=lambda)
    subll <- function(b, th) hl$LL(th,b)
    sg <- function(b, th) hl$grad(th,b, penalize=TRUE)
    grad <- hl$gradAll
    ## use to test if gradient at 0 lies within lambda
    sg0 <- function(b){
        hl$grad(rep(0, 4), b, penalize=FALSE)
    }
    p <- (length(theta)+1)/4   #actually p + 1
    blocks <- sapply(seq_len(p), function(x) which(x==coordmap(p)))
    round <- 0
    converged <- FALSE
    while(!converged){
        round <- round+1
        theta0 <- theta
        for(b in seq_len(p)){ #blockupdates
            ## test if setting block to 0 is minimum
            subtheta <- theta[blocks[[b]]]
            if(b != 1){
                if(sqrt(sum(sg0(b)^2)) < lambda){
                    theta[blocks[[b]]] <- 0
                    if(control$debug>2) print(noquote(paste0('Zeroed block ', b, ', grad0 = ', sqrt(sum(sg0()^2)))))
                    next
                } else if(sum(subtheta^2)<control$tol){
                    subtheta[] <- control$tol
                }
            }
            oo <- optim(subtheta, subll, sg, b=b, method='BFGS', control=control['maxit'])
            theta[blocks[[b]]] <- oo$par
        } # end blockupdates
        converged <- mean((theta-theta0)^2)<control$tol || round > control$maxrounds
        if(control$debug > 0){
            if(control$debug>1 || (round %% 10)==0) print(noquote(paste0('penll=', round(hl$LLall(theta), 4), ' theta= ', paste(round(theta, 2), collapse=','))))
        }
    } # end main loop
    structure(ifelse(abs(theta)<control$tol, 0, theta), flag=c(converged=converged, round=round))
}
