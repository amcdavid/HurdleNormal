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
cgpaths <- function(y.zif, this.model, nlambda=100, lambda, lambda.min.ratio=if(length(y.zif)<ncol(this.model)) .005 else .05, control=list(tol=1e-5, maxrounds=300, maxit=500, debug=1)){
    defaultControl <- list(tol=1e-6, maxrounds=300, maxit=500, debug=0)
    nocontrol <- setdiff(names(defaultControl), names(control))
    control[nocontrol] <- defaultControl[nocontrol]
    p <- ncol(this.model)+1
    n <- nrow(this.model)
    blocks <- sapply(seq_len(p), function(x) which(x==coordmap(p)))
    theta <- coordmap(p)
    theta[theta!=1] <- 0

    ## get likelihood object
    hl <- HurdleLikelihood(y.zif, this.model, theta=theta, lambda=0)
    subll <- function(b, th) hl$LL(th,b)
    sg <- function(b, th) hl$grad(th,b, penalize=TRUE)
    ## use to test if gradient at 0 lies within lambda
    sg0 <- function(b){
        hl$grad(rep(0, 4), b, penalize=FALSE)
    }
    LLall <- function(th) hl$LLall(th)
    ## value for empty sol
    l0 <- getLambda0(theta, subll, sg, sg0, control, p, blocks)
    if(missing(lambda)) lambda <- 10^seq(log10(1.05*l0), log10(l0*lambda.min.ratio), length.out=nlambda)
 #                                            
  #                                           if(p<n) 3 else (n/p), length.out=nlambda)
    out <- matrix(NA, nrow=length(lambda), ncol=length(theta))
    kktout <- matrix(NA, nrow=length(lambda), ncol=p)
    flout <- matrix(NA, nrow=length(lambda), ncol=4)
    colnames(out) <- parmap(p)
    rownames(kktout) <- rownames(flout) <- rownames(out) <- lambda
    jerr <- rep(0, length(lambda))

    
    for(l in seq_along(lambda)){
        hl$setLambda(lambda[l])
        sp <-solvePen(theta, lambda[l], control, p, blocks, subll, sg, sg0, LLall)      
        theta <- out[l,] <- sp
        kktout[l,] <- attr(sp, 'kkt')
        flout[l,] <- attr(sp, 'flag')
        if(control$debug>0) message('Lambda = ', lambda[l], 'rounds = ', attr(sp, 'flag')['round'])
        activeset <- sum(kktout[l,]>0)
        if(activeset >= n/2) break
        ## update jerr?
    }
    list(path=out, kktout=kktout, flout=flout)
}

getLambda0 <- function(theta, subll, sg, sg0, control, p, blocks){
    subtheta <- theta[blocks[[1]]]
    o <- optim(subtheta, subll, sg, b=1, method='BFGS')
    o <- optim(subtheta, subll, sg, b=1, method='BFGS')
    theta[blocks[[1]]] <- o$par
    subll(1, theta[blocks[[1]]])
    subgrad <- sapply(2:p, function(b) sqrt(sum(sg0(b)^2)))
    max(subgrad)
}

solvePen <- function(theta, lambda, control, p, blocks,  subll, sg, sg0, LLall){
    feval <- geval <- 0
    round <- 0
    converged <- FALSE
    while(!converged){
        round <- round+1
        theta0 <- theta
        for(b in seq_len(p)){ #blockupdates
            ## test if setting block to 0 is minimum
            subtheta <- theta[blocks[[b]]]
            if(b != 1){
                geval <- geval +1
                if(sqrt(sum(sg0(b)^2)) < lambda){
                    theta[blocks[[b]]] <- 0
                    if(control$debug>2) print(noquote(paste0('Zeroed block ', b, ', grad0 = ', sqrt(sum(sg0(b)^2)))))
                    next
                } else if(sum(subtheta^2)<control$tol){
                    subtheta[] <- control$tol
                }
            }
            oo <- optim(subtheta, subll, sg, b=b, method='BFGS', control=control['maxit'])
            feval <- feval+oo$counts['function']
            geval <- geval+oo$counts['gradient']
            theta[blocks[[b]]] <- oo$par
        } # end blockupdates
        converged <- mean((theta-theta0)^2)<(control$tol)^2 || round > control$maxrounds
        if(control$debug > 0){
            if(control$debug>1 || (round %% 10)==0) print(noquote(paste0('penll=', round(LLall(theta), 4), ' theta= ', paste(round(theta, 2), collapse=','))))
        }
    } # end main loop

    kkt <- sapply(seq_len(p), function(b){
        subtheta <- theta[blocks[[b]]]
        if(all(abs(subtheta)<control['tol'])){
            if(sqrt(sum(sg0(b)^2))<lambda) return(0L) else return(-99L)
        } 
        else return(sum(sg(b, subtheta)^2))
    })
        geval <- geval+p
    structure(ifelse(abs(theta)<control$tol, 0, theta), flag=c(converged=converged, round=round, geval=geval, feval=feval), kkt=kkt)
}
