##' Get a solution path for the CG model
##'
##' @param theta initial guess for parameter
##' @param y.zif (zero-inflated) response
##' @param this.model (zero-inflated) covariates
##' @param lambda penalty path desired
##' @param control optimization control parameters
##' @return matrix of parameters, one row per lambda
##' @useDynLib HurdleNormal
##' @importFrom Rcpp sourceCpp
cgpaths <- function(y.zif, this.model, nlambda=100, lambda, lambda.min.ratio=if(length(y.zif)<ncol(this.model)) .005 else .05, control=list(tol=1e-5, maxrounds=300, maxit=500, debug=1), standardize=FALSE){
    defaultControl <- list(tol=1e-6, maxrounds=300, maxit=500, debug=0, method='proximal', stepcontract=.5, stepsize=3, stepexpand=.1, precondition=TRUE)
    nocontrol <- setdiff(names(defaultControl), names(control))
    control[nocontrol] <- defaultControl[nocontrol]
    method <- match.arg(control$method, c('block', 'proximal'))
    p <- ncol(this.model)+1
    n <- nrow(this.model)
    ## mapping from blocks to indices in theta.  block[1] gives unpenalized parameters.
    ## block[p] gives last group (so actually p-1 penalized groups)
    blocks <- sapply(seq_len(p), function(x) which(x==coordmap(p)))
    theta <- coordmap(p)
    ## start at 0 (except kbb=1)
    theta[theta!=1] <- 0

    ## get likelihood object
    ## ultimately calls C++ code in src/hurdle_likelihood.cpp
    if(var(y.zif)<1e-5) stop("Response with no variance")
    hl <- HurdleLikelihood(y.zif, this.model, theta=theta, lambda=0)
    ## profile likelihood for block b.  other blocks held fixed
    subll <- function(b, th) hl$LL(th,b)
    ## profile gradient
    sg <- function(b, th) hl$grad(th,b, penalize=TRUE)
    ## use to test if gradient at 0 lies within lambda
    sg0 <- function(b){
        hl$grad(rep(0, 4), b, penalize=FALSE)
    }
    ## likelihood for all blocks
    LLall <- function(th, penalize=TRUE) hl$LLall(th, penalize)
    ## gradient for all blocks
    if(method=='proximal') gradAll <- hl$gradAll
    
    ## value for empty sol
    gl0 <- getLambda0(theta, subll, sg, sg0, control, p, blocks)
    ## penalty
    l0 <- gl0$lambda0
    ## value of unpenalized blocks
    theta <- gl0$theta0
    if(l0<0) stop('Negative lambda0 should not happen')
    if(missing(lambda)) lambda <- 10^seq(log10(1.05*l0), log10(l0*lambda.min.ratio), length.out=nlambda)
    ## solution path
    out <- matrix(NA, nrow=length(lambda), ncol=length(theta))
    ## kkt conditions
    kktout <- matrix(NA, nrow=length(lambda), ncol=p)
    ## diagnostic flags
    flout <- matrix(NA, nrow=length(lambda), ncol=5)
    colnames(out) <- parmap(p)
    rownames(kktout) <- rownames(flout) <- rownames(out) <- lambda
    jerr <- rep(0, length(lambda))


    ## loop over lambda
    for(l in seq_along(lambda)){
        hl$setLambda(lambda[l])
        if(method=='block'){
            sp <-solvePen(theta, lambda[l], control, p, blocks, subll, sg, sg0, LLall)
            control$stepsize <- attr(sp, 'flag')['gamma']
        } else{ #proximal
            sp <- solvePenProximal(theta, lambda[l], control, p, blocks, subll, sg0, sg, LLall, gradAll)
        }
        theta <- out[l,] <- as.numeric(sp)
        kktout[l,] <- attr(sp, 'kkt')
        flout[l,] <- attr(sp, 'flag')
        if(control$debug>0) message('Lambda = ', lambda[l], 'rounds = ', attr(sp, 'flag')['round'])
        activeset <- sum(kktout[l,]>0)
        if(activeset >= n/2) break
        ## update jerr?
    }
    list(path=out, kktout=kktout, flout=flout)
}

## Get profile MLE (holding all else at 0)
## and find max block 2-norm of gradient
getLambda0 <- function(theta, subll, sg, sg0, control, p, blocks){
    subtheta <- theta[blocks[[1]]]
    o <- optim(subtheta, subll, sg, b=1, method='BFGS')
    ## ensure theta update
    o <- optim(subtheta, subll, sg, b=1, method='BFGS')
    theta[blocks[[1]]] <- o$par
    subgrad <- sapply(2:p, function(b) sqrt(sum(sg0(b)^2)))
    list(lambda0=max(subgrad), theta0=theta)
}

blockHessian <- function(theta, sg, p){
    
}

## Solve the penalized function using proximal gradient descent
solvePenProximal <- function(theta, lambda, control, p, blocks, subll, sg0, sg, LLall, gradAll){
    feval <- geval <- 0
    ## total number of rounds
    round <- 0
    ## amount of consecutive successful (non-line search) prox's we've done
    step <- 1
    converged <- FALSE
    kkt <- Inf
    ## thetaPrime0 is previous iteration, thetaPrime1 is current and theta is proximal proposed point
    ## (used if we're extrapolating via FISTA)
    thetaPrime0 <- thetaPrime1 <- theta
    ## line search parameters
    gamma <- control$stepsize
    hess <- numDeriv::jacobian(gradAll, thetaPrime1, penalize=FALSE)
    ihess0 <- solve(hess[blocks[[1]],blocks[[1]]] +diag(3)*.1) ## use exact hessian plus fuzz for intercept
    hess0 <- solve(ihess0)
    thisll <- LLall(theta, penalize=TRUE)
    ## preconditioner
    if(control$precondition){
            pre <- c(0, sapply(2:p, function(b){
                1/max(c(diag(hess[blocks[[b]], blocks[[b]]]),.01))
            }))} else{
    pre <- rep(1, p)
}
    stopifnot(all(pre[2:p]>0))
    if(control$debug>2){
        debugval <- list(pre=matrix(NA, nrow=control$maxrounds, ncol=p), thisll=rep(NA, control$maxrounds),
                         kkt=matrix(NA, nrow=control$maxrounds, ncol=p),
                         theta=matrix(NA, nrow=control$maxrounds, ncol=length(theta)),
                         gamma=rep(NA, control$maxrounds))
    } else{
        debugval <- NULL
    }
    while(!converged){
        round <- round+1
        ##solve unpenalized block exactly (because of stepsize awkwardness)
        ## th0 <- thetaPrime[blocks[[1]]]        
        ## O <- optim(th0, subll, sg, method='BFGS', b=1)
        ## feval <- feval+O$counts['function']
        ## geval <- geval+O$counts['gradient']
        ## thetaPrime[blocks[[1]]] <- O$par
        ## thisll <- LLall(thetaPrime, penalize=FALSE)
        if(step>0){
            ##update gradient (no need if we haven't moved thetaPrime1)
            gr <- gradAll(thetaPrime1, penalize=FALSE)
        }
        ## update intercept
        theta[blocks[[1]]] <- thetaPrime1[blocks[[1]]] - gamma*crossprod(ihess0, gr[blocks[[1]]])
        #qapprox <- crossprod(theta[blocks[[1]]]-thetaPrime1[blocks[[1]]], hess0) %*% (theta[blocks[[1]]]-thetaPrime1[blocks[[1]]])
        ##message('grad=', paste(round(gr,2), collapse=','), '\ntheta=', paste(round(theta,2), collapse=','))
        geval <- geval+length(blocks)
        for(b in seq(2, p)){ ##project block onto group lasso penalty
            subtheta <- thetaPrime1[blocks[[b]]]-gamma*pre[b]*gr[blocks[[b]]]
            suml2 <- sqrt(sum(subtheta^2))
            contract <- (1-lambda*gamma*pre[b]/suml2)
            if(contract < 0){ # is 0 the minimum?
                theta[blocks[[b]]] <- 0
                if(control$debug>2) print(noquote(paste0('Zeroed block ', b, ', |subtheta|_2 = ', suml2)))
            } else{ #otherwise, just contract a bit
                theta[blocks[[b]]] <- subtheta*contract
            }
            #qapprox <- qapprox+sum((theta[blocks[[b]]]-subtheta)^2/pre[b])
        } # end blockupdates

        newll <- LLall(theta, penalize=TRUE)
        diffll <- newll-thisll
        
        ## boydCondition1 <- diffll -crossprod(gr, theta-thetaPrime1)-
        ##                   qapprox/(2* gamma)
        ## thispen <- sum(sapply(blocks[-1], function(b) sqrt(sum(thetaPrime1[b]^2))))
        ## newpen <- sum(sapply(blocks[-1], function(b) sqrt(sum(theta[b]^2))))
        
        boydCondition1 <- diffll
        move <- diffll<0
        if(!move){
            step <- 0
            gamma <- gamma*control$stepcontract
            if(gamma<.001){
                warning('Null step size')
                gamma <- .3
                move <- TRUE
            }
            if(control$debug>2) message('gamma= ', gamma)
        }
        if(move){
            thetaPrime1 <- theta #+ (round-2)/(round+3)*(theta-thetaPrime1)
            thisll <- newll
            step <- step+1
            if(step>4){
                #try increasing stepsize by shrinking towards control$stepsize
                gamma <- gamma*(1-control$stepexpand)+control$stepsize*control$stepexpand
            }
            
        }
        
        kkt <- sapply(seq_len(p), function(b){
            subtheta <- theta[blocks[[b]]]
            if(all(abs(subtheta)<control['tol'])){
                if(sqrt(sum(gr[blocks[[b]]]^2))<lambda) return(0L) else return(-99L)
            }
            else{
                return(sum(sg(b, subtheta)^2))
            }
        })

        
        ## message('kkt=', paste(round(kkt,2), collapse=','))
        if(control$debug>2){
            debugval$pre[round,] <- pre
            debugval$kkt[round,] <- kkt
            debugval$theta[round,] <- thetaPrime1
            debugval$gamma[round] <- gamma
        }
        converged <- all(abs(kkt)<control$tol) || (round >= control$maxrounds)
        feval <- feval+p+1
        if(control$debug > 0){
            if(control$debug>1 || (round %% 10)==0) print(noquote(paste0('penll=', round(LLall(thetaPrime1, penalize=TRUE), 5), ' theta= ', paste(round(thetaPrime1, 2), collapse=','), 'gamma= ', gamma )))
        }
    } # end main loop

    structure(ifelse(abs(theta)<control$tol, 0, theta), flag=c(converged=converged, round=round, geval=geval, feval=feval, gamma=gamma), kkt=kkt, debugval=debugval)        
}

solvePen <- function(theta, lambda, control, p, blocks,  subll, sg, sg0, LLall, loc, scal){
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

        kkt <- sapply(seq_len(p), function(b){
            subtheta <- theta[blocks[[b]]]
            if(all(abs(subtheta)<control['tol'])){
                if(sqrt(sum(sg0(b)^2))<lambda) return(0L) else return(-99L)
            } 
            else return(sum(sg(b, subtheta)^2))
        })

        converged <- all(abs(kkt)<control$tol) || (round > control$maxrounds)
        if(control$debug > 0){
            if(control$debug>1 || (round %% 10)==0) print(noquote(paste0('penll=', round(LLall(theta), 4), ' theta= ', paste(round(theta, 2), collapse=','))))
        }
    } # end main loop

        geval <- geval+p
    structure(ifelse(abs(theta)<control$tol, 0, theta), flag=c(converged=converged, round=round, geval=geval, feval=feval, gamma=1), kkt=kkt)
}
