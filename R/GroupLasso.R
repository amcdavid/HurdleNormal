defaultBix <- function(p){
    seq(2, 2*p, by=2)
}

makeModel <- function(zif, scale=TRUE, tol=.001){
    M <- cbind(zif, abs(zif)>tol)
    ord <- c(seq(2, ncol(zif)+1),seq(2, ncol(zif)+1))
    M <- M[,order(ord)]
    M <- scale(M, scale=scale)
    cbind(1, M)
}

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
cgpaths <- function(y.zif, this.model, nlambda=100, lambda, lambda.min.ratio=if(length(y.zif)<ncol(this.model)) .005 else .05, control=list(tol=1e-5, maxrounds=300, maxit=500, debug=1), standardize=FALSE, bidx=seq(2, ncol(this.model), by=2)){
    defaultControl <- list(tol=1e-6, maxrounds=300, maxit=500, debug=0, method='proximal', stepcontract=.5, stepsize=3, stepexpand=.1, updatehess=floor(300*.66))
    nocontrol <- setdiff(names(defaultControl), names(control))
    control[nocontrol] <- defaultControl[nocontrol]
    method <- match.arg(control$method, c('block', 'proximal'))
    p <- ncol(this.model)
    n <- nrow(this.model)
    ## mapping from blocks to indices in theta.  block[1] gives unpenalized parameters.
## block[p] gives last group (so actually p-1 penalized groups)
    bidx <- c(bidx, ncol(this.model)+1)
    blocks <- lapply(seq_len(length(bidx)-1), function(i){
        d <- seq(bidx[i], bidx[i+1]-1)
        c(d, d+p)
    })
                     
    ## start at 0 (except kbb=1)
    theta <- c(rep(0, 2*p), 1)

    ## get likelihood object
    ## ultimately calls C++ code in src/hurdle_likelihood.cpp
    if(var(y.zif)<1e-5) stop("Response with no variance")
    hl <- HurdleLikelihood(y.zif, this.model, theta=theta, lambda=0)
    ## profile likelihood for block b.  other blocks held fixed
    ## getsubll <- function(b, th, penalize){
    ##     fun <- function(th2){
    ##         th[blocks[[b]]] <- th2
    ##         hl$LLall(th, penalize)
    ##     }
    ##     return(fun)
    ## }
    ## profile gradient
    pen <- function(theta, lambda){
        tot <- 0
        for(b in seq_along(blocks)){
            tot <- tot+lambda*sqrt(sum(theta[blocks[[b]]]^2))
        }
        tot
    }

    gradNorm <- function(gradEval, lambda){
        tot <- rep(0, length(blocks))
        for(b in seq_along(blocks)){
            tot[b] <- sqrt(sum(gradEval[blocks[[b]]]^2))
        }
        tot
    }

    thetaNorm <- function(theta, lambda){
        tot <- rep(0, length(blocks))
        for(b in seq_along(blocks)){
            tot[b] <- sqrt(sum(theta[blocks[[b]]]^2))
        }
        tot
    }
        
    ## getsg <- function(b, th, penalize=TRUE){
    ##     fun <- function(th2){
    ##         th[blocks[[b]]] <- th2
    ##         hl$gradAll(th, penalize)[blocks[[b]]]
    ##     }
    ##     return(fun)
    ## }
    
    ## ## use to test if gradient at 0 lies within lambda
    ## sg0 <- function(b){
    ##     hl$grad(rep(0, 4), b, penalize=FALSE)
    ## }
    ## likelihood for all blocks
    LLall <- function(th, penalize=TRUE) hl$LLall(th, penalize)
    ## gradient for all blocks
    if(method=='proximal'){
        gradAll <- hl$gradAll
        gamma <- control$stepsize
    }
    
    ## value for empty sol
    interceptCol <- setdiff(seq_len(p), seq(min(bidx), max(bidx)))
    interceptPar <- c(interceptCol, interceptCol+p, length(theta))
    theta0 <- theta[interceptPar]
    hl0 <- HurdleLikelihood(y.zif, this.model[,interceptCol,drop=FALSE], theta=theta0, lambda=0)
    o <- optim(theta0, hl0$LLall, hl0$gradAll, penalize=FALSE, method='BFGS', control=list(maxit=5e4))
    if(o$convergence !=0) warning('Empty solution failed to converge')
    theta[interceptPar] <- o$par
    l0 <- max(gradNorm(hl$gradAll(theta)))
    
    if(l0<0) stop('Negative lambda0 should not happen')
    if(missing(lambda)) lambda <- 10^seq(log10(1.05*l0), log10(l0*lambda.min.ratio), length.out=nlambda)
    ## solution path
    out <- matrix(NA, nrow=length(lambda), ncol=length(theta))
    ## kkt conditions
    kktout <- matrix(NA, nrow=length(lambda), ncol=length(blocks)+1)
    ## diagnostic flags
    flout <- matrix(NA, nrow=length(lambda), ncol=6)
    colnames(flout) <- c('converged', 'round', 'geval', 'feval', 'gamma', 'nnz')
    colnames(out) <- c(outer(1:p, c('D', 'C'),paste0), 'kbb')
    rownames(kktout) <- rownames(flout) <- rownames(out) <- lambda
    jerr <- rep(0, length(lambda))

    message('grad block1 ', paste(hl0$gradAll(theta[interceptPar]), collapse=','))
    ## if(control$updatehess){
    ##     blockHessian(theta, sg, blocks, onlyActive=FALSE, control)
    ## }


    ## loop over lambda
    for(l in seq_along(lambda)){
        hl$setLambda(lambda[l])
        if(method=='block'){
            sp <-solvePen(theta, lambda[l], control, p, blocks, subll, sg, sg0, LLall)
        } else{ #proximal
            sp <- solvePenProximal(theta, lambda[l], control, blocks, LLall, gradAll, thetaNorm, pen, gamma)#, pre, ihess1, gamma)
            gamma <- as.numeric(attr(sp, 'flag')['gamma'])
        }
        theta <- out[l,] <- as.numeric(sp)
        kktout[l,] <- attr(sp, 'kkt')
        flout[l,] <- attr(sp, 'flag')
        if(control$debug>0) message('Lambda = ', round(lambda[l], 3), ' rounds = ', attr(sp, 'flag')['round'])
        activeset <- sum(kktout[l,]>0)
        if(activeset >= n/2) break
        ## update jerr?
    }
    names(blocks) <- c("(Intercepts)", colnames(this.model))
    list(path=out, kktout=kktout, flout=flout, blocks=blocks, lambda=lambda)
}


## blockHessian <- function(theta, sg, blocks, onlyActive=TRUE, control, fuzz=.1, pre){
##     message('Update hessian')
##     frame <- sys.frame(-1)
##     hess1 <- numDeriv::jacobian(sg, theta[blocks[[1]]], b=1, penalize=FALSE)
##     hess1 <- hess1+diag(nrow(hess1))*fuzz
##     sg(1, theta[blocks[[1]]], penalize=FALSE)
##     bsize <- length(blocks[[2]])
##     hess <- array(0, c(bsize, bsize, length(blocks)))
##     null <- rep(TRUE, length(blocks))
##     if(missing(pre)){
##         pre <- rep(1/max(diag(hess1)), length(blocks))
##     }
##     for(b in seq(2, length(blocks))){
##         subtheta <- theta[blocks[[b]]]
##         if(!onlyActive || sum(subtheta^2)>control$tol){
##             hess[,,b] <- numDeriv::jacobian(sg, subtheta, b=b, penalize=FALSE)
##             sg(b=b, subtheta, penalize=FALSE)
##             null[b] <- FALSE
##             pre[b] <- 1/max(diag(hess[,,b]), fuzz)
##         } else{
##             hess[,,b] <- diag(bsize)*fuzz
##         }
##     }
##     if(any(!null)){
##         hess[,,null] <- rowMeans(hess[,,!null], 2)
##     }
##     assign('ihess1', solve(hess1), frame)
##     assign('pre', pre, frame)
## }

## Solve the penalized function using proximal gradient descent
solvePenProximal <- function(theta, lambda, control, blocks, LLall, gradAll, thetaNorm, pen, gamma){
    feval <- geval <- 0
    ## total number of rounds
    round <- 0
    ## amount of consecutive successful (non-line search) prox's we've done
    step <- 1
    p <- length(blocks)
    converged <- FALSE
    kkt <- Inf
    ## thetaPrime0 is previous iteration, thetaPrime1 is current and theta is proximal proposed point
    ## (used if we're extrapolating via FISTA)
    thetaPrime0 <- thetaPrime1 <- theta
    ## line search parameters
    if(gamma<=0)     gamma <- control$stepsize
    thisll <- LLall(theta, penalize=FALSE)
    gr <- gradAll(thetaPrime1, penalize=FALSE)
    nonPen <- setdiff(seq_along(theta), unlist(blocks))
    ## preconditioner
    ## stopifnot(all(pre[2:p]>0))
    if(control$debug>2){
        debugval <- list(pre=matrix(NA, nrow=control$maxrounds, ncol=p+1), thisll=rep(NA, control$maxrounds),
                         kkt=matrix(NA, nrow=control$maxrounds, ncol=p+1),
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
        ## update intercept
        #theta[blocks[[1]]] <- thetaPrime1[blocks[[1]]] - gamma*crossprod(ihess1, gr[blocks[[1]]])
        #qapprox <- crossprod(theta[blocks[[1]]]-thetaPrime1[blocks[[1]]], hess1) %*% (theta[blocks[[1]]]-thetaPrime1[blocks[[1]]])
        ##message('grad=', paste(round(gr,2), collapse=','), '\ntheta=', paste(round(theta,2), collapse=','))
        geval <- geval+length(blocks)
        theta <- thetaPrime1-gamma*gr
        tNorm <- thetaNorm(theta, lambda)
        for(b in seq_len(p)){ ##project block onto group lasso penalty
            contract <- (1-lambda*gamma/tNorm[[b]])
            if(contract < 0){ # is 0 the minimum?
                theta[blocks[[b]]] <- 0
                if(control$debug>3) print(noquote(paste0('Zeroed block ', b, ', |subtheta|_2 = ', tNorm[[b]])))
            } else{ #otherwise, just contract a bit
                theta[blocks[[b]]] <- theta[blocks[[b]]]*contract
            }
            #qapprox <- qapprox+sum((theta[blocks[[b]]]-subtheta)^2/pre[b])
        } # end blockupdates

        newll <- LLall(theta, penalize=FALSE)
        boydCondition1 <- thisll +crossprod(gr, theta-thetaPrime1)+ sum((thetaPrime1-theta)^2) /(2* gamma)
        ## thispen <- sum(sapply(blocks[-1], function(b) sqrt(sum(thetaPrime1[b]^2))))
        ## newpen <- sum(sapply(blocks[-1], function(b) sqrt(sum(theta[b]^2))))
        move <- newll < boydCondition1
        if(!move){
            step <- 0
            gamma <- gamma*control$stepcontract
            if(gamma<.001){
                warning('Null step size')
                ##blockHessian(theta, sg, blocks, onlyActive=FALSE, control)
                gamma <- control$stepsize*control$stepcontract
                move <- TRUE
            }
            if(control$debug>2) message('gamma= ', gamma, appendLF=FALSE)
        }
        if(control$debug>2){
            message('LL: ', newll, ifelse(move, ' <', ' >='), ' Majorization: ', boydCondition1)
        }
        if(move){
            thetaPrime1 <- theta #+ (round-2)/(round+3)*(theta-thetaPrime1)
            thisll <- newll
            step <- step+1
            if(step>4){
                #try increasing stepsize by shrinking towards control$stepsize
                gamma <- gamma*(1-control$stepexpand)+control$stepsize*control$stepexpand
            }

            ##update gradient (no need if we haven't moved thetaPrime1)
            gr <- gradAll(thetaPrime1, penalize=FALSE)
            kkt <- sapply(seq_len(p), function(b){
                subtheta <- theta[blocks[[b]]]
                if(all(abs(subtheta)<control['tol'])){
                    if(sqrt(sum(gr[blocks[[b]]]^2))<lambda) return(0L) else return(-99L)
                }
                else{
                    pen <- lambda*subtheta/sqrt(sum(subtheta^2))
                    sg1 <- gr[blocks[[b]]]+pen
                    return(sum(sg1^2))
                }
            })
            kkt <- c(sum(gr[nonPen]^2), kkt)

        ## try updating hessian
        ## if(control$updatehess > 0 & (round %% control$updatehess)==0)  blockHessian(theta, sg, blocks, onlyActive=TRUE, control, pre=pre)
        
        ## message('kkt=', paste(round(kkt,2), collapse=','))
        if(control$debug>2){
#            debugval$pre[round,] <- pre
            debugval$kkt[round,] <- kkt
            debugval$theta[round,] <- thetaPrime1
            debugval$gamma[round] <- gamma
        }
        converged <- all(abs(kkt)<control$tol) || (round >= control$maxrounds)
        feval <- feval+p+1
        if(control$debug > 0){
            if(control$debug>1 && (round %% 10)==0) print(noquote(paste0('penll=', round(LLall(thetaPrime1, penalize=FALSE) + pen(thetaPrime1, lambda) , 5), ' theta= ', paste(round(thetaPrime1, 2), collapse=','), 'gamma= ', gamma )))
        }
        }
    } # end main loop

    structure(ifelse(abs(theta)<control$tol, 0, theta), flag=c(converged=converged, round=round, geval=geval, feval=feval, gamma=gamma, nnz=sum(kkt>0)-1), kkt=kkt, debugval=debugval)        
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
                    if(control$debug>3) print(noquote(paste0('Zeroed block ', b, ', grad0 = ', sqrt(sum(sg0(b)^2)))))
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
    structure(ifelse(abs(theta)<control$tol, 0, theta), flag=c(converged=converged, round=round, geval=geval, feval=feval, gamma=1, nnz=sum(kkt>0)-1), kkt=kkt)
}

plotSolPath <- function(cgpaths){
    l2 <- sapply(cgpaths$blocks, function(b){
        rowSums(cgpaths$path[,b]^2)
    })
    
    l2[l2<=0] <- NA
    l2 <- l2[,setdiff(colnames(l2), '(Intercepts)')]
    mp <- within(reshape2::melt(l2), {
    lambda <- as.numeric(X1)
    X2 <- as.character(X2)
})

    pathPlot <- ggplot(mp, aes(x=lambda, y=value, col=X2))+geom_line() + scale_x_log10()
    direct.label(pathPlot, 'last.points')
    invisible(l2)
}
