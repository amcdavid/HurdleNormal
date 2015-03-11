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
cgpaths <- function(y.zif, this.model, nlambda=100, lambda, lambda.min.ratio=if(length(y.zif)<ncol(this.model)) .005 else .05, control=list(tol=1e-5, maxrounds=300, maxit=500, debug=1), standardize=FALSE){
    defaultControl <- list(tol=1e-6, maxrounds=300, maxit=500, debug=0, method='proximal', stepsize=.6)
    nocontrol <- setdiff(names(defaultControl), names(control))
    control[nocontrol] <- defaultControl[nocontrol]
    method <- match.arg(control$method, c('block', 'proximal'))
    p <- ncol(this.model)+1
    n <- nrow(this.model)
    blocks <- sapply(seq_len(p), function(x) which(x==coordmap(p)))
    theta <- coordmap(p)
    theta[theta!=1] <- 0

    scal <- rep(1, p)
    loc <- rep(0, p)
    cond.scale <- function(x){
        xI <- abs(x)>control$tol
        if(sum(xI)<2){
            loc <- 0
            scale <- 1             
        } else{
            loc <- mean(x[xI])
            scale <- sd(x[xI])
        }
        list(z=(x-loc)/scale, loc=loc, scale=scale)
    }
    
    if(standardize){
        s <- cond.scale(y.zif)
        y.zif <-s$z
        scal[1] <- s$scale
        loc[1] <- s$loc
        for(i in seq(2,p)){
            s <- cond.scale(this.model[,i-1])
            this.model[,i-1] <- s$z
            scal[i] <- s$scale
            loc[i] <- s$loc
        }
    }

    ## get likelihood object
    hl <- HurdleLikelihood(y.zif, this.model, theta=theta, lambda=0)
    subll <- function(b, th) hl$LL(th,b)
    sg <- function(b, th) hl$grad(th,b, penalize=TRUE)
    ## use to test if gradient at 0 lies within lambda
    sg0 <- function(b){
        hl$grad(rep(0, 4), b, penalize=FALSE)
    }
    LLall <- function(th, penalize=TRUE) hl$LLall(th, penalize)

    if(method=='proximal') gradAll <- hl$gradAll
    
    ## value for empty sol
    gl0 <- getLambda0(theta, subll, sg, sg0, control, p, blocks)
    l0 <- gl0$lambda0
    theta <- gl0$theta0
    if(l0<0) stop('Negative lambda0 should not happen')
    if(missing(lambda)) lambda <- 10^seq(log10(1.05*l0), log10(l0*lambda.min.ratio), length.out=nlambda)
 #                                            
  #                                           if(p<n) 3 else (n/p), length.out=nlambda)
    out <- matrix(NA, nrow=length(lambda), ncol=length(theta))
    kktout <- matrix(NA, nrow=length(lambda), ncol=p)
    flout <- matrix(NA, nrow=length(lambda), ncol=5)
    colnames(out) <- parmap(p)
    rownames(kktout) <- rownames(flout) <- rownames(out) <- lambda
    jerr <- rep(0, length(lambda))

    
    for(l in seq_along(lambda)){
        hl$setLambda(lambda[l])
        if(method=='block'){
            sp <-solvePen(theta, lambda[l], control, p, blocks, subll, sg, sg0, LLall, loc, scal)      } else{ #proximal
                sp <- solvePenProximal(theta, lambda[l], control, p, blocks, subll, sg0, sg, LLall, gradAll, loc, scal)
        }
        spnull <- sp
        attributes(spnull) <- NULL
        theta <- out[l,] <- spnull
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
    theta[blocks[[1]]] <- o$par
    ## ensure theta update
    subll(1, theta[blocks[[1]]])
    subgrad <- sapply(2:p, function(b) sqrt(sum(sg0(b)^2)))
    list(lambda0=max(subgrad), theta0=theta)
}

solvePenProximal <- function(theta, lambda, control, p, blocks, subll, sg0, sg, LLall, gradAll, loc, scal){
    feval <- geval <- 0
    step <- round <- 0
    converged <- FALSE
    kkt <- Inf
    ## thetaPrime0 is previous iteration, thetaPrime1 is current and theta is proximal proposed point
    thetaPrime0 <- thetaPrime1 <- theta
    ## line search parameters
    beta <- beta0 <- 1
    betaw <- .95
    while(!converged){
        round <- round+1
        thisll <- LLall(thetaPrime1, penalize=FALSE)
        ##theta0 <- theta
        ##solve unpenalized block exactly (because of stepsize awkwardness)
        ## th0 <- thetaPrime[blocks[[1]]]        
        ## O <- optim(th0, subll, sg, method='BFGS', b=1)
        ## feval <- feval+O$counts['function']
        ## geval <- geval+O$counts['gradient']
        ## thetaPrime[blocks[[1]]] <- O$par
        ## thisll <- LLall(thetaPrime, penalize=FALSE)
        gr <- gradAll(thetaPrime1, penalize=FALSE)
        ## stopifnot(mean(gr[blocks[[1]]]^2)<mean(gr^2))
        theta <- thetaPrime1-beta*gr
        ##message('grad=', paste(round(gr,2), collapse=','), '\ntheta=', paste(round(theta,2), collapse=','))
        geval <- geval+length(blocks)
        for(b in seq(2, p)){ ##project block onto group lasso penalty
            subtheta <- theta[blocks[[b]]]
            suml2 <- sqrt(sum(subtheta^2))
            if(suml2 < lambda*beta){ # is 0 the minimum?
                theta[blocks[[b]]] <- 0
                if(control$debug>2) print(noquote(paste0('Zeroed block ', b, ', |subtheta|_2 = ', suml2)))
            } else{ #otherwise, just contract a bit
                theta[blocks[[b]]] <- subtheta*(1-lambda*beta/suml2)
            }
            
        } # end blockupdates

        newll <- LLall(theta, penalize=FALSE)
        diffll <- newll-thisll
        boydCondition1 <- diffll -(
            crossprod(gr, theta-thetaPrime1) +
            crossprod(theta-thetaPrime1)/(2*beta))
        if(boydCondition1>control$tol/10){
            step <- 0
            beta <- beta*control$stepsize
            if(beta<.001){
                warning('Null step size')
                browser()
            }
            if(control$debug>2) message('beta= ', beta)
        } else{
            thetaPrime1 <- theta# + (round-2)/(round+1)*(theta-thetaPrime1)
            thisll <- newll
            step <- step+1
            if(step>4){
                beta <- beta*betaw+beta0*(1-betaw)
            }
            
        }
        
        kkt <- sapply(seq_len(p), function(b){
            subtheta <- theta[blocks[[b]]]
            if(all(abs(subtheta)<control['tol'])){
                if(sqrt(sum(gr[blocks[[b]]]^2))<lambda) return(0L) else return(-99L)
            }
            else return(sum(gr[blocks[[b]]]^2))
        })

        
        ## message('kkt=', paste(round(kkt,2), collapse=','))
        converged <- all(abs(kkt)<control$tol) || (round > control$maxrounds)
#        thetaExtrap <- theta + (theta-theta0)*(round-1)/(round+10)
        feval <- feval+p+1
        if(control$debug > 0){
            if(control$debug>1 || (round %% 10)==0) print(noquote(paste0('penll=', round(LLall(thetaPrime1, penalize=TRUE), 5), ' theta= ', paste(round(thetaPrime1, 2), collapse=','), 'beta= ', beta )))
        }
    } # end main loop

    structure(ifelse(abs(theta)<control$tol, 0, theta), flag=c(converged=converged, round=round, geval=geval, feval=feval, beta=beta), kkt=kkt)        
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
    structure(ifelse(abs(theta)<control$tol, 0, theta), flag=c(converged=converged, round=round, geval=geval, feval=feval, beta=1), kkt=kkt)
}
