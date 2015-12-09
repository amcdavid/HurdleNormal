## Mapping from parameter indices to blocks
## Only used in simulation so far...redundant with Block()?
nativeMap <- function(p){
    covariates <- rep(seq(2, p), each=2)
    block <- c(1,                       #intercept
               covariates)
    total <- c(block, block, 1) #kbb
    total
}


##' Generate a design matrix corresponding to the 2-1-0 hurdle
##'
##' .. content for \details{} ..
##' @param zif zero-inflated covariates
##' @param ... arguments passed to `scale`
##' @param tol tolerance for declaring a parameter equal to zero
##' @return design matrix
##' @export
makeModel <- function(zif, ..., tol=.001){
    M <- cbind(zif, abs(zif)>tol)
    ord <- c(seq(2, ncol(zif)+1),seq(2, ncol(zif)+1))
    M <- M[,order(ord)]
    M <- scale(M, ...)
    cbind(1, M)
}


##' Generate a design matrix corresponding to the 2-1-0 hurdle
##'
##' .. content for \details{} ..
##' @param zif zero-inflated covariates
##' @param ... arguments passed to `scale`
##' @param tol tolerance for declaring a parameter equal to zero
##' @return design matrix
##' @export
makeModelOrthog <- function(zif, ..., tol=.001){
    M <- cbind(zif, abs(zif)>tol)
    browser()
    ord <- c(seq(2, ncol(zif)+1),seq(2, ncol(zif)+1))
    M <- M[,order(ord)]
    for(i in seq_len(ncol(zif))){
        cols <- seq(2*i-1, 2*i)
        UDV <- svd(M[,cols])
        M[,cols] <- UDV$u
    }
    M <- scale(M, scale=FALSE)
    cbind(1, M)
}




##' Define parameter groups
##'
##' Define parameter groups and relative penalization.
##' Only one of \code{this.model}, \code{bidx} or \code{bset} should be provided.
##' When \code{this.model} is provided, it is assumed to have been created by with \code{makeModel},
##' and the first column is assumed to be an unpenalized intercept, followed by pairs of interaction terms.
##' The grouping corresponds to the parameter set for a node.
##' Otherwise, only one of bidx or bset should be set, but this is not implemented yet.
##' @param this.model (optional) a matrix created by
##' @param bidx NOT IMPLEMENTED starting index of each group
##' @param group \code{character}: one of components, or none.
##' @param penalty.scale \code{numeric} of length bidx/bset
##' @param bset NOT IMPLEMENTED
##' @return a list containing a data.table `map` giving the mapping between parameters, groups and penalty scales and some other components
Block <- function(this.model, bidx, group='components', penalty.scale=NULL, bset){
    ## only one of this.model, bidx and bset should be non-missing
    ## group={components, none}
    ## penalty scale should be provided as a list of indices and group={groups, components,none}

    ## Always need groups -> parameters
    ## parameters -> genes
    ## groups -> penalties

    group <- match.arg(group, c('components', 'none'))
    if(!is.null(penalty.scale)){
        if(!is.list(penalty.scale)) stop("`penalty.scale` should be list with components `scale` and `group`")
        penalty.scale.group <- match.arg(penalty.scale$group, c('components', 'none', 'groups'))
        if(!is.numeric(penalty.scale$scale)) stop("`penalty.scale$scale` should give integer indices into groups or parameters")
        penalty.scale.scale <- penalty.scale$scale
    }
    
    if(!missing(this.model)){
        p <- ncol(this.model)
        nc <- (p+1)/2
        bidxMM <- c(1, rep(2:nc, each=2))
        map <- data.table(paridx=seq_len(2*p+1), block=c(bidxMM, bidxMM, 1), mmidx=c(1:p, 1:p, NA))
        if(is.null(penalty.scale)){
            penalty.map <- data.table(block=sort(unique(map$block)), lambda=c(0, rep(1, nc-1)))
            penalty.scale.group <- 'groups'
        }
    } else{
        stop("Not implemented")
    }

    if(penalty.scale.group=='groups'){
        map <- merge(map, penalty.map, by='block')
    } else if(penalty.scale.group=='components'){
        stop('Not implemented')
    } else{
        penalty.map <- data.table(paridx=seq_along(penalty.scale.scale), lambda=penalty.scale.lambda)
        map <- merge(map, penalty.map, by='paridx')
    }

    out <- list(map=map, nonpenpar=map[lambda<=0, paridx], nonpenMM=unique(map[lambda<=0 & !is.na(mmidx), mmidx]), nonpengrp=unique(map[lambda<=0, block]), lambdablock=map[,list(lambda=lambda[1]), keyby=block]$lambda)
    class(out) <-  'Block'
    out
}

## Let K^{-1} = A^T A = U^T D U
## Solve
## argmin_x ||y-A^T x||^2 + lambda * ||x||
## Then translate back to get solution of
## argmin_x ||y-x||^2 + lambda||A^{-T} x||
projectEllipse <- function(v, lambda, d, u, control){
    r_solve_thresh <- 1e-9
    max_r_iters <- 400
    r=0
    error_r=sum(v^2/(d*r+lambda)^2)-1
    r_iter=0
    while(abs(error_r)>=r_solve_thresh && r_iter<max_r_iters){
        r_iter=r_iter+1
        if(r_iter==max_r_iters){ warning("Linesearch iteration exceeded") }
        slope=2*sum(v^2*d/(d*r+lambda)^3)
        r=error_r/slope+r
        error_r=sum(v^2/(d*r+lambda)^2)-1
    } # end while(error_r>=...) i.e. we have solved for r
    ##return(u%*%diag(r/(d*r+lambda))%*%v)
    return(u %*% diag(sqrt(d)) %*%diag(r/(d*r+lambda))%*%v)
}

##' Get a solution path for the CG model
##'
##' @param y.zif (zero-inflated) response
##' @param this.model model matrix used for both discrete and continuous linear predictors
##' @param Blocks output from \code{Block} giving the grouping/scaling for the penalization.
##' @param nlambda if `lambda` is not provided, then the number of lambda to interpolate between
##' @param lambda.min.ratio if `lambda` is not provided, then the left end of the solution path as a function of the lambda0, the lambda for the empty model
##' @param lambda penalty path desired
##' @param penaltyFactor one of `full`, `diagonal` or `identity` giving how the penalty should be scaled \emph{blockwise}
##' @param control optimization control parameters
##' @param theta (optional) initial guess for parameter
##' @param blocks object of Blocks giving blocking and block-specific penalization
##' @return matrix of parameters, one row per lambda
##' @useDynLib HurdleNormal
##' @importFrom Rcpp sourceCpp
cgpaths <- function(y.zif, this.model, Blocks=Block(this.model), nlambda=100, lambda.min.ratio=if(length(y.zif)<ncol(this.model)) .005 else .05, lambda, penaltyFactor='full', control=list(tol=1e-3, maxrounds=300, debug=1), theta){
    defaultControl <- list(tol=1e-3, maxrounds=300, debug=1, stepcontract=.5, stepsize=1, stepexpand=.1)
    nocontrol <- setdiff(names(defaultControl), names(control))
    control[nocontrol] <- defaultControl[nocontrol]
    penaltyFactor <- match.arg(penaltyFactor, c('full', 'diagonal', 'identity'))
    
    p <- ncol(this.model)
    n <- nrow(this.model)
    ## mapping from blocks to indices in theta.
    blocklist <- split(Blocks$map$paridx, Blocks$map$block)
    
    ## start at 0 (except kbb=1)
    if(missing(theta))     theta <- c(rep(0, 2*p), 1)

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


    ## Gradient of a block, holding other blocks fixed at th
    getsg <- function(th, penalize=FALSE){
        fun <- function(th2, b){
            th[blocklist[[b]]] <- th2
            hl$gradAll(th, penalize)[blocklist[[b]]]
        }
        return(fun)
    }
    
    ## likelihood for all blocks
    LLall <- function(th, penalize=TRUE) hl$LLall(th, penalize)
    ## gradient for all blocks
    gradAll <- hl$gradAll
    gamma <- control$stepsize
    
    ## value for empty sol
    interceptCol <- Blocks$nonpenMM
    interceptPar <- Blocks$nonpenpar
    theta0 <- theta[interceptPar]
    hl0 <- HurdleLikelihood(y.zif, this.model[,interceptCol,drop=FALSE], theta=theta0, lambda=0)
    o <- optim(theta0, hl0$LLall, hl0$gradAll, penalize=FALSE, method='L-BFGS-B', control=list(maxit=control$maxrounds))
    if(o$convergence !=0) stop('Empty solution failed to converge')
    theta[interceptPar] <- o$par

    ## Hessian at empty solution
    sg <- getsg(theta)
    hess <- blockHessian(theta, sg, blocklist, onlyActive=FALSE, control, fuzz=.1)
    ## np <- solve(hess[[Blocks$nonpengrp]])*min(eigen(hess[[Blocks$nonpengrp]])$values)
    sqrtPen <- eigval <- eigvec <- penMat <- hess
    for(b in seq_along(blocklist)){
        this.lambda <- Blocks$lambdablock[b]
        if(penaltyFactor=='identity'){
            hess[[b]] <- diag(1, nrow(hess[[b]]))
        } else if(penaltyFactor=='diagonal'){
            hess[[b]] <- diag(diag(hess[[b]]))
        } # else leave it alone
        
        if(this.lambda>0){
            penMat[[b]] <-  solve(hess[[b]])/this.lambda^2
            hess[[b]] <- hess[[b]]*this.lambda^2
        } else{
            penMat[[b]] <- diag(1, nrow=nrow(hess[[b]]))
        }
        eig <- eigen(penMat[[b]])
        eigval[[b]] <- eig$values
        eigvec[[b]] <- eig$vectors
        sqrtPen[[b]] <- eig$vectors %*% sqrt(diag(eig$values)) %*% t(eig$vectors)
    }

    ## Returns proximal subtheta
    proxfun <- function(b, theta, gamma, lambda){
        subtheta <- theta[blocklist[[b]]]
        if(Blocks$lambdablock[b]<=0 || lambda <=0){
            return(subtheta)
        } else {
            u <- eigvec[[b]]
            d <- eigval[[b]]
            ## pm <- penMat[[b]]
            sqrtPen <- sqrtPen[[b]]
            v <- (t(u) %*% sqrtPen %*% subtheta)
            tnorm <- sqrt(sum(v^2))
            ## stopifnot(sum(abs(tnorm - sqrt(crossprod(subtheta, pm) %*% subtheta)))<1e-7)
            if(tnorm>gamma*lambda){
                ## obj <- function(x) sum((x-subtheta)^2)/2 + gamma*lambda*sqrt(crossprod(x, hess[[b]])%*%x)
                ## gr <- function(x) (x-subtheta) + gamma*lambda*as.vector(crossprod(x, hess[[b]]))/sqrt(crossprod(x, hess[[b]])%*%x)
                ## oo <- optim(subtheta, obj, gr, method='BFGS')
                res <- projectEllipse(v, gamma*lambda, d, u)
                return(res)
            } else{
                return(subtheta * 0)
            }
        }
    }


    ## Returns kkt divergence
    kktfun <- function(b, theta, nonpengrad, lambda, flagBad0=TRUE){
        sgrad <- nonpengrad[blocklist[[b]]]
        if(Blocks$lambdablock[b]<=0){
            return(sqrt(sum(sgrad^2)))
        } else{
            subtheta <- theta[blocklist[[b]]]
            tnorm <- sqrt(crossprod(subtheta, hess[[b]]) %*% subtheta)
            if(tnorm > 0){
                pgrad <- sgrad + lambda * as.vector(hess[[b]] %*% subtheta)/tnorm
                return(sqrt(sum(pgrad^2)))
            } else{
                psubgrad <- sqrt(sum(crossprod(sgrad, sqrtPen[[b]])^2))
                if(psubgrad <lambda){
                    return(0)
                } else if(flagBad0){
                    return(-99L)
                } else{
                    return(psubgrad)
                }
            }
        }
    }

    ## estimate lambda0 given penalty
    gr <- hl$gradAll(theta)
    l0block <- sapply(seq_along(blocklist), kktfun, theta=theta, nonpengrad=hl$gradAll(theta), lambda=0, flagBad0=FALSE)
    l0 <- max(l0block)
    if(l0<0) stop('Negative lambda0 should not happen')
    if(missing(lambda)) lambda <- 10^seq(log10(1.05*l0), log10(l0*lambda.min.ratio), length.out=nlambda)
    
    ## solution path
    out <- matrix(NA, nrow=length(lambda), ncol=length(theta))
    ## kkt conditions
    kktout <- matrix(NA, nrow=length(lambda), ncol=length(blocklist))
    ## diagnostic flags
    flout <- matrix(NA, nrow=length(lambda), ncol=6)
    colnames(flout) <- c('converged', 'round', 'geval', 'feval', 'gamma', 'nnz')
    colnames(out) <- c(outer(1:p, c('D', 'C'),paste0), 'kbb')
    rownames(kktout) <- rownames(flout) <- rownames(out) <- lambda
    jerr <- rep(0, length(lambda))

    message('grad block1 ', paste(hl0$gradAll(theta[interceptPar]), collapse=','))
    

    ## loop over lambda
    for(l in seq_along(lambda)){
        hl$setLambda(lambda[l])
        sp <- solvePenProximal(theta, lambda[l], control, blocklist, LLall, gradAll, proxfun, kktfun, pre0=np, hess, gamma)
        gamma <- as.numeric(attr(sp, 'flag')['gamma'])
        theta <- out[l,] <- as.numeric(sp)
        kktout[l,] <- attr(sp, 'kkt')
        flout[l,] <- attr(sp, 'flag')
        if(control$debug>0) message('Lambda = ', round(lambda[l], 3), ' rounds = ', attr(sp, 'flag')['round'])
        activeset <- sum(kktout[l,]>0)
        if(activeset >= n/2) break
        ## update jerr?
    }
    names(blocklist) <- paste0('B', seq(2, length(blocklist)))
    list(path=out, kktout=kktout, flout=flout, blocks=blocklist, lambda=lambda)
}


blockHessian <- function(theta, sg, blocks, onlyActive=TRUE, control, fuzz=.1, pre){
    message('Update hessian')
    hess <- vector('list', length(blocks))
    for(b in seq_along(blocks)){
        hess[[b]] <- numDeriv::jacobian(sg, theta[blocks[[b]]], b=b)
        hess[[b]] <- hess[[b]]+diag(nrow(hess[[b]]))*fuzz
    }
    hess
}


## Solve the penalized function using proximal gradient descent
solvePenProximal <- function(theta, lambda, control, blocklist, LLall, gradAll, proxfun, kktfun,pre0, hess, gamma){
    feval <- geval <- 0
    ## total number of rounds
    round <- 0
    ## amount of consecutive successful (non-line search) prox's we've done
    step <- 1
    p <- length(theta)
    converged <- FALSE

    ## thetaPrime0 is previous iteration, thetaPrime1 is current and theta is proposed proximal point
    ## (used if we're extrapolating via FISTA)
    thetaPrime0 <- thetaPrime1 <- theta
    ## line search parameters
    if(gamma<=0)     gamma <- control$stepsize
    thisll <- LLall(theta, penalize=FALSE)
    gr <- gradAll(thetaPrime1, penalize=FALSE)

    kkt <- vapply(seq_along(blocklist), kktfun, NA_real_, theta=thetaPrime1, nonpengrad=gr, lambda=lambda, flagBad0=TRUE)
    
    if(control$debug>2){
        debugval <- list(pre=matrix(NA, nrow=control$maxrounds, ncol=p), thisll=rep(NA, control$maxrounds),
                         kkt=matrix(NA, nrow=control$maxrounds, ncol=length(blocklist)),
                         theta=matrix(NA, nrow=control$maxrounds, ncol=length(theta)),
                         gamma=rep(NA, control$maxrounds))
    } else{
        debugval <- NULL
    }
    
    while(round < control$maxrounds && any(abs(kkt)>control$tol)){
        round <- round+1
        ##theta[blocks[[1]]] <- thetaPrime1[blocks[[1]]] - gamma*crossprod(ihess1, gr[blocks[[1]]])
        theta <- thetaPrime1-gamma*gr
        
        for(b in seq_along(blocklist)){
            ##apply proximal operator to block b
            theta[blocklist[[b]]] <- proxfun(b, theta, gamma, lambda)
        }
        feval <- feval+1
        newll <- LLall(theta, penalize=FALSE)
        boydCondition1 <- thisll +crossprod(gr, theta-thetaPrime1)+ sum((thetaPrime1-theta)^2) /(2* gamma)

        ## Has the majorization minimum decreased the unpenalized objective?
        move <- newll <= boydCondition1

        if(control$debug>2){
            message('LL: ', newll, ifelse(move, ' <', ' >='), ' Majorization: ', boydCondition1)
        }
        if(!move){
            step <- 0
            gamma <- gamma*control$stepcontract
            if(control$debug>2) message('gamma= ', gamma, appendLF=FALSE)
            if(gamma<.0001){
                ## This seems to only happen when only the intercepts are present and we have started at a stationary point
                warning('Null step size')
                gamma <- control$stepsize
                move <- TRUE
            }
        }
        if(move) {
            thetaPrime1 <- theta #+ (round-2)/(round+3)*(theta-thetaPrime1)
            thisll <- newll
            step <- step+1
            if(step>4){
                ##try increasing stepsize by shrinking towards control$stepsize
                gamma <- gamma*(1-control$stepexpand)+control$stepsize*control$stepexpand
            }

            ##update gradient and test kkt
            ## (no need if we haven't moved thetaPrime1)
            gr <- gradAll(thetaPrime1, penalize=FALSE)
            geval <- geval+length(blocklist)
            kkt <- vapply(seq_along(blocklist), kktfun, NA_real_, theta=thetaPrime1, nonpengrad=gr, lambda=lambda, flagBad0=TRUE)

            ## Debugging/tracing
            if(control$debug>1 && (round %% 10)==0 || control$debug>2){
                pen <- sapply(1:length(blocklist), function(b){
                    st <- theta[blocklist[[b]]]
                    lambda*sqrt(crossprod(st, hess[[b]]) %*% st)
                })
                print(noquote(paste0('penll=', round(LLall(thetaPrime1, penalize=FALSE) + sum(pen), 5), ' theta= ', paste(round(thetaPrime1, 2), collapse=','), 'gamma= ', sprintf('%2.2e', gamma))))
                debugval$thisll[round] <- LLall(thetaPrime1, FALSE)+sum(pen)
           }
            if(control$debug>2){
                debugval$kkt[round,] <- kkt
                debugval$theta[round,] <- thetaPrime1
                debugval$gamma[round] <- gamma
            }

        } # end move
    } # end main loop
    structure(theta, flag=c(converged=converged, round=round, geval=geval, feval=feval, gamma=gamma, nnz=sum(kkt>0)-1), kkt=kkt, debugval=debugval)     
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
