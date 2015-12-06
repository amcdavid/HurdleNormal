nativeMap <- function(p){
    covariates <- rep(seq(2, p), each=2)
    block <- c(1,                       #intercept
               covariates)
    total <- c(block, block, 1) #kbb
    total
}



makeModel <- function(zif, ..., tol=.001){
    M <- cbind(zif, abs(zif)>tol)
    ord <- c(seq(2, ncol(zif)+1),seq(2, ncol(zif)+1))
    M <- M[,order(ord)]
    M <- scale(M, ...)
    cbind(1, M)
}

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
##' @param blocks object of Blocks giving blocking and block-specific penalization
##' @param nlambda if `lambda` is not provided, then the number of lambda to interpolate between
##' @param lambda.min.ratio if `lambda` is not provided, then the left end of the solution path as a function of the lambda0, the lambda for the empty model
##' @param lambda penalty path desired
##' @param control optimization control parameters
##' @param theta (optional) initial guess for parameter
##' @return matrix of parameters, one row per lambda
##' @useDynLib HurdleNormal
##' @importFrom Rcpp sourceCpp
cgpaths <- function(y.zif, this.model, Blocks=Block(this.model), nlambda=100, lambda.min.ratio=if(length(y.zif)<ncol(this.model)) .005 else .05, lambda, penaltyFactor='full', control=list(tol=1e-5, maxrounds=300, maxit=500, debug=1), theta){
    defaultControl <- list(tol=1e-6, maxrounds=300, maxit=500, debug=0, method='proximal', stepcontract=.5, stepsize=3, stepexpand=.1, updatehess=floor(300*.66))
    nocontrol <- setdiff(names(defaultControl), names(control))
    control[nocontrol] <- defaultControl[nocontrol]
    method <- match.arg(control$method, c('block', 'proximal'))
    p <- ncol(this.model)
    n <- nrow(this.model)
    ## mapping from blocks to indices in theta.  block[1] gives unpenalized parameters.
## block[p] gives last group (so actually p-1 penalized groups)
    blocks <- split(Blocks$map$paridx, Blocks$map$block)
                     
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

  
    getsg <- function(th, penalize=FALSE){
        fun <- function(th2, b){
            th[blocks[[b]]] <- th2
            hl$gradAll(th, penalize)[blocks[[b]]]
        }
        return(fun)
    }
    
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
    interceptCol <- Blocks$nonpenMM
    interceptPar <- Blocks$nonpenpar
    theta0 <- theta[interceptPar]
    hl0 <- HurdleLikelihood(y.zif, this.model[,interceptCol,drop=FALSE], theta=theta0, lambda=0)
    o <- optim(theta0, hl0$LLall, hl0$gradAll, penalize=FALSE, method='BFGS', control=list(maxit=5e4))
    if(o$convergence !=0) warning('Empty solution failed to converge')
    theta[interceptPar] <- o$par

    sg <- getsg(theta)
    hess <- blockHessian(theta, sg, blocks, onlyActive=FALSE, control, fuzz=.1)
    np <- solve(hess[[Blocks$nonpengrp]])*min(eigen(hess[[Blocks$nonpengrp]])$values)
    sqrtPen <- eigval <- eigvec <- penMat <- hess
    for(b in seq_along(blocks)){
        this.lambda <- Blocks$lambdablock[b]
        ##hess[[b]] <- hess[[b]]+diag(1, nrow(hess[[b]]))
        ## Diagonalize
        hb <- matrix(0, nrow=nrow(hess[[b]]), ncol=ncol(hess[[b]]))
        diag(hb) <- diag(hess[[b]])
        if(penaltyFactor=='identity'){
            hess[[b]] <- diag(1, nrow(hess[[b]]))
        } else if(penaltyFactor=='diagonal'){
            hess[[b]] <- hb
        }
        
        ## End Diagonalize
        if(this.lambda>0){
            penMat[[b]] <-  solve(hess[[b]])/this.lambda^2
        } else{
            penMat[[b]] <- diag(1, nrow=nrow(hess[[b]]))
        }
        hess[[b]] <- hess[[b]]*this.lambda^2  #switcheroo
        eig <- eigen(penMat[[b]])
        eigval[[b]] <- eig$values
        eigvec[[b]] <- eig$vectors
        sqrtPen[[b]] <- eig$vectors %*% sqrt(diag(eig$values)) %*% t(eig$vectors)
    }

        ## Returns proximal subtheta
    proxfun <- function(b, theta, gamma, lambda){
        subtheta <- theta[blocks[[b]]]
        if(Blocks$lambdablock[b]<=0 || lambda <=0){
            return(subtheta)
        } else {
            u <- eigvec[[b]]
            d <- eigval[[b]]
            pm <- penMat[[b]]
            ## chl <- chol(penMat[[b]])
            sqrtPen <- sqrtPen[[b]]
            v <- (t(u) %*% sqrtPen %*% subtheta)
            tnorm <- sqrt(sum(v^2))
            stopifnot(sum(abs(tnorm - sqrt(crossprod(subtheta, pm) %*% subtheta)))<1e-7)
            if(tnorm>gamma*lambda){
                ## obj <- function(x) sum((x-subtheta)^2)/2 + gamma*lambda*sqrt(crossprod(x, hess[[b]])%*%x)
                ## gr <- function(x) (x-subtheta) + gamma*lambda*as.vector(crossprod(x, hess[[b]]))/sqrt(crossprod(x, hess[[b]])%*%x)

                ## oo <- optim(subtheta, obj, gr, method='BFGS')
                res <- projectEllipse(v, gamma*lambda, d, u)
                ## stopifnot(all(solve(chl, oo2$par)==oo$par))
                return(res)
            } else{
                return(subtheta * 0)
            }
        }
    }


    ## Returns kkt divergence
    kktfun <- function(b, theta, nonpengrad, lambda, flagBad0=TRUE){
        sgrad <- nonpengrad[blocks[[b]]]
        if(Blocks$lambdablock[b]<=0){
            return(sqrt(sum(sgrad^2)))
        } else{
            subtheta <- theta[blocks[[b]]]
            tnorm <- sqrt(crossprod(subtheta, hess[[b]]) %*% subtheta)
            if(tnorm > 0){
                pgrad <- sgrad + lambda * as.vector(hess[[b]] %*% subtheta)/tnorm
                #browser()
                return(sqrt(sum(pgrad^2)))
            } else{
                psubgrad <- sqrt(crossprod(sgrad, penMat[[b]]) %*% sgrad)
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
    ## return penalty matrix and kkt function (how related??)
    gr <- hl$gradAll(theta)
    l0block <- sapply(seq_along(blocks), kktfun, theta=theta, nonpengrad=hl$gradAll(theta), lambda=0, flagBad0=FALSE)
    l0 <- max(l0block)
    if(l0<0) stop('Negative lambda0 should not happen')
    if(missing(lambda)) lambda <- 10^seq(log10(1.05*l0), log10(l0*lambda.min.ratio), length.out=nlambda)
    ## solution path
    out <- matrix(NA, nrow=length(lambda), ncol=length(theta))
    ## kkt conditions
    kktout <- matrix(NA, nrow=length(lambda), ncol=length(blocks))
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
        if(method=='block'){
            sp <-solvePen(theta, lambda[l], control, p, Blocks, subll, sg, sg0, LLall)
        } else{ #proximal
            sp <- solvePenProximal(theta, lambda[l], control, Blocks, LLall, gradAll, proxfun, kktfun, pre0=np, hess, gamma)
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
    blocks <- c(interceptCol, blocks)
    names(blocks) <- c("(Intercepts)", paste0('B', seq(2, length(blocks))))
    list(path=out, kktout=kktout, flout=flout, blocks=blocks, lambda=lambda)
}


## need to minimize 1/2 * tcrossprod(y %*% Linv) - p%*%y, subject to crossprod(y) <= 1
## Pope https://tcg.mae.cornell.edu/pubs/Pope_FDA_08.pdf page 14

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
solvePenProximal <- function(theta, lambda, control, Blocks, LLall, gradAll, proxfun, kktfun,pre0, hess, gamma){
    feval <- geval <- 0
    ## total number of rounds
    round <- 0
    ## amount of consecutive successful (non-line search) prox's we've done
    step <- 1
    p <- length(theta)
    converged <- FALSE
    kkt <- Inf
    ## thetaPrime0 is previous iteration, thetaPrime1 is current and theta is proximal proposed point
    ## (used if we're extrapolating via FISTA)
    thetaPrime0 <- thetaPrime1 <- theta
    ## line search parameters
    if(gamma<=0)     gamma <- control$stepsize
    thisll <- LLall(theta, penalize=FALSE)
    gr <- gradAll(thetaPrime1, penalize=FALSE)
    blocks <- split(Blocks$map$paridx, Blocks$map$block)
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
        ## np <- blocks[[Blocks$nonpengrp]]
        ##  theta[np] <- thetaPrime1[np]-tcrossprod(gamma*gr[np], pre0)
        ##  theta[-np] <- thetaPrime1[-np]-gamma*gr[-np]
        theta <- thetaPrime1-gamma*gr
        
        for(b in seq_along(blocks)){ ##project block onto group lasso penalt
            theta[blocks[[b]]] <- proxfun(b, theta, gamma, lambda)
        }

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
            kkt <- vapply(seq_along(blocks), kktfun, NA_real_, theta=thetaPrime1, nonpengrad=gr, lambda=lambda, flagBad0=TRUE)

                ## try updating hessian
                ## if(control$updatehess > 0 & (round %% control$updatehess)==0)  blockHessian(theta, sg, blocks, onlyActive=TRUE, control, pre=pre)
                
                ## message('kkt=', paste(round(kkt,2), collapse=','))
                if(control$debug>2){
                                        #            debugval$pre[round,] <- pre
                    debugval$kkt[round,] <- kkt
                    debugval$theta[round,] <- thetaPrime1
                    debugval$gamma[round] <- gamma
                }
                converged <- all(abs(kkt)<control$tol)
                feval <- feval+p+1
                if(control$debug > 0){
                    if(control$debug>1 && (round %% 10)==0 || control$debug>2){
                        pen <- sapply(1:length(blocks), function(b){
                            st <- theta[blocks[[b]]]
                            lambda*sqrt(crossprod(st, hess[[b]]) %*% st)
                        })
                        print(noquote(paste0('penll=', round(LLall(thetaPrime1, penalize=FALSE) + sum(pen), 5), ' theta= ', paste(round(thetaPrime1, 2), collapse=','), 'gamma= ', sprintf('%2.2e', gamma))))
                        debugval$thisll[round] <- LLall(thetaPrime1, FALSE)+sum(pen)
                    }
                }
        }
        converged <- converged || (round >= control$maxrounds)
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
