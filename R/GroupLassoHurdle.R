##' Generate a design matrix corresponding to the 2-1-0 hurdle
##'
##' @param zif zero-inflated covariates
##' @param nodeId names to be used for nodes.  `colnames` of `zif` will be used if provided.
##' @param fixed optional matrix of fixed predictors
##' @param center center design matrix
##' @param scale scale design matrix
##' @param conditionalCenter in positive predictors, center only the non-zero components
##' (design matrix will still be centered, and orthogonal to indicators)
##' @param tol tolerance for declaring a parameter equal to zero
##' @return design matrix, with attribute fixedCols giving indices of unpenalized intercept columns. Each gene appears in adjacent columns with its continuous component first, followed by its binarization.
##' @export
makeModel <- function(zif, nodeId, fixed=NULL, center=TRUE, scale=FALSE, conditionalCenter=TRUE, tol=.001){
    nonZ <- not_zero(zif, tol)
    if(conditionalCenter){
        for(i in seq_len(ncol(zif))){
            zifPos <- zif[nonZ[,i],i]
            zif[nonZ[,i],i] <- zifPos-mean(zifPos)
        }
    }
    if(is.null(fixed)) fixed <- matrix(1, nrow=nrow(zif), ncol=1)
    if(!is.matrix(fixed) && !is.numeric(fixed)) stop('`Fixed` must be a numeric matrix')

    ## for(i in seq_len(ncol(zif))){
    ##     cols <- seq(2*i-1, 2*i)
    ##     UDV <- svd(M[,cols])
    ##     M[,cols] <- UDV$u
    ## }
    
    M <- cbind(zif, nonZ*1)
    ord <- c(seq(2, ncol(zif)+1),seq(2, ncol(zif)+1))
    M <- M[,order(ord)]
    M <- scale(M, center=center, scale=scale)
    MM <- cbind(fixed, M)
    if(missing(nodeId)){
        nodeId <- colnames(zif)
    }
    if(length(nodeId) != ncol(zif)){
        stop('Length of `nodeId` must match `ncols(zif)`')
    }
    structure(MM, fixedCols=seq_len(ncol(fixed)), nodeId=nodeId)
}


not_zero <- function(x, tol=0){
    abs(x)>tol   
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
    while(not_zero(error_r, r_solve_thresh) && r_iter<max_r_iters){
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
##' @param nodeId optional labels for the nodes.  will be used to stitch to
##' @param nlambda if `lambda` is not provided, then the number of lambda to interpolate between
##' @param lambda.min.ratio if `lambda` is not provided, then the left end of the solution path as a function of the lambda0, the lambda for the empty model
##' @param lambda penalty path desired
##' @param penaltyFactor one of `full`, `diagonal` or `identity` giving how the penalty should be scaled \emph{blockwise}
##' @param control optimization control parameters
##' @param theta (optional) initial guess for parameter
##' @return matrix of parameters, one row per lambda
##' @useDynLib HurdleNormal
##' @import Rcpp
##' @export
cgpaths <- function(y.zif, this.model, Blocks=Block(this.model), nodeId=NA_character_, nlambda=50, lambda.min.ratio=if(length(y.zif)<ncol(this.model)) .005 else .05, lambda, penaltyFactor='full', control=list(tol=5e-3, maxrounds=300, debug=1), theta){
    defaultControl <- list(tol=5e-3, maxrounds=300, debug=1, stepcontract=.5, stepsize=1, stepexpand=.1, FISTA=FALSE, newton0=FALSE, safeRule=2, returnHessian=FALSE, refit=TRUE)
    nocontrol <- setdiff(names(defaultControl), names(control))
    control[nocontrol] <- defaultControl[nocontrol]
    penaltyFactor <- match.arg(penaltyFactor, c('full', 'diagonal', 'identity'))
    
    p <- ncol(this.model)
    n <- nrow(this.model)
    ## mapping from blocks to indices in theta.
    blocklist <- split(Blocks$map$paridx, Blocks$map$block)
    
    ## start at 0 (except kbb=1)
    if(missing(theta))     theta <- c(rep(0, 2*p), 1)



    emptySol <- function(lambda){
        if(missing(lambda)){
            lambda <- 1
        }

        out <- Matrix::Matrix(0, nrow=length(lambda), ncol=length(theta), sparse=TRUE)
        ## kkt conditions
        kktout <- Matrix::Matrix(0, nrow=length(lambda), ncol=length(blocklist), sparse=TRUE)
        ## diagnostic flags
        flout <- matrix(NA, nrow=length(lambda), ncol=6)
        colnames(flout) <- c('converged', 'round', 'geval', 'feval', 'gamma', 'nnz')
        colnames(out) <- c(outer(1:p, c('D', 'C'),paste0), 'kbb')
        rownames(kktout) <- rownames(flout) <- rownames(out) <- lambda
        return(list(path=out, kktout=kktout, flout=flout, blocks=Blocks, lambda=lambda, nodeId=nodeId, path_np=out, loglik_np=rep(999999999, length(lambda))))
}

    ## get likelihood object
    ## ultimately calls C++ code in src/hurdle_likelihood.cpp
    if(var(y.zif)<1e-5){
        warning("Response with no variance")
        return(emptySol(lambda))
    }
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
    LLall <- function(th, penalize=FALSE) hl$LLall(th, penalize)
    ## gradient for all blocks
    gradAll <- hl$gradAll
    gamma <- control$stepsize
    
    ## value for empty sol
    o <- refitModel(theta, this.model, y.zif, activetheta=Blocks$nonpenpar, Blocks)
    if(o$convergence !=0){
        warning('Empty solution failed to converge')
        return(emptySol(lambda))
    }
    theta[o$activetheta] <- o$par

    ## Hessian at empty solution
    sg <- getsg(theta)
    hess <- blockHessian(theta, sg, Blocks, this.model, onlyActive=FALSE, control, fuzz=.1, hl, exact=TRUE)
    np <- solve(hess[[Blocks$nonpengrp]])*min(eigen(hess[[Blocks$nonpengrp]])$values)
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
            penMat[[b]] <- solve(hess[[b]])
        }
        eig <- eigen(penMat[[b]])
        eigval[[b]] <- eig$values
        eigvec[[b]] <- eig$vectors
        if(any(eig$values<0)){
            warning('Negative eigenvalues in block ', b, " :'(")
            eig$values <- pmax(eig$values, 0)
            browser()
            }
        sqrtPen[[b]] <- eig$vectors %*% sqrt(diag(eig$values)) %*% t(eig$vectors)
    }

    pre0 <- (penMat[[1]])/(max(eigval[[1]]))

    ## Returns proximal subtheta
    proxfun <- function(b, subtheta, gamma, lambda){
        if(Blocks$lambdablock[b]<=0 || lambda <=0){
            return(subtheta)
        } else {
            u <- eigvec[[b]]
            d <- eigval[[b]]
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
    kktfun <- function(b, theta, nonpengrad, lambda, gnormOnly=FALSE){
        sgrad <- nonpengrad[blocklist[[b]]]
        if(Blocks$lambdablock[b]<=0){
            gnorm <- sqrt(sum(sgrad^2))
            tol <- gnorm
        } else{
            subtheta <- theta[blocklist[[b]]]
            tnorm <- as.vector(sqrt(crossprod(subtheta, hess[[b]]) %*% subtheta))
            if(tnorm > 0){
                pgrad <- sgrad + lambda * as.vector(hess[[b]] %*% subtheta)/tnorm
                gnorm <- sqrt(sum(pgrad^2))
                tol <- gnorm
            } else{
                gnorm <- sqrt(sum(crossprod(sgrad, sqrtPen[[b]])^2))
                tol <- gnorm>lambda
            }
        }
        if(gnormOnly) return(gnorm) else return(c(gnorm, tol))
    }

    ## estimate lambda0 given penalty
    l0block <- sapply(seq_along(blocklist), kktfun, theta=theta, nonpengrad=hl$gradAll(theta, penalize = TRUE), lambda=0, gnormOnly=TRUE)
    gradblock <- data.table(g0=l0block, block=seq_along(l0block), active=FALSE)
    gradblock[block==1, active:=TRUE]
    l0 <- max(l0block)
    #Blocks$map <- merge(Blocks$map, , by='block')
    if(l0<0){
        warning('Negative lambda0 should not happen')
        return(emptySol(lambda))
    }
    if(missing(lambda)) lambda <- 10^seq(log10(1.01*l0), log10(l0*lambda.min.ratio), length.out=nlambda)
    
    path_np <- out <- matrix(0, nrow=length(lambda), ncol=length(theta))
    ## kkt conditions
    kktout <- matrix(0, nrow=length(lambda), ncol=length(blocklist))
    ## diagnostic flags
    flout <- matrix(NA, nrow=length(lambda), ncol=6)
    ## non-penalized log-likelihood
    nloglik_np <- rep(NA_real_, length(lambda))
    colnames(flout) <- c('converged', 'round', 'geval', 'feval', 'gamma', 'nnz')
    colnames(path_np) <- colnames(out) <- c(outer(1:p, c('D', 'C'),paste0), 'kbb')
    names(nloglik_np) <- rownames(path_np) <- rownames(kktout) <- rownames(flout) <- rownames(out) <- lambda

    ## loop over lambda
    refitted <- FALSE
    for(l in seq_along(lambda)){
        hl$lambda <- lambda[l]
        if(l>1){
            gradblock <- safeRule(gradblock, lambda[l], lambda[l-1], control$safeRule)
        }
        sp <- solvePenProximal(theta, lambda[l], control, blocklist, hl, proxfun, kktfun, hess, pre0, gamma, gradblock)
        gamma <- as.numeric(attr(sp, 'flag')['gamma'])
        theta <- out[l,] <- as.numeric(sp)
        kktout[l,] <- attr(sp, 'kkt')[2,]
        flout[l,] <- attr(sp, 'flag')
        gradblock <- attr(sp, 'gradblock')
        if(control$debug>0) message('Lambda = ', round(lambda[l], 3), ' rounds = ', attr(sp, 'flag')['round'], ' NNZ = ', attr(sp, 'flag')['nnz'], ' gamma = ', round(gamma, 3))
        activeset <- sum(kktout[l,]>0)
        if(control$refit){
            if(is.logical(refitted) || !all( not_zero(out[l-1,]) == not_zero(theta)) ){
                refitted <- refitModel(theta, this.model, y.zif, blocks=Blocks, control=list(factr=1e9))
            }
            nloglik_np[l] = refitted$value
            path_np[l,refitted$activetheta] = refitted$par
        }
        if(activeset >= n/2){
            out <- out[1:l,,drop=FALSE]
            kktout <- kktout[1:l,,drop=FALSE]
            flout <- flout[1:l,,drop=FALSE]
            l <- lambda[1:l]
            break
        }
    }
    
    res <- list(path=Matrix::Matrix(out, sparse=TRUE), kktout=Matrix::Matrix(kktout, sparse=TRUE), flout=flout, blocks=Blocks, lambda=lambda, nodeId=nodeId, hessian=if(control$returnHessian) hess else NULL, path_np=Matrix::Matrix(path_np, sparse=TRUE), loglik_np=-nloglik_np*nrow(this.model), nobs=length(y.zif))
    class(res) <- 'SolPath'
    res
}

globalVariables(c('active'))

refitModel <- function(theta_, this.model, y.zif, activetheta, blocks, fuzz=0, control=list()){
    m <- blocks$map
    if(missing(activetheta)){
        m[,theta:=theta_]
        activeblocks <- m[,.(active=any(not_zero(theta, fuzz))), keyby=block]
        activetheta <- sort(m[activeblocks[active==TRUE,],paridx,on='block'])
        }
    activemm = sort(unique(na.omit(blocks$map[list(paridx=activetheta),mmidx, on='paridx'])))
    hl0 <- HurdleLikelihood(y.zif, this.model[,activemm,drop=FALSE], theta=theta_[activetheta], lambda=0)
    o <- optim(theta_[activetheta], hl0$LLall, hl0$gradAll, penalize=FALSE, method='L-BFGS-B', control=control)
    c(o, activetheta=list(activetheta))
}

## .makeParams <- function(lambda, nlambda, theta, blocklist){
##      ## solution path
   
##     assign('out', out, pos=sys.frame(1))
##     assign('kktout', kktout, pos=sys.frame(1))
##     assign('flout', flout, pos=sys.frame(1))
## }

blockHessian <- function(theta, sg, Blocks, X, onlyActive=TRUE, control, fuzz=.1, hl, exact=FALSE){
    hess <- vector('list', length(Blocks))
    hl$LLall(theta, penalize=FALSE)
    gpart <- hl$gpart
    gplusc <- hl$gplusc
    hpart <- as.vector(hl$hpart)
    w1 <- as.vector(hl$cumulant2)
    w2 <- as.vector(exp(-gplusc)*w1^2)
    kbb <- theta[length(theta)]
    for(b in seq_along(Blocks$lambdablock)){
        bb <- na.omit(unique(Blocks$map[block==b, mmidx]))
        if(Blocks$lambdablock[b]==0){
            hess[[b]] <- numDeriv::jacobian(sg, theta[Blocks$map[block==b, paridx]], b=b)
            hess[[b]] <- zapsmall((hess[[b]]+t(hess[[b]]))/2+diag(nrow(hess[[b]]))*fuzz, digits=6)
        } else{
            lb <- length(bb)
            Xb <- X[,bb]
            tmp <- matrix(NA, nrow=2*lb, ncol=2*lb)
            tmp[1:lb, 1:lb] <- crossprod(Xb, Xb*w2)
            tmp[1:lb, (lb+1):(2*lb)] <- tmp[(lb+1):(2*lb), 1:lb] <- crossprod(Xb, Xb*w2*hpart)/kbb
            tmp[(lb+1):(2*lb), (lb+1):(2*lb)] <- crossprod(Xb*w1*(hpart^2*(1-w1)/kbb+1), Xb)/kbb
            hess[[b]] <- tmp/nrow(X) + diag(nrow(tmp))*fuzz
        }
    }
    hess
}

safeRule <- function(gradblock, l1, l0, rule){
    gradblock[active==FALSE, active:=(rule*l1-l0)<g0]
    gradblock
}

## Solve the penalized function using proximal gradient descent
solvePenProximal <- function(theta, lambda, control, blocklist, hl, proxfun, kktfun, hess, pre0, gamma, gradblock){
    feval <- geval <- 0
    ## total number of rounds
    mround <- round <- 0
    ## amount of consecutive successful (non-line search) prox's we've done
    step <- 1
    p <- length(theta)
    converged <- FALSE
    ## FISTA extrapolation parameter (default: no extrapolation)
    omega <- 0

    ## thetaPrime0 is previous iteration, thetaPrime1 is current and theta is proposed proximal point
    ## (used if we're extrapolating via FISTA)
    thetaPrime0 <- thetaPrime1 <- theta
    ## line search parameters
    if(gamma<=0)     gamma <- control$stepsize
    thisll <- hl$LLall(theta, penalize = FALSE)
    gr <- hl$gradFixed(penalize = FALSE)

    kkt <- vapply(seq_along(blocklist), kktfun, c(NA_real_, 1), theta=thetaPrime1, nonpengrad=gr, lambda=lambda)
    if(control$debug>2){
        debugval <- list(pre=matrix(NA, nrow=control$maxrounds, ncol=p), thisll=rep(NA, control$maxrounds),
                         kkt=matrix(NA, nrow=control$maxrounds, ncol=length(blocklist)),
                         theta=matrix(NA, nrow=control$maxrounds, ncol=length(theta)),
                         gamma=rep(NA, control$maxrounds))
    } else{
        debugval <- NULL
    }
    ## apparently indexing into gradblock is slow.
    bActive <- gradblock[active==TRUE,block]
    while(round < control$maxrounds && any(not_zero(kkt[2,], control$tol))){
        round <- round+1
        for(b in bActive){
            ##apply proximal operator to block b
            if(control$newton0 & b==1){
                theta[blocklist[[b]]] <- thetaPrime1[blocklist[[b]]] - gamma*crossprod(pre0, gr[blocklist[[b]]])
            } else{
                theta[blocklist[[b]]] <- proxfun(b, thetaPrime1[blocklist[[b]]]-gamma*gr[blocklist[[b]]], gamma, lambda)
            }
        }
        feval <- feval+1
        newll <- hl$LLall(theta, penalize=FALSE)
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
            if(control$FISTA){
                omega <- mround/(mround+5)
                #browser(expr=round>100)
            ## accept proposed proximal point, and extrapolate
                thetaTmp <- theta + omega*(theta-thetaPrime0)
                thetaPrime0 <- thetaPrime1
                thetaPrime1 <- thetaTmp
                newll <- hl$LLall(thetaPrime1, penalize=F)
            } else{
                thetaPrime0 <- thetaPrime1
                thetaPrime1 <- theta
            }
            mround <- mround+1
            thisll <- newll
            step <- step+1
            if(step>4){
                ##try increasing stepsize by shrinking towards control$stepsize
                gamma <- gamma*(1-control$stepexpand)+control$stepsize*control$stepexpand
            }

            ##update gradient and test kkt
            ## (no need if we haven't moved thetaPrime1)
            gr <- hl$gradFixed(penalize = FALSE)
            geval <- geval+length(blocklist)
            kktA <- vapply(bActive, kktfun, c(NA_real_, 1), theta=thetaPrime1, nonpengrad=gr, lambda=lambda)
            if(all(kktA[2,]<control$tol)){
                kkt[,bActive] <- kktA
                bInactive <- gradblock[active==FALSE,block]
                kktI <- vapply(bInactive, kktfun, c(NA_real_, 1), theta=thetaPrime1, nonpengrad=gr, lambda=lambda)
                kkt[,bInactive] <- kktI
                gradblock[active==FALSE,active:=kktI[2,]>0]
                bActive <- gradblock[active==TRUE,block]
            }

            ## Debugging/tracing
            if(control$debug>1 && (round %% 10)==0 || control$debug>2){
                pen <- sapply(1:length(blocklist), function(b){
                    st <- theta[blocklist[[b]]]
                    lambda*sqrt(crossprod(st, hess[[b]]) %*% st)
                })
                message(noquote(paste0('penll=', round(hl$LLall(thetaPrime1, penalize=FALSE) + sum(pen), 5), ' theta= ', paste(round(thetaPrime1, 2), collapse=','), 'gamma= ', sprintf('%2.2e', gamma))))
                debugval$thisll[round] <- hl$LLall(thetaPrime1, FALSE)+sum(pen)
           }
            if(control$debug>2){
                debugval$kkt[round,] <- kkt
                debugval$theta[round,] <- thetaPrime1
                debugval$gamma[round] <- gamma
            }

        } # end move
    } # end main loop
    gradblock[,g0:=kkt[1,]]
    structure(theta, flag=c(converged=converged, round=round, geval=geval, feval=feval, gamma=gamma, nnz=sum(kkt[2,]>0)-1), kkt=kkt, debugval=debugval, gradblock=gradblock)
}




#' Plot a solution path
#'
#' @param cgpaths result of \code{cgpaths}
#' @importFrom ggplot2 scale_x_log10
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
    directlabels::direct.label(pathPlot, 'last.points')
    invisible(l2)
}
