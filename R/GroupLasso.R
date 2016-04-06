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
##' @param zif zero-inflated covariates
##' @param nodeId names to be used for nodes.  `colnames` of `zif` will be used if provided.
##' @param Fixed optional matrix of fixed predictors
##' @param center center design matrix
##' @param scale scale design matrix
##' @param conditionalCenter in positive predictors, center only the non-zero components
##' (design matrix will still be centered, and orthogonal to indicators)
##' @param ... ignored
##' @param tol tolerance for declaring a parameter equal to zero
##' @return design matrix, with attribute fixedCols giving indices of unpenalized intercept columns
##' @export
makeModel <- function(zif, nodeId, Fixed=matrix(1, nrow=nrow(zif), ncol=1), center=TRUE, scale=FALSE, conditionalCenter=TRUE, ..., tol=.001){
    nonZ <- abs(zif)>tol
    if(conditionalCenter){
        for(i in seq_len(ncol(zif))){
            zifPos <- zif[nonZ[,i],i]
            zif[nonZ[,i],i] <- zifPos-mean(zifPos)
        }
    }

    ## for(i in seq_len(ncol(zif))){
    ##     cols <- seq(2*i-1, 2*i)
    ##     UDV <- svd(M[,cols])
    ##     M[,cols] <- UDV$u
    ## }
    
    M <- cbind(zif, nonZ*1)
    ord <- c(seq(2, ncol(zif)+1),seq(2, ncol(zif)+1))
    M <- M[,order(ord)]
    M <- scale(M, center=center, scale=scale)
    MM <- cbind(Fixed, M)
    if(missing(nodeId)){
        nodeId <- colnames(zif)
    }
    if(length(nodeId) != ncol(zif)){
        stop('Length of `nodeId` must match `ncols(zif)`')
    }
    structure(MM, fixedCols=seq_len(ncol(Fixed)), nodeId=nodeId)
}



##' Define parameter groups
##'
##' Define parameter groups and relative penalization.
##' If \code{this.model} is specified, then, neither \code{blist} nor \code{nlist} should be provided.
##' If \code{this.model} is not specified then both of the preceeding must be provided.
##' When \code{this.model} is provided, it is assumed to have been created by with \code{makeModel}.
##'  Column indices indicated in the attribute `fixedCols` are unpenalized intercepts, then followed by pairs of interaction terms.
##' The grouping corresponds to the parameter set for a node.
##' Otherwise, only one of bidx or bset should be set, but this is not implemented yet.
##' @details
##' There are four components that all need to be mapped between each other.
##' In increasing abstraction, with variable prefixes in parenthesis:
##' 1. (p)arameter vector. The map is given in these terms. 
##' 2. (mm) model matrix--columns from the covariate matrix.
##' 3. (b)locks -- penalty groups
##' 4. (n)odes. -- the graph-theoretic structure.
##' 5. lambda -- penalties as a function of blocks.
##' All of these components are provided in the `map`
##' @param this.model (optional) a matrix created by \code{\link{makeModel}}
##' @param blist a list of parameter indices, one per block.  By default the first block is assumed to be unpenalized.
##' @param mlist a list of parameter indices, one per column of the model.matrix. If omitted, assumed to equal to the identity.
##' @param nlist a named list of block indices, one per node
##' @param group \code{character}: one of components, or none.
##' @param penalty.scale optional list containing elements `scale` and `group`.
##' `group` should be one of 'block' or 'none'.  `scale` should be \code{numeric} of length `blist` or the sum of the `blist` lengths.
##' @return a list containing a data.table `map` giving the mapping between parameters, groups and penalty scales and some other components
##' @export
Block <- function(this.model, blist, mlist, nlist, lambda, group='components', penalty.scale=NULL){
    ## only one of this.model, bidx and bset should be non-missing
    ## group={components, none}
    ## penalty scale should be provided as a list of indices and group={groups, components,none}


    group <- match.arg(group, c('components', 'none'))
    if(!is.null(penalty.scale)){
        if(!is.list(penalty.scale)) stop("`penalty.scale` should be list with components `scale` and `group`")
        penalty.scale.group <- match.arg(penalty.scale$group, c('components', 'none', 'groups'))
        if(!is.numeric(penalty.scale$scale)) stop("`penalty.scale$scale` should give integer indices into groups or parameters")
        penalty.scale.scale <- penalty.scale$scale
    }
    ## this.model provided
    if(!missing(this.model)){
        nfixed <- length(attr(this.model, 'fixedCols'))
        nodeId <- c('(Fixed)', attr(this.model, 'nodeId'))

        ## used to assign nodeIds and blocks when group=='components'
        p <- ncol(this.model)
        nc <- (p-nfixed)/2 #penalized blocks
        ## block index into model matrix
        bidxMM <- c(rep(1, nfixed), (rep(1:nc, each=2)
                        +1)) #offset from non-pen block
        ## node index into parameter vector
        nidxPar <- c(bidxMM, bidxMM, 1) #variance

        ## group lasso
        if(group=='components'){
            ## block index into parameter
            bidxPar <- nidxPar
            ## node -> block map
        } else{
            nc <- (p-nfixed)               #regular lasso
            bidxMM <- c(rep(1, nfixed), ((1:nc)+1))
            bidxMM2 <- c(rep(1, nfixed), ((1:nc)+nc+1))
            bidxPar <- c(bidxMM, bidxMM2, 1)
        }
        nodeMap <- data.table(nodeId, nidx=seq_along(nodeId))
        map <- data.table(paridx=seq_len(2*p+1), block=bidxPar, mmidx=c(1:length(bidxMM), 1:length(bidxMM), NA), nidx=nidxPar)
        map <- merge(map, nodeMap, by='nidx')
    }
    ## end this.model provided
    else{
        bvec <- data.table(paridx=unlist(blist), block=rep(seq_along(blist), times=sapply(blist, length)))
        if(is.list(nlist)){
            nvec <- setNames(reshape2:::melt.list(nlist), c('block', 'nodeId'))
        } else{
            nvec <- data.table(block=seq_along(nlist), nodeId=nlist)
        }
        if(missing(mlist)) mlist <- as.list(bvec$paridx)
        mvec <- data.table(paridx=unlist(mlist), mmidx=rep(seq_along(mlist), times=sapply(mlist, length)))
        map <- merge(bvec, nvec, by='block')
        map <- merge(map, mvec, by='paridx')        
    }
    
    if(is.null(penalty.scale)){
        blocks <- sort(unique(map$block))
        penalty.map <- data.table(block=blocks, lambda=c(0, rep(1, length(blocks)-1)))
        penalty.scale.group <- 'groups'
    }    
    if(penalty.scale.group=='groups'){
        map <- merge(map, penalty.map, by='block')
    } else if(penalty.scale.group=='components'){
        stop('Not implemented')
    } else{
        penalty.map <- data.table(paridx=seq_along(penalty.scale.scale), lambda=penalty.scale.lambda)
        map <- merge(map, penalty.map, by='paridx')
    }

    setkey(map, paridx)
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
##' @export
cgpaths <- function(y.zif, this.model, Blocks=Block(this.model), nodeId=NA_character_, nlambda=100, lambda.min.ratio=if(length(y.zif)<ncol(this.model)) .005 else .05, lambda, penaltyFactor='full', control=list(tol=1e-3, maxrounds=300, debug=1), theta){
    defaultControl <- list(tol=1e-3, maxrounds=300, debug=1, stepcontract=.5, stepsize=1, stepexpand=.1, FISTA=FALSE, newton0=FALSE)
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
        return(list(path=out, kktout=kktout, flout=flout, blocks=Blocks, lambda=lambda, nodeId=nodeId))
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
    interceptCol <- Blocks$nonpenMM
    interceptPar <- Blocks$nonpenpar
    theta0 <- theta[interceptPar]
    hl0 <- HurdleLikelihood(y.zif, this.model[,interceptCol,drop=FALSE], theta=theta0, lambda=0)
    o <- optim(theta0, hl0$LLall, hl0$gradAll, penalize=FALSE, method='L-BFGS-B')
    if(o$convergence !=0){
        warning('Empty solution failed to converge')
        return(emptySol(lambda))
    }
    theta[interceptPar] <- o$par

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
        sqrtPen[[b]] <- eig$vectors %*% sqrt(diag(eig$values)) %*% t(eig$vectors)
    }

    pre0 <- (penMat[[1]])/(max(eigval[[1]]))

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
    if(l0<0){
        warning('Negative lambda0 should not happen')
        return(emptySol(lambda))
    }
    if(missing(lambda)) lambda <- 10^seq(log10(1.01*l0), log10(l0*lambda.min.ratio), length.out=nlambda)

     out <- matrix(0, nrow=length(lambda), ncol=length(theta))
    ## kkt conditions
    kktout <- matrix(0, nrow=length(lambda), ncol=length(blocklist))
    ## diagnostic flags
    flout <- matrix(NA, nrow=length(lambda), ncol=6)
    colnames(flout) <- c('converged', 'round', 'geval', 'feval', 'gamma', 'nnz')
    colnames(out) <- c(outer(1:p, c('D', 'C'),paste0), 'kbb')
    rownames(kktout) <- rownames(flout) <- rownames(out) <- lambda

    if(control$debug>1) message('grad block1 ', paste(hl0$gradAll(theta[interceptPar]), collapse=','))
    

    ## loop over lambda
    for(l in seq_along(lambda)){
        hl$setLambda(lambda[l])
        sp <- solvePenProximal(theta, lambda[l], control, blocklist, LLall, gradAll, proxfun, kktfun, hess, pre0, gamma)
        gamma <- as.numeric(attr(sp, 'flag')['gamma'])
        theta <- out[l,] <- as.numeric(sp)
        kktout[l,] <- attr(sp, 'kkt')
        flout[l,] <- attr(sp, 'flag')
        if(control$debug>0) message('Lambda = ', round(lambda[l], 3), ' rounds = ', attr(sp, 'flag')['round'], ' NNZ = ', attr(sp, 'flag')['nnz'], ' gamma = ', round(gamma, 3))
        activeset <- sum(kktout[l,]>0)
        if(activeset >= n/2){
            out <- out[1:l,,drop=FALSE]
            kktout <- kktout[1:l,,drop=FALSE]
            flout <- flout[1:l,,drop=FALSE]
            l <- lambda[1:l]
            break
        }
        ## update jerr?
    }
    
    list(path=Matrix::Matrix(out, sparse=TRUE), kktout=Matrix::Matrix(kktout, sparse=TRUE), flout=flout, blocks=Blocks, lambda=lambda, nodeId=nodeId)
}

## .makeParams <- function(lambda, nlambda, theta, blocklist){
##      ## solution path
   
##     assign('out', out, pos=sys.frame(1))
##     assign('kktout', kktout, pos=sys.frame(1))
##     assign('flout', flout, pos=sys.frame(1))
## }

blockHessian <- function(theta, sg, Blocks, X, onlyActive=TRUE, control, fuzz=.1, hl, exact=FALSE){
    hess <- vector('list', length(Blocks))
    #browser()
    hl$LLall(theta, penalize=FALSE)
    gpart <- hl$gpart()
    gplusc <- hl$gplusc()
    hpart <- as.vector(hl$hpart())
    w1 <- as.vector(hl$cumulant2())
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


## Solve the penalized function using proximal gradient descent
solvePenProximal <- function(theta, lambda, control, blocklist, LLall, gradAll, proxfun, kktfun, hess, pre0, gamma){
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
    thisll <- LLall(theta, penalize=FALSE)
    gr <- gradAll(penalize=FALSE)

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
        theta <- thetaPrime1-gamma*gr
        if(control$newton0){
            theta[blocklist[[1]]] <- thetaPrime1[blocklist[[1]]] - gamma*crossprod(pre0, gr[blocklist[[1]]])
        }
        
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
            if(control$FISTA){
                omega <- mround/(mround+5)
                browser(expr=round>100)
            ## accept proposed proximal point, and extrapolate
                thetaTmp <- theta + omega*(theta-thetaPrime0)
                thetaPrime0 <- thetaPrime1
                thetaPrime1 <- thetaTmp
                newll <- LLall(thetaPrime1, penalize=F)
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
            gr <- gradAll(penalize=FALSE)
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
