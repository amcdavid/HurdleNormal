## names of parameters in vector
parmap <- function(p){
     rep(c('gbb', 'gba', 'hbb', 'hab', 'hba', 'kba', 'kbb'), times=c(1, p-1, 1, p-1, p-1, p-1, 1))
}

## block groups
## block 1 is unpenalized intercepts
coordmap <- function(p){
    ## Gbb Gab  Hbb  Hab  Hba   Kba  Kbb
    c(1,   2:p, 1,   2:p, 2:p,  2:p, 1)
}

## the negative of the log-likelihood!
generatelogLik <- function(y, x, lambda=0, debug=FALSE, returnGrad=FALSE, xNames=NULL){
    TOL <- .001
    LARGE <- 500
    p <- ncol(x)+1
    xI <- abs(x)>TOL
    yI <- abs(y)>TOL
    gidx <- seq(2, p)
    par <- parmap(p)
    
    ## gradient of pieces of log-likelihood
    ## constant in theta
    dg <- list(gbb=1,
               gba=2*xI,
               hba=0,
               hbb=0,
               hab=x,
               kba=0,
               kbb=0)
    dh <- list(hba=xI,
               hbb=1,
               hab=0,
               kba=-x,
               gbb=0,
               gba=0,
               kbb=0)
    dk <- list(hab=0,
               hbb=0,
               hba=0,
               kba=0,
               gbb=0,
               gba=0,
               kbb=1)

    ## map from coordinate to the set of parameters for that coordinate
    ## own coordinate -> cmap==1 
    cmap <- coordmap(p)

    bookkeeping <- function(theta){
        frame <- sys.frame(-1)

        gBA <- as.matrix(theta[par=='gba'])
        hAB <- as.matrix(theta[par=='hab'])
        hBA <- as.matrix(theta[par=='hba'])
        kBA <- as.matrix(theta[par=='kba'])
        g0 <- theta[par=='gbb']
        h0 <- theta[par=='hbb']
        k0 <- theta[par=='kbb']
        stopifnot(names(g0)=='gbb' & names(h0)=='hbb' & names(k0)=='kbb')

        gpart <- g0+2*xI%*%gBA+x%*%hAB
        hpart <- h0+xI%*%hBA-x%*%kBA
        cpart <- -.5*log(k0/(2*pi)) + hpart^2/(2*k0)

        lgroup <- rep(0, p)
        for(i in seq(2, p)){
            lgroup[i] <- sqrt(sum(theta[cmap==i]^2))
        }
        stopifnot(lgroup[1]==0)
        

        ## prevent scaling problems
        cumulant <-ifelse(gpart + cpart>LARGE, gpart+cpart, log(1+ exp(gpart+cpart)))
        assign('gBA', gBA, pos=frame)
        assign('hAB', hAB, pos=frame)
        assign('hBA', hBA, pos=frame)
        assign('kBA', kBA, pos=frame)
        assign('g0', g0, pos=frame)
        assign('h0', h0, pos=frame)
        assign('k0', k0, pos=frame)
        assign('gpart', gpart, pos=frame)
        assign('hpart', hpart, pos=frame)
        assign('cpart', cpart, pos=frame)
        assign('cumulant', cumulant, pos=frame)
        assign('lgroup', lgroup, pos=frame)
    }

    ll <- function(theta){
        bookkeeping(theta)
        ## Note--could be a decent way to speed this up for sparse y
        loglik <- yI*gpart+y*hpart-.5*y^2*k0-cumulant
        sloglik <- mean(loglik)
        
        if(debug){
        cat(theta, ': ')
        cat(sloglik, '\n')
    }
        negpenll <- -(sloglik-lambda*sum(lgroup))
        ifelse(is.finite(negpenll), negpenll, Inf)
    }

    grad <- function(theta){
        bookkeeping(theta)
        ## sub parts to gradient
        hpart.expand <- do.call(cbind, rep(list(hpart), p-1))
        dC <- list(hbb=hpart/k0,
                   hba=xI*hpart.expand/k0,
                   kba=-x*hpart.expand/k0,
                   gbb=0,
                   gba=0,
                   hab=0,
                   kbb=-.5*(hpart^2+k0)/k0^2)
        ecumulant <- exp(gpart+cpart)
        rcumulant <- ifelse(gpart+cpart>LARGE, 1, ecumulant/(1+ecumulant))

        rcumulant.expand <- do.call(cbind, rep(list(rcumulant), p-1))

        grad <- matrix(NA, nrow=length(y), ncol=length(theta))

        ## loop through coordinates and apply chain rule
        for(p in unique(par)){
            this.rc <- if(sum(p==par)>1) rcumulant.expand else rcumulant
            grad[,p==par] <- dg[[p]]*yI + dh[[p]]*y - .5*y^2*dk[[p]] - (dg[[p]]+dC[[p]])*this.rc # ugh, dimensions varying
        }

        gradPen <- rep(0, length(par))
        if(lambda>0){
            for(i in seq_along(par)){
                if(cmap[i]==1) next
                gradPen[i] <- lambda*theta[i]/lgroup[cmap[i]]
            }
        }

        #                                browser()
        sgrad <- colSums(grad)/nrow(grad) - gradPen
        if(debug){
            cat('grad:', sgrad, '\n')
        }
        -sgrad
    }

    
    
    if(returnGrad) grad else ll
}

## ## the negative of the log-likelihood!
## generatelogLik <- function(y, x, lambda=0, debug=FALSE, xNames=NULL){
##     TOL <- .001
##     LARGE <- 500
##     p <- ncol(x)+1
##     xI <- abs(x)>TOL
##     yI <- abs(y)>TOL
##     gidx <- seq(2, p)
##     par <- parmap(p)
    
##     ## gradient of pieces of log-likelihood
##     ## constant in theta
##     dg <- list(gbb=1,
##                gba=2*xI,
##                hba=0,
##                hbb=0,
##                hab=x,
##                kba=0,
##                kbb=0)
##     dh <- list(hba=xI,
##                hbb=1,
##                hab=0,
##                kba=-x,
##                gbb=0,
##                gba=0,
##                kbb=0)
##     dk <- list(hab=0,
##                hbb=0,
##                hba=0,
##                kba=0,
##                gbb=0,
##                gba=0,
##                kbb=1)

##     ## map from coordinate to the set of parameters for that coordinate
##     ## own coordinate -> cmap==1 
##     cmap <- coordmap(p)

##     llAndGrad <- function(theta){
##         ## Common stuff and bookkeeping
##         gBA <- as.matrix(theta[par=='gba'])
##         hAB <- as.matrix(theta[par=='hab'])
##         hBA <- as.matrix(theta[par=='hba'])
##         kBA <- as.matrix(theta[par=='kba'])
##         g0 <- theta[par=='gbb']
##         h0 <- theta[par=='hbb']
##         k0 <- theta[par=='kbb']
##         stopifnot(names(g0)=='gbb' & names(h0)=='hbb' & names(k0)=='kbb')

##         gpart <- g0+2*xI%*%gBA+x%*%hAB
##         hpart <- h0+xI%*%hBA-x%*%kBA
##         cpart <- -.5*log(k0/(2*pi)) + hpart^2/(2*k0)

##         lgroup <- rep(0, p)
##         for(i in seq(2, p)){
##             lgroup[i] <- sqrt(sum(theta[cmap==i]^2))
##         }
##         lgroup <- lgroup*lambda
##         stopifnot(lgroup[1]==0)
        

##         ## prevent scaling problems
##         cumulant <-ifelse(gpart + cpart>LARGE, gpart+cpart, log(1+ exp(gpart+cpart)))

##         ## Log-likelihood
##         ## Note--could be a decent way to speed this up for sparse y
##         loglik <- yI*gpart+y*hpart-.5*y^2*k0-cumulant
##         sloglik <- sum(loglik)
        
##         if(debug){
##         cat(theta, ': ')
##         cat(sloglik, '\n')
##     }
##         negpenll <- -sloglik-sum(lgroup)
##         negpenll <- ifelse(is.finite(negpenll), negpenll, Inf)

##         ## Gradient:
##         ## sub parts
##         hpart.expand <- do.call(cbind, rep(list(hpart), p-1))
##         dC <- list(hbb=hpart/k0,
##                    hba=xI*hpart.expand/k0,
##                    kba=-x*hpart.expand/k0,
##                    gbb=0,
##                    gba=0,
##                    hab=0,
##                    kbb=-.5*(hpart^2+k0)/k0^2)
##         ecumulant <- exp(gpart+cpart)
##         rcumulant <- ifelse(gpart+cpart>LARGE, 1, ecumulant/(1+ecumulant))

##         rcumulant.expand <- do.call(cbind, rep(list(rcumulant), p-1))

##         grad <- matrix(NA, nrow=length(y), ncol=length(theta))

##         ## loop through coordinates and apply chain rule
##         for(p in unique(par)){
##             this.rc <- if(sum(p==par)>1) rcumulant.expand else rcumulant
##             grad[,p==par] <- dg[[p]]*yI + dh[[p]]*y - .5*y^2*dk[[p]] - (dg[[p]]+dC[[p]])*this.rc # ugh, dimensions varying
##         }

##         gradPen <- rep(0, length(par))
##         if(lambda>0){
##             for(i in seq_along(par)){
##                 if(cmap[i]==1) next
##                 theta[i]/lgroup[cmap[i]]
##             }
##     }

##         #                                browser()
##         sgrad <- colSums(grad)
##         if(debug){
##             cat('grad:', sgrad, '\n')
##         }
##         list(objective=negpenll, gradient=-sgrad)
##     }

##     llAndGrad
## }

LAYERMAP <- c('gbb'='G', 'gba'='G',
              'hbb'='Hdisc', 'hab'='Hdisc', 'hba'='Hcont',
              'kba'='K', 'kbb'='K')
OLAYER <- c(G=1, Hdisc=2, Hcont=3, K=4)

getJoint210 <- function(zmFits){
    s2 <- attr(zmFits, 'sigma2')
    genes <- which(!is.na(s2[,'dichotomous'])) #genes with errors will be NA                                    
    ng <- length(genes)
    ## Off-diagonal entries are all estimated twice:
    ## Once for A|B and once for B|A
    ## H matrix has more complicated relationship:
    ## Upper tri of A|B is lower tri of B|A
    ## So we'll keep the redundant estimates in two separate slices
    ## (Hdisc=upper tri of B|A, Hcont=lower tri of B|A)
    ## Then average them
    theta <- array(NA, dim=c(ng, ng, 4), #G, Hdisc, Hcont, K
                   dimnames=list(names(genes), names(genes), c('G', 'Hdisc', 'Hcont', 'K')))
    layer <-  OLAYER[LAYERMAP[parmap(ng)]]
    for(goodidx in seq_len(ng)){
        g <- genes[goodidx]             #original index
        par <- zmFits[[g,1]]
        this.genes <- attr(par, 'genes')
        good.genes.idx <- which(this.genes %in% names(genes))
        good.genes <- this.genes[good.genes.idx]
        browser(expr=length(intersect(names(genes), good.genes))!=ng)
        ##stopifnot(length(intersect(names(genes), good.genes))==ng)
        idx <- cbind(goodidx, match(good.genes, names(genes)), layer )
        theta[idx] <- par$coef[good.genes.idx]
    }
    G <- (theta[,,'G']+t(theta[,,'G']))/2
    Hdisc <- theta[,,'Hdisc']
    Hcont <- theta[,,'Hcont']
    diag(Hcont) <- diag(Hdisc)
    H <- (Hdisc+Hcont)/2
    K <- (theta[,,'K']+t(theta[,,'K']))/2
    list(G=G, H=H, K=K)
}

fix <- function(fun, fixedCoords, fixed){
    function(y){
        fixed[-fixedCoords] <- y
        fun(fixed)
    }
}
