## Generate a list of lambda (not necessarily overlapping) for each node
## Make some paths
## 

## anArray <- array(c(1, -1, 0,
##                    0, 1, 0,
##                    0, 1, 1), 

## toPL <- function(ar){
##     NULL
## }


context('Test path combining')
## Mock up some data
G <- diag(c(-14, -14, -15, -13))
H <- diag(c(5, 5, 8, 5))
K <- diag(c(1, 1, 2, 1))
G[1,3] <- G[3,1] <- -1
G[1,2] <- G[2,1] <- 1
K[3,4] <- K[4,3] <- .2
set.seed(1234)
rgh <-rGibbsHurdle(G, H, K, 8e3, 1e3, thin=.5)
colnames(rgh) <- LETTERS[1:ncol(rgh)]

lambdaList <- list(seq(.4, .1, by=-.1), seq(.4, .1, by=-.1), seq(.4, .1, by=-.1), seq(1, .2, by=-.2))


## Signed adjacency
adj <- G+H-K
## Neighborhoods
ngh <- lapply(1:nrow(adj), function(i) c(1, adj[i,-i]))
nlists <- list(list('(Fixed)'=1, B=2, C=3, D=4),
               list('(Fixed)'=1, A=2, C=3, D=4),
               list('(Fixed)'=1, A=2, B=3, D=4),
               list('(Fixed)'=1, A=2, B=3, C=4))

## Mocked-up paths
pathList <- lapply(seq_along(ngh), function (i){
    ll <- lambdaList[[i]]
    nl <- length(ll)
    ## scaling factor, repeated in each column for nl
    gamma <- matrix(rep(c(0, seq(0, 1, length.out=nl-1)), times=length(ngh)), ncol=nl, byrow=T)
    ## path
    blk <- Block(blist=list(1, 2, 3, 4), nlist=nlists[[i]])
    res <- list(path=Matrix::Matrix(t(ngh[[i]]*gamma), sparse=T), lambda=ll, blocks=blk, nodeId=colnames(rgh)[i])
    class(res) <- "SolPath"
    res
})



pathArray <- neighborhoodToArray(pathList)
test_that('pathList lambda range spans its members',{
    expect_equal(range(pathArray$lambda), range(do.call(c, lambdaList)))

})


test_that('pathList coincides with pathArray', {    
    pal <- pathArray$lambda
    testSlice <- function(iarray, ipath, nidx){
        if(length(iarray)>0){
            iarray <- iarray[length(iarray)]
            nzarray <- which(abs(pathArray$adjMat[[iarray]][nidx,])>0)
            nzpath <- which(abs(pathList[[nidx]]$path[ipath,])>0)
            return(list(nzarray, setdiff(nzpath, 1)))#no intercept (if provided)
        }
        list(numeric(0), numeric(0))
    }
    
    for(i in seq_along(pathList)){
        pl <- pathList[[i]]
        for(ipath in seq_along(pl$lambda)){
            l <- pl$lambda[ipath]            
            llo <- which(pal<l)
            arrayLess <- testSlice(llo, ipath, i)
            expect_gt(length(setdiff(arrayLess[[1]], arrayLess[[2]])), -0.1)

            lhi <- which(pal>l)
            arrayMore <- testSlice(llo, ipath, i)
            expect_gt(length(setdiff(arrayMore[[2]], arrayMore[[1]])), -0.1)
        }
    }

})

fit <- fitHurdle(rgh, parallel=FALSE, makeModelArgs=list(scale=FALSE, conditionalCenter=TRUE, center=TRUE), keepNodePaths=TRUE, nlambda=10, penalty='full', control=list(tol=5e-2, newton0=TRUE, debug=0))
al <- autoLogistic(rgh, nlambda=50, lambda.min.ratio=.01)
## irrelevant fixed predictor
al2 <- autoLogistic(rgh, fixed=cbind(1, fixedeff=rnorm(nrow(rgh))), family='gaussian',  nlambda=5, lambda.min.ratio=.1, keepNodePaths=TRUE)
## no fixed
al3 <- autoLogistic(rgh, family='gaussian',  nlambda=5, lambda.min.ratio=.1)
## relevant fixed
al4 <- autoLogistic(rgh, fixed=cbind(1, fixedeff=rgh[,4]+rnorm(nrow(rgh))/5), family='gaussian',  nlambda=5, lambda.min.ratio=.1)

test_that('Inherit from SolPath', {
        expect_true(inherits(attr(al2, 'nodePaths')[[1]], 'SolPath'))

})

test_that('Fixed columns appear in paths', {
    fixed <- attr(al2, 'nodePaths')[[1]]$path[, 'fixedeff']
    expect_true(all(abs(fixed)>0))
})


test_that("Including irrelevant fixed columns doesn't change adjacency (much)", {
    for(i in seq_along(al2$adjMat)){
        expect_lt(mean((al2$adjMat[[i]]-al3$adjMat[[i]])^2), .01)
    }
})


test_that("Including relevant fixed columns changes adjacency", {
    for(i in seq_along(al4$adjMat)){
        expect_true(all(abs(al4$adjMat[[i]][,4])==0))
    }
})



context('Edge interpolation')
test_that('True edges are monotone increasing', {
    ie <- interpolateEdges(fit$adjMat, fit$lambda, nknot=9)
    expect_equal(ie$trueEdges, sort(ie$trueEdges, dec=TRUE))

    ie <- interpolateEdges(al$adjMat, al$lambda, nknot=20)
    expect_equal(ie$trueEdges, sort(ie$trueEdges, dec=TRUE))
    
    ie <- interpolateEdges(al2$adjMat, al2$lambda, nknot=5)
    expect_equal(ie$trueEdges, sort(ie$trueEdges, dec=TRUE))

})
