## sbatch -n1 -c16 -t1-01 R --vanilla -f fitAlex.R
source('common.R')
load(system.file('data-raw', 'shalek.RData', package='HurdleNormal'))
sca_alex <- with(shalek, FromMatrix(exprs, cData, fData))
lps_pam_freq <- freq(MAST::subset(sca_alex, Time =='2h' & Stim %in% c('LPS', 'PAM')))
astim <- MAST::subset(sca_alex, Stim=='LPS' & Time=='2h')
astim <- astim[lps_pam_freq>.3,]

## Fit networks
alex_fits <- list()
parallel=TRUE

alex_fits[['regLPS']] <- fitSomeModels(astim,  lambda.min.ratio=.4, nlambda=70, control=list(tol=5e-3, maxrounds=80, FISTA=FALSE, debug=0, newton0=TRUE), parallel=parallel, checkpointDir='alex2hlps_chk', keepNodePaths=FALSE)
gc()

alex_fits[['ngoLPS']] <- fitSomeModels(astim,  lambda.min.ratio=.4, nlambda=70, control=list(tol=5e-3, maxrounds=80, FISTA=FALSE, debug=0, newton0=TRUE), parallel=parallel, checkpointDir='alex2hlps_ngo_chk', fixed=cbind(1, scale(colData(astim)$ngeneson)), keepNodePaths=FALSE)
gc()


## Set up GSEA db
## Calculate p values for marginal changes
testbystim <- lapply(list('LPS'),function(s){
    zz <- zlm(~Time+ngeneson, data=sca_alex[rownames(astim),] , method='bayesglm')
    wt <- lrTest(zz, 'Time')
    data.table(ALIAS=sapply(rownames(astim), simpleCap), pval=wt[,'hurdle', 'Pr(>Chisq)'], sign=sign(coef(zz, 'D')), stim=s)
})

testbystim2 <- rbindlist(testbystim, use.names=TRUE)
testbystim2[, fdr:=p.adjust(pval, 'bonferroni'),keyby=stim]
testbystim2[, cutfdr:=as.character(cut(fdr, c(0, 1e-3, 1e-2, 1e-1, 1), labels=c('***', '**', '*', "-"))),keyby=stim]
testbystim2[,GOID:=stringr::str_c(stim,':',  cutfdr, ':', sign.Time2h)]

## Count null populations for test
k <-sapply(rownames(astim), simpleCap)
godb<- select(Mus.musculus, keys=k, columns='GOID', key='ALIAS')
setDT(godb)
godb$GOID <- as.character(godb$GOID)
setkey(godb, GOID, ALIAS)
godb<- unique(godb)
godb<- rbind(godb, cbind(testbystim2[, .(GOID, ALIAS)], EVIDENCE='MAR', ONTOLOGY='BP'), fill=TRUE)
godb[,N:=.N, key=GOID]
godb<- godb[N>2 &!is.na(GOID)]
setkey(godb, GOID)
sp <- godb[,.(.N), key=GOID]
goids <- unique(sp$GOID)
## this overestimates the number of pairs by intersection(gene(i), genes(j)), but probably only matters for self-pairs
nullgo <- CJ(GOID.i=goids, GOID.j=goids)[, Nij:=as.vector(outer(sp$N, sp$N))]
nullgo[GOID.i==GOID.j, Nij:=Nij-sqrt(Nij)]

## Save go terms as data.table
goTerm <- select(Mus.musculus, k=sp$GOID, columns=c('GOID', 'TERM', 'DEFINITION'), keytype='GOID')
setDT(goTerm)
goTerm$GOID <- as.character(goTerm$GOID)

## Do permutation test
graph <- alex_fits[[1]][[1]]$edgeInterp[[1]]
N <- 700
np <- 4000
nedge1400perm <- mclapply(seq_len(np), function(ip){
        graph[] <- 0
        P <- ncol(graph)
        i <- sample(P*(P-1), N)
        graph[i] <- 1
        graph <- graph+Matrix::t(graph)
        out <-  edgeGoAnno(graph=graph, goAlias=godb, goTerm=goTerm, background=nullgo, annotate=FALSE)
        out[1:floor(np/20), phyper]
})
nedge1400perm <- do.call(cbind, nedge1400perm)
dimnames(nedge1400perm) <- list(rank=1:floor(np/20), rep=1:np)

save(alex_fits, astim, nullgo, godb, goTerm, nedge1400perm, file='alexNetworks.RData')
