source('simulation_library.R')
writeModels <- TRUE
fitModels <- TRUE
Pseq <- c(16, 32, 64, 128)
modelArgs <- expand.grid(type=c('sparse', 'dense'), P=Pseq, N=100, sparsity=.015, kcell=c(1, 10), contam=FALSE, stringsAsFactors=FALSE)
modelArgs <- rbind.fill(modelArgs, expand.grid(type=c('sparse', 'dense'), P=16, N=c(100, 500, 5000, 5e4), sparsity=.015, kcell=1, contam=c(TRUE, FALSE), stringsAsFactors=FALSE))
modelArgs = rbind.fill(modelArgs, expand.grid(type='ecoli', N=c(100, 500, 1000), kcell=c(1, 10), stringsAsFactors = FALSE))
modelList <- list()
for(i in seq_len(nrow(modelArgs))){
    thismodel <- modelArgs[i,'type']
    thisP <- modelArgs[i, 'P']
    thisN <- modelArgs[i, 'N']
    kcell <- modelArgs[i,'kcell']
    thisSparsity <- modelArgs[i, 'sparsity']
    contam = if(isTRUE(modelArgs[i,'contam'])) contaminate else NA
    if(thismodel=='sparse'){
        modelList[[i]] <- HurdleNormal:::simulateHurdle210(N=thisN, thisP, c('G'), structure='chain', structureArgs=list(sparsity=thisSparsity), intensityArgs=list(G=1), Gdiag=-4.5, Hdiag=2, gibbs_args=list(kcell=kcell, post_fun=contam))
    } else if(thismodel=='dense'){
        modelList[[i]] <- HurdleNormal:::simulateHurdle210(N=thisN, thisP, c('G', 'Hupper', 'Hlower', 'K'), structure='chain', structureArgs=list(sparsity=thisSparsity), intensityArgs=list(G=-.25,K=-.4, Hupper=-.75, Hlower=-.75), Gdiag=-.5, Hdiag=1, gibbs_args=list(kcell=kcell, post_fun=contam))
    } else if(thismodel=='ecoli'){
        
    } else{
        stop('bad model name')
    }
    }
}


if(writeModels){
    modelArgs <- rbind(modelArgs, addnArgs)
    save(modelList, modelArgs, file='simulation_graphs.RData')
}

if(fitModels){
    fittedmodels <- lapply(modelList, function(m){
        mclapply(1:30, function(k) fitAndSummarize(model=m, n=nrow(m$gibbs), parallel=FALSE))
    })
    save(fittedmodels, file='simulations.RData')
}
