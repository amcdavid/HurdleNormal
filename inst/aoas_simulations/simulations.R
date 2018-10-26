MC_CORES <- floor(as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE')))
if(is.na(MC_CORES)) MC_CORES <- 4

this_job <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
print(this_job)
if(is.na(this_job)) this_job <- 1:2
MODELS_PER_ARRAY <- 3
options(mc.cores=floor(MC_CORES/MODELS_PER_ARRAY))

devtools::load_all('../..')
source('common.R')
source('simulation_library.R')
writeModels <- !file.exists("simulation_graphs.RData")
fitModels <- TRUE

if(writeModels){
    ecoli_gnw <- readRDS('ecoli_network.rds')
    Pseq <- c(16, 32, 64, 128)
    Nseq <- c(100, 500, 2500, 12500)
    modelArgs <- expand.grid(type=c('sparse', 'dense'), P=Pseq, N=100, sparsity=.015, kcell=c(1, 10), contam='none', stringsAsFactors=FALSE)
    modelArgs <- rbind.fill(modelArgs, expand.grid(type=c('sparse', 'dense'), P=24, N=Nseq, sparsity=.05, kcell=1, contam=c('none', 't'), stringsAsFactors=FALSE))
    modelArgs <- rbind.fill(modelArgs, expand.grid(type=c('dense', 'ecoli'), P=24, N=Nseq[1:3], sparsity=.05, kcell=1, contam='selection', stringsAsFactors=FALSE))
    modelArgs = rbind.fill(modelArgs, expand.grid(type='ecoli', N=Nseq[1:3],
                                                   kcell=c(1, 10), contam='none', stringsAsFactors = FALSE))
    modelArgs <- modelArgs[nrow(modelArgs):1,]
    modelArgs$i <- 1:nrow(modelArgs)
    ## modelArgs <- modelArgs[c(8, 23, 24, 39),]
    modelList <- list()
    for(i in seq_len(nrow(modelArgs))){
        thismodel <- modelArgs[i,'type']
        thisP <- modelArgs[i, 'P']
        thisN <- modelArgs[i, 'N']
        kcell <- modelArgs[i,'kcell']
        thisSparsity <- modelArgs[i, 'sparsity']
        contam = switch(modelArgs[i,'contam'], t=contaminate, none=NA, selection=selectionModel)
        mvnorm <- modelArgs[i, 'contam']=='selection'
        if(thismodel=='sparse'){
            modelList[[i]] <- HurdleNormal:::simulateHurdle210(N=thisN, thisP, c('G'), structure='chain', structureArgs=list(sparsity=thisSparsity), intensityArgs=list(G=1), Gdiag=-4.5, Hdiag=2, gibbs_args=list(mvnorm=mvnorm, kcell=kcell, post_fun=contam))
        } else if(thismodel=='dense'){
            modelList[[i]] <- HurdleNormal:::simulateHurdle210(N=thisN, thisP, c('G', 'Hupper', 'Hlower', 'K'), structure='chain', structureArgs=list(sparsity=thisSparsity), intensityArgs=list(G=-.25,K=-.4, Hupper=-.75, Hlower=-.75), Gdiag=-.5, Hdiag=1, gibbs_args=list(mvnorm=mvnorm, kcell=kcell, post_fun=contam))
        } else if(thismodel=='ecoli'){
            ecoli_gnw$gibbs_args <- list(mvnorm=mvnorm, kcell=kcell, post_fun=contam)
            modelList[[i]] <- getGibbs(ecoli_gnw, Nkeep=thisN, thin=1)
        } else{
            stop('bad model name')
        }
        if(mvnorm){
            modelList[[i]]$true$G[] <-  modelList[[i]]$true$H[] <- 0
        }
    }
    save(modelList, modelArgs, file='simulation_graphs.RData')
}

if(fitModels){
    load('simulation_graphs.RData')
    lapply(this_job, function(k){
        set.seed(k+42450194)
        fittedmodels <- mclapply(modelList, function(m){
            message('Model dim= ', paste0(dim(m$gibbs), collapse='x'))
            fitAndSummarize(model=m, n=nrow(m$gibbs), parallel=TRUE)
        }, mc.cores=MODELS_PER_ARRAY, mc.preschedule=FALSE)
        saveRDS(fittedmodels, file.path('sim_chkpoint', sprintf('sim_%02d.rds', k)))
    })
    
     ## fittedmodels <- lapply(modelList, function(m){
     ##    fitAndSummarize(model=m, n=nrow(m$gibbs), parallel=TRUE)
     ## })
    
}
