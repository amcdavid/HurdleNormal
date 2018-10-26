library(MAST)
library(stringr)
source('common.R')
loadLausanne()
##load('lausanneStability.RData')
Subset <- c('HIV', 'Healthy')#,'TFH016', 'TFH017', 'TFH023', 'TFH041')
fsubset <- expand.grid(model=c('ngo', 'patient'), subset=Subset, stringsAsFactors=FALSE)
fsubset <- subset(fsubset, !((subset %like% 'TFH') & model %in% c('patient', 'ngo')))

samp1 <- subset(samp1[setdiff(names(which(freq(samp1)>.05)), c('CXCR5', 'PDCD1')),], SafeCS%in% c('pp'))
samp10 <- subset(samp10[rownames(samp1),], SafeCS %in% c('pp'))
manyfits <- list()
manystab <- list()
for(i in seq_len(nrow(fsubset))){
    sub <- fsubset[i, 'subset']
    mdl <- fsubset[i, 'model']
    if(sub %in% c('HIV', 'Healthy')){
        thiss1 <- subset(samp1, run==sub)
        thiss10 <- subset(samp10, run==sub)
    } else if(sub =='both'){
        thiss1 <- subset(samp1, run=='Healthy')
        thiss10 <- subset(samp10, run=='Healthy')
    } else{
        thiss1 <- subset(samp1, PatientID==sub)
        thiss10 <- subset(samp10, PatientID==sub)
    }
    if(mdl=='SafeCS'){
        form <- if(sub %like% 'TFH') ~scale(ngeneson) + SafeCS else ~scale(ngeneson) + PatientID + SafeCS
        fixed <- model.matrix(form, colData(thiss1), contrasts=list(PatientID='contr.sum', SafeCS='contr.sum'))
        fixed10 <- model.matrix(form, colData(thiss10), contrasts=list(PatientID='contr.sum', SafeCS='contr.sum'))
    } else if(mdl=='ngo'){
                fixed <- model.matrix(~scale(ngeneson), colData(thiss1), contrasts=list(PatientID='contr.sum'))
                fixed10 <- model.matrix(~scale(ngeneson), colData(thiss10), contrasts=list(PatientID='contr.sum'))

    } else if(mdl=='patient'){
        fixed <- model.matrix(~scale(ngeneson)+PatientID, colData(thiss1), contrasts=list(PatientID='contr.sum'))
        fixed10 <- model.matrix(~scale(ngeneson)+PatientID, colData(thiss10), contrasts=list(PatientID='contr.sum'))
    }
    message(paste(fsubset[i,], collapse=','))
    ee1 <- precenter(exprs(thiss1))
    control <- list(tol=3e-3, maxrounds=300, newton0=T, fista=FALSE, debug=0, returnHessian=FALSE)
    out <- fitSomeModels(thiss1, fixed=fixed, control=control, parallel=T, nlambda=50, lambda.min.ratio=.3)
    out10 <- setNames(fitSomeModels(thiss10, 'gaussianRaw', fixed=fixed10, parallel=T, nlambda=50, lambda.min.ratio=.3), 'gaussian10')
    ee10 <- exprs(thiss10)
    
    manyfits[[i]] <- c(out, out10)

    si <- setupStabilityIndex(ee1, B=30)
    si10 <- setupStabilityIndex(ee10, B=30)
    ss <- list()
    chkpt <- paste0('stab_scenario_', i, '/')
    if(!dir.exists(paste0('g10_', chkpt))){
    ss$hurdle <- stability(ee1, fixed, si, method=fitHurdle, stabilityCheckpointDir=paste0('h_', chkpt), parallel=T, nlambda=50, lambda.min.ratio=.3, control=control, keepNodePaths=FALSE)
    ss$gaussian <- stability(ee1, fixed, si, method=autoLogistic,  stabilityCheckpointDir=paste0('g_', chkpt), parallel=T, nlambda=50, lambda.min.ratio=.2, family='gaussian', keepNodePaths=FALSE)
    ss$logistic <- stability(ee1, fixed, si, method=autoLogistic, stabilityCheckpointDir=paste0('l_', chkpt), parallel=T, nlambda=50, lambda.min.ratio=.2, family='binomial', keepNodePaths=FALSE)
    ss$gaus10 <- stability(ee10, fixed10, si10, method=autoLogistic, parallel=T, nlambda=50, lambda.min.ratio=.3, stabilityCheckpointDir=paste0('g10_', chkpt), family='gaussian', keepNodePaths=FALSE)
}
    cs <- mclapply(str_c(c('h_', 'g_', 'l_', 'g10_'), chkpt), function(x) collectStability(theta=.04, stabilityCheckpointDir=x)[[3]], mc.cores=4)
    manystab[[i]] <- cs
}


mf2 <- do.call(c, manyfits)
dim(mf2) <- c(6, nrow(fsubset))
dimnames(mf2) <- list(names(manyfits[[1]]), do.call(paste, c(fsubset, sep=':')))
outFlat <- sparseCbind(lapply(mf2, '[[', 'edgeInterp'))
names(outFlat) <- outer(rownames(mf2), colnames(mf2), FUN=paste, sep=':')

## manyfits: fit on all data, using 10 cell too
## mf2: concatenation of manyfits
## outFlat: matricized mf2
## manystab 13 x 4 x P xP stabilities
## fsubset: data subset considered
save(manyfits, mf2, outFlat, manystab, fsubset, file='lausanneStability.RData')
