##
## export  OPENBLAS_NUM_THREADS=1
MC_CORES <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
if(is.na(MC_CORES)) MC_CORES <- 4

options(mc.cores=MC_CORES)

precenter <- function(samp){
    apply(samp, 2, function(x){
        xI <- abs(x)>0
        x[xI] <- scale(x[xI], scale=FALSE)
        x
    })
}

library(devtools)
library(MAST)
library(parallel)
library(Matrix)
library(data.table)
library(RColorBrewer)
library(Mus.musculus)
load_all('../../')

process_aracne <- function(mat, n=50){
    if(any(is.na(mat))) warning(sum(is.na(mat)), " NAs found in aracne estimate")
    mat[is.na(mat)] <- 0
    path <- seq(from=min(mat), to=max(mat), length.out=n)
    adjMat <- lapply(path, function(x) Matrix::Matrix((mat>x)*1,sparse=TRUE))
    list(adjMat=adjMat, lambda=path, BIC=rep(NA, length(path)))
}

fitSomeModels <- function(samp, models=c('hurdle', 'logistic', 'gaussian', 'gaussianRaw', 'aracne'), fixed=NULL, parallel=TRUE, control=list(debug=1, maxrounds=1000, tol=1e-3), lambda.min.ratio=.2, nlambda=100, checkpointDir=NULL){
    samp1c <- t(assay(samp))
    sampCenter <- precenter(t(assay(samp)))
    out <- setNames(as.list(models), models)
    if('hurdle' %in% models){
        hCNgo <- fitHurdle(sampCenter, parallel=parallel, control=control, makeModelArgs=list(scale=FALSE, fixed=fixed, center=TRUE, conditionalCenter=TRUE), lambda.min.ratio=lambda.min.ratio, nlambda=nlambda, checkpointDir=checkpointDir, keepNodePaths=TRUE)
        mCNgo <- interpolateEdges(hCNgo)
        out$hurdle <- mCNgo
        out$hurdle$paths <- attr(hCNgo, 'nodePaths')
    }


    if('logistic' %in% models){
        logfit <- autoGLM(samp1c, fixed=fixed, lambda.min.ratio=lambda.min.ratio, nlambda=nlambda, parallel=parallel)
        out$logistic <- interpolateEdges(logfit)
    }

    if('gaussian' %in% models){        
        hugefit <- autoGLM(sampCenter, fixed=fixed, family='gaussian', nlambda=floor(nlambda*1.2), lambda.min.ratio=lambda.min.ratio/2, parallel=parallel)
        hugemati <- interpolateEdges(hugefit)
        out$gaussian <- hugemati
    }

    if('gaussianRaw' %in% models){        
        hugefitRaw <- autoGLM(samp1c, fixed=fixed, family='gaussian', nlambda=floor(nlambda*1.2), lambda.min.ratio=lambda.min.ratio/2, parallel=parallel)
        hugemati <- interpolateEdges(hugefitRaw)
        out$gaussianRaw <- hugemati
    }

    if('aracne' %in% models){
        resids <- if(!is.null(fixed)) resid(lm(samp1c ~ fixed)) else samp1c
        hugefitRaw <- process_aracne(netbenchmark::aracne.wrap(resids))
        hugemati <- interpolateEdges(hugefitRaw)
        out$aracne <- hugemati
    }
    
    structure(out, colnames=colnames(samp1c), rownames=rownames(samp1c))
}


loadLausanne <- function(freq=.05){
    raw <- readRDS('data_private/AllCells.rds')
    rexprs <- raw$exprs
    colnames(rexprs) <- raw$fdata$shortName
    row.names(raw$fdata) <- raw$fdata$primerid <- raw$fdata$shortname
    fl <- FromMatrix("FluidigmAssay", exprsArray = t(rexprs), cData = raw$cdata, fData = raw$fdata)
    frq1 <- freq(subset(fl, ncells==1))
    fl <<- fl[frq1>freq,]
    samp10 <<- subset(fl, ncells==10)
    samp1 <<- subset(fl, ncells==1)
}

loadMait <- function(){
tsca_mait <- thresholdSCRNACountMatrix(2^exprs(sca_mait)-1, nbins=20)
sca_mait <- addlayer(sca_mait, 'et')
layer(sca_mait) <- 'et'
exprs(sca_mait) <- tsca_mait$counts_threshold
cData(sca_mait) <- data.frame(cData(sca_mait), ngeneson=cData(sca_mait)$nGeneOn)
colnames(sca_mait) <- make.unique(fData(sca_mait)$symbolid)
sca_mait
}

namedrbindlist <- function(li, nm=NULL, fill){
    if(!is.null(nm)) names(li) <- nm
    li <- lapply(names(li), function(linm){
               x <- li[[linm]]
               x <- data.table(x)
               x[,L1:=linm]
               x
               })
    rbindlist(li, fill=fill)

}


genNetwork <- function(hfit, nedge=40, min.components=3,min.degree.label=2,vnames, select, plot=TRUE, ...){
    ##interp <- length(hfit[[2]])-findInterval(nedge, rev(unlist(hfit[[2]])))
    interp <- which.min(abs(nedge-hfit$trueEdges))

    thisMat <- hfit[[1]][[interp]]
    if(!missing(vnames)) dimnames(thisMat) <- list(vnames, vnames)
    gadj <- graph.adjacency(abs(thisMat)>0, mode='max')
    sizes <- Matrix::summary(thisMat)
    sizes <- subset(sizes, i>j)[,'x']
    E(gadj)$width <- (abs(sizes)/max(abs(sizes)))^(1/4)
    E(gadj)$sign <- sign(sizes)
    clu <- components(gadj)
    if(!missing(select)){
        goodclu <- clu$membership[select]
        badclu <- setdiff(unique(clu$membership), goodclu)

    } else{
            badclu = which(clu$csize<min.components)
        }

    badv = clu$membership %in% badclu
    gadj <- delete.vertices(gadj, badv)
    
    #gadj <- igraph::simplify(gadj, edge.attr.comb='mean')
    V(gadj)$size <- igraph::degree(gadj)/sqrt(max(igraph::degree(gadj)))

    Npal <- 50
    pal <- colorRampPalette(brewer.pal(9, 'BrBG'))(Npal)
    wt <- E(gadj)$width*E(gadj)$sign
    minNeg <- min(wt[wt<0], -1e-3)*1.05
    maxPos <- max(wt[wt>0], 1e-3)*1.05
    wtkey <- c(seq(minNeg, 0, length.out=Npal/2), seq(0, maxPos, length.out=Npal/2))
    maxAll <- max(-minNeg, maxPos)
    wtkey <- seq(-maxAll, maxAll, length.out=Npal)
    paltable <- data.table(color=pal,wtkey=wtkey)
    setkey(paltable, wtkey)
    E(gadj)$color <- paltable[list(wt),color,roll=TRUE]
    V(gadj)$frame.color='blue'
    V(gadj)$label <- ifelse(igraph::degree(gadj)>min.degree.label, V(gadj)$name, '')
    V(gadj)$label.color='black'
    
    if(plot) plotNetwork(gadj, ...)
    
    structure(gadj, scale=paltable)
}

color.bar <- function(col, val, nticks=11) {
    min <- min(val)
    max <- max(val)
    scale = (length(col)-1)/(max-min)
    #dev.new(width=1.75, height=5)
    plot(c(0,2), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', bg='white')
    ticks=sort(c(0, seq(min, max, len=nticks-1)))
    #labels <- findInterval(ticks, val)
    axis(2, round(ticks, 1), las=1)
    for (i in 1:(length(col)-1)) {
     y = (i-1)/scale + min
     rect(0,y,2,y+1/scale, col=col[i], border=NA)
    }
}

plotNetwork <- function(gadj,  layout=layout_with_kk,  width.scale=2, colorbar=FALSE, colorbarx, colorbary, ...){
   
    
    plot(gadj,				#the graph to be plotted
         layout=layout,	# the layout method. see the igraph documentation for details
         vertex.label.dist=0.1,			#puts the name labels slightly off the dots
         #vertex.frame.color='blue', 		#the color of the border of the dots 
         vertex.label.font=1,			#the font of the name labels
         #vertex.label=,		#specifies the lables of the vertices. in this case the 'name' attribute is used
         vertex.label.cex=1,			#specifies the size of the font of the labels. can also be made to vary,
         edge.width=E(gadj)$width*width.scale,
         #edge.color=E(gadj)$color,
         ...)

    if(colorbar){
        scale <- attr(gadj, 'scale')
        if(missing(colorbarx)){
            colorbarx <- 1.2
            colorbary <- .9
            }
        #fields::colorbar.plot(x=colorbarx, y=colorbary, strip=scale[,wtkey], col=scale[,color], strip.width=.02, strip.length=.5)
        Hmisc::subplot(color.bar(scale[,color], scale[,wtkey]), colorbarx, colorbary, size=c(.1, 1.5))
        
        }

}

    
goCriteria <- c('GO:0003677', 'GO:0006355', #transcription factors
                'GO:0006351', #promoters
                'GO:0016070', #mRNA metabolic process
                'GO:0001190', #TF binding
                'GO:0017053') #transcription repressor
                


## Change a phrase to Capital letters
## Not vectorized.
simpleCap <- function(x) {
    x <- tolower(x)
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1, 1)), substring(s, 2),
          sep = "", collapse = " ")
     }

## Select a graph from the interpolation object
getNedgeInterp <- function(ginterp, n){
    gidx <- which.min(abs(ginterp$trueEdges-n))
    ginterp$edgeInterp[[gidx]]
}

## Return a data.table crosswalk of vertex numbers and names
getGraphAlias <- function(graph){
    alias <-  sapply(row.names(graph), simpleCap)
    desc <- data.table(ALIAS=alias)
    desc$i <- seq_along(desc$ALIAS)
    desc
}

##' Join vertex indices, index-ALIAS crosswalk and various databases
##'
##' Merge the godb with a list of vertex indices.  Multiple GOIDs can be assigned to each ALIAS
##' (it is a bipartite matching)
##' @param idx data.table
##' @param idxcol for joining idx and alias
##' @param alias idxcol-ALIAS table
##' @param godb ALIAS-GOID table.  Object of flavor Orgdb
##' @param evidence evidence codes to subset godb
##' @param otherdb (optional) data.table with other ALIAS-GOID records. Joined and unioned.
##' @return data.table of outer join
mergeWithAliasAndGo <- function(idx, idxcol, alias, go, evidence=c('EXP', 'IDA', 'PIP', 'IMP', 'IGI', 'ISS', 'ISO', 'ISA', 'ISM', 'IGC', 'IKR', 'IRD', 'IBA', 'IBD', 'MAR')){
    idx <- merge(idx, alias, by=idxcol)
    go <- go[!is.na(GOID) & EVIDENCE %in% evidence] #excludes some computational evidences.  see http://geneontology.org/page/guide-go-evidence-codes
    setkey(go, ALIAS, GOID)
    go <- unique(go)
    setkey(go, NULL)
    idxgo <- merge(idx, go, by='ALIAS', allow.cartesian=TRUE)
    idxgo[,N_GOID:=.N, keyby=GOID]
    idxgo
}

logodds <- function(tab){
    tab[1,1]*tab[2,2]/(tab[1,2]*tab[2,1])
    }


## Also check degree and clique completion
## Need a general way to merge GO with graph
## graphstatlike -> alias -> go -> stat
edgeGoAnno <- function(ginterp, n, graph, goAlias, goTerm, background, annotate=TRUE, nulldist, ...){
    ## get edge set
    if(missing(graph)){
        graph <- getNedgeInterp(ginterp, n)
    }
    eset <- Matrix::summary(graph)
    setDT(eset)
    
    ## Number of edges in complete graph vs sampled graph
    ncomplete <- nrow(graph)^2-nrow(graph)
    nedge <- nrow(eset)

    
    ## join vertex i to GO
    alias <- getGraphAlias(graph)
    eseti <- mergeWithAliasAndGo(eset, 'i', alias, goAlias, ...)

    ## merge again  vertex j
    setnames(eseti, c('EVIDENCE', 'GOID', 'N_GOID', 'ALIAS', 'ONTOLOGY'), c('EVIDENCE.i', 'GOID.i', 'N_GOID.i', 'ALIAS.i', 'ONTOLOGY.i'))
    setnames(alias, 'i', 'j')
    esetij <- mergeWithAliasAndGo(eseti, 'j', alias, goAlias, ...)
    setnames(esetij, c('EVIDENCE', 'GOID', 'N_GOID', 'ALIAS', 'ONTOLOGY'), c('EVIDENCE.j', 'GOID.j', 'N_GOID.j', 'ALIAS.j', 'ONTOLOGY.j'))

    ## Count go pairs
    gocount <- esetij[,.(N_GOID.ij=.N),keyby=list(GOID.i, GOID.j)]
    
    ## compare to this sample to number of edges in complete graph
    ## calculate tail prob.
    gcbg <- background[gocount]
    gcbg[,phyper:=phyper(N_GOID.ij, Nij, ncomplete-Nij, nedge, lower.tail=FALSE) ,keyby=Nij] #+ dhyper(N_GOID.ij, Nij, ncomplete-Nij, nedge)
    setkey(gcbg, phyper)

    if(annotate){
        gcbgSig <- gcbg[1:nrow(nulldist)][,rank:=1:nrow(nulldist)]
        gcbgSig[,pperm:=mean(c(1, nulldist[rank,]<phyper)), keyby=rank]
        gcbgSig[,phyperAdj:=p.adjust(phyper, 'BY', n=nrow(background))]

        ## join to go terms
        setkey(goTerm, GOID)
        goTerm <- goTerm[unique(gcbgSig$GOID.i)]
        gcbgSig <- merge(gcbgSig, goTerm[,.(GOID, TERM)], by.x='GOID.i', by.y='GOID', all.x=TRUE)
        gcbgSig <- merge(gcbgSig, goTerm[,.(GOID, TERM)], by.x='GOID.j', by.y='GOID', all.x=TRUE)

        ## get rid of duplicates (both pairs were concerned before)
        gcbgSig[, GOIDmin:=pmin(as.character(GOID.i), as.character(GOID.j))]
        gcbgSig[, GOIDmax:=pmax(as.character(GOID.i), as.character(GOID.j))]
        setkey(gcbgSig, GOIDmin, GOIDmax)
        gcbgSig <- unique(gcbgSig)
        gcbgSig[,':='(GOIDmin=NULL, GOIDmax=NULL)]

        gcbgSig[,':='(TERM.x=ifelse(is.na(TERM.x), GOID.i, TERM.x),
             TERM.y=ifelse(is.na(TERM.y), GOID.j, TERM.y))]
        gcbgSig[,':='(TERM.x=stringr::str_wrap(TERM.x, width=35),
                      TERM.y=stringr::str_wrap(TERM.y, width=35))]
        
        setkey(gcbgSig, phyperAdj)
        return(list(gcbgSig, goTerm, esetij[GOID.i%in% goTerm$GOID]))
    } else{
        return(gcbg)
    }
}

plotgenesee <- function(x, godb, goTerm,network, goids, additionalGenes=NULL){
    L1 <- x$L1[1]
    geneseeG <- igraph::graph.edgelist(as.matrix(x[,.(GOID.i, GOID.j)]), directed=FALSE)
    E(geneseeG)$weights <- sqrt(-log10(x$phyper)-6)
    E(geneseeG)$group <- x$L1
    if(missing(goids)){
        goids <- unique(x[Nij<7e4][1:9, c(GOID.i, GOID.j)])
    } else{
        goids <- unique(goids[GOID  %in% x$GOID.j |  GOID %in% x$GOID.i])
    }
    ## assign colors
    ncategories <- length(unique(goids$category))
    catColor <- goids[,color:=brewer.pal(ncategories, 'Set2')[.GRP],keyby=category]
    goTerm <- catColor[goTerm,,on="GOID"]
    goTerm[is.na(category), color:='white']
    goTerm[is.na(category), category:='Uncategorized']
    
    ## join genesee GOIDs to godb (and the genes)
    setkey(godb, GOID)
    sgo <- godb[goids[,.(GOID)],,on='GOID']
    sgo <- goTerm[sgo,,on='GOID']
     if(!'category' %in% names(sgo)){
        sgo[,category:=GOID]
    }
    sgo[,V:=toupper(ALIAS)]
    #sgo[,N:=.N, key=GOID]
    #sgo <- sgo[N<300,]

    ## top 10 terms no equal to unassociated set
    gselect <- unique(c(sgo$V, additionalGenes))
    networkG <- genNetwork(network, 1400, plot=F, select=gselect)
    networkG <- delete.vertices(networkG, igraph::degree(networkG)<1)
    V(networkG)$size <- 3
    E(networkG)$width <- 1
    ## select vertices in graph
    setkey(sgo, V)
    sgoGraph <- sgo[names(V(networkG)),,nomatch=NA, mult='first']
    sgoGraph[is.na(color), color:='white']
    ## ## add colors and plot
    ## setkey(sgoGraph, category)
    
    ## paltable <- unique(sgoGraph[, .(category, TERM)])
    ## paltable <- paltable[sgoGraph,,on='category']
    ## #paltable <- paltable[catColor,,on='GOID']
    ## #paltable[is.na(category), ':='(color='white', category='Uncategorized')]
    ## setkey(paltable, V)
    ## paltable <- promoteWarning(paltable[names(V(networkG)),,mult='first'], 'not a multiple of replacement length')
    paltable <- sgoGraph
    V(networkG)$color <- paltable[,color]
    V(networkG)$label <- paltable[,ifelse(!is.na(category) | V %in% additionalGenes, V, '')]
    V(networkG)$label.color <- paltable[,ifelse(V %in% additionalGenes, 'red', 'black')]
    plotNetwork(networkG, layout=layout_with_fr, vertex.label.cex=1)
    setkey(goTerm, color)
    leg <- unique(goTerm)
    legend('bottomright', fill=leg[,color], legend=leg[,category])

    setkey(goTerm, GOID)
    geneseeGtable <- goTerm[names(V(geneseeG)),,nomatch=NA, mult='first']
    V(geneseeG)$label <- ifelse(is.na(geneseeGtable$TERM), '', stringr::str_wrap(geneseeGtable$TERM, 40))
    plot(geneseeG, vertex.size=3, vertex.label.dist=.2, edge.width=E(geneseeG)$weights, layout=layout_with_kk, vertex.color=geneseeGtable$color, vertex.label.cex=1)
    legend('bottomright', fill=leg[,color], legend=leg[,category])

}
