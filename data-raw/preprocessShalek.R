## Source of .rds file below
writeAlex <- function(){
    install_github("RGLab/MAST@MASTClassic")
    install_github("RGLab/MASTdata")
    library(MAST)
    library(MASTDataPackage)
    data(MASTDataPackage)
    tsca_alex <- thresholdSCRNACountMatrix(2^exprs(sca_alex)-1, nbins=20)
    sca_alex <- addlayer(sca_alex, 'et')
    layer(sca_alex) <- 'et'
    exprs(sca_alex) <- tsca_alex$counts_threshold
    sca_alex
}

library(MAST)
u <- 'https://rochester.box.com/shared/static/v06gcp6x88ub8r4f8x1nopqoa9674cj4.rds'
tmp <- tempfile(fileext='.rds')
download.file(u, dest=tmp)
shalekli = readRDS(tmp)
shalek <- with(shalekli, FromMatrix(exprs, cdat, fdat))
subshalek = subset(shalek, Time %in% c('1h', '2h') & Stim == 'LPS')
subshalek = subshalek[freq(subshalek)>.3,]
listified = list(fData=as.data.frame(mcols(subshalek)), cData=as.data.frame(colData(subshalek)), exprs=assay(subshalek))
save(listified, file='data/shalek.RData', compress='xz')
