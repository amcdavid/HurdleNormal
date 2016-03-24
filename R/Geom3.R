stat_hurdle2d <- function(mapping = NULL, data = NULL, geom = "line",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE,tol=.01, ...) {
  layer(
    stat = StatHurdle2d, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, tol=tol, ...)
  )
}

StatHurdle2d <- ggplot2::ggproto('StatHurdle2d', ggplot2::Stat, compute_group=function(data, scales, tol){
vars=c('x', 'y')
    dv <- data[,vars,drop=FALSE]
    zero <- rowSums(abs(dv)<tol)>0
    data <- data[!zero,]
    Yprime <- lm(y ~ x, data=data)$fitted
    data[,'y'] <- Yprime
    data
})

StatHurdle1d <- ggplot2::ggproto("StatHurdle1d", ggplot2::Stat,
                        compute_group = function(data, scales, vars, tol){
                            maybevars <- c('x', 'y')
                            othervar <- setdiff(maybevars, vars)
                            dv <- data[,vars,drop=FALSE]
                            notdv <- data[1,setdiff(names(data), maybevars)]
                            good <- abs(data[,vars])>tol & abs(data[,othervar])<tol
                            m <- mean(dv[good,vars])
                            s2 <- mad(dv[good,vars])/sqrt(sum(good))
                            segloc <- data.frame(y=m, yend=m)
                            names(segloc) <- str_replace_all(names(segloc), '^y', vars)
                            segfix <- data.frame(x=-.2, xend=.2)
                            names(segfix) <- str_replace_all(names(segfix), '^x', othervar)
                            data.frame(segfix, segloc, notdv)#, linetype=1)
                        },
                        required_aes = c("x", "y")
                        )


stat_hurdle1d <- function(mapping = NULL, data = NULL, geom = "segment",
                       position = "identity", na.rm = FALSE, show.legend = NA, 
                       inherit.aes = TRUE,vars='x', tol=.01, ...) {
  layer(
    stat = StatHurdle1d, data = data, mapping = mapping, geom = geom, 
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, vars=vars, tol=tol, ...)
  )
}

##' Wrapper so ggally can use \code{geom_hurdle}
##'
##' .. content for \details{} ..
##' @param data 
##' @param mapping 
##' @param ... 
##' @return plot
##' @export
ggally_hurdle <- function(data, mapping, lwd.regression=1, lwd.axis=2, size.point=1, ...){
    p <- ggplot(data, mapping)+geom_point(size=size.point, ...)+stat_hurdle1d(vars='x', lwd=lwd.axis, color='red', ...)+stat_hurdle1d(vars='y', lwd=lwd.axis, color='red', ...)+stat_hurdle2d(lwd=lwd.regression, color='blue', ...)
    p$type <- 'continuous'
    p$subType <- 'hurdle'
    p
                  }

ggally_hmosaic <- function(data, mapping, ...){
    var1 <- (data[,deparse(mapping$x)]>0)*1
    var2 <- (data[,deparse(mapping$y)]>0)*1
    nm1 <- names(data)[2]
    nm2 <- names(data)[1]
    jointTable <- prop.table(table(var1, var2))
  ##  ## # mx <- colSums(jointTable)
  ##  ##  my <- rowSums(jointTable)
  ##   plotData <- as.data.frame(jointTable)
  ## plotData$marginVar1 <- prop.table(table(var1))
  ## plotData$var2Height <- plotData$Freq / plotData$marginVar1
  ## plotData$var1Center <- c(0, cumsum(plotData$marginVar1)[1:levVar1 -1]) +
  ##   plotData$marginVar1 / 2
  ##   browser()

  ## ggplot(plotData, aes(x=var1Center, y=var2Height)) +
  ##   geom_bar(stat = "identity", aes(width = marginVar1, fill = var2), col = "Black") + coord_flip()

    widths <- c(0, cumsum(apply(jointTable, 1, sum)))
    heights <- apply(jointTable, 1, function(x){c(0, cumsum(x/sum(x)))})

  alldata <- data.frame()
  allnames <- data.frame()
  for(i in 1:nrow(jointTable)){
    for(j in 1:ncol(jointTable)){
      alldata <- rbind(alldata, c(widths[i], widths[i+1], heights[j, i], heights[j+1, i]))
    }
  }
  colnames(alldata) <- c("xmin", "xmax", "ymin", "ymax")

  alldata[[nm1]] <- rep(dimnames(jointTable)[[1]],rep(ncol(jointTable), nrow(jointTable)))
  alldata[[nm2]] <- rep(dimnames(jointTable)[[2]],nrow(jointTable))

  p <- ggplot(alldata, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)) + 
    geom_rect(color="black", aes_string(fill=nm1)) + scale_x_continuous(breaks=NULL)+scale_y_continuous(breaks=NULL)# +
    #xlab(paste(xvar, "(count)")) + ylab(paste(yvar, "(proportion)"))
    p$type <- 'continuous'
    p$subType <- 'hmosaic'
    p
    
    
}
