##' Plot conditional regression line and conditional means
##'
##' @param mapping 
##' @param data 
##' @param stat 
##' @param position 
##' @param na.rm 
##' @param anno.colors 
##' @param tol 
##' @param ... 
##' @return plot
##' @import proto
##' @export
geom_hurdle <- function(mapping = NULL, data = NULL, stat = "identity", position = "identity", na.rm=FALSE, anno.colors=c('blue', 'red'), tol = .01, ...) {
  GeomHurdle$new(mapping = mapping, data = data, stat=stat, position = position, na.rm=FALSE, anno.colors=anno.colors, tol = tol, ...)
}

GeomHurdle <- proto::proto(ggplot2:::GeomPoint, {
  objname <- "hurdle"

  required_aes <- c("x", "y")

  default_aes <- function(.) c(.super$default_aes(), aes(linetype=1, linewidth=1))

  draw <- function(., data, ..., anno.colors=c('blue', 'red'), tol=.01){
      c2d <- calc_hurdle2d(data, tol=tol, color=anno.colors[1])
      cx <- calc_hurdle1d(data, 'x', tol, color=anno.colors[2])
      cy <- calc_hurdle1d(data, 'y', tol, color=anno.colors[2])
      cy <- cbind(x=-.2, xstart=-.2, xend=.2, cy)
      cx <- cbind(y=-.2, ystart=-.2, yend=.2, cx)
#      cbx <- GeomCrossbar$new(
      
      ggname(.$my_name(), grobTree(
          .super$draw(., data, ...),
          GeomLine$draw(c2d, ...),
          GeomSegment$draw(cx, ...),
          GeomSegment$draw(cy, ...)
         ))
  }

  
})

calc_hurdle2d <- function(data, vars=c('x', 'y'), tol, color){
    dv <- data[,vars,drop=FALSE]
    zero <- rowSums(abs(dv)<tol)>0
    data <- data[!zero,]
    Yprime <- lm(y ~ x, data=data)$fitted
    data[,'y'] <- Yprime
    data[,'colour'] <- color
    data[,'linetype'] <- 1
    data[,'size'] <- 1
    data
}

calc_hurdle1d <- function(data, vars, tol, color){
    maybevars <- c('x', 'y')
    othervar <- setdiff(maybevars, vars)
    dv <- data[,vars,drop=FALSE]
    notdv <- data[1,setdiff(names(data), c(maybevars, 'colour'))]
    good <- abs(data[,vars])>tol & abs(data[,othervar])<tol
    m <- mean(dv[good,vars])
    s2 <- mad(dv[good,vars])/sqrt(sum(good))
    #ms2 <- data.frame(x=0, ystart=m-1.96*s2, y=m, ymax=m+1.96*s2)
    ms2 <- data.frame(y=m, ystart=m, yend=m)
    names(ms2) <- str_replace_all(names(ms2), '^y', vars)
#    names(ms2)[1] <-  othervar
    notdv[,'size'] <- 1
    data.frame(notdv, ms2, colour=color)#, linetype=1)
}

##' Wrapper so ggally can use \code{geom_hurdle}
##'
##' .. content for \details{} ..
##' @param data 
##' @param mapping 
##' @param ... 
##' @return plot
##' @export
ggally_hurdle <- function(data, mapping, ...){
                      p <- ggplot(data, mapping)+geom_hurdle(...)
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
