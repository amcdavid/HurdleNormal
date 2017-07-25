##' Plot an igraph object
##'
##' @param gadj igraph object
##' @param layout layout function
##' @param width.scale scale factor for edges
##' @param colorbar \code{logical}, should a legend providing a scale for the edges be displayed? 
##' @param colorbarx coordinates for legend
##' @param colorbary coordinates for legend
##' @param ... passed to plot.igraph
##' @return plots the network 
##' @export
hn_plot <- function(gadj,  layout=layout_with_kk,  width.scale=2, colorbar=FALSE, colorbarx, colorbary, ...){
   
    
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
