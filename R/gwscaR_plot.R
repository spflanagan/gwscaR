#Author: Sarah P. Flanagan
#Purpose: Plot genome-wide statistics and perform other population genomics visualizations
#Package: gwscaR


#' Plot genome-wide statistics from a data frame.
#' @param fst.dat The data.frame containing at least three columns: the statistic to be plotted, the chromosome ID, and the BP ID for each locus. Each row is a locus.
#' @param pt.cols The color of the points. The default is c("darkgrey","lightgrey").
#' @param fst.name The name of the column containing the statistic to be plotted. Default is "Fst".
#' @param chrom.name The name of the column containing the chromosome information for each locus. Default is "Chrom".
#' @param bp.name The name of the column containing the basepair information for each locus. Default is "BP".
#' @param axis.size The value of cex.axis passed to the plotting of the y-axis. Default is 0.5. Set axis.size = 0 to suppress plotting of the y-axis.
#' @param scaffs.to.plot A vector or list containing the names of the chromosomes/scaffolds you want to include in the order in which you would like them to be plotted.
#' @param y.lim The limits for the y-axis.
#' @param pt.cex The size of the points
#' @param xlabels Either a boolean (TRUE) or a vector of labels for the x-axis
#' @param xlab.indices If the xlabels are not the same as the names of the chromosomes in fst.dat, then provide xlab.indices, which is a vector of integers that are the indices in scaffs.to.plot for which chromosome sets you would like to label.
#' @param scaffold.widths A data.frame with first column naming the chromosomes/linkage groups, second column having the last number of loci on the chromosome for indexing purposes.
#' e.g.:
#'  Chrom   End
#'  LG1     3945850
#'  LG2     435215
#' This parameter is used to make multiple plots from the same overall dataset have the same widths.
#' @return new.dat The fst.dat data frame with new values in BP, scaled to be the sequence in which points are plotted based on their position on the chromosome and the order the chromosomes are plotted in.
#' @seealso Flanagan & Jones 2017
#' @export
fst.plot<-function(fst.dat,scaffold.widths=NULL,scaffs.to.plot=NULL,
                   fst.name="Fst", chrom.name="Chrom", bp.name="BP",
                   y.lim=NULL,xlabels=NULL,xlab.indices=NULL,
                   axis.size=0.5,pt.cols=c("darkgrey","lightgrey"),pt.cex=0.5,...){
  #make sure fsts are numeric
  fst.dat[,fst.name]<-as.numeric(as.character(fst.dat[,fst.name]))
  #if scaffold widths aren't provided, we add them ourselves
  if(is.null(scaffold.widths)){
    scaffold.widths<-data.frame(Chrom=levels(as.factor(fst.dat[,chrom.name])),
                                Max=tapply(as.numeric(as.character(fst.dat[,bp.name])),
                                           fst.dat[,chrom.name],max))
  }
  colnames(scaffold.widths)<-c("Chrom","Max")
  #if no specific scaffolds are provided, we plot them all
  if(is.null(scaffs.to.plot)){
    scaffs.to.plot<-scaffold.widths[,1]#unique(fst.dat[,chrom.name])
  }
  
  #set up new dataframe
  new.dat<-data.frame(stringsAsFactors = F)
  last.max<-0
  for(i in 1:length(scaffs.to.plot)){
    #pull out the data for this scaffold
    if(nrow(scaffold.widths[scaffold.widths$Chrom %in% scaffs.to.plot[i],])>0){ #sanity check
    chrom.dat<-fst.dat[fst.dat[,chrom.name] %in% scaffs.to.plot[i],]
    if(nrow(chrom.dat)>0){
      chrom.dat$plot.pos<-as.numeric(as.character(chrom.dat[,bp.name]))+last.max
      new.dat<-rbind(new.dat,chrom.dat)
      #last.max<-max(chrom.dat$plot.pos)+
      #               as.numeric(scaffold.widths[scaffold.widths[,1] %in% scaffs.to.plot[i],2])
    }
    last.max<-last.max+
      as.numeric(scaffold.widths[scaffold.widths$Chrom %in% scaffs.to.plot[i],2])
    }else{
      print(paste(scaffs.to.plot[i], "has no designated width."))
    }
  }
  #make sure everything is the correct class
  new.dat$plot.pos<-as.numeric(as.character(new.dat$plot.pos))
  
  new.dat[,chrom.name]<-as.factor(as.character(new.dat[,chrom.name]))
  scaffs.to.plot<-as.factor(as.character(scaffs.to.plot))
  #determine the axis limits
  x.max<-max(new.dat$plot.pos,na.rm=T)
  x.min<-min(new.dat$plot.pos,na.rm=T)
  if(is.null(y.lim)){
    y.max<-max(fst.dat[,fst.name])+0.1*max(fst.dat[,fst.name])
    y.min<-min(fst.dat[,fst.name])-0.1*min(fst.dat[,fst.name])
    if(min(fst.dat[,fst.name]) < 0) {
      y.min<-min(fst.dat[,fst.name]) - 0.1*min(fst.dat[,fst.name])
    } else {
      y.min<-0
    }
    
    y.lim<-c(y.min,y.max)
  }
  
  #plot
  colors<-data.frame(lg=as.character(new.dat[,chrom.name]),col=rep(pt.cols[1],nrow(new.dat)),
                     stringsAsFactors = F)
  colors[as.numeric(factor(colors$lg,levels=scaffs.to.plot))%%2==0,"col"]<-pt.cols[2] #defining the levels maintains the order
  plot(new.dat$plot.pos,new.dat[,fst.name],
       xlim=c(x.min,x.max),ylim=y.lim,...,
       cex=pt.cex, col=colors$col,
       bty="n",axes=F, xlab="", ylab="")
  #optionally add the y-axis
  if(axis.size>0){
    axis(2, at = seq(round(y.lim[1],2),round(y.lim[2],2),
                     round((y.lim[2]-y.lim[1])/2, digits=2)),
         ylim =y.lim, pos=0,
         labels=seq(round(y.lim[1],2),round(y.lim[2],2),
                    round((y.lim[2]-y.lim[1])/2, digits=2)),
         las=1, xlab="", ylab="", cex.axis=axis.size)
  }
  #optionally add the x-axis
  if(!is.null(xlabels)){
    #is this a boolean or a list?
    if(class(xlabels)=="logical"){
      xlabels<-as.character(scaffs.to.plot)
    }
    #if there are not indices, they should just be the ones that match scaffolds
    if(is.null(xlab.indices)){
      xlab.indices<-which(scaffs.to.plot %in% xlabels)
    }
    for(i in 1:length(xlabels)){
      text(x=median(new.dat[new.dat[,chrom.name] %in% scaffs.to.plot[xlab.indices[i]],"plot.pos"]),
           y=(y.lim[1]-((y.lim[2]-y.lim[1])/20)),
           labels=xlabels[i],adj=1,xpd=TRUE,cex=as.numeric(axis.size))
    }
  }
  return(new.dat)
}

#' Plot genome-wide statistics from a data frame with rectangles.
#' @note This fst.plot function is more efficient but this will work for most cases (might be buggy in some cases)
#' @param fst.dat The data.frame containing at least three columns: the statistic to be plotted, the chromosome ID, and the BP ID for each locus. Each row is a locus.
#' @param ci.dat A vector containing two values, upper and lower cutoff values (in that order).
#' Points above or below the cutoffs will be colored based on sig.col. Default is NULL, which turns this option off.
#' @param sig.col A vector containing two color values, the first for points above the upper cutoff value and the second for points below the lower cutoff value.
#' The defaults are red and yellow.
#' @param pt.col The color of the points. The default is grey7.
#' @param fst.name The name of the column containing the statistic to be plotted. Default is "Fst".
#' @param chrom.name The name of the column containing the chromosome information for each locus. Default is "Chrom".
#' @param bp.name The name of the column containing the basepair information for each locus. Default is "BP".
#' @param axis.size The value of cex.axis passed to the plotting of the y-axis. Default is 0.5. Set axis.size = 0 to suppress plotting of the y-axis.
#' @param scaffold.order A vector or list containing the names of the chromosomes/scaffolds in the order in which you would like them to be plotted.
#' @param groups A vector indicating which chromosomes to plot (generally to exclude those scaffolds not found in this particular set of statistics due to pruning/filters)
#' @param print.names A TRUE/FALSE value indicating whether chromosome/scaffold IDs should be printed beneath the x-axis.
#' @param y.lim The limits for the y-axis.
#' @param group.boundaries A data.frame with first column naming the chromosomes/linkage groups, second column having the last number of loci on the chromosome for indexing purposes.
#' e.g.:
#'  Chrom   End
#'  LG1     3945850
#'  LG2     435215
#' This parameter is used to make multiple plots from the same overall dataset have the same widths.
#' @return xes The fst.dat data frame with new values in BP, scaled to be the sequence in which points are plotted based on their position on the chromosome and the order the chromosomes are plotted in.
#' @seealso Flanagan & Jones 2017
#' @export
fst.plot.rect<-function(fst.dat,ci.dat=NULL, sig.col=c("red","yellow"),pt.col="grey7",
                        fst.name="Fst", chrom.name="Chrom", bp.name="BP",axis.size=0.5,
                        scaffold.order=NULL,groups=NULL,print.names=FALSE,y.lim=NULL,
                        group.boundaries=NULL,pt.cex=0.5,...){
  
  if(!is.null(scaffold.order)){
    scaff.ord<-scaffold.order$component_id
    lgs<-scaffold.order$object
  } else{
    if(!is.null(group.boundaries)){
      scaff.ord<-unique(group.boundaries[,1])
      lgs<-scaff.ord
    }else{
      scaff.ord<-unique(fst.dat[,chrom.name])
      lgs<-scaff.ord
    }
  }
  if(!is.null(groups)){
    lgs<-groups
    scaff.ord<-groups
  }
  fst.dat[,fst.name]<-as.numeric(as.character(fst.dat[,fst.name]))
  these.scaffs<-scaff.ord[scaff.ord %in% unique(fst.dat[,chrom.name])]
  all.scaff<-split(fst.dat, factor(fst.dat[,chrom.name]))
  
  #keep only the scaffolds you're plotting
  if(!is.null(groups))
  {
    all.scaff<-all.scaff[names(all.scaff) %in% groups]
  }
  
  group.boundaries<-group.boundaries[scaff.ord,]
  last.max<-0
  rect.xs<-NULL
  addition.values<-0
  xlist<-NULL
  #set up the rectangles
  for(i in 1:length(scaff.ord)){
    if(is.null(group.boundaries)){
      new.max<-last.max+nrow(all.scaff[[as.character(scaff.ord[i])]])
    }else{
      new.max<-as.numeric(group.boundaries[(group.boundaries[1,]) %in% scaff.ord[i],2])+
        as.numeric(as.character(last.max))
    }
    rect.xs<-rbind(rect.xs,c(last.max, new.max))
    rownames(rect.xs)[i]<-scaff.ord[i]
    addition.values<-c(addition.values, new.max)
    last.max<-new.max
  }
  rownames(rect.xs)<-scaff.ord
  for(i in 1:length(all.scaff)){
    if(is.null(group.boundaries)){
      #then this is done by relative position only
      all.scaff[[i]]<-all.scaff[[i]][order(all.scaff[[i]][,bp.name]),]
      all.scaff[[i]][,bp.name]<-
        seq(rect.xs[rownames(rect.xs) %in% names(all.scaff)[i],1]+1,
            rect.xs[rownames(rect.xs) %in% names(all.scaff)[i],1]+nrow(all.scaff[[i]]),1)
    }else{
      #then it takes genome space into account
      
      all.scaff[[i]]<-
        all.scaff[[i]][order(all.scaff[[i]][,bp.name]),]
      all.scaff[[i]][,bp.name]<-
        as.numeric(rect.xs[rownames(rect.xs) %in% names(all.scaff)[i],1])+
        as.numeric(as.character(all.scaff[[i]][,bp.name]))
      
    }
  }
  #change BP to plot
  x.max<-max(rect.xs,na.rm=T)
  x.min<-min(rect.xs,na.rm=T)
  if(is.null(y.lim)){
    y.max<-max(fst.dat[,fst.name])+0.1*max(fst.dat[,fst.name])
    y.min<-min(fst.dat[,fst.name])-0.1*min(fst.dat[,fst.name])
    if(min(fst.dat[,fst.name]) < 0) {
      y.min<-min(fst.dat[,fst.name]) - 0.1*min(fst.dat[,fst.name])
    } else {
      y.min<-0
    }
    
    y.lim<-c(y.min,y.max)
  }
  displacement<-y.lim[1]-((y.lim[2]-y.lim[1])/30)
  plot(c(x.min,x.max),y.lim,xlim=c(x.min,x.max),
       ylim=y.lim, col=pt.col,...,
       bty="n",type="n",	axes=F, xlab="", ylab="")
  for(i in 1:nrow(rect.xs)){
    if(i%%2 == 0) {
      rect.color<-"white"
    } else {
      rect.color<-"gray75"
    }
    rect(rect.xs[i,1],y.lim[1],rect.xs[i,2],y.lim[2],
         col=rect.color, border=NA)
    if(print.names==T){
      text(x=mean(rect.xs[i,]),
           y=displacement,labels=rownames(rect.xs)[i],
           adj=1,xpd=T,srt=45)
    }
  }
  for(i in 1:length(all.scaff)){
    points(all.scaff[[i]][,bp.name],
           all.scaff[[i]][,fst.name],
           pch=19, cex=pt.cex,col=pt.col,
           xlim=c(x.min,x.max),ylim=y.lim)
    if(!is.null(ci.dat)){
      temp.sig<-all.scaff[[i]][all.scaff[[i]][,fst.name] >= ci.dat[1],]
      points(temp.sig[,bp.name], temp.sig[,fst.name],
             col=sig.col[1], pch=19, cex=0.5)
      temp.sig<-all.scaff[[i]][all.scaff[[i]][,fst.name] <= ci.dat[2],]
      points(temp.sig[,bp.name], temp.sig[,fst.name],
             col=sig.col[2], pch=19, cex=0.5)
    }
  }
  if(axis.size>0){
    axis(2, at = seq(round(y.lim[1],2),round(y.lim[2],2),
                     round((y.lim[2]-y.lim[1])/2, digits=2)),
         ylim =y.lim, pos=0,
         labels=seq(round(y.lim[1],2),round(y.lim[2],2),
                    round((y.lim[2]-y.lim[1])/2, digits=2)),
         las=1,tck = -0.01, xlab="", ylab="", cex.axis=axis.size)
  }
  xes<-do.call("rbind",all.scaff)
  return(xes)
}


#' This plots the output of structure
#' @param structure.out A data.frame containing the structure output.
#' @param k The value of K used in Structure
#' @param pop.order The order in which the populations are plotted.
#' @param filename The name of the jpeg for the output file.
#' @param make.file A boolean (TRUE/FALSE) indicating whether the jpeg should be created
#' @param plot.new A boolean (TRUE/FALSE) indicating whether this is being added to an existing plot
#' @param colors A list of colors for the structure colors (if not provided, defaults to rainbow palette)
#' @param xlabel A boolean (TRUE/FALSE) indicating whether x labels should be plotted (default is TRUE)
#' @param ylabel An optional label for the y-axis.
#' @param xlabcol An optional vector of colors for x-axis labels
#' @param lab.cex Size of the labels
#' @export

plotting.structure<-function(structure.out, k, pop.order,
                             filename=paste("str.k",k,".jpeg",sep=""),make.file=TRUE,lab.cex=1,
                             plot.new=TRUE,colors=NULL,xlabel=TRUE,ylabel=NULL,xlabcol=NULL){
  str.split<-split(structure.out,structure.out[,1])
  if(is.null(colors)){
    bar.colors<-rainbow(k,s=0.5)
  } else {
    bar.colors<-colors
  }
  if(is.null(xlabcol)){
    xlabcol<-rep("black",length(str.split))
  }
  if(make.file==TRUE){
    jpeg(filename,width=7, height=1.25, units="in", res=300)
    par(mfrow=c(1,length(str.split)),mar=c(1,0,0,0), oma=c(1,0,0,0),cex=0.5)
  } else {
    if(plot.new==TRUE){
      par(mfrow=c(1,length(str.split)),mar=c(1,0,0,0), oma=c(1,0,0,0),cex=0.5)
    }
  }
  for(i in 1:length(str.split)){
    pop.index<-pop.order[i]
    barplot(height=as.matrix(t(str.split[[pop.index]][,-1])),
            beside=FALSE, space=0,	border=NA, col=bar.colors,
            xlab="", ylab="", xaxt='n', yaxt='n')#, new=plot.new)
    if(xlabel==TRUE){
      mtext(pop.index, 1, line=0.5, cex=lab.cex, outer=F,col=xlabcol[i])}
    if(!is.null(ylabel)){
      if(i == 1) { mtext(ylabel,2,cex=lab.cex) }
    }
  }
  if(make.file==TRUE) {dev.off()}
}

#' Plot covariance matrix from treemix
#' @param stem The filename basic stem for that run
#' @param poporder The list of populations in the order to plot them.
#' @return cp The covariance matrix in the plotting order.
treemix.cov.plot<-function(stem,poporder,colors=c("blue","yellow","red"),...){
  if("package:lattice" %in% search() == FALSE || "package:grid" %in% search()==FALSE){
    stop("ERROR: packages lattice and grid must be loaded")
  }
  cov<-read.table(gzfile(paste(stem, ".cov.gz", sep = "")), as.is = T, head = T, quote = "", comment.char = "")
  #reorder
  covplot = data.frame(matrix(nrow = nrow(cov), ncol = ncol(cov)))
  for(i in 1:length(poporder)){
    for( j in 1:length(poporder)){
      
      covplot[i, j] = cov[which(names(cov)==poporder[i]), which(names(cov)==poporder[j])]
      rownames(covplot)[i]<-poporder[i]
      colnames(covplot)[j]<-poporder[j]
    }
  }
  cp<-as.matrix(covplot)
  cp[lower.tri(cp)]<-NA
  cp[upper.tri(cp)]<-covplot[upper.tri(covplot)]
  #set up the colors
  pal<-colorRampPalette(colors)
  ncol=80
  cols<-pal(ncol)
  cpl<-levelplot(cp,col.regions=cols,alpha.regions=0.7,
                 scales = list(x=list(rot=90),tck = 0),xlab="",ylab="")
  print(cpl,...)
  trellis.focus("legend", side="right", clipp.off=TRUE, highlight=FALSE)
  grid.text("covariance", 0.2, 0, hjust=0.5, vjust=1.2)
  trellis.unfocus()
  return(cp)
}

#' Add a legend to the margins of a multiplot figure
#' @param ... legend parameters
#' @return Nothing is returned
#' @notes Modified from https://stackoverflow.com/questions/3932038/plot-a-legend-outside-of-the-plotting-area-in-base-graphics
#' @example 
#' plot(rnorm(50), rnorm(50), col=c("steelblue", "indianred"), pch=c(15,19))
#' outer.legend("top",legend=c("A","B"),pch=c(15,19),col=c("steelblue", "indianred"),ncol=2)
#' @export
outer.legend <- function(...) {
  opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
              mar=c(0, 0, 0, 0), new=TRUE)
  on.exit(par(opar))
  plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
  legend(...)
}

#' Create violin plots with 
#' @param ... legend parameters
#' @return Optionally returns the summary statistics of the data.
#' @notes Modified from vioplot() in library(vioplot)
#' @example 
#' mu<-2
#' si<-0.6
#' bimodal<-c(rnorm(1000,-mu,si),rnorm(1000,mu,si))
#' uniform<-runif(2000,-4,4)
#' normal<-rnorm(2000,0,3)
#' gwsca.vioplot(bimodal,uniform,normal,col=c("red","blue","green"))
#' @include sm
#' @export
gwsca.vioplot <- function(x,...,range=1.5,h=NULL,ylim=NULL,names=NULL, horizontal=FALSE,
                        col="magenta", border="black", lty=1, lwd=1, rectCol="black", colMed="white", pchMed=19, at, add=FALSE, wex=1,
                        drawRect=TRUE,plot.axes=TRUE,axis.box=FALSE,plot.ann=TRUE)
{
  # process multiple datas
  datas <- list(x,...)
  n <- length(datas)
  
  if(missing(at)) at <- 1:n
  
  # pass 1
  #
  # - calculate base range
  # - estimate density
  #
  
  # setup parameters for density estimation
  upper  <- vector(mode="numeric",length=n)
  lower  <- vector(mode="numeric",length=n)
  q1     <- vector(mode="numeric",length=n)
  q3     <- vector(mode="numeric",length=n)
  med    <- vector(mode="numeric",length=n)
  base   <- vector(mode="list",length=n)
  height <- vector(mode="list",length=n)
  baserange <- c(Inf,-Inf)
  
  # global args for sm.density function-call
  args <- list(display="none")
  
  if (!(is.null(h)))
    args <- c(args, h=h)
  
  for(i in 1:n) {
    data<-datas[[i]]
    
    # calculate plot parameters
    #   1- and 3-quantile, median, IQR, upper- and lower-adjacent
    data.min <- min(data)
    data.max <- max(data)
    q1[i]<-quantile(data,0.25)
    q3[i]<-quantile(data,0.75)
    med[i]<-median(data)
    iqd <- q3[i]-q1[i]
    upper[i] <- min( q3[i] + range*iqd, data.max )
    lower[i] <- max( q1[i] - range*iqd, data.min )
    
    #   strategy:
    #       xmin = min(lower, data.min))
    #       ymax = max(upper, data.max))
    #
    
    est.xlim <- c( min(lower[i], data.min), max(upper[i], data.max) )
    
    # estimate density curve
    smout <- do.call("sm.density", c( list(data, xlim=est.xlim), args ) )
    
    # calculate stretch factor
    #
    #  the plots density heights is defined in range 0.0 ... 0.5
    #  we scale maximum estimated point to 0.4 per data
    #
    hscale <- 0.4/max(smout$estimate) * wex
    
    # add density curve x,y pair to lists
    base[[i]]   <- smout$eval.points
    height[[i]] <- smout$estimate * hscale
    
    # calculate min,max base ranges
    t <- range(base[[i]])
    baserange[1] <- min(baserange[1],t[1])
    baserange[2] <- max(baserange[2],t[2])
    
  }
  
  # pass 2
  #
  # - plot graphics
  
  # setup parameters for plot
  if(!add){
    xlim <- if(n==1)
      at + c(-.5, .5)
    else
      range(at) + min(diff(at))/2 * c(-1,1)
    
    if (is.null(ylim)) {
      ylim <- baserange
    }
  }
  if (is.null(names)) {
    label <- 1:n
  } else {
    label <- names
  }
  # setup colors and borders
  if(length(col) < n){
    col<-c(col,rep(col,n-length(col)))
  }
  if(length(border)<n){
    border<-c(border,rep(border,n-length(border)))
  }
  if(length(colMed)<n){
    colMed<-c(colMed,rep(colMed,n-length(colMed)))
  }
  boxwidth <- 0.05 * wex
  
  # setup plot
  if(!add)
    plot.new()
  if(!horizontal) {
    if(!add){
      plot.window(xlim = xlim, ylim = ylim)
      if(plot.axes){
        if(plot.ann){
          axis(2)
          axis(1,at = at, label=label )
        }else{
          axis(2,labels=F)
          axis(1,at = at, labels=F )
        }
      }
    }
    
    if(axis.box){ box() }
    for(i in 1:n) {
      # plot left/right density curve
      polygon( c(at[i]-height[[i]], rev(at[i]+height[[i]])),
               c(base[[i]], rev(base[[i]])),
               col = col[i], border=border[i], lty=lty, lwd=lwd)
      
      if(drawRect){
        # plot IQR
        lines( at[c( i, i)], c(lower[i], upper[i]) ,lwd=lwd, lty=lty)
        
        # plot 50% KI box
        rect( at[i]-boxwidth/2, q1[i], at[i]+boxwidth/2, q3[i], col=rectCol)
        
        # plot median point
        points( at[i], med[i], pch=pchMed, col=colMed )
      }
    }
    
  }
  else {
    if(!add){
      plot.window(xlim = ylim, ylim = xlim,bty=axis.bty)
      if(plot.axes){
        if(plot.ann){
          axis(2)
          axis(1,at = at, label=label )
        }else{
          axis(2,labels=F)
          axis(1,at = at, labels=F )
        }
      }
    }
    
    if(axis.box){ box() }
    for(i in 1:n) {
      # plot left/right density curve
      polygon( c(base[[i]], rev(base[[i]])),
               c(at[i]-height[[i]], rev(at[i]+height[[i]])),
               col = col[i], border=border[i], lty=lty, lwd=lwd)
      
      if(drawRect){
        # plot IQR
        lines( c(lower[i], upper[i]), at[c(i,i)] ,lwd=lwd, lty=lty)
        
        # plot 50% KI box
        rect( q1[i], at[i]-boxwidth/2, q3[i], at[i]+boxwidth/2,  col=rectCol)
        
        # plot median point
        points( med[i], at[i], pch=pchMed, col=colMed )
      }
    }
  }
  invisible (list( upper=upper, lower=lower, median=med, q1=q1, q3=q3))
}