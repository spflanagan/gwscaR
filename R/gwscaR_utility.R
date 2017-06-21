#' Standardize a variable by its mean
#' @note Useful for generate Bayenv2 environmental variables
#' @param x A vector of numeric values
#' @return newx A vector of the input values after standardizing to the mean.
#' @export
std.by.mean<-function(x){
  m<-mean(x)
  s<-sd(x)
  newx<-(x-m)/s
  return(newx)
}

#' Calculate the standard error of the mean
#' @param x A list of values
#' @return sem A value of the standard error of the mean
#' @export
sem<-function(x){
  sem<-sd(x)/sqrt(length(x))
  return(sem)
}


#' A function to reorder a data.frame
#' @param dat A data.frame. Group ID should be in column 1
#' @param order.list A list of factors for reordering
#' @return dat.new The reordered dataframe
#' @export
changeorder.df<-function(dat,order.list){
  #dat has to have the grouping IDs in row 1
  #those grouping ids must match the factors in order.list
  dat.sep<-split(dat, dat[,1])
  dat.new<-dat.sep[[order.list[1]]]
  for(i in 2:length(order.list)){
    dat.new<-rbind(dat.new, dat.sep[[order.list[i]]])
  }
  return(dat.new)
}

#' Choose one SNP per RAD locus from a vcf
#' @param vcf A data.frame in vcf format
#' @return new.vcf A data.frame with only one SNP per RAD locus
#' @export
choose.one.snp<-function(vcf){
  keep.col<-colnames(vcf)
  vcf$id.pos<-paste(vcf$ID,vcf$POS,sep=".")
  sub.vcf<-tapply(vcf$id.pos,vcf$ID, sample,size=1)
  new.vcf<-vcf[vcf$id.pos %in% sub.vcf,keep.col]
  return(new.vcf)
}
