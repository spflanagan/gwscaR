#' Standardize a variable by its mean
#' @note Useful for generate Bayenv2 environmental variables
#' @param x A vector of numeric values
#' @return newx A vector of the input values after standardizing to the mean.
#' @export
std.by.mean<-function(x){
  m<-mean(x)
  s<-stats::sd(x)
  newx<-(x-m)/s
  return(newx)
}

#' Calculate the standard error of the mean
#' @param x A list of values
#' @return sem A value of the standard error of the mean
#' @export
sem<-function(x){
  sem<-stats::sd(x)/sqrt(length(x))
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

#' Convert a vcf df to a genepop file
#' @param vcf A data.frame in vcf format
#' @param pop.list A list of population names (the individuals must have names containing the population names)
#' @param pop.map A data frame with col 1 = ID names and col 2 = pop names
#' @param gpop.name A name for the output genepop file
#' @return gpop A dataframe in genepop format
#' @export
#' @examples
#' vcf.file<-system.file("extdata", "example.vcf.txt",package = "gwscaR")
#' vcf<-parse.vcf(vcf.file)
#' gpop<-vcf2gpop(vcf,
#'     pop.list=c("FEM","PRM","OFF"),
#'     gpop.name="out.genepop")
vcf2gpop<-function(vcf,pop.list=NULL,pop.map=NULL,gpop.name){#without the SNP column
  locusids<-paste(vcf$`#CHROM`,as.character(vcf$POS),sep="_")
  locusids<-gsub("\\.","_",locusids)
  indids<-colnames(vcf)[10:ncol(vcf)]
  gpop.mat<-extract.gt.vcf(vcf[,colnames(vcf)!="SNP"])
  gpop<-t(gpop.mat[,10:ncol(gpop.mat)])
  gpop[gpop=="0/0" | gpop=="0|0"]<-"0101"
  gpop[gpop=="0/1" | gpop=="0|1"]<-"0102"
  gpop[gpop=="1/0" | gpop=="1|0"]<-"0201"
  gpop[gpop=="1/1" | gpop=="1|1"]<-"0202"
  gpop[gpop=="./." | gpop==".|."]<-"0000"
  #write to file
  utils::write.table(locusids,gpop.name,sep='\n',quote=FALSE,
              col.names = paste("Title line: ",gpop.name,sep=""),row.names=FALSE)
  if(!is.null(pop.list)){
    for(i in 1:length(pop.list)){
      pop<-gpop[grep(pop.list[i],rownames(gpop)),,drop=FALSE]
      utils::write.table(paste("POP",pop.list[i],sep=" "),gpop.name,quote=FALSE,col.names = FALSE,row.names=FALSE,append=TRUE)
      rownames(pop)<-paste(rownames(pop),",",sep="")
      utils::write.table(pop,gpop.name,quote=FALSE,col.names=FALSE,row.names=TRUE,sep=" ",append=TRUE)
    }
  } else if(!is.null(pop.map)){
    for(i in 1:length(unique(pop.map[,2]))){

      pname<-unique(pop.map[,2])[i]
      pop<-gpop[rownames(gpop) %in% pop.map[pop.map[,2]%in%pname,1],,drop=FALSE]
      if(nrow(pop)>0){
        utils::write.table(paste("POP",pname,sep=" "),gpop.name,quote=FALSE,col.names = FALSE,row.names=FALSE,append=TRUE)
        rownames(pop)<-paste(rownames(pop),",",sep="")
        utils::write.table(pop,gpop.name,quote=FALSE,col.names=FALSE,row.names=TRUE,sep=" ",append=TRUE)
      }
    }
  } else {
    stop("Must provide a pop.list or pop.map")
  }
  colnames(gpop)<-locusids
  return(gpop)
}


#' Convert a vcf df to a coancestry input file
#' @param vcf A data.frame in vcf format
#' @param out.name A name for a an output file (default is coancestry_afs.txt)
#' @return co.afs A dataframe with all of the allele frequencies
#' @export
#' @examples
#' vcf.file<-system.file("extdata", "example.vcf.txt",package = "gwscaR")
#' vcf<-parse.vcf(vcf.file)
#' co.afs<-vcf2coanAF(vcf,"coancestry_afs.txt")
vcf2coanAF<-function(vcf,out.name="coancestry_afs.txt"){
  co.afs<-do.call(rbind,apply(vcf,1,function(vcf.row,out.name){
    af<-calc.afs.vcf(vcf.row)
    bases<-c("A","C","G","T")
    rname<-which(bases==vcf.row[["REF"]])
    aname<-which(bases==vcf.row[["ALT"]])
    co.af<-cbind(round(af$RefFreq,4),round(af$AltFreq,4))
    colnames(co.af)<-c(rname,aname)
    suppressWarnings(utils::write.table(co.af,out.name,sep='\t',append = TRUE,quote=FALSE,row.names = FALSE,col.names = TRUE))
    row.names(co.af)<-vcf.row["ID"]
    colnames(co.af)<-NULL
    return(data.frame(co.af))
  },out.name=out.name))
  return(co.afs)
}

#' Convert a vcf df to a coancestry input file
#' @param vcf A data.frame in vcf format
#' @param out.name A name for a an output file (default is coancestry_gty.txt)
#' @return co.gt A dataframe with all of the genotypes
#' @export
#' @examples
#' vcf.file<-system.file("extdata", "example.vcf.txt",package = "gwscaR")
#' vcf<-parse.vcf(vcf.file)
#' co.gts<-vcf2coanGT(vcf,"coancestry_gty.txt")
vcf2coanGT<-function(vcf,out.name="coancestry_gty.txt"){
  gts<-extract.gt.vcf(vcf)
  co.gt<-do.call(cbind,apply(gts,1,function(gt){
    bases<-c("A","C","G","T")
    rname<-which(bases==gt[["REF"]])
    aname<-which(bases==gt[["ALT"]])
    g<-do.call(rbind,lapply(gt[10:length(gt)],function(x) {
      strsplit(as.character(x),"/")[[1]]}))
    #print(rname)
    g[g=="0"]<-rname
    g[g=="1"]<-aname
    g[g=="."]<-0
    return(as.data.frame(g))
  }))

  utils::write.table(co.gt,out.name,row.names=TRUE,col.names=FALSE,quote=FALSE,sep='\t')
  return(co.gt)
}
