#Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
#Purpose: Calculate useful population genetics statistics 
#Package: gwscaR


#' sliding window pi and rho - across all snps
#' @note pi = 1-sum((ni choose 2)/(n choose i)); ni is number of alleles i in sample, n = sum(ni)
#' @note Jones et al. (2012) used 2500bp sliding windows with a step size 500bp<-more than just SNPs, but I'll just focus on SNPs
#' @note Hohenlohe did a similar thing and weighted pi by all nt sites (not just SNPs) but rho by SNPs only
#' @param vcf.row A row of a vcf file; use this in conjunction with apply
#' @return pi The nucleotide diversity at that site
#' @example all.pi<-apply(vcf,1,pi)
calc.pi<-function(vcf.row){
  alleles<-vcf.alleles(vcf.row)
  af.num<-table(alleles)
  n<-sum(af.num)
  pi<-1-sum(choose(af.num,2))/choose(n,2)
  return(pi)
}

#' Calculate rho (private alleles)
#' @note rho=1 if allele in pop j is only found in that pop and at least one ind was genotyped at that site in each pop; rho = 0 otherwise
calc.rho<-function(vcf.row,pop.list){
  pop.alleles<-lapply(pop.list,function(pop){
    pop.vcf<-cbind(vcf.row[1:9],vcf.row[grep(pop,colnames(vcf.row))])
    alleles<-vcf.alleles(pop.vcf)
    af.num<-table(alleles)
    return(names(af.num))
  })
  rho<-0
  pop.counts<-lapply(pop.alleles,length)
  if(length(pop.counts[pop.counts==0])>0){
    unique.counts<-table(unlist(pop.alleles[pop.counts==1]))
    uni<-unique.counts[unique.counts==1]
    if(length(uni)>1){
      uni.matches<-grep(names(uni),pop.alleles)
      if(length(uni.matches)==1){ #it should only match itself
        rho<-1}
    }}
  return(rho)
}

#' Calculate a sliding average
sliding.avg<-function(dat,win.start,width){
  if((win.start+width)>nrow(dat)){
    win.end<-nrow(dat)
  }else {
    win.end<-(win.start+width)
  }
  win.dat<-dat[win.start:win.end,]
  avg.dat<-data.frame(Avg.Pos=mean(dat[win.start:win.end,1]),
                      Avg.Pi=mean(win.dat[,2]))
  return(avg.dat)
}

#' Calculate a sliding window
sliding.window<-function(vcf,chr,stat="pi",width=250,pop.list=NULL){
  avg.dat<-lapply(chr,function(chr){
    chr.vcf<-vcf[vcf[,1] %in% chr,]
    if(stat=="pi"){
      dat<-data.frame(Pos=chr.vcf$POS,Pi=unlist(apply(chr.vcf,1,calc.pi))) }
    if(stat=="rho"){
      dat<-data.frame(Pos=chr.vcf$POS,Rho=unlist(apply(chr.vcf,1,calc.rho,pop.list=pop.list))) }
    steps<-seq(1,nrow(chr.vcf),50)
    if(stat=="pi"){
      avg.stat<-do.call("rbind",lapply(steps,sliding.avg,dat=dat,width=width)) }
    if(stat=="rho"){
      avg.stat<-do.call("rbind",lapply(steps,sliding.avg,dat=dat,width=width)) }
    avg.stat$Chr<-rep(chr,nrow(avg.stat))
    head(avg.stat)
    return(avg.stat)
  })
  return(avg.dat)
}

#' Generate a distance matrix for a SNP using fsts
#' @param vcf.row A row of a vcf file
#' @param pop.list A list of populations in the order you want them to appear in the matrix
#' @return fst.tree A neighbor-joining tree from the ape package
#' @example fst.trees<-list()
#' for(vcf.row in 1: nrow(vcf)){
#'   fst.trees<-c(fst.trees,get.dist(vcf[vcf.row,],pop.list)) #getting an error
#' }
#' @export
get.dist<-function(vcf.row,pop.list){
  if("package:ape" %in% search() == FALSE){
    stop("ERROR: You must load the package ape for get.dist() to run.")
  }
  fst.matrix<-matrix(nrow=length(pop.list),ncol=length(pop.list))
  for(i in 1:(length(pop.list)-1)){
    for(j in (i+1):length(pop.list)){
      pop1<-colnames(vcf.row)[grep(pop.list[i],colnames(vcf.row))]
      pop2<-colnames(vcf.row)[grep(pop.list[j],colnames(vcf.row))]
      fst<-fst.one.vcf(vcf.row, pop1,pop2,maf=0,cov.thresh = 0)
      fst.matrix[i,j]<-fst$Fst
    }
  }
  colnames(fst.matrix)<-pop.list
  rownames(fst.matrix)<-pop.list
  fst.nj<-ape::njs(as.dist(t(fst.matrix)))
  fst.tree<-write.tree(fst.nj,digits=0)
  return(fst.tree)
}

#' Generate a treemix file from vcf
#' 
treemix.from.vcf<-function(vcf,pop.list){
  tm.df<-matrix(nrow=nrow(vcf),ncol=length(pop.list))
  for(i in 1:nrow(vcf)){
    vcf.row<-vcf[i,]
    #tm.df<-as.matrix(t(apply(vcf,1,function(vcf.row){ #freaking apply was giving me weird results
    all.alleles<-names(table(vcf.alleles(vcf.row))) 
    if(length(all.alleles)==2){
      this.loc<-do.call("cbind",  
                        lapply(pop.list,function(pop){
                          pop.vcf.row<-cbind(vcf.row[1:9],vcf.row[grep(pop,colnames(vcf.row))])
                          pop.alleles<-vcf.alleles(pop.vcf.row)
                          pop.counts<-table(pop.alleles)
                          tm.pop<-"0,0"
                          if(length(pop.counts)==1){
                            allele<-names(pop.counts)
                            if(grep(allele,all.alleles)==1){ #maintain order
                              tm.pop<-paste(pop.counts[[1]],0,sep=",") 
                            }else{
                              tm.pop<-paste(0,pop.counts[[1]],sep=",") }
                          }
                          if(length(pop.counts)==2){
                            tm.pop<-paste(pop.counts[[all.alleles[1]]],pop.counts[[all.alleles[2]]],sep=",") 
                          }
                          return(tm.pop)
                        }))
    }
    #return(this.loc)
    tm.df[i,]<-as.vector(this.loc)
  }#)))
  colnames(tm.df)<-pop.list
  return(tm.df)
}

#' Calculate isolation by distance per locus
#' @param ped.file A data.frame in ped format
#' @param dist.mat A distance matrix containing the geeographic distances between populations
#' @param pop.order A list with the order for the populations
#' @return results.mantel A data.frame with two columns containing the mantel test results
#' @export
fst.ibd.byloc<-function(ped.file,dist.mat,pop.order){
  results.mantel<-data.frame()
  for(i in seq(7,ncol(ped.file),2)){
    res<-ade4::mantel.rtest(
      as.dist(t(pairwise.fst(ped.file,i,i+1,pop.order))),
      as.dist(t(dist.mat)), nrepet=9999)
    results.mantel<-rbind(results.mantel,cbind(res$obs,res$pvalue))
  }
  results.mantel<-as.data.frame(results.mantel)
  colnames(results.mantel)<-c("Obs","P")
  return(results.mantel)
}

#' Calculate pairwise Pst between population pairs
#' @param dat A dataframe with the trait values, first column must be the pop ID
#' @param pop.order A list of the order of the populations
#' @return dat.var A data.frame with the pairwise Pst values
#' @export
pairwise.pst<-function(dat, pop.order){
  #first column must be pop id/grouping factor
  dat.split<-split(dat, factor(dat[,1]))
  dat.var<-as.data.frame(setNames(
    replicate(length(pop.order),numeric(0), simplify = F), pop.order))
  for(i in 1:(length(pop.order)-1)){
    for(j in (i+1):length(pop.order)){
      temp.data<-rbind(as.data.frame(dat.split[[pop.order[i]]]),
                       as.data.frame(dat.split[[pop.order[j]]]))
      colnames(temp.data)<-c("PopID","Var")
      temp.data$PopID<-factor(temp.data$PopID)
      anv <- nlme::lme(fixed=Var ~ 1, random=~1|PopID,data=temp.data)
      varcomp <- nlme::VarCorr(anv)
      v.btwn<- as.numeric(varcomp[1])
      v.wthn <- as.numeric(varcomp[2])
      pst <- v.btwn/(v.btwn+2*v.wthn)
      dat.var[pop.order[i],pop.order[j]]<-pst
      #aov.var<-summary.aov(
      #	aov(temp.data[,2]~temp.data[,1]))[[1]]$`Sum Sq`
      #aov.df<-summary.aov(
      #	aov(temp.data[,2]~temp.data[,1]))[[1]]$`Df`
      #dat.var[pop.order[i],pop.order[j]]<-aov.var[2]/(aov.var[2]+
      #	(2*(aov.var[1]/(aov.df[2]-1))))
    }
  }
  dat.var<-rbind(dat.var,rep(NA, ncol(dat.var)))
  rownames(dat.var)<-colnames(dat.var)
  return(dat.var)
}

#' Calulate pairwise Psts and test for IBD
#' @param trait.df A data frame with trait values for all traits and all individuals
#' @param comp.df A dataframe with the distance values
#' @param id.index The index value for the trait being calculated
#' @param pop.order A vector of population names in the order you want them.
#' @return results.mantel A data.frame with the results of the mantel test
#' @export
pst.mantel<-function(trait.df,comp.df,id.index,pop.order){
  results.mantel<-data.frame()
  for(i in 3:ncol(trait.df)){
    res<-ade4::mantel.rtest(
      as.dist(t(pairwise.pst(trait.df[,c(id.index,i)],pop.order))),
      as.dist(t(comp.df)), nrepet=9999)
    results.mantel<-rbind(results.mantel,cbind(res$obs,res$pvalue))
  }
  results.mantel<-as.data.frame(results.mantel)
  rownames(results.mantel)<-colnames(trait.df)[3:ncol(trait.df)]
  colnames(results.mantel)<-c("Obs","P")
  return(results.mantel)
}


#' Compare Pst and Fst locus by locus from a ped file
#' @param ped.file A data.frame with data in ped format
#' @param trait.df A data.frame containing all of the traits data
#' @param pop.order An order in which the populations should be analyzed
#' @param trait.ind The trait index in the trait.df
#' @return results.list A dataframe containing the restuls of the mantel test (one column for observations and one for p-values)
#' @export
fst.pst.byloc<-function(ped.file,trait.df,pop.order,trait.ind){
  results.list<-list()
  for(j in 3:ncol(trait.df)){
    results.mantel<-data.frame()
    for(i in seq(7,ncol(ped.file),2)){
      res<-ade4::mantel.rtest(
        as.dist(t(pairwise.fst(ped.file,i,i+1,pop.order))),
        as.dist(t(pairwise.pst(trait.df[,c(trait.ind,j)],pop.order))),
        nrepet=9999)
      results.mantel<-rbind(results.mantel,cbind(res$obs,res$pvalue))
    }
    results.mantel<-as.data.frame(results.mantel)
    colnames(results.mantel)<-c("Obs","P")
    results.list<-append(results.list,data.frame(results.mantel))
  }
  #names(results.list)<-colnames(trait.df)[3:ncol(trait.df)]
  return(results.list)
}

