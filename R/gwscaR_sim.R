


#' Subset a ped file and run pcadapt on it
#' @param ped The ped data.frame or file
#' @param num_inds The number of individuals per population to simulate
#' @param subname The name for the subset image file
#' @return A data.frame with ped data
#' @export
pcadapt_subset<-function(ped, num_inds,subname){
  #install and load pcadapt if you haven't already
  if("pcadapt" %in% rownames(utils::installed.packages())){
    do.call('library',list("pcadapt"))
  }else{
    utils::install.packages("pcadapt",dependencies = TRUE)
    do.call("library",list("pcadapt"))
  }
  if(ped %in% ls()){
    ped<-ped
  }else if(ped %in% list.files()){
    ped<-utils::read.delim(ped)
  }else{
    stop("ped does not exist as an object or a file")
  }
  if(is.null(subname)){
    subname<-paste("subset",num_inds,"inds1.ped",sep="")
  }
  if(num_inds==0){
    keepinds<-ped[,2]
  }else{
    keepinds<-unlist(tapply(ped[,2],ped[,1],sample,size=num_inds,replace=FALSE))
  }
  pedsub<-ped[ped[,2] %in% keepinds,]
  utils::write.table(pedsub,subname,col.names = FALSE,row.names = FALSE,quote=FALSE,sep=" ")
  filename<-pcadapt::read.pcadapt(subname,type="ped")
  x<-pcadapt::pcadapt(filename,K=10, min.maf=0.001)
  grDevices::png(gsub("ped","png",subname),height=7,width=7,units="in",res=300)
  graphics::plot(x,option="scores",pop=pedsub[,1])#K=6
  grDevices::dev.off()
  return(x)
}


#' Simulate genetic data for populations with different allele frequencies
#' @param npops The number of populations to simulate
#' @param inds The number of individuals per population to simulate
#' @param nloc The number of loci to simulate
#' @param p A vector of average minor allele frequencies per locus to simulate
#' @param outname The output name for the file
#' @param analyze A logical value for whether or not to run the analysis (defaults to TRUE)
#' @param ns A vector of numbers of individuals to simulate
#' @return A data.frame with ped data
#' @export
popgen.sim<-function(npops=8,inds=10,nloc=10000,p=rep(0.01,8),outname="simulated.ped",analyze=TRUE,ns=c(12,24,36,48,96)){

  #install and load pcadapt if you haven't already
  if("pcadapt" %in% rownames(utils::installed.packages())){
    do.call('library',list("pcadapt"))
  }else{
    utils::install.packages("pcadapt",dependencies = TRUE)
    do.call("library",list("pcadapt"))
  }
  #create ped first columns
  simped<-data.frame(cbind(rep(1:npops,each=inds),rep(paste("Ind",rep(1:inds),sep=""),npops),
                           rep(0,inds*npops),rep(0,inds*npops),rep(0,inds*npops),rep(-9,inds*npops)),
                     stringsAsFactors = FALSE)

  #simulate allele frequencies for each locus
  if(length(p)>1)
  {
    print("simulating different allele frequencies for each population")
    frqs<-lapply(p,function(f){
      x<-stats::rnorm(f,n=nloc,sd=f)
      z<-unlist(lapply(x,function(y){ #keeps it between 0 and 1
        while(y < 0 | y > 1){
          y<-stats::rnorm(f,n=1,sd=f)
        }
        return(y)
      }))
      return(z)
    })
  }else{
    x<-stats::rnorm(p,n=nloc,sd=p)
    frqs<-unlist(lapply(x,function(y){ #keeps it between 0 and 1
      while(y < 0 | y > 1){
        y<-stats::rnorm(p,n=1,sd=p)
      }
      return(y)
    }))

  }



  #for each locus, assign genotypes for each population
  for(n in seq(1,nloc*2,by=2)){
    #designate columns
    col1<-n+6
    col2<-n+7

    #simulate alleles
    if(length(p)>1){
      al1<-unlist(lapply(frqs,function(f){ stats::rbinom(f[(n+1)/2],n=inds,size=1) } ))
      al2<-unlist(lapply(frqs,function(f){ stats::rbinom(f[(n+1)/2],n=inds,size=1) } ))
    }else{
      al1<-al2<-NULL
      for(i in 1:npops){
        al1<-c(al1,stats::rbinom(frqs[(n+1)/2],n=inds,size=1))
        al2<-c(al2,stats::rbinom(frqs[(n+1)/2],n=inds,size=1))
      }

    }

    if(stats::rbinom(1,1,0.5)==0){ #randomly decide if it's going to be G/C or A/T locus
      if(stats::rbinom(1,1,0.5)==0){ #randomly decide if major allele is G
        al1[al1==0]<-"G"
        al1[al1==1]<-"C"
        al2[al2==0]<-"G"
        al2[al2==1]<-"C"
      }else{ #or if it's C
        al1[al1==0]<-"C"
        al1[al1==1]<-"G"
        al2[al2==0]<-"C"
        al2[al2==1]<-"G"
      }
    }else{
      if(stats::rbinom(1,1,0.5)==0){ #randomly decide if major allele is A
        al1[al1==0]<-"A"
        al1[al1==1]<-"T"
        al2[al2==0]<-"A"
        al2[al2==1]<-"T"
      }else{ #or if it's T
        al1[al1==0]<-"T"
        al1[al1==1]<-"A"
        al2[al2==0]<-"T"
        al2[al2==1]<-"A"
      }
    }
    #put this in the ped object
    simped[,col1]<-al1
    simped[,col2]<-al2
  }
  utils::write.table(simped,outname,col.names = FALSE,row.names = FALSE,quote=FALSE,sep=" ")

  if(isTRUE(analyze)){
    # now run pcadapt
    filename<-pcadapt::read.pcadapt(outname,type="ped")
    x<-pcadapt::pcadapt(filename,K=20,min.maf=0.001)
    grDevices::png(paste(outname,".png",sep=""),height = 7, width = 7, units="in",res=300)
    graphics::plot(x,option="scores",pop=simped[,1])
    grDevices::dev.off()
  }

  #subset it for further analysis
  pcasubs<-lapply(ns,function(n){
    sub<-pcadapt_subset(ped = simped, num_inds = n, subname=paste("sub",n,outname,sep=""))
  })
  return(simped)
}
