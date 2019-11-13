

#' Convert a vcf file to a SNPs input file for dadi
#' @param vcf The vcf dataframe (from parse.vcf) or the name of the vcf file
#' @param filename The name fo the file you want the data written to. If not provided, no file is written
#' @param pop.list Optional list of population names. If this is provided, individuals should have the population names in their IDs. If neither pop.map or pop.list are provided, all individuals are placed in one population.
#' @param pop.map Optional map of individual IDs in column 1 and their population names in column 2 (same as a Stacks population map). Can either be a file or an R object. If neither pop.map or pop.list are provided, all individuals are placed in one population.
#' @param whitelist Optional list of SNPs to include in the output file. If not provided, all SNPs in the vcf will be output.
#' @param outgroup Optional 2-column dataframe containing ingroup (column 1) and outgroup (column 2) sequences. If not provided, both will be "---".
#' @return a dataframe containing the contents of the vcf file, including headers.
#' @examples
#' vcf.file<-system.file("extdata", "example.vcf.txt",package = "gwscaR")
#' dadi<-vcf2dadiSNPs(vcf.file)
#' @seealso Gutenkunst et al. 2013
#' @export
vcf2dadiSNPs<-function(vcf, filename=NULL,pop.list=NA,pop.map=NULL,whitelist=NULL,outgroup=NA){
  #check to see if vcf needs to be read in from file.
  if(!exists(vcf)){
    print("Parsing vcf file")
    vcf<-parse.vcf(vcf)
  }
  #check to see if there is a pop.list or pop.map
  if(is.na(pop.list)){
    if(is.null(pop.map)){
      npops<-1
      pop.list<-"Pop1"
      pop.map<-data.frame(Individuals=colnames(vcf[,10:ncol(vcf)]),
                          Pop=rep("Pop1",length(colnames(vcf[,10:ncol(vcf)]))))
    }else{
      if(!exists(pop.map)){
        pop.map<-read.delim(pop.map)
      }
      npops<-length(levels(as.factor(pop.map[,2])))
      pop.list<-levels(as.factor(pop.map[,2]))
    }
  }else{
    npops<-length(pop.list)
    pop.map<-data.frame()
    for(i in 1:npops){
      inds<-grep(pop.list[i],colnames(vcf),value=TRUE)
      pop.map<-rbind(Individuals=pop.map,
                     Pop=cbind(inds, rep(pop.list[i],length(inds))))
    }

  }

  #generate other pieces of info for the output
  allele1<-vcf$REF
  allele2<-vcf$ALT

  if(!is.null(whitelist)){
    gts<-extract.gt.vcf(vcf[vcf$ID %in% whitelist,])
  }else{
    gts<-extract.gt.vcf(vcf)
  }

  if(is.na(outgroup)){
    ingroup<-rep("---",nrow(vcf))
    outgroup<-rep("---",nrow(vcf))
  } else{
    ingroup<-outgroup[,1]
    outgroup<-outgroup[,2]
  }

  #create the new dataframe
  dadi.snps<-data.frame(Ingroup=ingroup,Outgroup=outgroup,Allele1=allele1,Allele2=allele2,GeneID=gts$ID,Position=gts$POS)

  #count alleles per population
  for(i in 1:length(pop.list)){
    inds<-unlist(pop.map[pop.map[,2] %in% pop.list[i],1])
    pop.dat<-data.frame(pop_al1=rep(0,nrow(gts)),pop_al2=rep(0,nrow(gts)))
    for(ii in 1:nrow(gts)){
      alleles<-vcf.alleles(gts[ii,c(colnames(gts)[1:9],
                      colnames(gts)[colnames(gts)%in%inds])])
      al1<-length(alleles[alleles%in%gts[ii,"REF"]])
      al2<-length(alleles[alleles%in%gts[ii,"ALT"]])
      pop.dat[ii,]<-c(al1,al2)
    }
    dadi.snps<-cbind(dadi.snps,pop.dat)
    colnames(dadi.snps)[c(ncol(dadi.snps)-1,ncol(dadi.snps))]<-
      c(paste(pop.list[i],"_1",sep=""),paste(pop.list[i],"_2",sep = ""))
  }

  #reorder the columns
  dadi.snps<-dadi.snps[,c("Ingroup","Outgroup","Allele1",grep("_1",colnames(dadi.snps),value=TRUE),"Allele2",
                          grep("_2",colnames(dadi.snps),value=TRUE),"GeneID","Position")]
  colnames(dadi.snps)<-gsub("_1","",colnames(dadi.snps))
  colnames(dadi.snps)<-gsub("_2"," ",colnames(dadi.snps))
  #write to file
  if(!is.null(filename)){
    utils::write.table(dadi.snps,filename,col.names = TRUE,row.names = FALSE,sep='\t',quote=FALSE)
  }
  return(dadi.snps)
}
