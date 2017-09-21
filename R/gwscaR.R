#Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
#Purpose: Run Fst-based selection components analysis

#' Read in a vcf file
#' @param filename The name of the vcf file
#' @return a dataframe containing the contents of the vcf file, including headers.
#' @examples
#' vcf<-parse.vcf(system.file("extdata/example.vcf",package = "gwscaR"))
#' @seealso Flanagan & Jones 2017
#' @export
parse.vcf<-function(filename){
  #if(substr(filename,nchar(filename)-3,nchar(filename)) != ".vcf") { filename<-paste(filename,"vcf",sep=".") }
  vcf<-read.delim(filename,comment.char="#",sep='\t',header=F,stringsAsFactors = F,strip.white = T)
  header.start<-grep("#CHROM",scan(filename,what="character"))
  header<-scan(filename,what="character")[header.start:(header.start+ncol(vcf)-1)]
  colnames(vcf)<-header
  return(vcf)
}

#' Calculate per-locus coverage from a vcf file
#' @param vcf A data.frame containing data in vcf format.
#' @param subset A list of the column names of the individuals to be used (optional)
#' @return cov.dat A data.frame with one row for each locus and 14 columns:
#'    Chrom: The chromosome
#'    Pos: The position/BP on the chromosome
#'    Locus: The Locus ID
#'    NumMissing: The number of individuals not gentoyped at this locus
#'    NumPresent: The number of individuals genotyped at this locus
#'    PropMissing: The proportion of individuals not genotyped at this locus
#'    AvgCovRef: The average coverage of the reference allele in genotyped individuals
#'    AvgCovAlt: The average coverage of the alternative allele in genotyped individuals
#'    AvgCovRatio: The average ration of Reference/Alternative allele coverage
#'    AvgCovTotal: The average of the number of reference + alternative reads per individual
#'    CovVariance: The variance in coverage among individuals
#'    NumHet: The number of individuals genotyped as heterozygotes
#'    PropHet: The proportion of individuals genotyped as heterozygotes
#'    TotalNumReads: The total number of reads at this locus
#' @export
vcf.cov.loc<-function(vcf,subset=NULL){
  if(is.null(subset)){
    subset<-colnames(vcf)[10:ncol(vcf)]
  }
  cov.dat<-do.call("rbind",apply(vcf,1,function(vcf.row){
    cov<-unlist(lapply(vcf.row[subset],function(x){
      c<-strsplit(as.character(x),split=":")[[1]][3]
      return(c)
    }))
    miss<-length(cov[cov==".,."])
    pres<-length(cov[cov!=".,."])
    ref<-sum(as.numeric(unlist(lapply(cov[cov!=".,."],
                                      function(x){
                                        strsplit(as.character(x),",")[[1]][1]
                                      }))))/pres
    alt<-sum(as.numeric(unlist(lapply(cov[cov!=".,."],
                                      function(x){
                                        strsplit(as.character(x),",")[[1]][2]
                                      }))))/pres
    tot<-sum(as.numeric(unlist(lapply(cov[cov!=".,."],
                                      function(x){
                                        as.numeric(strsplit(as.character(x),",")[[1]][1]) +
                                          as.numeric(strsplit(as.character(x),",")[[1]][2])
                                      }))))
    var.cov<-var(as.numeric(unlist(lapply(cov[cov!=".,."],
                                          function(x){
                                            as.numeric(strsplit(as.character(x),",")[[1]][1]) +
                                              as.numeric(strsplit(as.character(x),",")[[1]][2])
                                          }))))
    het<-unlist(lapply(vcf.row[subset],function(x){
      strsplit(as.character(x),split=":")[[1]][1]
    }))
    het<-length(het[het=="0/1" | het=="1/0"])
    return(data.frame(Chrom=vcf.row[1],Pos=vcf.row[2],Locus=vcf.row[3],
                      NumMissing=miss, NumPresent=pres,PropMissing=miss/(miss+pres),
                      AvgCovRef=ref,AvgCovAlt=alt, AvgCovRatio=ref/alt,AvgCovTotal=tot/pres, CovVariance=var.cov,
                      NumHet=het,PropHet=het/pres,TotalNumReads = tot,stringsAsFactors = F))
  }))
  return(data.frame(cov.dat))
}

#' Calculate per-individual coverage from a vcf file
#' @param vcf A data.frame containing data in vcf format.
#' @return icov A data.frame with one row for each individual and columns:
#'    NumMissing: The number of missing loci
#'    NumPresent: The number of loci genotyped in this individual
#'    AvgCovRef: The average coverage for reference alleles across all genotyped loci in this individual
#'    AvgCovAlt: The average coverage for alternative alleles across all genotyped loci in this individual
#'    AvgCovTot: The average total coverage (ref + alt coverage) across all genotyped loci in this individual
#'    PropHet: Proportion of loci at which this individual is heterozygous
#'    NumReads: The total number of reads for this individual
#' @export
vcf.cov.ind<-function(vcf){
  icov<-as.data.frame(do.call("rbind",apply(vcf,2,function(vcf.col){
    cov<-unlist(lapply(vcf.col,function(x){
      c<-strsplit(as.character(x),split=":")[[1]][3]
      return(c)
    }))
    miss<-length(cov[cov==".,."])
    pres<-length(cov[cov!=".,."])
    ref<-sum(as.numeric(unlist(lapply(cov[cov!=".,."],
                                      function(x){
                                        strsplit(as.character(x),",")[[1]][1]
                                      }))))/pres
    alt<-sum(as.numeric(unlist(lapply(cov[cov!=".,."],
                                      function(x){
                                        strsplit(as.character(x),",")[[1]][2]
                                      }))))/pres
    tot<-sum(as.numeric(unlist(lapply(cov[cov!=".,."],
                                      function(x){
                                        as.numeric(strsplit(as.character(x),",")[[1]][1]) +
                                          as.numeric(strsplit(as.character(x),",")[[1]][2])
                                      }))))/pres
    het<-unlist(lapply(vcf.col,function(x){
      strsplit(as.character(x),split=":")[[1]][1]
    }))
    het<-length(het[het=="0/1" | het=="1/0"])
    return(list(NumMissing=miss,NumPresent=pres,AvgCovRef=ref,AvgCovAlt=alt,AvgCovTot=tot,PropHet=het/pres, NumReads=tot*pres))
  })))
  return(as.data.frame(icov))
}

#' Calculate allele frequencies from a vcf file
#' @param vcf.row A single row of a vcf (used in conjunctioin with apply)
#' @return A data.frame with the columns: Chrom,Pos,Ref,RefFreq,Alt,AltFreq
#'  
#' @export
calc.afs.vcf<-function(vcf.row){
  #use in conjunction with apply
  #e.g. apply(vcf,1,afs.vcf)
  gt1<-unlist(lapply(vcf.row,function(x){
    c<-strsplit(as.character(x),split=":")[[1]][1]
    return(c)
  }))
  gt1<-gt1[gt1 %in% c("0/0","1/0","0/1","1/1")]
  gt1[gt1=="1/0"]<-"0/1"
  gt1<-gsub(pattern = "0",replacement = vcf.row["REF"],gt1)
  gt1<-gsub(pattern = "1",replacement = vcf.row["ALT"],gt1)
  al1<-unlist(strsplit(as.character(gt1),split = "/"))
  #calculate frequencies
  freq1<-summary(factor(al1))/sum(summary(factor(al1)))
  if(length(freq1)==1)
  {
    if(names(freq1)==vcf.row["REF"])
    {
      freq1<-c(freq1,0)
      names(freq1)<-unlist(c(vcf.row["REF"],vcf.row["ALT"]))
    }
    else
    {
      freq1<-c(freq1,0)
      names(freq1)<-unlist(c(vcf.row["ALT"],vcf.row["REF"]))
    }
  }
  if(length(freq1)==0){
    freq1<-c(0,0)
    names(freq1)<-unlist(c(vcf.row["ALT"],vcf.row["REF"]))
  }
  return(data.frame(Chrom=vcf.row["#CHROM"], Pos=vcf.row["POS"], Ref=vcf.row["REF"],
                    RefFreq=freq1[names(freq1) %in% vcf.row["REF"]],
                    Alt=vcf.row["ALT"],AltFreq=freq1[names(freq1) %in% vcf.row["ALT"]]))
}


#' Extract alleles from a vcf row
#' @param vcf.row A vcf row (containing only the individuals and locus info)
#' @return al1 A list of alleles
#' @export
vcf.alleles<-function(vcf.row){
  gt1<-unlist(lapply(vcf.row,function(x){
    c<-strsplit(as.character(x),split=":")[[1]][1]
    return(c)
  }))
  gt1<-gt1[gt1 %in% c("0/0","1/0","0/1","1/1")]
  gt1[gt1%in%"1/0"]<-"0/1"
  gt1<-gsub(pattern = "0",replacement = vcf.row["REF"],gt1)
  gt1<-gsub(pattern = "1",replacement = vcf.row["ALT"],gt1)
  al1<-unlist(strsplit(as.character(gt1),split = "/"))
  return(al1)
}




#' A function to pull out only the genotype fields from vcf-format data
#' @param vcf A data.frame with genotype data in vcf format
#' @return vcf A new vcf with only the GT data
#' @export
extract.gt.vcf<-function(vcf){
  if(length(strsplit(as.character(vcf[1,10]),":")[[1]])>1){
    new<-vcf[,1:9]
    for(i in 10:ncol(vcf)){
      new<-cbind(new,
                 sapply(vcf[,i],function(x) {
                   strsplit(as.character(x),":")[[1]][1]})
      )
    }
    colnames(new)<-colnames(vcf[,c(1:9,10:ncol(vcf))])
    vcf<-new
  }
  return(vcf)
}

#' A function to infer maternal alleles from a vcf file
#' @param dad.kid A data.frame with two columns, each one with matching father and offspring IDs
#' @param vcf A data.frame with genomic data in vcf format
#' @return mat A data.frame of inferred maternal alleles
#' @export
infer.mat.alleles<-function(dad.kid, vcf){
  #dad.kid<-read.table("both.dad.kid.pairs.txt")
  vcf<-extract.gt.vcf(vcf)
  mat<-apply(dad.kid,1,function(x){
    mom<-apply(vcf,1,function(y){
      d1<-strsplit(y[x[1]],"/")[[1]][1]
      d2<-strsplit(y[x[1]],"/")[[1]][2]
      k1<-strsplit(y[x[2]],"/")[[1]][1]
      k2<-strsplit(y[x[2]],"/")[[1]][2]
      mom_allele <- "."
      if (d1 == d2 && k1 == k2) #the case where both are homozygous
      {
        if (d1 == k1)
          mom_allele <- d1
      }
      if (d1 == d2 && k1 != k2) #the case where dad is homozygous but kid is het
      {
        if (d1 == k1)
          mom_allele <- k2
        if (d1 == k2)
          mom_allele <- k1
      }
      if (d1 != d2 && k1 == k2) #the case where dad is het but off is hom
      {
        if (d1 == k1 || d2 == k1)
          mom_allele <- k1;
      }
      if (d1 != d2 && k1 != k2) #if they're both hets you can't do anything with it.
      mom_allele <- "."
      if (d1 == "." || k1 == ".")
        mom_allele <- "."
      return(c(y[1:9],mom_allele))
    })

  })
  return(mat)
}

#' A function to merge two vcf data frames
#' @param vcf1 A dataframe with genomic data in vcf format
#' @param vcf2 A dataframe with genomic data in vcf format
#' @param vcf.name A name for a file with the new vcf file (default is merge.vcf)
#' @return vcf A new data.frame in vcf format
#' @export
combine.vcfs<-function(vcf1,vcf2, vcf.name="merge.vcf"){
  vcf1<-extract.gt.vcf(vcf1)
  vcf2<-extract.gt.vcf(vcf2)
  col.id<-c(colnames(vcf1)[1:9],colnames(vcf1)[!(colnames(vcf1) %in%
                                                   colnames(vcf2))])
  vcf1a<-vcf1[,col.id]
  vcf1a$index<-paste(vcf1a$`#CHROM`,vcf1a$ID,vcf1a$POS,sep=".")
  vcf2$index<-paste(vcf2$`#CHROM`,vcf2$ID,vcf2$POS,sep=".")
  vcf<-merge(vcf1a,vcf2, by="index")
  addedon<-vcf[duplicated(vcf$index),"index"]
  if(!is.null(dim(addedon))) vcf<-vcf[!(vcf$index %in% addedon),]
  drops<-c("index","#CHROM.y","POS.y","ID.y","REF.y","ALT.y","QUAL.y",
           "FILTER.y","INFO.y","FORMAT.y")
  vcf<-vcf[,!(colnames(vcf) %in% drops)]
  colnames(vcf)[1:9]<-colnames(vcf1)[1:9]
  if(vcf.name!=""){
    write.table(vcf,vcf.name,col.names=T,row.names=F,
              quote=F,sep='\t')
  }
  return(vcf)
}

#' Find significant loci using chi-square test
#' @param fst.df A data.frame containing at least the following columns: Fst, Num1 (number of individuals genotyped in pop 1), Num2 (number of individuals genotyped in pop2)
#' @return fst.df The same dataframe but with additional columns
#' @export
fst.sig<-function(fst.df){
  fst.df<-fst.df[!is.na(fst.df$Fst),]
  fst.df$Chi<-2*((fst.df$Num1+fst.df$Num2)/2)*fst.df$Fst
  fst.df$Chi.p<-1-pchisq(fst.df$Chi,1)
  fst.df$Chi.p.adj<-p.adjust(fst.df$Chi.p,method="BH")
  return(fst.df)
}

#' Conduct selection components analysis
#' @param vcf A data.frme with genotype data in vcf format
#' @param locus.info A list of column names with the locus info (e.g. c(#CHROM,POS))
#' @param group1 A list of column names with individuals from the first group
#' @param group2 A list of column names with individuals from the second group
#' @param prop.ind.thresh A proportion of indidivudals requried in each population (default is 0.5)
#' @param maf.cutoff A minimum allele frequency cutoff (default is 0.05)
#' @return sel A dataframe with the Fst values and Chi-squared values
#' @export
gwsca<-function(vcf,locus.info,group1,group2,prop.ind.thresh=0.5,maf.cutoff=0.05){
  sel<-do.call(rbind,apply(vcf,1,fst.one.vcf,group1=c(locus.info,group1),group2=c(locus.info,group2),
    cov.thresh=prop.ind.thresh,maf=maf.cutoff))
  sel<-sel[!is.na(sel$Fst),]
  sel$Chi<-2*((sel$Num1+sel$Num2)/2)*sel$Fst
  sel$Chi.p<-1-pchisq(sel$Chi,1)
  sel$Chi.p.adj<-p.adjust(sel$Chi.p,method="BH")
  return(sel)
}





