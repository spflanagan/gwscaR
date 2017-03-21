#Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
#Purpose: Calculate useful population genetics statistics, plot genome-wide statistics,
#and run Fst-based selection components analysis

#' Plot genome-wide statistics from a data frame.
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
fst.plot<-function(fst.dat,ci.dat=NULL, sig.col=c("red","yellow"),pt.col="grey7",
                   fst.name="Fst", chrom.name="Chrom", bp.name="BP",axis.size=0.5,
                   scaffold.order=NULL,groups=NULL,print.names=FALSE,y.lim=NULL,
                   group.boundaries=NULL,pt.cex=0.5,...){

  if(!is.null(scaffold.order)){
    scaff.ord<-scaffold.order$component_id
    lgs<-scaffold.order$object
  } else{
    if(!is.null(group.boundaries)){
      scaff.ord<-levels(factor(group.boundaries[,1]))
      lgs<-scaff.ord
    }else{
      scaff.ord<-levels(factor(fst.dat[,chrom.name]))
      lgs<-scaff.ord
    }
  }
  if(!is.null(groups)){
    lgs<-groups
    scaff.ord<-groups
  }
  fst.dat[,fst.name]<-as.numeric(as.character(fst.dat[,fst.name]))
  these.scaffs<-scaff.ord[scaff.ord %in% factor(fst.dat[,chrom.name])]
  all.scaff<-split(fst.dat, factor(fst.dat[,chrom.name]))

  group.boundaries<-group.boundaries[scaff.ord,]
  last.max<-0
  rect.xs<-NULL
  addition.values<-0
  xlist<-NULL
  #set up the rectangles
  for(i in 1:length(scaff.ord)){
    if(is.null(group.boundaries)){
      new.max<-last.max+nrow(all.scaff[[scaff.ord[i]]])
    }else{
      new.max<-as.numeric(group.boundaries[scaff.ord[i],2])+as.numeric(last.max)
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
        seq(rect.xs[rownames(rect.xs) %in% names(all.scaff)[i],1]+1,rect.xs[rownames(rect.xs) %in% names(all.scaff)[i],1]+nrow(all.scaff[[i]]),1)
    }else{
      #then it takes genome space into account

      all.scaff[[i]]<-
        all.scaff[[i]][order(all.scaff[[i]][,bp.name]),]
      all.scaff[[i]][,bp.name]<-
        as.numeric(rect.xs[rownames(rect.xs) %in% names(all.scaff)[i],1])+as.numeric(all.scaff[[i]][,bp.name])

    }
  }
  #change BP to plot
  x.max<-max(rect.xs)
  x.min<-min(rect.xs)
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


#' Read in a vcf file
#' @param filename The name of the vcf file
#' @return a dataframe containing the contents of the vcf file, including headers.
#' @examples
#' vcf<-parse.vcf(system.file("extdata/example.vcf",package = "gwscaR"))
#' @seealso Flanagan & Jones 2017
#' @export
parse.vcf<-function(filename){
  #if(substr(filename,nchar(filename)-3,nchar(filename)) != ".vcf") { filename<-paste(filename,"vcf",sep=".") }
  vcf<-read.delim(filename,comment.char="#",sep='\t',header=F,stringsAsFactors = F)
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
  cov.dat<-do.call("rbind",apply(vcf[,subset],1,function(vcf.row){
    cov<-unlist(lapply(vcf.row,function(x){
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
    return(data.frame(Chrom=vcf.row[1],Pos=vcf.row["POS"],Locus=vcf.row["ID"],
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

#' Calculate pairwise fst between two separate vcf files
#' @param vcf1 A data.frame containing genotype information in vcf format
#' @param vcf2 A data.frame containing genotype information in vcf format (for another set of individuals/samples)
#' @param match.index The name of the column that should be used to pair the two vcfs (I usually create a new column in which Chrom and Pos are concatenated)
#' @param cov.thresh A coverage threshold specifying the proportion of individuals that should be included in each population. The default is 0.2
#' @return A data.frame with the columns:
#'  Chrom
#'  Pos
#'  Hs1 = expected heterozygosity in population 1 (vcf1)
#'  Hs2 = expected heterozygosity in population 2 (vcf2)
#'  Hs = weighted average expected heterozygosity within populations
#'  Ht = expected heterozygosity among populations
#'  Fst= Fst
#'  NumAlleles=The number of alleles (2 for biallelic loci)
#'  Num1=number of individuals genotyped & included in Fst calculations in population 1 (vcf1)
#'  Num2=number of individuals genotyped & included in Fst calculations in population 2 (vcf2)
#'  @export

fst.two.vcf<-function(vcf1.row,vcf2,match.index, cov.thresh=0.2){
  #match.index is the column used to match the two
  #use in conjunction with apply
  #e.g. apply(vcf,1,fst.two.vcf,vcf2=vcf.2,match.index="SNP")
  hs1<-hs2<-hs<-ht<-0
  freqall<-gt1<-gt2<-NULL
  vcf2.row<-vcf2[vcf2[,match.index]%in%vcf1.row[match.index],]
  if(nrow(vcf2.row)>1)#first make sure we have one reading per locus
  {
    print("Multiple instances in vcf2.")
    fst<-NA
  }
  else{
    if(nrow(vcf2.row)==0)
    {
      print("No instances in vcf2.")
      fst<-NA
    }else #we're good to go
    {
      gt1<-unlist(lapply(vcf1.row,function(x){
        c<-strsplit(as.character(x),split=":")[[1]][1]
        return(c)
      }))
      num.ind<-length(gt1)-10
      gt1<-gt1[gt1 %in% c("0/0","1/0","0/1","1/1")]
      gt1[gt1=="1/0"]<-"0/1"
      gt1<-gsub(pattern = "0",replacement = vcf1.row["REF"],gt1)
      gt1<-gsub(pattern = "1",replacement = vcf1.row["ALT"],gt1)
      if(length(gt1)/num.ind>=cov.thresh){
        al1<-unlist(strsplit(as.character(gt1),split = "/"))
        gt2<-unlist(lapply(vcf2.row,function(x){
          c<-strsplit(as.character(x),split=":")[[1]][1]
          return(c)
        }))
        num.ind<-length(gt2)
        gt2<-gt2[gt2 %in% c("0/0","1/0","0/1","1/1")]
        gt2[gt2=="1/0"]<-"0/1"
        gt2<-gsub(pattern = "0",replacement = vcf2.row["REF"],gt2)
        gt2<-gsub(pattern = "1",replacement = vcf2.row["ALT"],gt2)
        if(length(gt2)/num.ind>=cov.thresh){
          al2<-unlist(strsplit(as.character(gt2),split="/"))
          #calculate frequencies
          freq1<-summary(factor(al1))/sum(summary(factor(al1)))
          freq2<-summary(factor(al2))/sum(summary(factor(al2)))
          freqall<-summary(as.factor(c(al1,al2)))/
            sum(summary(as.factor(c(al1,al2))))
          hets<-c(names(freq1)[2],names(freq2)[2])
          if(length(freq1)>1 & length(freq2)>1){ #both must be polymorphic
            hs1<-2*freq1[1]*freq1[2]
            hs2<-2*freq2[1]*freq2[2]
            hs<-mean(c(hs1,hs2))
            ht<-2*freqall[1]*freqall[2]
            fst<-(ht-hs)/ht
          } else {
            hs1<-1-sum(freq1*freq1)
            hs2<-1-sum(freq2*freq2)
            if(length(freqall)<=1){ fst<-0 }
            else{
              ht<-2*freqall[1]*freqall[2]
              fst<-NA
            }
          }
        }
        else {
          fst<-NA #gt2 doesn't pass coverage threshold
        }
      }else {
        fst<-NA #it doesn't pass the coverage threshold
      }
    }#end else good to go
  }#end else vcf2

  return(data.frame(Chrom=vcf1.row["#CHROM"],Pos=vcf1.row["POS"],
                    Hs1=hs1,Hs2=hs2,Hs=hs,Ht=ht,Fst=fst,NumAlleles=length(factor(freqall)),
                    Num1=length(gt1),Num2=(length(gt2))))
}#end function

#' Calculate allele frequencies from a vcf file
#' @param vcf.row A single row of a vcf (used in conjunctioin with apply)
#' @return A data.frame with the columns:
#'  Chrom
#'  Pos
#'  Ref
#'  RefFreq 
#'  Alt 
#'  AltFreq 
#'  @export
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

#' Calculate Fst using Nei's formulation (1-(Hw/2Hb))
#' @param al1 A table of allele frequencies for pop 1
#' @param al2 A table of allele frequencies for pop 2
#' @return data.frame containing the columns:
#'  Hs1 = expected heterozygosity in population 1 (vcf1)
#'  Hs2 = expected heterozygosity in population 2 (vcf2)
#'  Hs = weighted average expected heterozygosity within populations
#'  Ht = expected heterozygosity among populations
#'  Fst= Fst
#'  NumAlleles=The number of alleles (2 for biallelic loci)
#'  Num1=number of individuals genotyped & included in Fst calculations in population 1 (vcf1)
#'  Num2=number of individuals genotyped & included in Fst calculations in population 2 (vcf2)
#'  @export
calc.fst.nei<-function(al1,al2){
  hw<-hb<-fst<-0
  freq1<-summary(factor(al1))/sum(summary(factor(al1)))
  freq2<-summary(factor(al2))/sum(summary(factor(al2)))
  al12<-c(al1,al2)
  freqall<-summary(as.factor(al12))/
    sum(summary(as.factor(al12)))
  if(length(freq1)>1 & length(freq2)>1){ #both must be polymorphic
    hs1<-1-sum(freq1*freq1)
    hs2<-1-sum(freq2*freq2)
    hw<-sum(hs1,hs2)
    hb<-1-sum(freqall*freqall)
    fst<-1-(hw/(2*hb))
  } else {
    hs1<-1-sum(freq1*freq1)
    hs2<-1-sum(freq2*freq2)
    if(length(freqall)<=1){
      hw<-1
      fst<-NA
    } else{
      hw<-1-sum(freqall*freqall)
      fst<-NA
    }
  }
  return(data.frame(Hs1=hs1,Hs2=hs2,Hs=hb,Ht=hw,Fst=as.numeric(fst),NumAlleles=length(factor(freqall)),
                    Num1=length(al1),Num2=length(al2)))
}

#' Calculate Fst using Wright (1943)'s formulation ((ht-hs)/ht)
#' @param al1 A table of allele frequencies for pop 1
#' @param al2 A table of allele frequencies for pop 2
#' @return data.frame containing the columns:
#'  Hs1 = expected heterozygosity in population 1 (vcf1)
#'  Hs2 = expected heterozygosity in population 2 (vcf2)
#'  Hs = weighted average expected heterozygosity within populations
#'  Ht = expected heterozygosity among populations
#'  Fst= Fst
#'  NumAlleles=The number of alleles (2 for biallelic loci)
#'  Num1=number of individuals genotyped & included in Fst calculations in population 1 (vcf1)
#'  Num2=number of individuals genotyped & included in Fst calculations in population 2 (vcf2)
#'  @export
calc.fst.wright<-function(al1,al2){
  hs<-ht<-fst<-0
  freq1<-summary(factor(al1))/sum(summary(factor(al1)))
  freq2<-summary(factor(al2))/sum(summary(factor(al2)))
  al12<-c(al1,al2)
  freqall<-summary(as.factor(al12))/
    sum(summary(as.factor(al12)))
  if(length(freq1)>1 & length(freq2)>1){ #both must be polymorphic
    hs1<-1-sum(freq1*freq1)
    hs2<-1-sum(freq2*freq2)
    hs<-(hs1*length(al1)+hs2*length(al2))/(length(al1)+length(al2))
    ht<-1-sum(freqall*freqall)
    fst<-(ht-hs)/ht
  } else {
    hs1<-1-sum(freq1*freq1)
    hs2<-1-sum(freq2*freq2)
    if(length(freqall)<=1){
      ht<-1
      fst<-NA
    } else{
      ht<-1-sum(freqall*freqall)
      fst<-NA
    }
  }
  return(data.frame(Hs1=hs1,Hs2=hs2,Hs=hs,Ht=ht,Fst=as.numeric(fst),NumAlleles=length(factor(freqall)),
                    Num1=length(al1),Num2=length(al2)))
}

#' Calculate fsts from a single vcf
#' @param vcf A data.frame in vcf format
#' @param group1 A list of individual names for pop1
#' @param group2 A list of individual names for pop2
#' @param cov.thresh A proportion of individuals required for genotyping individuals (default is 0.2)
#' @param maf A minimum allele frequency (default is 0.05)
#' @return out A data.frame with columns:
#' Hs1 (expected heterozygosity in pop 1), Hs2 (expected heterozygosity in pop2),
#' Hs (expected heterozygosity within pops), Ht (expected heterozygosity between pops),
#' Fst, NumAlleles, Num1 (number of individuals in pop1), and Num2 (number of individuals in pop2)
#' @export
fst.one.vcf<-function(vcf.row,group1,group2, cov.thresh=0.2, maf=0.05){

  hs1<-hs2<-hs<-ht<-0
  freqall<-gt1<-gt2<-NULL
  gt1<-unlist(lapply(vcf.row[group1],function(x){
    c<-strsplit(as.character(x),split=":")[[1]][1]
    return(c)
  }))
  num.ind<-length(gt1)
  gt1<-gt1[gt1 %in% c("0/0","1/0","0/1","1/1")]
  gt1[gt1=="1/0"]<-"0/1"
  gt1<-gsub(pattern = "0",replacement = vcf.row["REF"],gt1)
  gt1<-gsub(pattern = "1",replacement = vcf.row["ALT"],gt1)
  if(length(gt1)/num.ind>=cov.thresh){
    al1<-unlist(strsplit(as.character(gt1),split = "/"))
    gt2<-unlist(lapply(vcf.row[group2],function(x){
      c<-strsplit(as.character(x),split=":")[[1]][1]
      return(c)
    }))
    num.ind<-length(gt2)
    gt2<-gt2[gt2 %in% c("0/0","1/0","0/1","1/1")]
    gt2[gt2=="1/0"]<-"0/1"
    gt2<-gsub(pattern = "0",replacement = vcf.row["REF"],gt2)
    gt2<-gsub(pattern = "1",replacement = vcf.row["ALT"],gt2)
    if(length(gt2)/num.ind>=cov.thresh){
      al2<-unlist(strsplit(as.character(gt2),split="/"))
      #calculate frequencies
      freq1<-summary(factor(al1))/sum(summary(factor(al1)))
      freq2<-summary(factor(al2))/sum(summary(factor(al2)))
      freqall<-summary(as.factor(c(al1,al2)))/
        sum(summary(as.factor(c(al1,al2))))
      hets<-c(names(freq1)[2],names(freq2)[2])
      if(length(freq1)>1 & length(freq2)>1 #both must be polymorphic
         & min(freq1,freq2) >= maf){ #and match the maf
        hs1<-2*freq1[1]*freq1[2]
        hs2<-2*freq2[1]*freq2[2]
        hs<-mean(c(hs1,hs2))
        ht<-2*freqall[1]*freqall[2]
        fst<-(ht-hs)/ht
      } else {
        hs1<-1-sum(freq1*freq1)
        hs2<-1-sum(freq2*freq2)
        if(length(freqall)<=1){ fst<-0 }
        else{
          ht<-2*freqall[1]*freqall[2]
          fst<-NA
        }
      }
    }
    else {
      fst<-NA #gt2 doesn't pass coverage threshold
    }
  }else {
    fst<-NA #it doesn't pass the coverage threshold
  }

  return(data.frame(Chrom=vcf.row["#CHROM"],Pos=vcf.row["POS"],
                    Hs1=hs1,Hs2=hs2,Hs=hs,Ht=ht,Fst=fst,NumAlleles=length(factor(freqall)),
                    Num1=length(gt1),Num2=length(gt2)))
}#end function fst.one.vcf

#' A function that calculates the allele frequencies from a vcf
#' @param vcf A data.frame with genotype data in vcf format
#' @return out A data.frame with columns Chrom, Pos, Ref (the reference allele),
#' RefFreq (frequency of reference allele), Alt (alternative allele), and AltFreq (frequency of the alternative allele)
#' @export
calc.afs.vcf<-function(vcf){
  #use in conjunction with apply
  out<-apply(vcf,1,function(vcf.row){
    al1<-vcf.alleles(vcf.row)
    #calculate frequencies
    if(length(al1)>0){
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
      return(data.frame(Chrom=vcf.row["#CHROM"], Pos=vcf.row["POS"], Ref=vcf.row["REF"],
                        RefFreq=freq1[names(freq1) %in% vcf.row["REF"]],
                        Alt=vcf.row["ALT"],AltFreq=freq1[names(freq1) %in% vcf.row["ALT"]]))
    }else{
      return(data.frame(Chrom=vcf.row["#CHROM"], Pos=vcf.row["POS"], Ref=vcf.row["REF"],
                        RefFreq=0,Alt=vcf.row["ALT"],AltFreq=0))
    }
  })
  return(out)

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


#' Calculate pairwise fsts from a dataset in ped format
#' @param raw A data.frame with data in ped format
#' @param group1 A list with the individuals in group 1
#' @param group2 A list with the individuals in group 2
#' @param cov.thresh A threshold for the number of individuals in the populations (default is 0.2)
#' @return fst.dat A data.frame with columns:
#'  Locus Name = the locus name
#'  Hs1 = expected heterozygosity in pop 1
#'  Hs2 = expected heterozygosity in pop 2
#'  Hs = expected heterozygosity within populations
#'  Ht = expected heterozygosity among populations
#'  Fst = Fst
#'  NumAlleles = number of alleles at the locus
#'  Num1 = the number of individuals in pop 1
#'  Num2 = the number of individuals in pop 2
#' @export
fst.one.plink<-function(raw,group1, group2, cov.thresh=0.2){
  fst.dat<-data.frame(Locus=character(),
                      Hs1=numeric(),Hs2=numeric(),Hs=numeric(),Ht=numeric(),Fst=numeric(),NumAlleles=numeric(),
                      Num1=numeric(),Num2=numeric(),stringsAsFactors=F)
  grp1<-raw[raw$IID %in% group1,]
  grp2<-raw[raw$IID %in% group2,]
  for(i in 7:ncol(raw)){
    na1<-length(grp1[is.na(grp1[,i]),i])/nrow(grp1)
    na2<-length(grp2[is.na(grp2[,i]),i])/nrow(grp2)
    gt1<-grp1[!is.na(grp1[,i]),i]
    gt2<-grp2[!is.na(grp2[,i]),i]
    gt1[gt1=="1"]<-"1/2"
    gt1[gt1=="2"]<-"2/2"
    gt1[gt1=="0"]<-"1/1"
    gt2[gt2=="1"]<-"1/2"
    gt2[gt2=="2"]<-"2/2"
    gt2[gt2=="0"]<-"1/1"

    if(na1<=(1-cov.thresh)){
      al1<-unlist(strsplit(as.character(gt1),split = "/"))
      if(na2<=(1-cov.thresh)){
        al2<-unlist(strsplit(as.character(gt2),split="/"))
        #calculate frequencies
        fst<-calc.fst.wright(al1,al2)
      }
      else {
        fst<-data.frame(Hs1=NA,Hs2=NA,Hs=NA,Ht=NA,Fst=NA,NumAlleles=length(factor(c(al1,al2))),
                        Num1=length(al1),Num2=length(al2)) #gt2 doesn't pass coverage threshold
      }
    }else {
      fst<-data.frame(Hs1=NA,Hs2=NA,Hs=NA,Ht=NA,Fst=NA,NumAlleles=length(factor(c(al1,al2))),
                      Num1=length(al1),Num2=length(al2)) #it doesn't pass the coverage threshold
    }
    fst.dat[(i-6),]<-cbind(as.character(colnames(raw)[i]),fst["Hs1"],fst["Hs2"],as.numeric(fst["Hs"]),fst["Ht"],
                           as.numeric(fst["Fst"]),fst["NumAlleles"],fst["Num1"],fst["Num2"])

  }
  return(fst.dat)
}#end fst.one.plink

#' A function to pull out only the genotype fields from vcf-format data
#' @param vcf A data.frame with genotype data in vcf format
#' @return vcf A new vcf with only the GT data
#' @export
extract.gt.vcf<-function(vcf){
  if(length(strsplit(as.character(vcf[1,10]),":")[[1]])>1){
    new<-vcf[,1:3]
    for(i in 10:ncol(vcf)){
      new<-cbind(new,
                 sapply(vcf[,i],function(x) {
                   strsplit(as.character(x),":")[[1]][1]})
      )
    }
    colnames(new)<-colnames(vcf[,c(1:3,10:ncol(vcf))])
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
  col.id<-c(colnames(vcf1)[1:3],colnames(vcf1)[!(colnames(vcf1) %in%
                                                   colnames(vcf2))])
  vcf1a<-vcf1[,col.id]
  vcf1a$index<-paste(vcf1a$`#CHROM`,vcf1a$ID,vcf1a$POS,sep=".")
  vcf2$index<-paste(vcf2$`#CHROM`,vcf2$ID,vcf2$POS,sep=".")
  vcf<-merge(vcf1a,vcf2, by="index")
  addedon<-vcf[duplicated(vcf$index),"index"]
  if(!is.null(dim(addedon))) vcf<-vcf[!(vcf$index %in% addedon),]
  drops<-c("index","#CHROM.y","POS.y","ID.y")
  vcf<-vcf[,!(colnames(vcf) %in% drops)]
  colnames(vcf)[1:3]<-c("CHROM","POS","ID")
  write.table(vcf,vcf.name,col.names=T,row.names=F,
              quote=F,sep='\t')
  return(vcf)
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
  sel<-fst.one.vcf(vcf,c(locus.info,group1),c(locus.info,group2),
    cov.thresh=prop.ind.thresh,maf=maf.cutoff)
  sel<-sel[!is.na(sel$Fst),]
  sel$Chi<-2*((sel$Num1+sel$Num2)/2)*sel$Fst
  sel$Chi.p<-1-pchisq(sel$Chi,1)
  sel$Chi.p.adj<-p.adjust(sel$Chi.p,method="BH")
  return(sel)
}

#' Calculate pairwise Fst values
#' @param ped A dataframe with data in ped format, where each allele for each locus is in its own column
#' @param allele1 An index for the first allele at the locus
#' @param allele2 An index for the second allele at the locus
#' @param pop.order A list with the Pop IDs in the correct order
#' @return dat.var A matrix of pairwise fst values (calculated as (ht-hs)/ht)
#' @export
pairwise.fst<-function(ped,allele1,allele2,pop.order){
  #V1 of ped should be pop index
  ped.split<-split(ped[,c(allele1,allele2)], factor(ped[,1]))
  dat.var<-as.data.frame(setNames(
    replicate(length(pop.order),numeric(0), simplify = F), pop.order))
  for(i in 1:(length(pop.order)-1)){
    for(j in (i+1):length(pop.order)){
      pop1<-factor(ped.split[[pop.order[i]]][ped.split[[pop.order[i]]]!="0"])
      pop2<-factor(ped.split[[pop.order[j]]][ped.split[[pop.order[j]]]!="0"])
      freq1<-summary(pop1)/sum(summary(pop1))
      freq2<-summary(pop2)/sum(summary(pop2))
      freqall<-summary(as.factor(c(pop1,pop2)))/
        sum(summary(as.factor(c(pop1,pop2))))
      if(length(freq1)>1){ hs1<-2*freq1[1]*freq1[2]
      } else {
        hs1<-0
      }
      if(length(freq2)>1){ hs2<-2*freq2[1]*freq2[2]
      } else {
        hs2<-0
      }
      if(length(freqall)>1){
        hs<-(hs1*length(pop1)+hs2*length(pop2))/(length(pop1)+length(pop2))
        ht<-2*freqall[1]*freqall[2]
        fst<-(ht-hs)/ht
      }
      if(length(freqall)<=1){ fst<-1 }
      dat.var[pop.order[i],pop.order[j]]<-fst
    }
  }
  dat.var<-rbind(dat.var,rep(NA, ncol(dat.var)))
  rownames(dat.var)<-colnames(dat.var)
  return(as.matrix(dat.var))
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



