


#' Calculate pairwise fst between two separate vcf files
#' @param vcf1.row A row of a data.frame containing genotype information in vcf format
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

  return(data.frame(Chrom=vcf1.row[1],Pos=vcf1.row["POS"],
                    Hs1=hs1,Hs2=hs2,Hs=hs,Ht=ht,Fst=fst,NumAlleles=length(factor(freqall)),
                    Num1=length(gt1),Num2=(length(gt2)),stringsAsFactors=F))
}#end function


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
  freqall<-mean(c(freq1[1],freq2[1]))
  if(length(freq1)>1 & length(freq2)>1){ #both must be polymorphic
    hs1<-1-sum(freq1*freq1)
    hs2<-1-sum(freq2*freq2)
    hw<-sum(c(hs1,hs2))
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
  freqall<-mean(c(freq1[1],freq2[1]))
  if(length(freq1)>1 & length(freq2)>1){ #both must be polymorphic
    hs1<-1-sum(freq1*freq1)
    hs2<-1-sum(freq2*freq2)
    # hs<-(hs1*length(al1)+hs2*length(al2))/(length(al1)+length(al2))
    hs<-mean(c(hs1,hs2))
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
  return(data.frame(Hs1=hs1,Hs2=hs2,Hs=hs,Ht=ht,Fst=as.numeric(fst),NumAlleles=length(unique(al12)),
                    Num1=length(al1),Num2=length(al2)))
}

#' Calculate fsts from a single vcf
#' @param vcf.row A row of a data.frame in vcf format
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
  num.ind1<-length(gt1[gt1 %in% c("0/0","1/0","0/1","1/1","./.")])
  gt1<-gt1[gt1 %in% c("0/0","1/0","0/1","1/1")]
  gt1[gt1=="1/0"]<-"0/1"
  gt1<-gsub(pattern = "0",replacement = vcf.row["REF"],gt1)
  gt1<-gsub(pattern = "1",replacement = vcf.row["ALT"],gt1)
  al1<-unlist(strsplit(as.character(gt1),split = "/"))
  if(length(gt1)/num.ind1>=cov.thresh){
    gt2<-unlist(lapply(vcf.row[group2],function(x){
      c<-strsplit(as.character(x),split=":")[[1]][1]
      return(c)
    }))
    num.ind2<-length(gt2[gt2 %in% c("0/0","1/0","0/1","1/1","./.")])
    gt2<-gt2[gt2 %in% c("0/0","1/0","0/1","1/1")]
    gt2[gt2=="1/0"]<-"0/1"
    gt2<-gsub(pattern = "0",replacement = vcf.row["REF"],gt2)
    gt2<-gsub(pattern = "1",replacement = vcf.row["ALT"],gt2)
    al2<-unlist(strsplit(as.character(gt2),split="/"))
    if(length(gt2)/num.ind2>=cov.thresh){
      #calculate frequencies
      freq1<-summary(factor(al1))/sum(summary(factor(al1)))
      freq2<-summary(factor(al2))/sum(summary(factor(al2)))
      freqall<-mean(c(freq1[1],freq2[1]))
      # freqall<-summary(as.factor(c(al1,al2)))/
      #  sum(summary(as.factor(c(al1,al2))))
      hets<-c(names(freq1)[2],names(freq2)[2])
      if(length(freq1)>1 & length(freq2)>1 #both must be polymorphic
         & min(freq1,freq2) >= maf){ #and match the maf
        hs1<-2*freq1[1]*freq1[2]
        hs2<-2*freq2[1]*freq2[2]
        #hs<-((hs1*length(gt1))+(hs2*length(gt2)))/(length(gt1)+length(gt2)) #weighted avg
        hs<-mean(c(hs1,hs2))
        ht<-2*freqall*(1-freqall)
        fst<-(ht-hs)/ht
        num.al<-length(unique(c(al1,al2)))
      } else {

        hs1<-1-sum(freq1[1]*freq1[1])
        hs2<-1-sum(freq2[1]*freq2[1])
        hs<-mean(c(hs1,hs2))
        ht<-2*freqall*(1-freqall)
        if(length(freqall)<=1){ fst<-0 }
        else{ fst<-NA }
        num.al<-length(unique(c(al1,al2)))
      }
    }
    else {
      fst<-NA #gt2 doesn't pass coverage threshold
      num.al<-length(unique(c(al1,al2)))
    }
  }else {
    fst<-NA #it doesn't pass the coverage threshold
    num.al<-length(unique(al1))
  }

  return(data.frame(Chrom=vcf.row[1],Pos=vcf.row["POS"],
                    Hs1=hs1,Hs2=hs2,Hs=hs,Ht=ht,Fst=fst,NumAlleles=num.al,
                    Num1=length(gt1),Num2=length(gt2),stringsAsFactors = FALSE,row.names=NULL))
}#end function fst.one.vcf


#' Calculate pairwise fsts from a dataset in ped format
#' @param raw A data.frame with data in ped format
#' @param group1 A list with the individuals in group 1
#' @param group2 A list with the individuals in group 2
#' @param cov.thresh A threshold for the number of individuals in the populations (default is 0.2)
#' @param loc_names Locus names to be included in the output (optional)
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
fst.one.plink<-function(raw,group1, group2, cov.thresh=0.2,loc_names=NULL){
  fst.dat<-data.frame(Locus=character(),
                      Hs1=numeric(),Hs2=numeric(),Hs=numeric(),Ht=numeric(),Fst=numeric(),NumAlleles=numeric(),
                      Num1=numeric(),Num2=numeric(),stringsAsFactors=F)
  grp1<-raw[raw[,2] %in% group1,]
  grp2<-raw[raw[,2] %in% group2,]
  n<-1
  for(i in seq(7,ncol(raw),by = 2)){
    # get all the alleles
    alleles<-unique(unlist(grp1[,i:(i+1)],grp2[,i:(i+1)]))
    alleles<-alleles[alleles!=0]
    # figure out how many are NA in each group
    na1<-length(grp1[grp1[,i]==0,i])/nrow(grp1)
    na2<-length(grp2[grp2[,i]==0,i])/nrow(grp2)
    # get the alleles for each group
    al1<-unlist(c(as.character(grp1[grp1[,i]%in% alleles,i]),as.character(grp1[grp1[,i+1]%in% alleles,i+1])))
    al2<-unlist(c(as.character(grp2[grp2[,i]%in% alleles,i]),as.character(grp2[grp2[,i+1]%in% alleles,i+1])))

    if(na1<=(1-cov.thresh) & na2<=(1-cov.thresh)){
      #calculate fst
      fst<-calc.fst.wright(al1,al2)
    }
    else {
      fst<-data.frame(Hs1=NA,Hs2=NA,Hs=NA,Ht=NA,Fst=NA,NumAlleles=length(factor(c(al1,al2))),
                      Num1=length(al1),Num2=length(al2)) #gt2 doesn't pass coverage threshold
    }

    if(is.null(loc_names)){
      fst.dat[n,]<-cbind(n,fst["Hs1"],fst["Hs2"],as.numeric(fst["Hs"]),fst["Ht"],
                         as.numeric(fst["Fst"]),fst["NumAlleles"],fst["Num1"]/2,fst["Num2"]/2)
    }else{
      fst.dat[n,]<-cbind(loc_names[n],fst["Hs1"],fst["Hs2"],as.numeric(fst["Hs"]),fst["Ht"],
                         as.numeric(fst["Fst"]),fst["NumAlleles"],fst["Num1"]/2,fst["Num2"]/2)
    }

    n<-n+1
  }
  return(fst.dat)
}#end fst.one.plink

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
      freqall<-mean(c(freq1[1],freq2[1]))
      if(length(freq1)>1){ hs1<-2*freq1[1]*freq1[2]
      } else {
        hs1<-0
      }
      if(length(freq2)>1){ hs2<-2*freq2[1]*freq2[2]
      } else {
        hs2<-0
      }
      if(length(freqall)>1){
        # hs<-(hs1*length(pop1)+hs2*length(pop2))/(length(pop1)+length(pop2))
        hs<-mean(c(hs1,hs2))
        ht<-2*freqall*(1-freqall)
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

#' Calculate mean pairwise Fsts
#' @param vcf A data.frame containing the contents of a vcf file
#' @param pop.list1 A list of populations. This will be used to subset the vcf file to extract individuals for each pairwise comparison so it must contain a string found in individual IDs.
#' @param pop.list2 A list of populations. Same requirements as pop.list1. Individuals from populations in pop.list1 will be compared to individuals in pop.list2.
#' @param maf.cutoff An optional value for the minimum allele frequency. Must be between 0 and 1 (default is 0.05)
#' @param cov.thresh Coverage threshold between 0 and 1 (default is 0.2)
#' @return mu A data.frame containing the columns: Chrom, Pos, SNP, Sum.Fst, Count, and Mean.Fst. Sum.Fst is the total sum of Fst values, which are divided by Count (the number of comparisons the locus was present in) to generate Mean.Fst.
#' @export
calc.mean.fst <- function(vcf,pop.list1,pop.list2, maf.cutoff = 0.05,cov.thresh=0.2) {
  mu<-data.frame(Chrom=vcf$`#CHROM`,Pos=vcf$POS,SNP=vcf$SNP,
                 Sum.Fst=rep(0,nrow(vcf)),Count=rep(0,nrow(vcf)),
                 Mean.Fst=rep(0,nrow(vcf)),stringsAsFactors=FALSE)
  for(i in 1:length(pop.list1)){
    for(j in 1:length(pop.list2)){
      fsts<-gwsca(vcf=vcf,locus.info=locus.info,
                  group1=colnames(vcf)[grep(pop.list1[i],colnames(vcf))],
                  group2=colnames(vcf)[grep(pop.list2[j],colnames(vcf))],
                  maf.cutoff = maf.cutoff,prop.ind.thresh = cov.thresh)
      fsts$Fst[fsts$Fst==0]<-NA #replace ones that weren't calc'd with NA
      fsts$SNP<-paste(fsts$Chrom,as.numeric(as.character(fsts$Pos)),sep=".")
      new.mu<-do.call("rbind",apply(mu,1,function(x){
        new.x<-data.frame(Chrom=x["Chrom"],Pos=x["Pos"],SNP=x["SNP"],
                          Sum.Fst=as.numeric(x["Sum.Fst"]),Count=as.numeric(x["Count"]),
                          Mean.Fst=as.numeric(x["Mean.Fst"]),stringsAsFactors=FALSE)
        this.fst<-fsts[fsts$SNP %in% x["SNP"],]

        if(nrow(this.fst)>0){
          for(k in 1:nrow(this.fst)){
            if(!is.na(this.fst[k,"Fst"])){
              if(!is.na(new.x["Sum.Fst"])){
                new.x["Sum.Fst"]<-new.x["Sum.Fst"]+this.fst[k,"Fst"]
                new.x["Count"]<-new.x["Count"]+1
              }else{
                new.x["Sum.Fst"]<-this.fst[k,"Fst"]
                new.x["Count"]<-new.x["Count"]+1
              }}
          }
        }
        if(nrow(this.fst)>1){ #if there are multiple records, keep the one with the maximum sample size
          this.fst<-this.fst[this.fst$Num1+this.fst$Num2 == max(this.fst$Num1+this.fst$Num2),]
          if(nrow(this.fst)>1){#if they have the same sample size, randomly choose one.
            this.fst<-this.fst[sample(nrow(this.fst),size = 1,replace = FALSE),]
          }
        }
        if(nrow(this.fst)==1){
          if(!is.na(this.fst["Fst"])){
            if(!is.na(new.x["Sum.Fst"])){
              new.x["Sum.Fst"]<-new.x["Sum.Fst"]+this.fst["Fst"]
              new.x["Count"]<-new.x["Count"]+1
            }else{
              new.x["Sum.Fst"]<-this.fst["Fst"]
              new.x["Count"]<-new.x["Count"]+1
            }
          }}

        return(new.x)
      }))
      mu<-new.mu
    }
  }
  mu$Mean.Fst<-mu$Sum.Fst/mu$Count
  return(mu)
}
