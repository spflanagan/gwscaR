#Author: Sarah P. Flanagan (spflanagan.phd@gmail.com)
#Purpose: Calculate useful population genetics statistics, plot genome-wide statistics, 
#and run Fst-based selection components analysis

plot.genome.wide<-function(bp,var,y.max,x.max, rect.xs=NULL,y.min=0,x.min=0, 
	plot.new=FALSE, plot.axis=TRUE, rect.color="white",plot.rect=TRUE, 
	pt.cex=1, pt.col="black"){
	#********************************************
	#this function plots a variable without scaffold info. 
	#feed it the basepair (x) values and variable (y) values 
	#*********************************************
	if(plot.new==TRUE){ par(new=new) }
	plot(bp, var,xlab="",ylab="", new=plot.new,
		type="n", bg="transparent", axes=F, bty="n", 
		xlim=c(x.min,x.max),ylim=c(y.min, y.max))
	if(plot.rect==TRUE){
		num.rect<-nrow(rect.xs)
		if(is.null(num.rect)) {
			rect(rect.xs[1],y.min,rect.xs[2],y.max, 
				col=rect.color, border=NA)
		} else {
			for(i in 1:nrow(rect.xs)){
				rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
					col=rect.color, border=NA)
			}
		}
	}
	if(plot.axis){
	axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
		ylim = c(y.min, y.max), pos=0,
		las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)}
	points(bp, var, pch=19, cex=pt.cex,col=pt.col,
		xlim=c(x.min,x.max),ylim=c(y.min, y.max))
}

fst.plot<-function(fst.dat,ci.dat=NULL, sig.col=c("red","yellow"),pt.col="grey7",
	fst.name="Fst", chrom.name="Chrom", bp.name="BP",axis.size=0.5,
	scaffold.order=NULL,groups=NULL,print.names=FALSE,y.lim=NULL){
	if(!is.null(scaffold.order)){
		scaff.ord<-scaffold.order$component_id
		lgs<-scaffold.order$object
	} else{
		scaff.ord<-levels(factor(fst.dat[,chrom.name]))
		lgs<-scaff.ord
	}
	if(!is.null(groups)){
		lgs<-groups
		scaff.ord<-groups
	}
  fst.dat[,fst.name]<-as.numeric(as.character(fst.dat[,fst.name]))
  
	all.scaff<-split(fst.dat, factor(fst.dat[,chrom.name]))
	last.max<-0
	rect.xs<-NULL
	addition.values<-0
	xlist<-NULL
	xs<-NULL
	for(i in 1:length(scaff.ord)){
		all.scaff[[scaff.ord[i]]]<-
			all.scaff[[scaff.ord[i]]][order(all.scaff[[scaff.ord[i]]][,bp.name]),]	
		all.scaff[[scaff.ord[i]]][,bp.name]<-
			seq(last.max+1,last.max+nrow(all.scaff[[scaff.ord[i]]]),1)
		xs<-c(xs, seq(last.max+1,last.max+nrow(all.scaff[[scaff.ord[i]]]),1))
		new.max<-max(xs)
		#scaffold.order[i,"new_start"]<-last.max
		#scaffold.order[i,"new_end"]<-new.max
		rect.xs<-rbind(rect.xs,c(last.max, new.max))
		rownames(rect.xs)[i]<-scaff.ord[i]
		addition.values<-c(addition.values, new.max)
		last.max<-new.max
	}
	#change BP to plot
	x.max<-max(xs)
	x.min<-min(xs)
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
		ylim=y.lim, col=pt.col,
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
			text(x=mean(all.scaff[[scaff.ord[i]]][
				all.scaff[[scaff.ord[i]]]$Chrom==rownames(rect.xs)[i],
				bp.name]),
				y=displacement,labels=rownames(rect.xs)[i],
				adj=1,xpd=T,srt=45)
		}
	}
	for(i in 1:length(scaff.ord)){
		points(all.scaff[[scaff.ord[i]]][,bp.name], 
			all.scaff[[scaff.ord[i]]][,fst.name], 
			pch=19, cex=0.5,col=pt.col,
			xlim=c(x.min,x.max),ylim=y.lim)
	  if(!is.null(ci.dat)){
  		temp.sig<-all.scaff[[scaff.ord[i]]][all.scaff[[scaff.ord[i]]][,fst.name] >= ci.dat[1],]
  		points(temp.sig[,bp.name], temp.sig[,fst.name], 
  			col=sig.col[1], pch=19, cex=0.5)
  		temp.sig<-all.scaff[[scaff.ord[i]]][all.scaff[[scaff.ord[i]]][,fst.name] <= ci.dat[2],]
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


parse.vcf<-function(filename){
  vcf<-read.delim(filename,comment.char="#",sep='\t',header=F,stringsAsFactors = F)
  header.start<-grep("#CHROM",scan(filename,what="character"))
  header<-scan(filename,what="character")[header.start:(header.start+ncol(vcf)-1)]
  colnames(vcf)<-header
  return(vcf)
}

vcf.cov.loc<-function(vcf.row,subset){
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
  return(data.frame(Chrom=vcf.row[1],Pos=vcf.row["POS"],Locus=vcf.row["ID"],
                    NumMissing=miss, NumPresent=pres,PropMissing=miss/(miss+pres),
                    AvgCovRef=ref,AvgCovAlt=alt, AvgCovRatio=ref/alt,AvgCovTotal=tot/pres, CovVariance=var.cov,
                    NumHet=het,PropHet=het/pres,TotalNumReads = tot,stringsAsFactors = F))
}

vcf.cov.ind<-function(vcf.col){
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
}

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
      al1<-vcf.alleles(vcf1.row)
      al2<-vcf.alleles(vcf2.row)
      if((length(al1)/2)/(length(vcf1.row)-10)>=cov.thresh & (length(al2)/2)/(ncol(vcf2.row)-10)>=cov.thresh){
        #calculate frequencies
        fst<-calc.fst.wright(al1,al2)
      }else {
        # print(paste(vcf.row["#CHROM"],vcf.row["POS"],"fails cov thresh"),sep=" ")
        fst<-data.frame(Hs1=NA,Hs2=NA,Hs=NA,Ht=NA,Fst=NA,NumAlleles=length(factor(c(al1,al2))),
                        Num1=length(al1),Num2=length(al2)) #it doesn't pass the coverage threshold
      }
    }#good to go
  }
  return(data.frame(Chrom=vcf1.row["#CHROM"],Pos=vcf1.row["POS"],
                    Hs1=fst["Hs1"],Hs2=fst["Hs2"],Hs=fst["Hs"],Ht=fst["Ht"],Fst=as.numeric(fst["Fst"]),NumAlleles=fst["NumAlleles"],
                    Num1=fst["Num1"],Num2=fst["Num2"],stringsAsFactors=FALSE))
}#end function

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

fst.one.vcf<-function(vcf.row,group1,group2, cov.thresh=0.2, maf=0.05){
  al1<-vcf.alleles(vcf.row[group1])
  al2<-vcf.alleles(vcf.row[group2])
  f1<-table(al1)/sum(table(al1))
  f2<-table(al2)/sum(table(al2))
  if(((length(al2)/2)/(length(group2)-10))>=cov.thresh & ((length(al1)/2)/(length(group1)-10))>=cov.thresh & min(f1,f2)>=maf){
        fst<-calc.fst.wright(al1,al2)
  }else {
    # print(paste(vcf.row["#CHROM"],vcf.row["POS"],"fails cov thresh"),sep=" ")
    fst<-data.frame(Hs1=NA,Hs2=NA,Hs=NA,Ht=NA,Fst=NA,NumAlleles=length(summary(factor(c(al1,al2)))),
                    Num1=length(al1),Num2=length(al2)) #it doesn't pass the coverage threshold
  }
  
  return(data.frame(Chrom=vcf.row[1],Pos=vcf.row[2],
                    Hs1=fst["Hs1"],Hs2=fst["Hs2"],Hs=fst["Hs"],Ht=fst["Ht"],Fst=as.numeric(fst["Fst"]),NumAlleles=fst["NumAlleles"],
                    Num1=fst["Num1"],Num2=fst["Num2"],stringsAsFactors=FALSE))
  
  return(fst)
}
#fsts.both<-do.call("rbind",apply(both.sub,1,fst.one.vcf,group1=c(locus.info,o.ind),group2=c(locus.info,d.ind),cov.thresh=0.5))

calc.afs.vcf<-function(vcf.row){
  #use in conjunction with apply
  #e.g. apply(vcf,1,afs.vcf)
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
  
}


choose.one.snp<-function(vcf){
  keep.col<-colnames(vcf)
  vcf$id.pos<-paste(vcf$ID,vcf$POS,sep=".")
  sub.vcf<-tapply(vcf$id.pos,vcf$ID, sample,size=1)
  new.vcf<-vcf[vcf$id.pos %in% sub.vcf,keep.col]
  return(new.vcf)
}



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
      return(c(y[1:9],mom.allele))
    })
    
  })
}

merge.vcfs<-function(vcf1,vcf2, vcf.name="merge.vcf"){
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

gwsca<-function(vcf,locus.info,group1,group2,prop.ind.thresh=0.5,maf.cutoff=0.05){
  sel<-do.call("rbind",apply(vcf,1,fst.one.vcf,c(locus.info,group1),c(locus.info,group2), 
    cov.thresh=prop.ind.thresh,maf=maf.cutoff))
  sel<-sel[!is.na(sel$Fst),]
  sel$Chi<-2*((sel$Num1+sel$Num2)/2)*sel$Fst
  sel$Chi.p<-1-pchisq(sel$Chi,1)
  sel$Chi.p.adj<-p.adjust(sel$Chi.p,method="BH")
  return(sel)
}

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
        hs<-mean(c(hs1,hs2))
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

#***************************************************************************#
#CALCULATE ISOLATION BY DISTANCE PER LOCUS
#***************************************************************************#
fst.ibd.byloc<-function(ped.file,dist.mat,pop.order){
  results.mantel<-data.frame()
  for(i in seq(7,ncol(ped.file),2)){
    res<-mantel.rtest(
      as.dist(t(pairwise.fst(ped.file,i,i+1,pop.order))),
      as.dist(t(dist.mat)), nrepet=9999)
    results.mantel<-rbind(results.mantel,cbind(res$obs,res$pvalue))
  }
  results.mantel<-as.data.frame(results.mantel)
  colnames(results.mantel)<-c("Obs","P")
  return(results.mantel)
}

#***************************************************************************#
#CALCULATE PAIRWISE PST BETWEEN POPULATION PAIRS
#***************************************************************************#
pairwise.pst<-function(dat, pop.order){
  #first column must be pop id/grouping factor
  library(nlme)
  dat.split<-split(dat, factor(dat[,1]))
  dat.var<-as.data.frame(setNames(
    replicate(length(pop.order),numeric(0), simplify = F), pop.order))
  for(i in 1:(length(pop.order)-1)){
    for(j in (i+1):length(pop.order)){
      temp.data<-rbind(as.data.frame(dat.split[[pop.order[i]]]),
                       as.data.frame(dat.split[[pop.order[j]]]))
      colnames(temp.data)<-c("PopID","Var")
      temp.data$PopID<-factor(temp.data$PopID)
      anv <- lme(fixed=Var ~ 1, random=~1|PopID,data=temp.data)
      varcomp <- VarCorr(anv)
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

#***************************************************************************#
#CALCULATE PAIRWISE PSTS FOR ALL TRAITS AND TEST FOR IBD
#***************************************************************************#
all.traits.pst.mantel<-function(trait.df,comp.df,id.index){
  results.mantel<-data.frame()
  for(i in 3:ncol(trait.df)){
    res<-mantel.rtest(
      as.dist(t(pairwise.pst(trait.df[,c(id.index,i)],pop.order))),
      as.dist(t(comp.df)), nrepet=9999)
    results.mantel<-rbind(results.mantel,cbind(res$obs,res$pvalue))
  }
  results.mantel<-as.data.frame(results.mantel)
  rownames(results.mantel)<-colnames(trait.df)[3:ncol(trait.df)]
  colnames(results.mantel)<-c("Obs","P")
  return(results.mantel)
}


#***************************************************************************#
#COMPARE FST AND PST PER LOCUS
#***************************************************************************#
fst.pst.byloc<-function(ped.file,trait.df,pop.order,trait.ind){
  results.list<-list()
  for(j in 3:ncol(trait.df)){
    results.mantel<-data.frame()
    for(i in seq(7,ncol(ped.file),2)){
      res<-mantel.rtest(
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

#***************************************************************************#
#CALCULATE STANDARD ERROR OF THE MEAN
#***************************************************************************#
sem<-function(x){
  sem<-sd(x)/sqrt(length(x))
  return(sem)
}

#***************************************************************************#
#PLOT ANY GENOME-WIDE STATISTIC
#***************************************************************************#
plotting.genome.wide<-function(bp,var,y.max,x.max, rect.xs=NULL,y.min=0,x.min=0, 
                               plot.new=FALSE, plot.axis=TRUE, rect.color="white",plot.rect=TRUE, 
                               pt.cex=1, pt.col="black"){
  #********************************************
  #this function plots a variable without scaffold info. 
  #feed it the basepair (x) values and variable (y) values 
  #*********************************************
  if(plot.new==TRUE){ par(new=new) }
  plot(bp, var,xlab="",ylab="", new=plot.new,
       type="n", bg="transparent", axes=F, bty="n", 
       xlim=c(x.min,x.max),ylim=c(y.min, y.max))
  if(plot.rect==TRUE){
    num.rect<-nrow(rect.xs)
    if(is.null(num.rect)) {
      rect(rect.xs[1],y.min,rect.xs[2],y.max, 
           col=rect.color, border=NA)
    } else {
      for(i in 1:nrow(rect.xs)){
        rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
             col=rect.color, border=NA)
      }
    }
  }
  if(plot.axis){
    axis(2, at = seq(y.min,y.max,round((y.max-y.min)/2, digits=2)),
         ylim = c(y.min, y.max), pos=0,
         las=1,tck = -0.01, xlab="", ylab="", cex.axis=0.75)}
  points(bp, var, pch=19, cex=pt.cex,col=pt.col,
         xlim=c(x.min,x.max),ylim=c(y.min, y.max))
}

#***************************************************************************#
##AUTOMATED FST PLOTTING
#***************************************************************************#
plotting.fsts.scaffs<-function(fst.dat, fst.name="Fst",chrom.name="Chr",
                               bp.name="BP", pt.lty=0,pt.col="grey7",new=T,
                               ci.dat=NULL, pt.cex=1,y.lab=NULL,axis.size=0.5,scaffold.order=NULL,
                               groups=NULL,print.names=FALSE,pt.pch=19,
                               sig.col="dark green", col.pt.pch=8, col.pt.cex=2){
  if(!is.null(scaffold.order)){
    scaff.ord<-scaffold.order$component_id
    lgs<-scaffold.order$object
  } else{
    scaff.ord<-levels(factor(fst.dat[,chrom.name]))
    lgs<-scaff.ord
  }
  if(!is.null(groups)){
    lgs<-groups
    scaff.ord<-groups
  }
  all.scaff<-split(fst.dat, factor(fst.dat[,chrom.name]))
  last.max<-0
  rect.xs<-NULL
  addition.values<-0
  xlist<-NULL
  xs<-NULL
  for(i in 1:length(scaff.ord)){
    all.scaff[[scaff.ord[i]]]<-
      all.scaff[[scaff.ord[i]]][order(all.scaff[[scaff.ord[i]]][,bp.name]),]	
    all.scaff[[scaff.ord[i]]][,bp.name]<-
      seq(last.max+1,last.max+nrow(all.scaff[[scaff.ord[i]]]),1)
    xs<-c(xs, seq(last.max+1,last.max+nrow(all.scaff[[scaff.ord[i]]]),1))
    new.max<-max(xs)
    #scaffold.order[i,"new_start"]<-last.max
    #scaffold.order[i,"new_end"]<-new.max
    rect.xs<-rbind(rect.xs,c(last.max, new.max))
    rownames(rect.xs)[i]<-scaff.ord[i]
    addition.values<-c(addition.values, new.max)
    last.max<-new.max
  }
  #change BP to plot
  x.max<-max(xs)
  x.min<-min(xs)
  y.max<-max(fst.dat[,fst.name])+0.1*max(fst.dat[,fst.name])
  y.min<-min(fst.dat[,fst.name])-0.1*min(fst.dat[,fst.name])
  if(min(fst.dat[,fst.name]) < 0) {
    y.min<-min(fst.dat[,fst.name]) - 0.1*min(fst.dat[,fst.name])
  } else {
    y.min<-0
  }
  displacement<-y.min-((y.max-y.min)/30)
  if(new==T){
    plot(c(x.min,x.max),c(y.min,y.max),xlim=c(x.min,x.max), 
         ylim=c(y.min, y.max), 
         bty="n",type="n",	axes=F, xlab="", ylab="")
    for(i in 1:nrow(rect.xs)){
      if(i%%2 == 0) {
        rect.color<-"white"
      } else {
        rect.color<-"gray75"
      }
      rect(rect.xs[i,1],y.min,rect.xs[i,2],y.max, 
           col=rect.color, border=NA)
      if(print.names==T){
        text(x=mean(all.scaff[[scaff.ord[i]]][
          all.scaff[[scaff.ord[i]]]$Chrom==rownames(rect.xs)[i],
          bp.name]),
          y=displacement,labels=rownames(rect.xs)[i],
          adj=1,xpd=T,srt=45)
      }
    }
  }
  for(i in 1:length(scaff.ord)){
    if(pt.lty==0){
      points(all.scaff[[scaff.ord[i]]][,bp.name], 
             all.scaff[[scaff.ord[i]]][,fst.name], 
             pch=pt.pch, cex=0.5,col=pt.col,
             xlim=c(x.min,x.max),ylim=c(y.min, y.max))
    } else {
      lines(all.scaff[[scaff.ord[i]]][,bp.name], 
            all.scaff[[scaff.ord[i]]][,fst.name], 
            lty=pt.lty, cex=0.5,col=pt.col,
            xlim=c(x.min,x.max),ylim=c(y.min, y.max))
      
    }
    temp.sig<-all.scaff[[scaff.ord[i]]][all.scaff[[scaff.ord[i]]][,fst.name] >= ci.dat[1],]
    points(temp.sig[,bp.name], temp.sig[,fst.name], 
           col=sig.col[1], pch=col.pt.pch, cex=col.pt.pch)
    temp.sig<-all.scaff[[scaff.ord[i]]][all.scaff[[scaff.ord[i]]][,fst.name] <= ci.dat[2],]
    points(temp.sig[,bp.name], temp.sig[,fst.name], 
           col=sig.col[2], pch=col.pt.pch, cex=col.pt.pch)
  }
  if(axis.size>0){
    axis(2, at = seq(round(y.min,2),round(y.max,2),
                     round((y.max-y.min)/2, digits=2)),
         ylim = c(y.min, y.max), pos=0,
         labels=seq(round(y.min,2),round(y.max,2),
                    round((y.max-y.min)/2, digits=2)),
         las=1,tck = -0.01, xlab="", ylab="", cex.axis=axis.size)
  }
  xes<-do.call("rbind",all.scaff)
  return(xes)
}

#***************************************************************************#
##REORDER A DATAFRAME
#***************************************************************************#
reorder.df<-function(dat,order.list){
  #dat has to have the grouping IDs in row 1
  #those grouping ids must match the factors in order.list
  dat.sep<-split(dat, dat[,1])
  dat.new<-dat.sep[[order.list[1]]]
  for(i in 2:length(order.list)){
    dat.new<-rbind(dat.new, dat.sep[[order.list[i]]])
  }
  return(dat.new)
}




#***************************************************************************#
#PLOT A STRUCTURE BARPLOT
#***************************************************************************#

plotting.structure<-function(structure.out, k, pop.order, 
                             filename=paste("str.k",k,".jpeg",sep=""),make.file=TRUE,
                             plot.new=TRUE,colors=NULL,xlabel=TRUE,ylabel=NULL){
  str.split<-split(structure.out,structure.out[,1])
  if(is.null(colors)){
    bar.colors<-rainbow(k,s=0.5)
  } else {
    bar.colors<-colors
  }
  if(make.file==TRUE){
    jpeg(filename,width=7, height=1.25, units="in", res=300)
    par(mfrow=c(1,length(str.split)))
  } 
  #par(mfrow=c(1,length(str.split)),mar=c(1,0,0,0), oma=c(1,0,0,0),cex=0.5)
  for(i in 1:length(str.split)){
    pop.index<-pop.order[i]
    barplot(height=as.matrix(t(str.split[[pop.index]][,-1])),
            beside=FALSE, space=0,	border=NA, col=bar.colors,
            xlab="", ylab="", xaxt='n', yaxt='n')#, new=plot.new)
    if(xlabel==TRUE){
      mtext(pop.index, 1, line=0.5, cex=1, outer=F)}
    if(!is.null(ylabel)){
      if(i == 1) { mtext(ylabel,2,cex=1) }
    }	
  }
  if(make.file==TRUE) {dev.off()}
}



