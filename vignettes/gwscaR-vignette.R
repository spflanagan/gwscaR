## ------------------------------------------------------------------------
library(gwscaR)
vcf.file<-system.file("extdata", "example.vcf.txt",package = "gwscaR")
vcf<-parse.vcf(vcf.file)

## ------------------------------------------------------------------------
vcf$SNP<-paste(vcf$`#CHROM`,vcf$POS,sep=".")
head(vcf$SNP)
locus.info<-c("#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","SNP")

## ------------------------------------------------------------------------
grp1<-grep("FEM",colnames(vcf),value=T)
grp2<-c(grep("PRM",colnames(vcf),value=T),grep("NPM",colnames(vcf),value=T))

## ------------------------------------------------------------------------
sel<-do.call(rbind,apply(vcf,1,fst.one.vcf,group1=c(locus.info,grp1),group2=c(locus.info,grp2),
    cov.thresh=0.5,maf=0.05))
sel<-sel[!is.na(sel$Fst),] #Remove the ones that weren't polymorphic
head(sel)

## ------------------------------------------------------------------------
sel$Chi<-2*((sel$Num1+sel$Num2)/2)*sel$Fst
sel$Chi.p<-1-pchisq(sel$Chi,1)
sel$Chi.p.adj<-p.adjust(sel$Chi.p,method="BH")

## ----eval=FALSE----------------------------------------------------------
#  sel<-gwsca(vcf,locus.info,grp1,grp2,prop.ind.thresh=0.5,maf.cutoff=0.05)

## ----fig.show='hold'-----------------------------------------------------
lgs<-seq(1,22)
par(mar=c(2,2,2,0),oma=c(1,1,1,1))
sel.plot<-fst.plot(sel, fst.name="Fst",axis.size=0.6,
                   xlabels=lgs,xlab.indices = seq(1,22),
                   chrom.name="Chrom",bp.name="Pos")

#plot without the scaffolds
sel$Chrom <-as.numeric(as.character(sel$Chrom)) #the lgs and chrom have to be the same class
sel.plot<-fst.plot(sel, fst.name="Fst",axis.size = 0.6,
             chrom.name="Chrom",bp.name="Pos",xlabels=TRUE,
             scaffs.to.plot=lgs)


## ---- fig.show='hold'----------------------------------------------------
#create a data.frame with the starts and stops for the chromosomes
bounds<-data.frame(levels(as.factor(vcf$`#CHROM`)),tapply(as.numeric(as.character(vcf$POS)),vcf$`#CHROM`,max))

par(mar=c(2,2,2,0),oma=c(1,1,1,1))
sel.plot<-fst.plot(sel, fst.name="Fst",axis.size=0.6,
             chrom.name="Chrom",bp.name="Pos",xlabels=lgs,
             scaffold.widths = bounds)

#Or if you just want to plot the linkage groups
lgs<-seq(1,22)
sel.plot<-fst.plot(sel, fst.name="Fst",axis.size = 0.6,
             chrom.name="Chrom",bp.name="Pos",xlabels=T,
             scaffold.widths = bounds, scaffs.to.plot=lgs)

#And you can highlight the points that have adjusted p-values < 0.05
points(sel.plot$plot.pos[sel.plot$Chi.p.adj<0.05],sel.plot$Fst[sel.plot$Chi.p.adj<0.05],
       col="cornflowerblue",pch=19,cex=0.75) #only one point is significant in this example


## ----eval=F--------------------------------------------------------------
#  ibd.by.loc<-fst.ibd.byloc(sub.ped,dist,pop.list)

## ----eval=F--------------------------------------------------------------
#  fem.psts<-apply(fem.phenotype.data[,3:10],2,function(x){
#  	pst<-pairwise.pst(data.frame(fem.phenotype.data[,1],x),pop.list)
#  	return(pst)
#  })
#  
#  fem.pst.fst.loc<-fst.pst.byloc(sub.ped,fem.phenotype.data,pop.list,1)
#  
#  fem.dist<-.pst.mantel(fem.phenotype.data,dist,1)
#  

