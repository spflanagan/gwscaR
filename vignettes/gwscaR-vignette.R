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
par(mar=c(2,2,2,0),oma=c(1,1,1,1))
sel.plot<-fst.plot(sel, fst.name="Fst",axis.size=1,
             chrom.name="Chrom",bp.name="Pos")

#plot without the scaffolds
lgs<-seq(1,22)
sel$Chrom <-as.numeric(as.character(sel$Chrom)) #the lgs and chrom have to be the same class
sel.plot<-fst.plot(sel, fst.name="Fst",axis.size = 1,
             chrom.name="Chrom",bp.name="Pos",
             groups=lgs)


## ---- fig.show='hold'----------------------------------------------------
#create a data.frame with the starts and stops for the chromosomes
bounds<-data.frame(levels(as.factor(vcf$`#CHROM`)),tapply(as.numeric(as.character(vcf$POS)),vcf$`#CHROM`,max))

par(mar=c(2,2,2,0),oma=c(1,1,1,1))
sel.plot<-fst.plot(sel, fst.name="Fst",axis.size=1,
             chrom.name="Chrom",bp.name="Pos",
             group.boundaries = bounds)
#highlight points that have adjusted p-values < 0.05
#points(sel.plot$Pos[sel.plot$Chi.p.adj<0.05],sel.plot$Fst[sel.plot$Chi.p.adj<0.05],col="cornflowerblue",pch=19)

#Or if you just want to plot the linkage groups
lgs<-seq(1,22)
sel.plot<-fst.plot(sel, fst.name="Fst",axis.size = 1,
             chrom.name="Chrom",bp.name="Pos",
             group.boundaries = bounds, groups=lgs)


## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)
points(1:10,col="blue")

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

