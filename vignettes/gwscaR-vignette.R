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
#  head(sel)

## ------------------------------------------------------------------------
head(sel)
colnames(sel)
vcf[1:5,1:20]
locus.info

## ------------------------------------------------------------------------
fst.plot(sel, fst.name="Fst", 
             chrom.name="Chrom", axis.size=0,bp.name="Pos",
             sig.col=c("black","black"))

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

