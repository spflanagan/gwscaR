## ------------------------------------------------------------------------
library(gwscaR)
sessionInfo()
vcf.file<-system.file("extdata", "example.vcf.txt",package = "gwscaR")
print(vcf.file)
#vcf<-parse.vcf(vcf.file)

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

