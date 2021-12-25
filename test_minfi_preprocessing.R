library("devtools")
library("roxygen2")
install("../DMRscaler", quick=T)

library(doParallel)
registerDoParallel()


results_dir <-paste("./results/")


library("GEOquery")

gse <- getGEO("GSE74432", GSEMatrix = TRUE)
phen <- gse$GSE74432_series_matrix.txt.gz@phenoData@data
phen <- phen[intersect( grep("[Cc]ontrol", phen$characteristics_ch1.1),grep("whole blood",phen$characteristics_ch1.2) ),]
rm(gse)

## get methylation data as idat files (NOTE: this saves files locally in working directory, unpacked size is 2.01 Gb
if(length(list.files("GSE74432/idat", pattern = "idat$"))==0){
  getGEOSuppFiles("GSE74432")
  untar("GSE74432/GSE74432_RAW.tar", exdir = "GSE74432/idat")
  file.remove("GSE74432/GSE74432_RAW.tar")


  idat_files <- list.files("GSE74432/idat", pattern = "idat.gz$", full = TRUE)
  sapply(idat_files, gunzip, overwrite = TRUE); rm(idat_files)
}

##  Preprocessing
library("minfi")
###  Reading of idat files done with minfi library ###
idats_dir<-"GSE74432/idat"
targets <- data.frame("Basename"= stringr::str_split_fixed(basename(phen$supplementary_file), "_Grn", 2)[,1] )

RGSet <- read.metharray.exp(base = idats_dir, targets = targets)
GRset.funnorm <- preprocessFunnorm(RGSet);rm(RGSet)
snps <- getSnpInfo(object = GRset.funnorm)
GRset.funnorm <- dropLociWithSnps(GRset.funnorm, snps=c("SBE", "CpG"), maf=0);rm(snps)
rm(idats_dir)
