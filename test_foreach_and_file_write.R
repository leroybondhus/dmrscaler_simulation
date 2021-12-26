
  ## install DMRscaler

library("devtools")
library("roxygen2")
install("../DMRscaler", quick=T)

library(doParallel)
registerDoParallel()

results_dir <-paste("./results/")

library("GEOquery")


# Pre-processing
###  Reading of idat files done with minfi library ###

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
```


controls <- grep("[Cc]ontrol",phen$title)
locs <- getLocations(GRset.funnorm)
locs <- data.frame("names"=locs@ranges@NAMES, "pos"=locs@ranges@start, "chr" = rep(locs@seqnames@values, locs@seqnames@lengths))
B <- getBeta(GRset.funnorm)



# Generate Simulation Sets

simul_sc_tbl <- data.frame("min_cgs"=c(3,6,9,12), "size"=c(1e3,1e4,1e5,1e6), "count"=c(50,50,50,50) )

dmr_buffer <- 10
simul_pars <- merge(merge(merge(data.frame("num_samples"=c(8)),
                                data.frame("delta_beta"=c(0.1, 0.4))),
                          data.frame("noise" = c(0,  0.5))),
                    data.frame("rep"=c(1:2))
)
simul_pars <- simul_pars[order(simul_pars$num_samples, simul_pars$delta_beta, simul_pars$noise, simul_pars$rep),]
simul_pars_list <- split(simul_pars, 1:nrow(simul_pars))
for(i in 1:length(simul_pars_list)){
  names(simul_pars_list)[i] <-paste(colnames(simul_pars),simul_pars[i,], sep = "__", collapse = "___")
}
rm(simul_pars)

simul_constructor_list <- foreach(simul_pars = simul_pars_list, .final = function(x) setNames(x, names(simul_pars_list)), .errorhandling = "pass") %dopar% {
  locs$in_dmr <- FALSE
  locs$dmr_size <- NA
  locs$dmr_name <- NA
  for(i in 1:nrow(simul_sc_tbl)){
    size <- simul_sc_tbl[i,]$size
    min_cgs <- simul_sc_tbl[i,]$min_cgs
    count <- simul_sc_tbl[i,]$count
    while(count > 0){
      size_sampled <-sample(floor(size/10)+1:size,1)
      x <- sample(1:nrow(locs),1)
      x_loc <-locs[x,]
      x_which <- which(locs$chr == x_loc$chr & locs$pos >= x_loc$pos & locs$pos <= x_loc$pos+size_sampled )
      x_which_left <- min(x_which)-dmr_buffer
      x_which_right <- max(x_which)+dmr_buffer
      if(x_which_left <  min(which(locs$chr == x_loc$chr)) ){next;} ## test if valid left bound
      if(x_which_right > max(which(locs$chr == x_loc$chr)) ){next;} ## test if valid right bound
      if(length(x_which) >= min_cgs  ## at least min_cgs
         & sum(locs$in_dmr[x_which_left:x_which_right])==0   ## not overlapping a dmr that already exists
         & (locs$pos[max(x_which)]-locs$pos[min(x_which)] >= size/10)  ##  size is between size/10 and size
      ){ ## min number cgs achieved and none in dmr already...
        x_which <- unique(c(min(x_which), sample(x_which, ceiling((1- simul_pars$noise) * length(x_which) )  ) ,max(x_which))) ## sample left and right most cg and introduce noise
        locs$in_dmr[x_which] <- TRUE
        locs$dmr_size[x_which] <- locs$pos[max(x_which)]-locs$pos[min(x_which)]
        locs$dmr_name[x_which] <- paste("dmr", size, count, sep = "_")
        count<-count-1
      }
    }
  }
  simul_list <- list(pars = simul_pars)
  simul_list$dmr_locs <- locs[which(locs$in_dmr),]

  ## identify g1 and g2
  g12 <- sample(colnames(B), 2*simul_pars$num_samples, replace = FALSE)
  simul_list$g1 <- g12[1:simul_pars$num_samples]
  simul_list$g2 <- g12[(simul_pars$num_samples+1):length(g12)]
  simul_list
}

## NOTE: with simul_constructor_list, simulations are deterministic i.e. the random sampling happens within definition of the simul_constructor_list

for(i in 1:length(simul_constructor_list)){
  write.csv(simul_constructor_list[[i]]$g1, paste(results_dir,"test_",i,".csv", sep = ""))
}
