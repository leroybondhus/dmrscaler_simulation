library(IRanges)
library(GenomicRanges)
library(data.table)
library(stringr)
library(foreach)
library(doParallel)
registerDoParallel(detectCores()-1)
library(ggplot2)
library(ggExtra)
library(scales)


load("fdr_simul_setup.Rdata")

## output dir
output_dir <- "./results/tables/opstats_tables/"

##  input filenames
filenames <- list.files(path = "./results/fdr_testing", pattern = "*result.csv", full.names = T)


args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1){
  stop("must supply arguments for SIMUL_RESULT_ID")
}

SIMUL_RESULT_ID <- as.numeric(args[1])

### add granges for locs and simulated dmrs
locs_gr <- GRanges(seqnames=locs$chr, ranges = IRanges(start=locs$pos, width=1))

chr_seqlengths <- numeric(length=length(unique(locs$chr)))
names(chr_seqlengths) <- unique(locs$chr)
for(chr in unique(locs$chr)){
  chr_seqlengths[chr] <- max(locs$pos[which(locs$chr==chr)])+1
}
chr_seqinfo <- Seqinfo(names(chr_seqlengths), seqlengths = chr_seqlengths)

for(i in 1:length(simul_constructor_list)){
  simul <- simul_constructor_list[[i]]
  dmr_ids <- unique(simul$dmr_locs$dmr_name)
  temp <- data.frame(chr=character(length = length(dmr_ids)),
                     start=numeric(length = length(dmr_ids)),
                     end=numeric(length = length(dmr_ids)),
                     names=character(length = length(dmr_ids)))
  for(j in 1:length(dmr_ids)){
    which <- which(simul$dmr_locs$dmr_name==dmr_ids[j])
    temp$chr[j] <- as.character(unique(simul$dmr_locs$chr[which]))
    temp$start[j] <- min(simul$dmr_locs$pos[which])
    temp$end[j] <- max(simul$dmr_locs$pos[which])
    temp$names[j] <- dmr_ids[j]
  }

  simul_constructor_list[[i]]$grange <- GRanges(seqnames = temp$chr,
                                                ranges = IRanges(start=temp$start,
                                                                 end=temp$end,
                                                                 names=temp$names),
                                                strand = rep("*",nrow(temp)),
                                                seqinfo = chr_seqinfo)
}

### add granges for called dmrs

simul_results <- list()

for(simul_id in 1:length(simul_constructor_list)){
  for(method_id in 1:length(method_set_list)){
    filename <- filenames[grep(paste("simul_set_",simul_id,
                                     "__method_set_",method_id,
                                     sep=""), filenames)]
    result <- fread(filename)
    if(method_set_list[[method_id]]$method=="dmrscaler"){
      result <- result[grep("64",result$layer),]
      result$pval_region <- result$pval_region_adj
    }
    simul_results[[basename(filename)]]$grange <- GRanges(
      seqnames = result$chr,
      ranges = IRanges(start = result$start,
                       end = result$stop,
                       pval_region = result$pval_region),
      strand = rep("*",nrow(result)),
      seqinfo = chr_seqinfo
    )
    ### need to set a cutoff to region p-value for TP,FP, FN, etc
    which <- which(simul_results[[basename(filename)]]$grange$pval_region < 0.05)
    simul_results[[basename(filename)]]$grange <- simul_results[[basename(filename)]]$grange[which]
  }
}




i <- SIMUL_RESULT_ID

simul_id <- as.numeric(str_match(names(simul_results)[i],
                                 "simul_set_([0-9]*)")[2])
simul <- simul_constructor_list[[simul_id]]
method_id <- as.numeric(str_match(names(simul_results)[i],
                                  "method_set_([0-9]*)")[2])

called_grs <- simul_results[[i]]$grange
simul_grs <- simul$grange
### need to set a cutoff to region p-value for TP,FP, FN, etc
called_grs <- called_grs[which(called_grs$pval_region < 0.05)]
## order called_grs
called_grs <- called_grs[order(called_grs$pval_region)]
called_grs$TP <- -1 ##
cs_overlap <-  findOverlaps( called_grs, simul_grs)
for(j in 1:length(called_grs)){
  if(!is.element(j, cs_overlap@from)){ called_grs$TP[j] <- 0;  next;}
  temp_intersect <- intersect(simul_grs,called_grs[j] )
  called_grs$TP[j] <- sum(temp_intersect@ranges@width) / called_grs[j]@ranges@width
}

temp_seq <- floor(seq(from=1,to=length(called_grs),
                      length.out = min(length(called_grs),20)))
temp_colnames <- c("method","delta_beta","noise","rep",
                   "Precision", "Recall", "FDR", "Power" )
temp_df <- data.frame(matrix(NA, nrow=length(temp_seq),
                             ncol=length(temp_colnames)) )
colnames(temp_df) <- temp_colnames
temp_df$method <- names(method_set_list)[method_id]
temp_df$delta_beta <- simul$pars$delta_beta
temp_df$noise <- simul$pars$noise
temp_df$rep <- simul$pars$rep

OP_STATS <- list()
OP_STATS[["feature"]] <- temp_df
OP_STATS[["basepair"]] <- temp_df
OP_STATS[["cpg_probe"]] <- temp_df


for(j in 1:length(temp_seq)){
  #print(j)
  called_grs_sub <- called_grs[1:temp_seq[j]]
  cs_overlap <-  findOverlaps( called_grs_sub, simul_grs)
  simul_grs$TP <- 0
  for(k in 1:length(simul_grs)){
    if(!is.element(k, cs_overlap@to)){ simul_grs$TP[k] <- 0;  next;}
    temp_intersect <- intersect(simul_grs[k], called_grs_sub )
    simul_grs$TP[k] <- sum(temp_intersect@ranges@width) / simul_grs[k]@ranges@width
  }

  ## RECORD OP_STATS HERE
  ### feature
  OP_STATS$feature$Precision[j] <- mean(called_grs_sub$TP)
  OP_STATS$feature$Recall[j] <- mean(simul_grs$TP)
  OP_STATS$feature$FDR[j] <- 1 - OP_STATS$feature$Precision[j]
  OP_STATS$feature$Power[j] <- OP_STATS$feature$Recall[j]
  ### basepair
  OP_STATS$basepair$Precision[j] <- sum((width(called_grs_sub) * called_grs_sub$TP)) /
    sum(width(called_grs_sub))
  OP_STATS$basepair$Recall[j] <- sum((width(simul_grs) * simul_grs$TP )) /
    sum(width(simul_grs))
  OP_STATS$basepair$FDR[j] <- 1 - OP_STATS$basepair$Precision[j]
  OP_STATS$basepair$Power[j] <- OP_STATS$basepair$Recall[j]
  ### cpg_probe
  cs_intersect <-  intersect( called_grs_sub, simul_grs)
  OP_STATS$cpg_probe$Precision[j] <- length(intersect(cs_intersect,locs_gr)) /
    length(intersect(called_grs_sub, locs_gr))
  OP_STATS$cpg_probe$Recall[j] <- length(intersect(cs_intersect,locs_gr)) /
    length(intersect(simul_grs, locs_gr))
  OP_STATS$cpg_probe$FDR[j] <- 1 - OP_STATS$cpg_probe$Precision[j]
  OP_STATS$cpg_probe$Power[j] <- OP_STATS$cpg_probe$Recall[j]

}

for(f in names(OP_STATS)){
  out_filename <- str_replace(names(simul_results)[i], "_result.csv","_OS.csv")
  out_filename <- paste(output_dir,"fdr_",f,"_",out_filename ,sep="")
  write.csv(OP_STATS[[f]], out_filename)
}



