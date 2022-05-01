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
output_dir <- "./results/tables/global_opstats_tables/"

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

temp_colnames <- c("method","delta_beta","noise","rep",
                   "cg_fdr_cutoff","region_cutoff",
                   "P", "N", "Called_P", "Called_N", ## positive count, negative count
                   "TP_called", "TP_simul",
                   "TN_called" ,"TN_simul",
                   "CG_PROP_DIFF"
)
OP_STATS <- list()
OP_STATS[["feature"]] <- data.frame(matrix(NA, nrow=1,
                                           ncol=length(temp_colnames)) )
OP_STATS[["basepair"]] <- data.frame(matrix(NA, nrow=1,
                                            ncol=length(temp_colnames)) )
OP_STATS[["cpg_probe"]] <- data.frame(matrix(NA, nrow=1,
                                             ncol=length(temp_colnames)) )

colnames(OP_STATS[["feature"]]) <- temp_colnames
colnames(OP_STATS[["basepair"]]) <- temp_colnames
colnames(OP_STATS[["cpg_probe"]]) <- temp_colnames



simul_id <- as.numeric(str_match(names(simul_results)[i],
                                 "simul_set_([0-9]*)")[2])
simul <- simul_constructor_list[[simul_id]]
method_id <- as.numeric(str_match(names(simul_results)[i],
                                  "method_set_([0-9]*)")[2])

called_grs <- simul_results[[i]]$grange
simul_grs <- simul$grange
simul_cg_grs <- GRanges(seqnames=simul$dmr_locs$chr, IRanges(start=simul$dmr_locs$pos, width=1))
### need to set a cutoff to region p-value for TP,FP, FN, etc
called_grs <- called_grs[which(called_grs$pval_region < 0.05)]



for(j in 1:length(OP_STATS)){
  op_stat_type <- names(OP_STATS)
  OP_STATS[[j]]$method <-names(method_set_list)[method_id]
  OP_STATS[[j]]$delta_beta <- simul$pars$delta_beta
  OP_STATS[[j]]$noise <- simul$pars$noise
  OP_STATS[[j]]$rep <- simul$pars$rep
  # OP_STATS[[j]]$cg_order_rand <- simul$pars$cg_order_rand
  OP_STATS[[j]]$cg_fdr_cutoff <- simul$pars$cg_fdr_cutoff
  OP_STATS[[j]]$region_cutoff <- simul$pars$region_cutoff
}

## start : add TP rate to grs features ###
cs_intersect <-  intersect( called_grs, simul_grs)
cs_overlap <-  findOverlaps( called_grs, simul_grs)
called_inverse_grs <- gaps(called_grs)
called_inverse_grs <- called_inverse_grs[which(as.character(strand(called_inverse_grs))=="*" )]
simul_inverse_grs <- gaps(simul_grs)
simul_inverse_grs <- simul_inverse_grs[which(as.character(strand(simul_inverse_grs))=="*" )]
cs_inverse_intersect <-  intersect( called_grs, simul_inverse_grs)
cs_inverse_overlap <-  findOverlaps( called_grs, simul_inverse_grs)
c_inverse_s_inverse_intersect <-  intersect( called_inverse_grs, simul_inverse_grs)
c_inverse_s_inverse_overlap <-  findOverlaps( called_inverse_grs, simul_inverse_grs)

if(length(called_grs)==0){
  for(j in 1:length(OP_STATS)){
    OP_STATS[[j]]$TP_called <- 0
    OP_STATS[[j]]$TP_simul <- 0
    OP_STATS[[j]]$Called_P <- 0
  }
  OP_STATS$feature$P <- length(simul_grs)
  OP_STATS$feature$N <- length(simul_inverse_grs)
  OP_STATS$feature$Called_N <- length(called_inverse_grs)
  OP_STATS$feature$TN_called <- length(called_inverse_grs)
  OP_STATS$feature$TN_simul <- length(simul_inverse_grs)
  OP_STATS$feature$CG_PROP_DIFF <- NA

  OP_STATS$basepair$P <- sum(simul_grs@ranges@width)
  OP_STATS$basepair$N <- sum(simul_inverse_grs@ranges@width)
  OP_STATS$basepair$TN_called <- sum(called_inverse_grs@ranges@width)
  OP_STATS$basepair$TN_simul <- sum(simul_inverse_grs@ranges@width)
  OP_STATS$basepair$CG_PROP_DIFF <- NA

  OP_STATS$cpg_probe$P <- sum( countOverlaps(simul_grs, locs_gr) )
  OP_STATS$cpg_probe$N <- sum( countOverlaps(simul_inverse_grs, locs_gr) )
  OP_STATS$cpg_probe$TN_called <- sum( countOverlaps(called_inverse_grs, locs_gr) )
  OP_STATS$cpg_probe$TN_simul <- sum( countOverlaps(simul_inverse_grs, locs_gr) )
  OP_STATS$cpg_probe$CG_PROP_DIFF <- NA

}else{
  simul_grs$TP <- -1
  for(j in 1:length(simul_grs)){
    if(!is.element(j, cs_overlap@to)){ simul_grs[j]$TP <- 0; next;}
    temp_gr <- intersect(cs_intersect, simul_grs[j])
    simul_grs[j]$TP <- sum(temp_gr@ranges@width) / simul_grs[j]@ranges@width
  }
  simul_inverse_grs$TP <- -1
  for(j in 1:length(simul_inverse_grs)){
    if(!is.element(j, cs_inverse_overlap@to)){ simul_inverse_grs[j]$TP <- 1; next;}
    temp_gr <- intersect(cs_inverse_intersect, simul_inverse_grs[j])
    simul_inverse_grs[j]$TP <- 1-(sum(temp_gr@ranges@width) / simul_inverse_grs[j]@ranges@width)
  }
  called_grs$TP <- -1
  for(j in 1:length(called_grs)){
    if(!is.element(j, cs_overlap@from)){ called_grs[j]$TP <- 0; next;}
    temp_gr <- intersect(cs_intersect, called_grs[j])
    called_grs[j]$TP <- sum(temp_gr@ranges@width) / called_grs[j]@ranges@width
  }
  called_inverse_grs$TP <- -1
  for(j in 1:length(called_inverse_grs)){
    if(!is.element(j, c_inverse_s_inverse_overlap@from)){ called_inverse_grs[j]$TP <- 0; next;}
    temp_gr <- intersect(c_inverse_s_inverse_intersect, called_inverse_grs[j])
    called_inverse_grs[j]$TP <- sum(temp_gr@ranges@width) / called_inverse_grs[j]@ranges@width
  }
  ## end : add TP rate to grs features ###

  ## start: feature operating characteristics
  OP_STATS$feature$P <- length(simul_grs)
  OP_STATS$feature$N <- length(simul_inverse_grs)
  OP_STATS$feature$Called_P <- length(called_grs)
  OP_STATS$feature$Called_N <- length(called_inverse_grs)
  OP_STATS$feature$TP_called <- mean(called_grs$TP)*length(called_grs)
  OP_STATS$feature$TP_simul <- mean(simul_grs$TP)*length(simul_grs)
  OP_STATS$feature$TN_called <- mean(called_inverse_grs$TP)*length(called_inverse_grs)
  OP_STATS$feature$TN_simul <- mean(simul_inverse_grs$TP)*length(simul_inverse_grs)
  OP_STATS$feature$CG_PROP_DIFF <- mean(countOverlaps(called_grs,simul_cg_grs)/countOverlaps(called_grs,locs_gr))
  ## end: feature operating characteristics
  ## start: basepair operating characteristics
  OP_STATS$basepair$P <- sum(simul_grs@ranges@width)
  OP_STATS$basepair$N <- sum(simul_inverse_grs@ranges@width)
  OP_STATS$basepair$Called_P <- sum(called_grs@ranges@width)
  OP_STATS$basepair$Called_N <- sum(called_inverse_grs@ranges@width)
  OP_STATS$basepair$TP_called <- sum(called_grs$TP * called_grs@ranges@width)
  OP_STATS$basepair$TP_simul <- sum(simul_grs$TP * simul_grs@ranges@width)
  OP_STATS$basepair$TN_called <- sum(simul_inverse_grs@ranges@width) - sum(cs_inverse_intersect@ranges@width)
  OP_STATS$basepair$TN_simul <- sum(simul_inverse_grs@ranges@width) - sum(cs_inverse_intersect@ranges@width)
  OP_STATS$basepair$CG_PROP_DIFF <- sum(countOverlaps(called_grs,simul_cg_grs))/sum(countOverlaps(called_grs,locs_gr))
  ## end: basepair operating characteristics
  ## start: cpg_probe operating characteristics
  OP_STATS$cpg_probe$P <- sum( countOverlaps(simul_grs, locs_gr) )
  OP_STATS$cpg_probe$N <- sum( countOverlaps(simul_inverse_grs, locs_gr) )
  OP_STATS$cpg_probe$Called_P <- sum( countOverlaps(called_grs, locs_gr) )
  OP_STATS$cpg_probe$Called_N <- sum( countOverlaps(called_inverse_grs, locs_gr) )
  OP_STATS$cpg_probe$TP_called <- length(findOverlaps(locs_gr, cs_intersect)@from)
  OP_STATS$cpg_probe$TP_simul <- length(findOverlaps(locs_gr, cs_intersect)@from)
  temp <- length(findOverlaps(locs_gr, cs_inverse_intersect)@from)
  OP_STATS$cpg_probe$TN_called <- (OP_STATS$cpg_probe$N -
                                     length(findOverlaps(locs_gr, cs_inverse_intersect)@from))
  OP_STATS$cpg_probe$TN_simul <- (OP_STATS$cpg_probe$N -
                                    length(findOverlaps(locs_gr, cs_inverse_intersect)@from))
  OP_STATS$cpg_probe$CG_PROP_DIFF <- sum(countOverlaps(called_grs,simul_cg_grs))/sum(countOverlaps(called_grs,locs_gr))
}


OP_STATS
## end: cpg_probe operating characteristics

for(f in names(OP_STATS)){
  out_filename <- str_replace(names(simul_results)[i], "_result.csv","_GOS.csv")
  out_filename <- paste(output_dir,"fdr_",f,"_",out_filename ,sep="")
  write.csv(OP_STATS[[f]], out_filename)
}



