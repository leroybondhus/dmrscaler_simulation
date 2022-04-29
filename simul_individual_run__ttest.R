library("devtools")
library("roxygen2")
library("DMRscaler")
library("bumphunter")
library("DMRcate")
library("doParallel")
registerDoParallel()


output_dir <-paste("./results/ttest_out/")

load("simul_setup.Rdata")

#################
run_TT <- function (control_indices, case_indices, Beta, mu=0){
  p_val <- vector(mode = "numeric", length = nrow(Beta))
  for (i in 1:nrow(Beta)) {
    temp <- t.test(Beta[i, case_indices], y = Beta[i,
                                                   control_indices], paired = FALSE, mu=mu)
    p_val[i] <- temp$p.value
  }
  tt <- data.frame(p_val)
  rownames(tt) <- rownames(Beta)
  return(tt)
}
#################

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("must supply arguments for simulation parameters (SIMUL_SET) and method parameters (METHOD_SET)")
}

SIMUL_SET_ID <- as.numeric(args[1])
METHOD_SET_ID <- as.numeric(args[2])
# filename <- paste(output_dir,"args_test_",SIMUL_SET_ID,"__method_set_",METHOD_SET_ID,".csv", sep="" )
# TEST <- data.frame("SIMUL_SET_ID"=SIMUL_SET_ID, "METHOD_SET_ID"=METHOD_SET_ID )
# write.table(TEST, filename, row.names = F)


simul_set <- simul_constructor_list[[SIMUL_SET_ID]]
simul_set_name <- names(simul_constructor_list)[SIMUL_SET_ID]

method_set <- method_set_list[[METHOD_SET_ID]]
method_set_name <- names(method_set_list)[METHOD_SET_ID]


if(grepl("bumphunter", method_name, ignore.case = TRUE)){stop("stopping")}
if(grepl("dmrcate", method_name, ignore.case = TRUE)){stop("stopping")}

### Start: Simulate DMRs ####
g12 <- c(simul_set$g1,simul_set$g2)
g1 <- simul_set$g1
g2 <- simul_set$g2
B_mod <- as.matrix(B[,g12])

dmr_names <- unique(simul_constructor_list[[1]]$dmr_locs$dmr_name)
for(i in 1:length(dmr_names)){
  dmr_i <- dmr_names[i]
  locs_i <- simul_set$dmr_locs$names[which(simul_set$dmr_locs$dmr_name==dmr_i)]
  mean_g1_i <- mean(as.matrix(B_mod[locs_i,g1]))
  mean_g2_i <- mean(as.matrix(B_mod[locs_i,g2]))
  if(mean_g1_i > mean_g2_i){
    B_mod[locs_i, g1] <- B_mod[locs_i, g1]+simul_set$pars$delta_beta
  }else{
    B_mod[locs_i, g2] <- B_mod[locs_i, g2]+simul_set$pars$delta_beta
  }
}
B_mod[which(B_mod <= 0)] <- min(B_mod[which(B_mod>0)])
B_mod[which(B_mod >= 1)] <- max(B_mod[which(B_mod<1)])

### End: Simulate DMRs ####





method_name <- method_set$method

if(grepl("dmrscaler", method_name, ignore.case = TRUE)){
  mwr <- run_TT(g1,g2,B_mod)
  locs$pval <- mwr$p_val
  pval_cutoff_1 <- DMRscaler::get_loc_fdr_pval(B_mod, g1,g2, wilcox.test, fdr=0.1)
  pval_cutoff_2 <- DMRscaler::get_loc_fdr_pval(B_mod, g1,g2, wilcox.test, fdr=0.01)
  region_pval_cutoff <- 0.01
} else if(grepl("bumphunter", method_name, ignore.case = TRUE)){
  design <- rep(-1,length(colnames(B_mod)))
  design[which(is.element(colnames(B_mod),g1))] <- 1
  design <- cbind(rep(1,length(colnames(B_mod) ) ), design )
  colnames(design) <- c("(Intercept)","(Intercept)")

} else if(grepl("dmrcate", method_name, ignore.case = TRUE)){
  design <- rep(-1,length(colnames(B_mod)))
  design[which(is.element(colnames(B_mod),g1))] <- 1
  design <- cbind(rep(1,length(colnames(B_mod) ) ), design )
  colnames(design)<- c("(Intercept)","(Intercept)")
  M <- log2(B_mod / (1-(B_mod)) )
  myannotation <- cpg.annotate("array", object=M, what="M", arraytype = "450K", analysis.type = "differential", design = design,  coef = 2)

} else if(grepl("comb", method_name, ignore.case = TRUE)){
  mwr <- run_TT(g1,g2,B_mod)
  locs$pval <- mwr$p_val
  combp_input_bed <- data.frame(chrom=locs$chr,start=locs$pos,end=locs$pos+1,pval=mwr$p_val)
  combp_input_bed <- combp_input_bed[order(combp_input_bed$chrom),]
  colnames(combp_input_bed)[1] <- "chrom"
  combp_input_bed <- combp_input_bed[order(as.character(combp_input_bed$chrom)),]
  combp_temp_file_prefix <- paste(output_dir,"simul_set_",SIMUL_SET_ID,"__method_set_",METHOD_SET_ID,"_combp",sep="")
  filename <- paste(combp_temp_file_prefix,"_input.bed",sep="")
  data.table::fwrite(combp_input_bed, file = filename, row.names = F,col.names = T, sep = "\t")
  combp_out_filename <- paste(combp_temp_file_prefix,".regions-t.bed",sep="")

  method_set$function_call <-  paste(method_set$function_call,
                                     combp_temp_file_prefix,
                                     filename, "\")",  sep=" " )

} else {
  stop("method_name not found")
}

t1 <- Sys.time()
method_set_result <- eval(parse(text=method_set$function_call))
t2 <- Sys.time()
filename <- paste(output_dir,"simul_set_",SIMUL_SET_ID,"__method_set_",METHOD_SET_ID,"_TIME.csv", sep="" )
TIME <- data.frame("run"=basename(filename), "time"=as.numeric(difftime(t2,t1, units="secs")))
write.table(TIME, filename, row.names = F)


## convert all method outputs to standard GRange objects with chr,start,stop,pval
if(grepl("dmrscaler", method_name, ignore.case = TRUE)){
  out_df <- method_set_result[[1]][0,]
  for(i in 1:length(method_set_result)){
    if(nrow(method_set_result[[i]])>0 ){
      out_df <- rbind(out_df, data.frame(method_set_result[[i]],layer=names(method_set_result)[i] ))
    }
  }

} else if(grepl("bumphunter", method_name, ignore.case = TRUE)){
  out_df <- method_set_result$table
  colnames(out_df)[which(colnames(out_df)=="end")] <- "stop"
  colnames(out_df)[which(colnames(out_df)=="p.value")] <- "pval_region"

} else if(grepl("dmrcate", method_name, ignore.case = TRUE)){
  out_df <- data.frame(coord=method_set_result@coord,
                       stouffer=method_set_result@Stouffer,
                       HMFDR=method_set_result@HMFDR,
                       Fisher=method_set_result@Fisher)
  out_df$chr <- stringr::str_split_fixed(out_df$coord,":",2)[,1]
  out_df$start <- stringr::str_split_fixed(stringr::str_split_fixed(out_df$coord,":",2)[,2],"-",2)[,1]
  out_df$stop <- stringr::str_split_fixed(stringr::str_split_fixed(out_df$coord,":",2)[,2],"-",2)[,2]
  out_df$pval_region <- out_df$stouffer

} else if(grepl("comb", method_name, ignore.case = TRUE)){
  out_df <- data.table::fread(combp_out_filename)
  colnames(out_df)[which(colnames(out_df)=="end")] <- "stop"
  colnames(out_df)[which(colnames(out_df)=="#chrom")] <- "chr"
  colnames(out_df)[which(colnames(out_df)=="z_p")] <- "pval_region"
  ## add line removing combp intermediary files
  system(paste("rm ", paste(combp_temp_file_prefix,"*",sep = "") ))

} else {
  stop("method_name not found")
}

filename <- paste(output_dir,"simul_set_",SIMUL_SET_ID,"__method_set_",METHOD_SET_ID,"_result.csv", sep="" )
write.table(out_df, file = filename, row.names = F)
