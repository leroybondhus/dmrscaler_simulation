library("devtools")
library("roxygen2")
install("../DMRscaler", quick=T)

library(doParallel)
registerDoParallel()

output_dir <-paste("./results/")

load("simul_setup.Rdata")


args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2){
  stop("must supply arguments for simulation parameters (SIMUL_SET) and method parameters (METHOD_SET)")
}

SIMUL_SET_INDEX <- args[1]
METHOD_SET_INDEX <- args[2]

simul_set <- simul_constructor_list[[SIMUL_SET_INDEX]]
simul_set_name <- names(simul_constructor_list)[SIMUL_SET_INDEX]

method_set <- method_set_list[[METHOD_SET_INDEX]]
method_set_name <- names(method_set_list)[METHOD_SET_INDEX]



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
  mwr <- DMRscaler::run_MWW(g1,g2,B_mod)
  locs$pval <- mwr$p_val
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
  mwr <- DMRscaler::run_MWW(g1,g2,B_mod)
  locs$pval <- mwr$p_val
  combp_input_bed <- data.frame(chrom=locs$chr,start=locs$pos,end=locs$pos+1,pval=10^-mwr$p_val)
  combp_input_bed <- combp_input_bed[order(combp_input_bed$chrom),]
  colnames(combp_input_bed)[1] <- "#chrom"
  filename <- paste(output_dir,"simul_set_",SIMUL_SET_INDEX,"__method_set_",METHOD_SET_INDEX,"_combp_input.bed",sep="")
  data.table::fwrite(combp_input_bed, file = filename, row.names = F,col.names = T, sep = "\t")
  filename <- paste(output_dir,"simul_set_",SIMUL_SET_INDEX,"__method_set_",METHOD_SET_INDEX,"_combp_output.bed",sep="")
  method_set$function_call <-

} else {
  stop("method_name not found")
}

method_set_result <- eval(parse(text=method_set$function_call))


if(grepl("dmrscaler", method_name, ignore.case = TRUE)){

} else if(grepl("bumphunter", method_name, ignore.case = TRUE)){

} else if(grepl("dmrcate", method_name, ignore.case = TRUE)){

} else if(grepl("comb", method_name, ignore.case = TRUE)){

} else {
  stop("method_name not found")
}

filename <- paste(output_dir,"simul_set_",SIMUL_SET_INDEX,"__method_set_",METHOD_SET_INDEX,"_result.csv" )
write.table()
