library("devtools")
library(roxygen2)

library(minfi)
library(doParallel)
library(rlang)

library("valr")
library(IRanges)
registerDoParallel()

source("3rd_draft_islands.R")

project_dir <- "/u/project/arboleda/leroybon/PROJECTS/DNAme_Simulation/"
results_dir <- paste(project_dir, "results/", sep = "")
data_dir <- paste(project_dir, "intermediate_data/", sep = "")


document("/u/project/arboleda/leroybon/PACKAGES/dmrscaler")
install("/u/project/arboleda/leroybon/PACKAGES/dmrscaler")

promoterInfo<-read.csv(paste(data_dir, "annotations_short.csv", sep=""),row.names = 1)
controlBetaMatrix<-read.csv(paste(data_dir, "Control_Beta_Matrix.csv", sep = ""),row.names = 1)
controlCGlocs<-read.csv(paste(data_dir,"Control_CG_locs.csv", sep = ""),row.names = 1)
controlPhenoData<-read.csv(paste(data_dir,"Control_Pheno_Data.csv", sep = ""),row.names=1)
par_table <- read.csv(paste(data_dir,"par_table.csv", sep = ""), header=FALSE)

#subset the control beta matrix for cpgs from chr5
controlCGlocs<-controlCGlocs[controlCGlocs$V3=="chr5",]
controlBetaMatrix<-controlBetaMatrix[rownames(controlBetaMatrix)%in%controlCGlocs$V1,]

############################      Set up data      ############################

#pars <- list(delta_beta = c(0.05 ) ) #, 0.1, 0.15, 0.2, 0.25, 0.3) ) #,
#             noise = c(0, 0.1, 0.2, 0.3, 0.4, 0.5),
#             num_samples = c(5, 6, 7, 8, 10, 12, 15),
#             dmr_count = c(10, 50, 100, 200))

args <- commandArgs(trailingOnly = TRUE)
arg <- as.numeric(args[1])

name <- par_table[arg,1]
par_var <- par_table[arg,2]

#for(name in names(pars)){
 # par_var <- pars[[name]]
 # try(
#for(element in 1:length(par_var)){
  ## defaults first
n_promoters <- 50 ## how many DMRs to introduce
delta_beta <- 0.1
noise <- 0.1
num_cases <- 5
num_controls <- num_cases
if(name == "delta_beta"){ delta_beta <- par_var }
if(name == "noise"){ noise <- par_var }
if(name == "num_samples"){
  num_cases <- par_var
  num_controls <- par_var
}
if(name == "dmr_count"){ n_promoters <- par_var }    

nrows = 10  ## number to run for each setting 
results_stats <- data.frame( seq(1,nrows), 
                             delta_beta=rep(delta_beta, nrows),
                             noise=rep(noise, nrows),
                             num_samples=rep(num_cases, nrows),
                             dmr_count=rep(n_promoters, nrows),
                             FP=numeric(length = nrows),
                             FN=numeric(length = nrows),
                             Jaccard=numeric(length = nrows))
results_stat <- data.frame(delta_beta=delta_beta, 
                           noise=noise, 
                           num_samples=num_cases,
                           dmr_count=n_promoters,
                           meanFP=0, stdevFP=0, meanFN=0, stdevFN=0 )
result_file <- paste("deltaB", delta_beta,"noise",noise, "num_samples",num_cases, "dmr_count", n_promoters, sep = "_")
result_file <- paste("/home/lbondhus/Desktop/PROJECTS/dmrscaler_simulation/", result_file, ".csv", sep = "")

for(test_num in 1:nrows) {
  
  ### Angela's Code
  twogroups <- randomizeTwoGroups(controlBetaMatrix)
  randprom <- randomizeSelectRegions(x=promoterInfo, n=n_promoters)
  randbed <- bedForRegions(randprom)
  randcgs <- cpgTableForRegions(randprom)
  randcgs_dropped <- randcgs[sample(1:nrow(randcgs), size = floor(noise*nrow(randcgs))),]
  randcgs  <- randcgs[setdiff(randcgs$Name,randcgs_dropped$Name),]
  
  simulatedbeta <- alterByMu(mu=delta_beta, betaMatrix = controlBetaMatrix, expDesign = twogroups, cpgTable = randcgs)
  
  # ### quick visual test 
  # temp <- rowMeans(simulatedbeta[,twogroups[twogroups$group=="A","sampleID"]])-rowMeans(simulatedbeta[,twogroups[twogroups$group=="Aprime","sampleID"]])
  # hist(temp, breaks = 50)
  # temp <- controlBetaMatrix-simulatedbeta
  # hist(as.vector(as.matrix(temp)), breaks = 50)
  # ###
  
  
  all_case <- which(is.element(colnames(simulatedbeta), twogroups[twogroups$group=="A","sampleID"]))
  all_control <- which(is.element(colnames(simulatedbeta), twogroups[twogroups$group=="Aprime","sampleID"]))
  
  case_index <- sample(all_case, num_cases)
  control_index <- sample(all_control, num_controls)
  
  
  ############################     run everything         #########################
  num_perm <- 10
  clt_reps <- 5e5
  
  rim <- dmrscaler::generate_rand_index_matrix(num_controls = length(control_index),
                                               num_cases = length(case_index),
                                               num_permutations = num_perm)
  mrp <- dmrscaler::run_MWW_rand_permutation(index_matrix = rim, 
                                             Beta = simulatedbeta,
                                             num_permutations = num_perm)
  mrp <- -log10(mrp)
  mwr <- dmrscaler::run_MWW(control_indices = control_index ,
                            case_indices = case_index,
                            Beta = simulatedbeta)
  mwr <- -log10(mwr)
  fdt <- dmrscaler::write_FDR_table(real_table = mwr,
                                    rand_table = mrp)
  
  fdrscaler <- dmrscaler::get_FDR_scalar(MWW_FDR_table = fdt,
                                         MWW_FDR_threshold = 0.25)
  if(is.na(fdrscaler)){fdrscaler <- 1}
  cltable <- dmrscaler::write_CLT_lookup_table(num_reps = clt_reps ,
                                               data_to_sample = mwr$p_val,
                                               FDR_scaler = fdrscaler,
                                               clt_numCGs = c(2, 5, 10, 25, 50))
  ## data <-  names, chr, pos,  scoring_value (-log10pval)
  data <- controlCGlocs
  colnames(data)<-c("names","pos","chr")
  data$scoring_values <- mwr$p_val
  data$chr <- droplevels(data$chr)
  layer_sizes <- c(2,4,8,16,32,64)
  layers<-list()
  for(i in 1:length(layer_sizes)){
    #print(paste("layer", i, sep="_"))
    nn_step_fraction = 2 
    nn_step_size <- floor(layer_sizes[i] / nn_step_fraction)
    nn <- dmrscaler::n_nearest_window_scoring_func(indat = data, n_nearest = layer_sizes[i], step_size = nn_step_size, FDR = fdrscaler)
    signn <- dmrscaler::determine_significant_windows(window_results=nn, indat=data, quants=cltable , quants_significance_cutoff = "0.99999" )
    signn <- dmrscaler::add_significance(three_window_list =  signn, lookup_table = cltable)
    layers[[i]]<-signn
  }
  names(layers)<-paste("layer", layer_sizes, sep="_")
  for(i in 1:length(layers)){
    temp<-as.data.frame(layers[[i]][[1]])
    # print(length(layers[[i]]))
    if(length(layers[[i]]) <= 1){ 
      layers[[i]]<-temp
      next
    }
    for(j in 2:length(layers[[i]])){
      temp<-rbind(temp, layers[[i]][[j]])
    }
    layers[[i]]<-temp[[1]]
  }
  atomic_layer <- data
  for(i in 1:length(layers)){
    for(k in 1:length(layers[[i]]$start_pos)){
      layers[[i]]$start_index[k]<-which(atomic_layer$pos==layers[[i]]$start_pos[k] & as.character(atomic_layer$chr) == as.character(layers[[i]]$chr[k]))
      layers[[i]]$stop_index[k]<-which(atomic_layer$pos==layers[[i]]$stop_pos[k] & as.character(atomic_layer$chr) == as.character(layers[[i]]$chr[k]))
    }
  }
  built_layers <- list()
  built_layers[[1]] <- in_layer_merge(dmrs = layers[[1]], CG_table = atomic_layer, FDR_scaler = fdrscaler, lookup_table = cltable)
  for(i in 2:length(layers)){
    #print(i)
    built_layers[[i]] <- build_next_layer(prev_layer = built_layers[[i-1]], 
                                          windows_to_layer = layers[[i]], 
                                          CG_table = atomic_layer,
                                          FDR_scaler=fdrscaler,
                                          lookup_table = cltable)
    # print("done with one")
  }
  
  ############### Report statistics
  ## 
  ## jaccard
  
  
  
  layer_result <- built_layers[[6]]
  layer_res_bed <- layer_result[,1:3]
  colnames(layer_res_bed) <- c("chrom","start", "end")
  layer_res_bed <- tbl_interval(layer_res_bed)
  layer_res_bed$chrom <- as.character(layer_res_bed$chrom)
  
  randbed_true_pos <- randbed[,2:4]
  colnames(randbed_true_pos) <- c("chrom","start", "end")
  randbed_true_pos<-tbl_interval(randbed_true_pos)
  randbed_true_pos$chrom <- as.character(randbed_true_pos$chrom)
  randbed_true_pos$start <- as.numeric(randbed_true_pos$start)
  randbed_true_pos$end <- as.numeric(randbed_true_pos$end)
  
  results_stats$Jaccard[test_num] <- bed_jaccard(x=layer_res_bed, y=randbed_true_pos)$jaccard
  
  ###
  
  irang_true_pos <- IRanges(start = randbed_true_pos$start ,end = randbed_true_pos$end, names = randbed_true_pos$chrom)
  irang_layers <- IRanges(start = layer_res_bed$start ,end = layer_res_bed$end, names = layer_res_bed$chrom)
  
  found_dms <- overlapsAny(irang_layers, irang_true_pos)  ## FP are FALSE here
  results_stats$FP[test_num] <- length(which(!found_dms)) / length(found_dms) ## FP proportion
  
  true_dms <- overlapsAny(irang_true_pos, irang_layers)  ## FN are FALSE here
  results_stats$FN[test_num] <- length(which(!true_dms)) / length(true_dms)  # FN rate
  
  ##
  
}
results_stat$meanFP <- mean(results_stats$FP)
results_stat$stdevFP <- sd(results_stats$FP)
results_stat$meanFN <- mean(results_stats$FN)
results_stat$stdevFN <- sd(results_stats$FN)


write.csv(results_stat, result_file)
#}  

