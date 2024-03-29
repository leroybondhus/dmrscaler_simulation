library("devtools")
library(roxygen2)

library(minfi)
library(doParallel)
library(rlang)
library(MESS)

library("valr")
library(IRanges)
registerDoParallel()

source("3rd_draft_islands.R")
promoterInfo<-read.csv("/home/lbondhus/Desktop/PROJECTS/dmrscaler_simulation/intermediate_data/annotations_short.csv",row.names = 1)

document("/home/lbondhus/Desktop/PROJECTS/dmrscaler")
install("/home/lbondhus/Desktop/PROJECTS/dmrscaler")

results_dir <-paste("/home/lbondhus/Desktop/PROJECTS/dmrscaler_simulation/results/")


controlBetaMatrix<-read.csv("/home/lbondhus/Desktop/PROJECTS/dmrscaler_simulation/intermediate_data/Control_Beta_Matrix.csv",row.names = 1)
controlCGlocs<-read.csv("/home/lbondhus/Desktop/PROJECTS/dmrscaler_simulation/intermediate_data/Control_CG_locs.csv",row.names = 1)
controlPhenoData<-read.csv("/home/lbondhus/Desktop/PROJECTS/dmrscaler_simulation/intermediate_data/Control_Pheno_Data.csv",row.names=1)
#subset the control beta matrix for cpgs from chr5
controlCGlocs<-controlCGlocs[controlCGlocs$V3=="chr5",]
controlBetaMatrix<-controlBetaMatrix[rownames(controlBetaMatrix)%in%controlCGlocs$V1,]

promoterInfo <- promoterInfo[which(is.element(promoterInfo$Name, rownames(controlBetaMatrix))),]

############################      Set up data      ############################
all_results <- list() 
pars <- list( delta_beta = c(0.05 , 0.1, 0.15, 0.2, 0.25, 0.3) 
              , noise = c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
              , num_samples = c(5, 6, 7, 8, 10, 12, 15)
#             , dmr_count = c(10, 50, 100, 200)
             )


for(name in names(pars)){
  par_var <- pars[[name]]
  try(
  for(element in 1:length(par_var)){
    ## defaults first
    n_promoters <- 50 ## how many DMRs to introduce
    delta_beta <- 0.1
    noise <- 0.1
    num_cases <- 5
    num_controls <- num_cases
    if(name == "delta_beta"){ delta_beta <- par_var[element] }
    if(name == "noise"){ noise <- par_var[element] }
    if(name == "num_samples"){
      num_cases <- par_var[element]
      num_controls <- par_var[element]
    }
    if(name == "dmr_count"){ n_promoters <- par_var[element] }    
    
    nrows = 5  ## number to run for each setting 
    results_stats <- data.frame( seq(1,nrows), 
                                 delta_beta=rep(delta_beta, nrows),
                                 noise=rep(noise, nrows),
                                 num_samples=rep(num_cases, nrows),
                                 dmr_count=rep(n_promoters, nrows),
                                 auc=rep(-1, nrows),
#                                 FP=numeric(length = nrows),
#                                 FN=numeric(length = nrows),
                                 jaccard=rep(-1, nrows))
  try(
    for(test_num in 1:nrows) {
      
      ### Angela's Code
      run_prom <- FALSE  ## FALSE means run rand
      twogroups <- randomizeTwoGroups(controlBetaMatrix)
      
      if(run_prom){
        randprom <- randomizeSelectRegions(x=promoterInfo, n=n_promoters)
        randbed <- bedForRegions(randprom)
        randcgs <- cpgTableForRegions(randprom)
        randcgs_dropped <- randcgs[sample(1:nrow(randcgs), size = floor(noise*nrow(randcgs))),]
        randcgs  <- randcgs[setdiff(randcgs$Name,randcgs_dropped$Name),]
      } else{
        
          # fixed width
        if(FALSE){ # fixed cg count
          a <- sort(sample(1:(length(promoterInfo$Name)-50), 50))
          b <- a + (1:50)   ## define dmrs of length 1+[1, 2, 3, 4, 5, ..., 50]
          rand_dmrs <- IRanges(start=a, end=b, names = rep("chr5", length(a)) )
          rand_dmrs <- reduce(rand_dmrs)
          dist_to_next <- rand_dmrs@start[2:length(rand_dmrs)] - (rand_dmrs@start[1:(length(rand_dmrs)-1)]+rand_dmrs@width[1:(length(rand_dmrs)-1)])
          which(dist_to_next < 50)
          rand_dmrs <- rand_dmrs[-which(dist_to_next < 50)]   ## drop dmrs that are too close to one another

          randcgs <- data.frame(matrix(ncol = 4, nrow = 0))
          colnames(randcgs) <- c("Name", "chr", "pos", "Islands_Name")
          ### give things same names as in angela's code...
          for(dm_index in 1:length(rand_dmrs)){
            temp <- promoterInfo[ start(rand_dmrs)[dm_index]:end(rand_dmrs)[dm_index], c( "Name", "chr", "pos") ]
            temp$Islands_Name <- as.character(rep(dm_index, nrow(temp)))
           # temp$Islands_Name <- paste("name_", temp$Islands_Name, sep = "")
            randcgs <- rbind(randcgs, temp)
          } 
          randbed <- data.frame(matrix(ncol = 4, nrow = 0))
          randbed <- data.frame(chr=rep("chr5", length(rand_dmrs)), start = promoterInfo[start(rand_dmrs), "pos"], stop = promoterInfo[end(rand_dmrs), "pos"])
          #colnames(randbed) <- c("Islands_Name", chr, start, stop)
        }
          
          if(TRUE){ ## fixed length
            width <- 1e5
            a <- sort(sample(promoterInfo$pos, 50))
            rand_dmrs <- IRanges(start=a, end=a+width-1, names = rep("chr5", length(a)) )
            rand_dmrs <- reduce(rand_dmrs)
            dist_to_next <- rand_dmrs@start[2:length(rand_dmrs)] - (rand_dmrs@start[1:(length(rand_dmrs)-1)]+rand_dmrs@width[1:(length(rand_dmrs)-1)])
            
            rand_dmrs_drop <- rep(FALSE,length(rand_dmrs))
            ### drop dmrs if <3 cgs, if width is changed by overlapping others, if distance between dmrs is not at least 1 dmr width
            for(dm_index in 1:length(rand_dmrs)){
              which <- which((promoterInfo$pos >= start(rand_dmrs)[dm_index]) & (promoterInfo$pos <= end(rand_dmrs)[dm_index]))
              if( (length(which) < 3) | (rand_dmrs@width[dm_index] != width) | (dm_index > length(dist_to_next)) | (dist_to_next[dm_index] < width) ){rand_dmrs_drop[dm_index] <- TRUE} # require 3 CGs for dmr
            } 
            rand_dmrs<-rand_dmrs[-which(rand_dmrs_drop),]
            
            randcgs <- data.frame(matrix(ncol = 4, nrow = 0))
            colnames(randcgs) <- c("Name", "chr", "pos", "Islands_Name")            
            for(dm_index in 1:length(rand_dmrs)){
              which <- which((promoterInfo$pos >= start(rand_dmrs)[dm_index]) & (promoterInfo$pos <= end(rand_dmrs)[dm_index]))
              temp <- promoterInfo[ which, c( "Name", "chr", "pos") ]
              temp$Islands_Name <- as.character(rep(dm_index, nrow(temp)))
              # temp$Islands_Name <- paste("name_", temp$Islands_Name, sep = "")
              randcgs <- rbind(randcgs, temp)
            }
            randcgs_dropped <- randcgs[sample(1:nrow(randcgs), size = floor(noise*nrow(randcgs))),]
            randcgs  <- randcgs[setdiff(randcgs$Name,randcgs_dropped$Name),]
            
            randbed <- data.frame(matrix(ncol = 4, nrow = 0))
            randbed <- data.frame(chr=rep("chr5", length(rand_dmrs)), start = start(rand_dmrs) , stop = end(rand_dmrs))
            
                        
          }
      }
      
      
      
      
      simulatedbeta <- alterByMu(mu=delta_beta, betaMatrix = controlBetaMatrix, expDesign = twogroups, cpgTable = randcgs)
      
      ### quick visual test
      temp <- rowMeans(simulatedbeta[,twogroups[twogroups$group=="A","sampleID"]])-rowMeans(simulatedbeta[,twogroups[twogroups$group=="Aprime","sampleID"]])
      hist(temp, breaks = 50)
      temp <- controlBetaMatrix-simulatedbeta
      hist(as.vector(as.matrix(temp)), breaks = 50)
      ###
      
      
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
        signn <- dmrscaler::determine_significant_windows(window_results=nn, indat=data, quants=cltable , quants_significance_cutoff = "0.99" )
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
      
  
      layer_num = 1
      layer_result <- built_layers[[layer_num]]
      
      

      irang_true_pos <- IRanges(start = as.numeric(randbed$start) ,end = as.numeric(randbed$stop), names = randbed$chr)
     # irang_true_pos <- IRanges(start = randbed_true_pos$start ,end = randbed_true_pos$end, names = randbed_true_pos$chrom)
      irang_layers <- IRanges(start = layer_result$start_pos ,end = layer_result$stop_pos, names = layer_result$chr)
      
      found_dms <- overlapsAny(irang_layers, irang_true_pos)  ## FP are FALSE here
      layer_result$true_pos <- found_dms
      true_dms <- overlapsAny(irang_true_pos, irang_layers)  ## FN are FALSE here
      randbed$found <- true_dms

      layer_res_bed <- layer_result[which(layer_result$true_pos),1:3]
      colnames(layer_res_bed) <- c("chrom","start", "end")
      layer_res_bed <- tbl_interval(layer_res_bed)
      layer_res_bed$chrom <- as.character(layer_res_bed$chrom)
      
      randbed_true_pos <- randbed[which(randbed$found),c("chr","start","stop")]
      colnames(randbed_true_pos) <- c("chrom","start", "end")
      
      randbed_true_pos <- tbl_interval(randbed_true_pos)
      randbed_true_pos$chrom <- as.character(randbed_true_pos$chrom)
      randbed_true_pos$start <- as.numeric(randbed_true_pos$start)
      randbed_true_pos$end <- as.numeric(randbed_true_pos$end)
      
      results_stats$jaccard[test_num] <- bed_jaccard(x=layer_res_bed, y=randbed_true_pos)$jaccard    
  
      #### PRECISION RECALL MEASUREMENT HERE 
      
      sorted_result <- layer_result[order(-layer_result$unsigned_bin_score),] ## sort descending
      precision_recall <- data.frame(
        precision = rep(-1, nrow(sorted_result)),
        recall = rep(-1, nrow(sorted_result))
      )
      for(i in 1:nrow(sorted_result)){
        ## precision = TP / TP+FP 
        ## recall = TP / TP+FN
        temp_res_bed <- sorted_result[1:i, c(1,2,3)]
        colnames(temp_res_bed) <- c("chrom","start", "end")
        temp_res_iran <- IRanges(start = temp_res_bed$start ,end = temp_res_bed$end, names = temp_res_bed$chrom)
        
        called <- overlapsAny(temp_res_iran, irang_true_pos)  ## FP are FALSE here
        tp <- length(which(called==TRUE))
        tp_p_fp <- length(called)
        precision_recall$precision[i] <- tp / tp_p_fp
        
        all_real <- overlapsAny(irang_true_pos, temp_res_iran)  ## FN are FALSE here
        precision_recall$recall[i] <- length(which(all_real == TRUE)) / length(all_real)
      }
      precision_recall
      plot(precision_recall$recall,precision_recall$precision)
      ##
    

      precision_recall<-rbind(precision_recall, c(0,max(precision_recall$recall)))
      auc <- auc(precision_recall$recall, precision_recall$precision)
      results_stats$auc[test_num] <- auc
      ### AUC 
    }
  )
    result_file <- paste("deltaB", delta_beta,"noise",noise, "num_samples",num_cases, "dmr_count", n_promoters, "layerCGcount", layer_sizes[layer_num] ,  sep = "_")
    result_file <- paste(results_dir, result_file, ".csv", sep = "")
    write.csv( results_stats, file = result_file)
    all_results[[paste(name,element, sep = "_")]] <- results_stats
  }

  )
}
