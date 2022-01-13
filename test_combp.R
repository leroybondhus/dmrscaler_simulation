library("devtools")
library("roxygen2")
library("DMRscaler")
library("bumphunter")
library("DMRcate")
library("doParallel")
registerDoParallel()


output_dir <-paste("./results/")

load("simul_setup.Rdata")


system("comb-p pipeline -c 4 --dist 1000 --step 5000 \
       --seed 1e-3 -p ./test_combp_out \
       --region-filter-p 0.1 ./test_combp_input.bed ")

