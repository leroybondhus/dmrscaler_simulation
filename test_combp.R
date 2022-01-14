library("devtools")
library("roxygen2")
library("doParallel")
registerDoParallel()


output_dir <-paste("./results/")

load("simul_setup.Rdata")




combp_1 = list(method="combp",
               function_call="system(\"comb-p pipeline -c 4 --dist 1000 --step 5000 \\
       --seed 1e-3 -p ./test_combp_out \\
       --region-filter-p 0.1 ./results/test_combp_input.bed \")")

method_set_result <- eval(parse(text=combp_1$function_call))


