library("devtools")
library("roxygen2")
library("doParallel")
registerDoParallel()


output_dir <-paste("./results/")

load("simul_setup.Rdata")


out_file_prefix <- "./results/test_combp_out"
in_file_name <- "./results/combp_temp_input.bed"
combp_1 = list(method="combp",
               function_call=
                 paste("system(\"comb-p pipeline -c 4 --dist 1000 --step 500 ",
                       " --seed 1e-3 --region-filter-p 0.1 -p ", out_file_prefix,
                       " ",in_file_name,  " \") ", sep="" )
               )

method_set_result <- eval(parse(text=combp_1$function_call))


