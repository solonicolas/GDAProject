
# let's upload some csv results from our analysis

setwd("/Users/solonicolas/DSSC/Genomic Data Analytics/Giulio/exam")
ttype = "Ovary-AdenoCA"
csv_results_dir = paste0(ttype,"/csv_results/")

metadata = read_csv(paste0(csv_results_dir,"original_metadata.csv"))
cnaqc_res = read_csv(paste0(csv_results_dir,"cnaqc_res.csv"))
ccf_res = read_csv(paste0(csv_results_dir,"ccf_res.csv"))

mobster_res_1_1 = read_csv(paste0(csv_results_dir,"mobster_res_1_1.csv"))
mobster_res_2_0 = read_csv(paste0(csv_results_dir,"mobster_res_2_0.csv"))
mobster_res_2_1 = read_csv(paste0(csv_results_dir,"mobster_res_2_1.csv"))
mobster_res_2_2 = read_csv(paste0(csv_results_dir,"mobster_res_2_2.csv"))
mobster_res_max = read_csv(paste0(csv_results_dir,"mobster_res_max.csv"))

evol_res_1_1 = read_csv(paste0(csv_results_dir,"evolutionary_parameters_1_1.csv"))
evol_res_2_0 = read_csv(paste0(csv_results_dir,"evolutionary_parameters_2_0.csv"))
evol_res_2_1 = read_csv(paste0(csv_results_dir,"evolutionary_parameters_2_1.csv"))
evol_res_2_2 = read_csv(paste0(csv_results_dir,"evolutionary_parameters_2_2.csv"))
evol_res_max = read_csv(paste0(csv_results_dir,"evolutionary_parameters_max.csv"))

