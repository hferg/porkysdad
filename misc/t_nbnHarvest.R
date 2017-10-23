setwd("/home/hfg/Documents/projects/pearson_richard/sdm_pipeline_tim_advent/nbn_atlas_harvest")
c4_sp <- read.csv("c4_names.csv")
ref_list <- c4_sp$NBN_Name
dr_nums <- read.csv("nbn_dataset_drnums.csv")

library(porkysdad)

# test
