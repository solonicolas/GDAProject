#### Prepare environment ####
require(dplyr)
require(CNAqc)
require(tidyverse)
require(mobster)

# set the working directory and the file path
setwd("/Users/solonicolas/DSSC/Genomic Data Analytics/Giulio/exam")
file_path = "./rds_giasone"

#### Load and prepare the data ####
pcawg_files = list.files(file_path, full.names = TRUE)

# for now let's take only 10 samples
pcawg_files = pcawg_files[1:10]
all_pcawg_rds = lapply(pcawg_files, readRDS)


# save some metadata of the samples in a df
metadata = Reduce(rbind, lapply(all_pcawg_rds,
                       function(x) {
                         x$metadata
                       }))

# write the metadata in a csv file
metadata %>% readr::write_csv("metadata.csv")

# put the value of 1 in the CFF columns if it is made of only NA values
# to avoid the error in the cnaqc package
all_pcawg_rds = lapply(all_pcawg_rds,
                       function(x) {
                         if(all(is.na(x$cna$CCF))==TRUE) { # check if there are all NA
                           x$cna = x$cna %>% 
                             mutate(CCF = 1)
                         }
                         return (x)
                       })


# create CNAqc objects
all_cnaqc = lapply(all_pcawg_rds,
                   function(x) {
                     CNAqc::init(
                       snvs = x$mutations,
                       cna = x$cna,
                       purity = x$metadata$purity,
                       ref = 'hg19'
                     )
                   })

#### Assessing the quality of the data with CNAqc ####

# peak analysis
# it's going to take some minutes (even for only 10 samples)
all_cnaqc = lapply(all_cnaqc, function(x) {CNAqc::analyze_peaks(x)})


# take out some results from the cnaqc object
cnaqc_res = data.frame(matrix(ncol = 4))
colnames(cnaqc_res) = c("sample",
                     "karyotype",
                     "mutation_multiplicity",
                     "cna_QC")

for(i in 1:length(all_cnaqc)){
  x = all_cnaqc[[i]]
  sample = x$snvs$sample[1]
  karyotypes = x$peaks_analysis$matches$karyotype
  mutation_multiplicity = x$peaks_analysis$matches$mutation_multiplicity
  QC = x$peaks_analysis$matches$QC
  
  for(i in 1:length(karyotypes)){
    cnaqc_res[nrow(cnaqc_res)+1,] = c(sample,
                                karyotypes[i],
                                mutation_multiplicity[i],
                                QC[i])
  }
} 

# the first line comes NA, so I remove it manually
cnaqc_res = cnaqc_res[2:dim(cnaqc_res)[1],]

# write the first cnaqc results in a csv files
cnaqc_res %>% readr::write_csv("cnaqc_res.csv")


# CCF analysis
# it's going to take less than before, but still a bit
all_cnaqc = lapply(all_cnaqc, function(x) {CNAqc::compute_CCF(x)})

# take out some results from the cnaqc object
ccf_res = data.frame(matrix(ncol = 3))
colnames(ccf_res) = c("sample",
                     "karyotype",
                     "ccf_QC")

for(i in 1:length(all_cnaqc)){
  
  x = all_cnaqc[[i]]
  sample = x$snvs$sample[1]
 
  for(ccf in x$CCF_estimates){
    ccf_res[nrow(ccf_res)+1,] = c(sample,
                                  ccf$QC_table$karyotype,
                                  ccf$QC_table$QC)
  }
} 

# the first line comes NA, so I remove it manually
ccf_res = ccf_res[2:dim(ccf_res)[1],]

# write the first cnaqc results in a csv files
ccf_res %>% readr::write_csv("ccf_res.csv")



# horrible double for loops sorry!!
# I tried the solution with the apply function and part of it is below. 
# But I miss the last part: putting together everything in a df

# sample <- sapply(all_cnaqc, function(x) {x$snvs$sample[1]})
# purity <- sapply(all_cnaqc, function(x) {x$purity})
# ploidy <- sapply(all_cnaqc, function(x) {x$ploidy})
# 
# mutation_multiplicity <- sapply(all_cnaqc, function(x) {x$peaks_analysis$matches$mutation_multiplicity})
# QC <- sapply(all_cnaqc, function(x) {x$peaks_analysis$matches$QC})
# karyotype <- sapply(all_cnaqc, function(x) {x$peaks_analysis$matches$karyotype})

  
