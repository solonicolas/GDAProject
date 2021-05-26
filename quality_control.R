#### Prepare environment ####
require(dplyr)
require(CNAqc)
require(tidyverse)
require(mobster)

# set the working directory and the file path
setwd("/Users/solonicolas/DSSC/Genomic Data Analytics/Giulio/exam")
file_path = "./rds_giasone"


#### Load the data ####
pcawg_files = list.files(file_path, full.names = TRUE)

# for now let's take only 20 samples
pcawg_files = pcawg_files[1:20]
all_pcawg_rds = lapply(pcawg_files, readRDS)


# save some metadata of the samples in a df
metadata = Reduce(rbind, lapply(all_pcawg_rds,
                                function(x) {
                                  # compute the avg coverage per sample
                                  x$metadata$avg_coverage = mean(x$mutations$DP, na.rm=T)
                                  return (x$metadata)
                                }))

# write the metadata in a csv file
metadata %>% readr::write_csv("metadata.csv")

# index each row with the sample name
names(all_pcawg_rds) = metadata$sample

#### Choose the filters to apply later ####
min_purity = 0.1
min_avg_coverage = 10
#ttypes = c(NA, 'Liver-HCC', 'Prost-AdenoCA') # take some tumor types
ttypes = unique(metadata$ttype) # take all the tumor types
min_karyotype_mut = 10
min_avg_coverage_karyotype = 30


# apply some primary filters on PURITY, COVERAGE and TUMOR TYPE
metadata = metadata %>% 
  filter(purity>min_purity & avg_coverage>min_avg_coverage & ttype %in% ttypes) %>% 
  select(sample,
         purity,
         ploidy,
         avg_coverage,
         ttype,
         mutation_drivers)

all_pcawg_rds = all_pcawg_rds[metadata$sample]


# put the value of 1 in the CCF columns if it is made of only NA values
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

# peak and CCF analysis
# it's going to take some minutes (even for only 20 samples)

all_cnaqc = lapply(all_cnaqc,  
                   function(x) {
                     x = x %>% 
                       CNAqc::analyze_peaks(
                         min_absolute_karyotype_mutations = min_karyotype_mut
                        ) %>% 
                       CNAqc::compute_CCF() 
                     return(x)
                   })


# take out some results from the cnaqc object
cnaqc_res = data.frame(matrix(ncol = 5))
colnames(cnaqc_res) = c("sample",
                        "karyotype",
                        "mutation_multiplicity",
                        "avg_coverage_karyotype",
                        "cna_QC")

ccf_res = data.frame(matrix(ncol = 3))
colnames(ccf_res) = c("sample",
                      "karyotype",
                      "ccf_QC")

for(x in all_cnaqc){
  
  # check if peak and CCF analysis went well 
  # they could be NULL because of the filter on the number of mutations
  if(!is.null(x$peaks_analysis) & !is.null(x$CCF_estimates)) {
    
    sample = x$snvs$sample[1]
    karyotypes = x$peaks_analysis$matches$karyotype
    mutation_multiplicity = x$peaks_analysis$matches$mutation_multiplicity
    QC = x$peaks_analysis$matches$QC
    
    for(i in 1:length(karyotypes)){
      
      avg_coverage_karyotype = x$snvs %>%
        filter(karyotype==karyotypes[i]) %>%
        summarise(mean = mean(DP)) %>% 
        select(mean)
      
      cnaqc_res[nrow(cnaqc_res)+1,] = c(sample,
                                        karyotypes[i],
                                        mutation_multiplicity[i],
                                        avg_coverage_karyotype,
                                        QC[i])
    }
    
    for(ccf in x$CCF_estimates){
      ccf_res[nrow(ccf_res)+1,] = c(sample,
                                    ccf$QC_table$karyotype,
                                    ccf$QC_table$QC)
    }
  }
} 

# the first line comes NA, so I remove it manually
cnaqc_res = cnaqc_res[2:dim(cnaqc_res)[1],]
ccf_res = ccf_res[2:dim(ccf_res)[1],]

# write the cnaqc results in a csv files
cnaqc_res %>% readr::write_csv("cnaqc_res.csv")
ccf_res %>% readr::write_csv("ccf_res.csv")


# apply last filters on the PASS/FAIL status and on the COVERAGE per karyotype

filtered_cnaqc_res = cnaqc_res %>% 
  filter(cna_QC=='PASS' & avg_coverage_karyotype>min_avg_coverage_karyotype)

filtered_ccf_res = ccf_res %>% 
  filter(ccf_QC=='PASS')

filtered_results = full_join(filtered_cnaqc_res, filtered_ccf_res) %>% 
  filter(cna_QC=='PASS' & ccf_QC=='PASS')

filtered_cnaqc = all_cnaqc[unique(filtered_results$sample)]

filtered_cnaqc = lapply(filtered_cnaqc,
                        function(x) {
                          
                          kar = filtered_results %>% 
                            filter(sample==x$snvs$sample[1]) %>%
                            select(karyotype) %>%
                            unique() %>%
                            unlist()
                          
                          x = x %>% subset_by_segment_karyotype(kar)
                          return(x)
                          })
  
#

# all_cnaqc$'00b9d0e6-69dc-4345-bffd-ce32880c8eef' %>% plot_peaks_analysis()
# all_cnaqc$'00b9d0e6-69dc-4345-bffd-ce32880c8eef' %>% plot_CCF()
# plot_qc(all_cnaqc$'00b9d0e6-69dc-4345-bffd-ce32880c8eef')
