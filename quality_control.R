#### Prepare environment ####
require(dplyr)
require(CNAqc)
require(tidyverse)
require(mobster)
library(cowplot)

# take the time
start_time <- Sys.time()

setwd("/Users/solonicolas/DSSC/Genomic Data Analytics/Giulio/exam")

file_path = "rds_PCAWG_selected.rds"

#### Load the data ####

# if the input data are in a directory
# pcawg_files = list.files('rds_giasone', full.names = TRUE)
# pcawg_files = pcawg_files[1:30]
# all_pcawg_rds = lapply(pcawg_files, readRDS)

# if the input data are in a unique RDS file
all_pcawg_rds = readRDS(file_path)

# save some metadata of the samples in a df
columns = c('sample','purity','ploidy','ttype','mutation_drivers')
metadata = Reduce(rbind, lapply(all_pcawg_rds,
                                function(x) {
                                  if(all(columns %in% names(x$metadata))==TRUE){
                                    # compute the avg coverage per sample
                                    x$metadata$avg_coverage = mean(x$mutations$DP, na.rm=T)
                                    return (x$metadata)
                                  } else{return (NA)}
                                })) %>%
  select(sample,purity,ploidy,avg_coverage,ttype,mutation_drivers) %>%
  filter(!across(everything(), is.na))

  
# choose the tumor to analyse
ttype = 'Ovary-AdenoCA'
ttypes = c(ttype) # take some tumor types
#ttypes = unique(metadata$ttype) # take all the tumor types

# create dirs in order to store the results

ttype_dir = paste0(ttype,'/')
csv_results_dir = paste0(ttype_dir,'csv_results/')
plots_dir = paste0(ttype_dir,'plots/')
first_statistics = paste0(ttype_dir,'first_statistics/')

if(dir.exists(ttype_dir)){unlink(ttype_dir, recursive = TRUE)}
dir.create(ttype_dir)
dir.create(csv_results_dir)
dir.create(plots_dir)
dir.create(first_statistics)

# write the metadata in a csv file
metadata %>% write_csv(paste0(csv_results_dir,"original_metadata.csv"))

# index each row with the sample name
# names(all_pcawg_rds) = metadata$sample

#### Choose the filters to apply later ####
min_purity = 0.6
min_avg_coverage = 40
min_absolute_karyotype_mutations = 30
min_karyotype_size = 0.05
tail_cut_off = 0.05

# apply some primary filters on PURITY, COVERAGE and TUMOR TYPE
metadata = metadata %>% 
  filter(purity>min_purity & avg_coverage>min_avg_coverage & ttype %in% ttypes)

all_pcawg_rds = all_pcawg_rds[metadata$sample]
         
# put the value of 1 in the CCF columns if it is made of only NA values
# to avoid the error in the cnaqc package
# moreover, keep only mutations with VAF > 0.05
all_pcawg_rds = lapply(all_pcawg_rds,
                       function(x) {
                         if(all(is.na(x$cna$CCF))==TRUE) { # check if there are all NA
                           x$cna = x$cna %>% 
                             mutate(CCF = 1)
                         }
                         
                         x$mutations = x$mutations %>% 
                           filter(VAF > tail_cut_off)
                         
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
# it's going to take some minutes
all_cnaqc = lapply(all_cnaqc,  
                   function(x) {
                     x = x %>% 
                       CNAqc::analyze_peaks(
                         min_absolute_karyotype_mutations = min_absolute_karyotype_mutations,
                         min_karyotype_size = min_karyotype_size
                       )
                     return(x)
                   })

# take out some results of the peaks analysis
cnaqc_res = data.frame(matrix(ncol =5))
colnames(cnaqc_res) = c("sample",
                        "karyotype",
                        "mutation_multiplicity",
                        "cna_QC",
                        "overall_cna_QC")

for(x in all_cnaqc){
  
  # check if peak analysis went well 
  # they could be NULL because of the filter on the number of mutations
  if(!is.null(x$peaks_analysis)) {
    
    karyotypes = x$peaks_analysis$matches$karyotype
    mutation_multiplicity = x$peaks_analysis$matches$mutation_multiplicity
    QC = x$peaks_analysis$matches$QC
    overall_QC = x$peaks_analysis$QC
    
    for(i in 1:length(karyotypes)){
      cnaqc_res[nrow(cnaqc_res)+1,] = c(x$snvs$sample[1],
                                        karyotypes[i],
                                        mutation_multiplicity[i],
                                        QC[i],
                                        overall_QC)
    }
  }
}

# the first line comes NA, so I remove it manually
cnaqc_res = cnaqc_res[2:dim(cnaqc_res)[1],]

# write the results
cnaqc_res %>% write_csv(paste0(csv_results_dir,"cnaqc_res.csv"))

# apply filters on the cnaqc PASS/FAIL status
filtered_cnaqc_res = cnaqc_res %>% 
  filter(cna_QC=='PASS' & overall_cna_QC=='PASS')

filtered_cnaqc = all_cnaqc[unique(filtered_cnaqc_res$sample)]

# plot peaks analysis
lapply(filtered_cnaqc,
       function(x) {
         
         # create dir for each sample
         path = paste0(plots_dir, x$snvs$sample[1], '/')
         dir.create(path)
         
         pdf(paste0(path, 'peaks_analysis.pdf'),width=10,height=9/16*10)
         print(x %>% plot_peaks_analysis(empty_plot=FALSE))
         dev.off()
       })

# subset only the PASS karyotypes
subset_cnaqc = lapply(filtered_cnaqc,
                        function(x) {
                          
                          kar = filtered_cnaqc_res %>%
                            filter(sample==x$snvs$sample[1]) %>%
                            select(karyotype) %>%
                            unique() %>%
                            unlist()
                          
                          x = x %>% subset_by_segment_karyotype(kar)
                          return(x)
                        })

# CCF calculation
# it's going to take some minutes
all_ccf = lapply(subset_cnaqc,
                 function(x) {
                   x = x %>% CNAqc::compute_CCF()
                   return(x)
                  })

# take out some results of the CCF analysis
ccf_res = data.frame(matrix(ncol = 3))
colnames(ccf_res) = c("sample","karyotype","ccf_QC")

for(x in all_ccf){
  
  # check if CCF analysis went well 
  # it could be NULL because of the filter on the number of mutations
  if(!is.null(x$CCF_estimates)) {
    
    for(ccf in x$CCF_estimates){
      ccf_res[nrow(ccf_res)+1,] = c(x$snvs$sample[1],
                                    ccf$QC_table$karyotype,
                                    ccf$QC_table$QC)
    }
  }
} 

# the first line comes NA, so I remove it manually
ccf_res = ccf_res[2:dim(ccf_res)[1],]

# write the results
ccf_res %>% write_csv(paste0(csv_results_dir,"ccf_res.csv"))

# apply last filters on the CCF PASS/FAIL status
filtered_ccf_res = ccf_res %>%
  filter(ccf_QC=='PASS')

filtered_ccf = all_ccf[unique(filtered_ccf_res$sample)]

# plot CCF results
lapply(filtered_ccf,
       function(x) {
         path = paste0(plots_dir, x$snvs$sample[1], '/')
         
         pdf(paste0(path, 'CCF_analysis.pdf'),width=10,height=9/16*10)
         print(x %>% plot_CCF(empty_plot=FALSE))
         dev.off()
       })

# subset only the PASS karyotypes
subset_ccf = lapply(filtered_ccf,
                    function(x) {
                      
                      kar = filtered_ccf_res %>%
                        filter(sample==x$snvs$sample[1]) %>%
                        select(karyotype) %>%
                        unique() %>%
                        unlist()
                      
                      x = x %>% subset_by_segment_karyotype(kar)
                      return(x)
                    })


#### A bit of statistics ####

# let's plot the distribution of the karyotypes

# situation after peaks analysis
peaks_kar_plots = lapply(subset_cnaqc,
                   function(x) {
                     p = plot_karyotypes(x)
                     p$data$call = ""
                     p$labels$title = ""
                     p$labels$x = ""
                     p$labels$y = ""
                     p$theme$plot.margin=margin(0,0,0,0)
                     p$theme$strip.text$margin=margin(0,0,0,0)
                     p$theme$panel.spacing=margin(0)
                     return(p)
                   })

pdf(paste0(first_statistics,'karyotypes_after_peaks'))
ggarrange(plotlist=peaks_kar_plots,
          common.legend = TRUE,
          legend = 'bottom')
dev.off()

# situation after both peaks and CCF analysis
ccf_kar_plots = lapply(subset_ccf,
                   function(x) {
                     p = plot_karyotypes(x)
                     p$data$call = ""
                     p$labels$title = ""
                     p$labels$x = ""
                     p$labels$y = ""
                     p$theme$plot.margin=margin(0,0,0,0)
                     p$theme$strip.text$margin=margin(0,0,0,0)
                     p$theme$panel.spacing=margin(0)
                     return(p)
                   })

pdf(paste0(first_statistics,'karyotypes_after_CCF'))
ggarrange(plotlist=ccf_kar_plots,
          common.legend = TRUE,
          legend = 'bottom')
dev.off()



#### Mobster deconvolution ####

# here we should decide whether to take diploid mutations using VAF
# or aneuploid mutations using CCF values
final_karyotype = 'max' #'2:1'
all_mobster = subset_ccf # subset_cnaqc
filtered_res = filtered_ccf_res # filtered_cnaqc_res

if(final_karyotype == 'max') {
    
  # subset only the karyotype with the maximum number of mutations
  all_mobster = lapply(all_mobster,
                       function(x) {
                         
                         # max karyotype
                         max_kar = names(x$n_karyotype[1])
                         
                         x = x %>% 
                           subset_snvs() %>% 
                           subset_by_segment_karyotype(c(max_kar))
                         return(x)
                         })
    
} else {
  # keep only sample with 'final_karyotype'
  all_mobster = all_mobster[unique(
    filtered_res$sample[filtered_res$karyotype==final_karyotype])]
    
  # subset only the final karyotype
  all_mobster = lapply(all_mobster,
                       function(x) {
                         x = x %>% 
                           subset_snvs() %>% 
                           subset_by_segment_karyotype(c(final_karyotype))
                         return(x)
                        })
}

# if we chose NON diploid mutations
# we should add the CCF column to the data
# we calculate again CCf values!
all_mobster = lapply(all_mobster,
                     function(x) {
                       x = x %>% CNAqc::compute_CCF()
                       x$snvs$VAF = x$CCF_estimates[[1]]$mutations$VAF
                       return(x)
                       })

# perform mobster deconvolution
mobster_fits = lapply(all_mobster,
                      function(x) {
                        x = x$snvs %>% mobster_fit(auto_setup = 'FAST')
                        return(x)
                      })

# take out some results of the mobster deconvolution
mobster_res = data.frame(matrix(ncol = 5))
colnames(mobster_res) = c("sample","kar","clusters","tail","n")

for(x in mobster_fits){
  
  mobster_res[nrow(mobster_res)+1,] = c(x$best$Call$X$sample[1],
                                        x$best$Call$X$karyotype[1],
                                        x$best$Kbeta,
                                        x$best$fit.tail,
                                        x$best$N)
}

# the first line comes NA, so I remove it manually
mobster_res = mobster_res[2:dim(mobster_res)[1],]

# write the results
kar = ifelse(final_karyotype=='max',
           final_karyotype,
           gsub(":","_",final_karyotype))

mobster_res %>% write_csv(paste0(csv_results_dir,
                                 "mobster_res_",kar,".csv"))


# plot some mobster decovolution results
lapply(mobster_fits,
       function(x) {
         
         path = paste0(plots_dir, x$best$Call$X$sample[1], '/')
         pdf(paste0(path,'deconvolution_',
                    ifelse(kar!='max',kar,
                           paste0(gsub(":",
                                       "_",
                                       x$best$Call$X$karyotype[1]),"_max")),'.pdf'),
             width=10,
             height=9/16*10)
         print(x$best %>% plot())
         dev.off()
       })


# filter only sample with the tail
tail_res = mobster_res %>% filter(tail==TRUE)

tail_samples = mobster_fits[unique(tail_res$sample)]

# calculation of the evolutionary parameters
evol_params = data.frame(matrix(ncol = 6))
colnames(evol_params) = c("sample","kar","mu","clusters","s","t")

for(x in tail_samples){
  
  params = x %>% evolutionary_parameters()
  k = x$best$Call$X$karyotype[1]
  
  # if there are subclones
  if((dim(params)[2]>2 & k=="1:1") | (dim(params)[1]>1 & k!="1:1")) {
    for(i in 1:dim(params)[1]) {
      evol_params[nrow(evol_params)+1,] = c(x$best$Call$X$sample[1],
                                            k,
                                            params$mu[i],
                                            params$cluster[i],
                                            params$s[i],
                                            params$time[i])
    }
  } else { # if there are NO subclones
    evol_params[nrow(evol_params)+1,] = c(x$best$Call$X$sample[1],
                                          k,
                                          params$mu,
                                          NA,
                                          NA,
                                          NA)
  }
}

# the first line comes NA, so I remove it manually
evol_params = evol_params[2:dim(evol_params)[1],]

# write the results
evol_params %>% write_csv(paste0(csv_results_dir,
                                 "evolutionary_parameters_",kar,".csv"))

# take the time
end_time <- Sys.time()
total_time = end_time - start_time
print(total_time)
