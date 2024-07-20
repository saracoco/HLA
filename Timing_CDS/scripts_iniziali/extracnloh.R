# load dependences
library(ggplot2)
library(tidyverse)
library(rstan)
library(cmdstanr)
source("scripts/fit.R")
source("scripts/utils.R")
source("scripts/simulate_mutations.R")
source("scripts/plot.R")


UPN01 <- readRDS("Data/extra_cnloh/alpha_beta/UPN01/mutations.rds")
UPN02 <- readRDS("Data/extra_cnloh/alpha_beta/UPN02/mutations.rds")
UPN03 <- readRDS("Data/extra_cnloh/alpha_beta/UPN03/mutations.rds")
UPN04 <- readRDS("Data/extra_cnloh/alpha_beta/UPN04/mutations.rds")
UPN05 <- readRDS("Data/extra_cnloh/alpha_beta/UPN05/mutations.rds")
UPN10 <- readRDS("Data/extra_cnloh/alpha_beta/UPN10/mutations.rds")
UPN11 <- readRDS("Data/extra_cnloh/alpha_beta/UPN11/mutations.rds")


# ALTERNATIVE FILTERING timing_classification 

UPN01_NV = UPN01 %>% filter(timing_classification %in% c("alpha private", "beta", "beta private"), #ok
                            PASS == TRUE)

UPN02_NV = UPN02 %>% filter(timing_classification %in% c("alpha private", "beta", "beta private"), #no
                            PASS == TRUE)

UPN03_NV = UPN03 %>% filter(timing_classification %in% c("alpha private", "beta", "beta private"), #no
                            PASS == TRUE)

UPN04_NV = UPN04 %>% filter(timing_classification %in% c("alpha private", "beta", "beta private"), #ok
                            PASS == TRUE)

UPN05_NV = UPN05 %>% filter(timing_classification %in% c("alpha private", "beta", "beta private"), #ok
                            PASS == TRUE)

UPN10_NV = UPN10 %>% filter(timing_classification %in% c("alpha private", "beta", "beta private"), #ok
                            PASS == TRUE)

UPN11_NV = UPN11 %>% filter(timing_classification %in% c("alpha private", "beta", "beta private"), #no
                            PASS == TRUE)

data <- list(UPN01_NV = UPN01_NV, UPN02_NV = UPN02_NV, UPN03_NV = UPN03_NV, UPN04_NV = UPN04_NV, UPN05_NV = UPN05_NV, UPN10_NV = UPN10_NV, UPN11_NV = UPN11_NV)
names <- c("UPN01_NV","UPN02_NV","UPN03_NV","UPN04_NV","UPN05_NV","UPN10_NV", "UPN11_NV")


data_cna <- dplyr::tibble()
inference_results <- dplyr::tibble()
summarized_results <- dplyr::tibble()
accepted_mutations <- dplyr::tibble()

purity = 1 # set purity
for(i in 1:length(names)){
  
  Input_data <- data[[names[i]]]
  
  #Input_data = UPN05_NV
  data_cna_single <- Input_data %>% mutate (Major = unlist(strsplit(Input_data$segment.REL[1], ":"))[1],
                                            minor = unlist(strsplit(Input_data$segment.REL[1], ":"))[2],
                                            purity = purity,
                                            from = min(Input_data$from), ### GIUSTO?
                                            to = max(Input_data$to),
                                            segment_name = names[i])
  
  #data_cna <- (seq_df_SC$seg_id %>% unique) aggiusta select only one segment 
  
  mutations <- Input_data %>% mutate (NV = NV.REL, DP = DP.REL) 
  
  fit <- fit_timing(segments = data_cna_single[1,], mutations = mutations, purity=purity)
  
  
  inference_results <- dplyr::bind_rows(inference_results, fit$inference_results)
  summarized_results <- dplyr::bind_rows(summarized_results, fit$summarized_results)
  
  
  data_cna_single <- data_cna_single[1,] %>% select(from, to, chr, segment_name) 
  data_cna <- dplyr::bind_rows(data_cna, data_cna_single)
  accepted_mutations <- dplyr::bind_rows(accepted_mutations, fit$accepted_mutations)
  
}


fit <- list(inference_results = inference_results, summarized_results = summarized_results, accepted_mutations=accepted_mutations)




p_istogram <- inference_results %>%
  ggplot(mapping = aes(x=tau, fill=segment_name, col = segment_name)) +
  geom_histogram(binwidth=0.02, alpha = 0.5, position = "identity") +
  theme(legend.position = 'top')


table(fit$accepted_mutations$segment)