#
#

library(patchwork)
library(ggplot2)
library(tidyverse)
library(rstan)
library(cmdstanr)
library(grid)
library(gridExtra)

setwd("E:/u/cdslab/scocomello/scratch/CDS_ORFEO/Timing_CDS")
setwd("C:/Users/sarac/HLA/Timing_CDS")


source("scripts/fit.R") #fit_no_filtering
source("scripts/utils.R")
#source("scripts/simulate_mutations.R")
source("scripts/plot.R")



# LOAD DATA #
UPN04_extra <- readRDS("../Data/extra_cnloh/alpha_beta/UPN04/mutations.rds")
UPN05_extra <- readRDS("../Data/extra_cnloh/alpha_beta/UPN05/mutations.rds")

UPN04_alpha_beta <- readRDS("../Data/alpha_beta/UPN04/mutations.rds")
UPN05_alpha_beta <- readRDS("../Data/alpha_beta/UPN05/mutations.rds")


# FILTERING #
UPN04_extra_NV = UPN04_extra %>% filter(timing_classification %in% c("alpha private", "beta"),PASS == TRUE)
UPN05_extra_NV = UPN05_extra %>% filter(timing_classification  %in% c("alpha private", "beta"),chr == "chr1",PASS == TRUE)

UPN04_alpha_beta_NV = UPN04_alpha_beta %>% ungroup %>% filter(timing_classification %in% c("alpha"),PASS == TRUE)
UPN05_alpha_beta_NV = UPN05_alpha_beta %>% filter(timing_classification %in% c("alpha", "beta"),PASS == TRUE)


data <- list(UPN04 = UPN04_extra_NV, UPN05 = UPN05_extra_NV, UPN04_LSH = UPN04_alpha_beta_NV, UPN05_LSH = UPN05_alpha_beta_NV)
names <- c("UPN04","UPN04_LSH", "UPN05", "UPN05_LSH")

# names <- c("UPN04","UPN04_LSH")
# names <- c("UPN05", "UPN05_LSH")


# for (i in 1:length(names)){
#   assign("x", data[[names[i]]])
#   hist_data <- x %>% 
#   ggplot(mapping = aes(x=(NV.REL/DP.REL), fill=timing_classification, col = timing_classification)) + #
#   geom_histogram(binwidth=0.01, alpha = 0.5, position = "identity") +
#   theme(legend.position = 'top') +
#   ggtitle( names[i], data[[names[i]]]$segment.REL)
# 
#   print(hist_data)
#   
#   ggsave(paste0("plots/data",names[i],".png"))
# }
# 
# 
# for (i in 1:length(names)){
#   assign("x", data[[names[i]]])
#   hist_data <- x %>% 
#   ggplot(mapping = aes(x=(NV.REL), fill=timing_classification, col = timing_classification)) + #
#   geom_histogram(binwidth=0.01, alpha = 0.5, position = "identity") +
#   theme(legend.position = 'top') +
#   ggtitle( names[i], data[[names[i]]]$segment.REL)
# 
#   print(hist_data)
#   
#   ggsave(paste0("plots/data",names[i],"_NV.png"))
# }



data_cna <- dplyr::tibble()
inference_results <- dplyr::tibble()
summarized_results <- dplyr::tibble()
accepted_mutations <- dplyr::tibble()
component_binomial <- dplyr::tibble()
y_rep <- dplyr::tibble()
omega_1 <- dplyr::tibble()
omega_2 <- dplyr::tibble()


purity = 0.98
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
  
  saveRDS(fit, paste0("results/fit", names[i],".rds"))
  
  
  inference_results <- dplyr::bind_rows(inference_results, fit$inference_results)
  summarized_results <- dplyr::bind_rows(summarized_results, fit$summarized_results)
  component_binomial <- dplyr::bind_rows(component_binomial, fit$component_binomial)
  y_rep <- dplyr::bind_rows(y_rep, fit$y_rep)
  omega_1 <- dplyr::bind_rows(omega_1, fit$omega_1)
  omega_2 <- dplyr::bind_rows(omega_2, fit$omega_2)


  
  data_cna_single <- data_cna_single[1,] %>% select(from, to, chr, segment_name) 
  data_cna <- dplyr::bind_rows(data_cna, data_cna_single)
  accepted_mutations <- dplyr::bind_rows(accepted_mutations, fit$accepted_mutations)
  
}



fit <- list(omega_1 = omega_1, omega_2 = omega_2, inference_results = inference_results, summarized_results = summarized_results, accepted_mutations=accepted_mutations, component_binomial = component_binomial, y_rep = y_rep)

saveRDS(fit, "fit_extra_alpha_beta_LOH.rds")
fit_extra_alpha_beta<-readRDS("fit_extra_alpha_beta_LOH.rds")

fit_extra_alpha_beta$summarized_results




library(wesanderson)

# hist_data_extra <- as.data.frame(data["UPN04"]$UPN04) %>% 
#   ggplot(mapping = aes(x=(NV.REL/DP.REL), fill=timing_classification, col = timing_classification)) +
#   geom_histogram(binwidth=0.01, alpha = 0.5, position = "identity") +
#   theme(legend.position = 'top') +
#   ggtitle( names[i], data[[names[i]]]$segment.REL)


# overlap on the same graph
hist_data_extra <- as.data.frame(data["UPN04"]$UPN04) %>% 
  ggplot(mapping = aes(x=(NV.REL/DP.REL),  fill = timing_classification)) +
  scale_color_manual(values = c("black"))+
  scale_fill_manual(values = wes_palette("Darjeeling1"))+
  geom_histogram(binwidth=0.01, alpha = 0.5, position = "identity") +
  theme(legend.position = 'top') +
  xlim(0, 1) +
  ggtitle( "UPN04", data["UPN04"]$UPN04$segment.REL)


hist_data_LSH <- as.data.frame(data["UPN04_LSH"]$UPN04_LSH) %>% 
  ggplot(mapping = aes(x=(NV.REL/DP.REL), fill=timing_classification)) +
  scale_color_manual(values = c("black"))+
  scale_fill_manual(values = wes_palette("Darjeeling1"))+
  geom_histogram(binwidth=0.01, alpha = 0.5, position = "identity") +
  theme(legend.position = 'top') +
  xlim(0, 1) +
  ggtitle("UPN04_LSH", data["UPN04_LSH"]$UPN04$segment.REL)


p_istogram <- fit_extra_alpha_beta$inference_results %>% filter(segment_name == "UPN04" | segment_name == "UPN04_LSH") %>%
  ggplot(mapping = aes(x=tau, fill=segment_name)) +
  scale_fill_manual(values = wes_palette("Cavalcanti1"))+
  geom_histogram(binwidth=0.02, alpha = 0.7, position = "identity") +
  theme(legend.position = 'top') +
  xlim(0, 1) 


plot_inference_04 <- (p_istogram + (hist_data_extra / hist_data_LSH)) + 
      plot_layout(widths = c(10,6), heights = c(10,1)) +
      plot_annotation(
        title = 'UPN04 LSH and extra event ',
        subtitle = "", #Extra event is on chr 
        caption = "" #caption
      ) & theme(text = element_text(size = 12), plot.title = element_text(size = 15), plot.subtitle = element_text(size = 12), axis.text = element_text(size = 10), plot.caption = element_text(size = 5))
plot_inference_04

ggsave(paste0("plots/inference_UPN04.png"), width = 18, height = 10)






hist_data_extra <- as.data.frame(data["UPN05"]$UPN05) %>% 
  ggplot(mapping = aes(x=(NV.REL/DP.REL), fill=timing_classification)) +
  geom_histogram(binwidth=0.01, alpha = 0.5, position = "identity") +
  scale_fill_manual(values = wes_palette("Darjeeling1"))+
  theme(legend.position = 'top') +
  xlim(0, 1) +
  ggtitle( "UPN05", data["UPN05"]$UPN05$segment.REL)


hist_data_LSH <- as.data.frame(data["UPN05_LSH"]$UPN05_LSH) %>% 
  ggplot(mapping = aes(x=(NV.REL/DP.REL), fill=timing_classification)) +
  geom_histogram(binwidth=0.01, alpha = 0.5, position = "identity") +
  scale_fill_manual(values = wes_palette("Darjeeling1"))+
  theme(legend.position = 'top') +
  xlim(0, 1) +
  ggtitle( "UPN05_LSH", data["UPN05_LSH"]$UPN05_LSH$segment.REL)


p_istogram <- fit_extra_alpha_beta$inference_results %>% filter(segment_name == "UPN05" | segment_name == "UPN05_LSH") %>%
  ggplot(mapping = aes(x=tau, fill=segment_name)) + #
  scale_fill_manual(values = wes_palette("Cavalcanti1"))+
  geom_histogram(binwidth=0.02, alpha = 0.7, position = "identity") +
  xlim(0, 1) +
  theme(legend.position = 'top')

  
plot_inference_05 <- (p_istogram + (hist_data_extra / hist_data_LSH)) + 
  plot_layout(widths = c(10,6), heights = c(10,1)) +
  plot_annotation(
        title = 'UPN05 LSH and extra event ',
        subtitle = "", #Extra event is on chr 
        caption = "" #caption
      ) & theme(text = element_text(size = 12), plot.title = element_text(size = 15), plot.subtitle = element_text(size = 12), axis.text = element_text(size = 10), plot.caption = element_text(size = 5))
plot_inference_05

ggsave(paste0("plots/inference_UPN05.png"), width = 18, height = 10)


# Facet on segments 
# p_istogram_2 <- fit_extra_alpha_beta$inference_results %>%
#   ggplot(mapping = aes(x=tau, fill=segment_name, col = segment_name)) + 
#   geom_histogram(binwidth=0.01, alpha = 0.5, position = "identity") +
#   theme(legend.position = 'none') +
#   facet_grid( vars(segment_name) , scales = "free") 
# p_istogram_2


