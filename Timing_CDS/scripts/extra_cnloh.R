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


UPN01_NV = UPN01 %>% filter(NV.DIA %in% c(0, 'NA'),
                            NV_PG.REL %in% c(0, 'NA'), 
                            NV_DG.REL %in% c(0, 'NA'),
                            NV_PG.REL2 %in% c(0, 'NA'),
                            NV_DG.REL2 %in% c(0, 'NA'),
                            PASS == TRUE)

# UPN02_NV = UPN02 %>% filter(NV_PG.REL %in% c(0, 'NA'),
#                             NV_DG.REL %in% c(0, 'NA'),
#                             NV.PRE %in% c(0, 'NA'), 
#                             NV_PG.PRE %in% c(0, 'NA'),
#                             PASS == TRUE)

UPN02_NV = UPN02 %>% filter(NV.PRE %in% c(0, 'NA'), 
                            PASS == TRUE)


# UPN03_NV = UPN03 %>% filter(NV_PG.REL %in% c(0, 'NA'), 
#                             NV_DG.REL %in% c(0, 'NA'),
#                             NV_PG.PRE %in% c(0, 'NA'),
#                             NV.PRE %in% c(0, 'NA'),
#                             PASS == TRUE)

UPN03_NV = UPN03 %>% filter(NV.PRE %in% c(0, 'NA'),
                            PASS == TRUE)

# UPN04_NV = UPN04 %>% filter(NV_PG.REL %in% c(0, 'NA'), 
#                             NV_DG.REL %in% c(0, 'NA'),
#                             NV_PG.PRE %in% c(0, 'NA'),
#                             NV.PRE %in% c(0, 'NA'),
#                             NV.DIA %in% c(0, 'NA'),
#                             NV_PG.DIA %in% c(0, 'NA'),
#                             PASS == TRUE)

UPN04_NV = UPN04 %>% filter(NV.PRE %in% c(0, 'NA'),        #ACCEPTED MUTATIONS > 0
                            PASS == TRUE)

# UPN05_NV = UPN05 %>% filter(NV_PG.REL %in% c(0, 'NA'), 
#                             NV_DG.REL %in% c(0, 'NA'),
#                             NV_PG.PRE %in% c(0, 'NA'),
#                             NV.PRE %in% c(0, 'NA'),
#                             PASS == TRUE)

UPN05_NV = UPN05 %>% filter(NV.PRE %in% c(0, 'NA'),        #ACCEPTED MUTATIONS > 0
                            PASS == TRUE)


# UPN10_NV = UPN10 %>% filter(NV_PG.REL %in% c(0, 'NA'), 
#                             NV_DG.REL %in% c(0, 'NA'),
#                             NV.DIA %in% c(0, 'NA'),
#                             NV_PG.DIA %in% c(0, 'NA'),
#                             PASS == TRUE)

UPN10_NV = UPN10 %>% filter(NV.DIA %in% c(0, 'NA'), #NO PRE! CHIEDERE   #ACCEPTED MUTATIONS > 0
                            PASS == TRUE)

# UPN11_NV = UPN11 %>% filter(NV_PG.REL %in% c(0, 'NA'), 
#                             NV_DG.REL %in% c(0, 'NA'),
#                             NV.DIA %in% c(0, 'NA'),
#                             NV_PG.DIA %in% c(0, 'NA'),
#                             PASS == TRUE)

UPN11_NV = UPN11 %>% filter(NV.DIA %in% c(0, 'NA'), #NO PRE! CHIEDERE  
                            PASS == TRUE)



data <- list(UPN01_NV, UPN02_NV, UPN03_NV, UPN04_NV, UPN05_NV, UPN10_NV, UPN11_NV)
names <- c("UPN01","UPN02","UPN03","UPN04","UPN05","UPN10","UPN11")









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

purity = 1
# i=4
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




# data_cna
# table(inference_results$segment_name)

fit_02 <- list(inference_results = fit$inference_results %>% filter(segment_name == "UPN02_NV"), summarized_results = fit$summarized_results %>% filter(segment == "UPN01_NV"))
data_cna_02 <- data_cna %>% filter(segment_name == "UPN02_NV")


# plot tau distributions for different copy number events
p <- plot_inference(fit_results = fit, segments = data_cna)

#fallo simile a quello di intervals bayesplot
p_istogram <- inference_results %>%
  ggplot(mapping = aes(x=tau, fill=segment_name, col = segment_name)) +
  geom_histogram(binwidth=0.02, alpha = 0.5, position = "identity") +
  theme(legend.position = 'top')


table(fit$accepted_mutations$segment)














#PROVO A FARE FIT SENZA FILTERING PHASE? as they alreadi categorize the mutations? 
Input_data <- data[[7]]

#Input_data = UPN05_NV
data_cna <- Input_data %>% mutate (Major = unlist(strsplit(Input_data$segment.REL[1], ":"))[1],
                                   minor = unlist(strsplit(Input_data$segment.REL[1], ":"))[2],
                                   purity = purity,
                                   from = min(Input_data$from),
                                   to = max(Input_data$to))

mutations <- Input_data %>% mutate (NV = NV.REL, DP = DP.REL) 


segments = data_cna[1,]
mutations = mutations 
purity = purity
possible_k = c("2:1", "2:2", "2:0")
alpha = .05
min_mutations_number = 2
beta_binomial = F
beta_binomial_disp = 0.01

#model <- cmdstanr::cmdstan_model("models/mixture_CNA_timing_binomial.stan")






segments <- segments %>%
  drop_na(Major, minor)

n_segments <- nrow(segments)
inference_results <- dplyr::tibble()
summarized_results <- dplyr::tibble()

for (segment_idx in 1:n_segments) {
  print(segment_idx)
  
  # Segments
  segment <- segments[segment_idx, ]
  chr <- segment$chr
  
  segment_id <- paste(chr, segment$from, segment$to, sep = "_")
  
  # Get karyotype
  Major <- segment$Major
  minor <- segment$minor
  
  k <- paste(Major, minor, sep=':')
  
  peaks <- get_clonal_peaks(k, purity)
  
  if (k %in% possible_k) {
    # Get info for mutations
    segment_mutations <- mutations %>%
      filter(chr == segment$chr,from > segment$from, to < segment$to) %>%
      drop_na(DP)
    
    accepted_mutations <- data.frame(segment_mutations)
    
    
    if (nrow(accepted_mutations) >= min_mutations_number) {
      
      cli::cli_alert_info("Fitting segment with index {.val {segment_idx}}")
      
      # model and input data
      input_data <- list(
        N = nrow(accepted_mutations),
        NV = accepted_mutations$NV,
        DP = accepted_mutations$DP,
        peaks = peaks,
        beta_dispersion = beta_binomial_disp
      )
      
      fit <- model$sample(data=input_data, iter_warmup=2000, iter_sampling=2000, chains=8, parallel_chains=8)
      # fit <- model$pathfinder(data=input_data, algorithm = 'single')
      
      # fit <- rstan::sampling(model, input_data, chains = 4, warmup = 1000, iter = 2000)
      
      # Compute tau posteriors
      tau_posteriors <- get_tau_posteriors(fit, k)$tau %>% unname() %>% as.numeric()
      
      q1 <- alpha / 2
      q2 <- 1 - alpha / 2
      tau_low <- quantile(tau_posteriors, q1) %>% unname()
      tau_high <- quantile(tau_posteriors, q2) %>% unname()
      tau_mean <- mean(tau_posteriors)
      
      inference_results <- dplyr::bind_rows(inference_results, dplyr::tibble(tau = tau_posteriors, segment = segment_idx, karyotype = k, chr = chr, segment_id = segment_id))
      
      summarized_results <- dplyr::bind_rows(summarized_results, dplyr::tibble(tau_low = tau_low, tau_mean = tau_mean, tau_high = tau_high, segment = segment_idx, karyotype = k, chr = chr, segment_id = segment_id))
    }
  }
}






#fallo simile a quello di intervals bayesplot
p_istogram <- inference_results %>%
  ggplot(mapping = aes(x=tau, fill=segment_id, col = segment_id)) +
  geom_histogram() +
  theme(legend.position = 'none')



