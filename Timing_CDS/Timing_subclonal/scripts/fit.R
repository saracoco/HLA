fit_timing = function(segments, mutations, purity,
                      possible_k = c("2:1", "2:2", "2:0"), alpha = .05, min_mutations_number = 2,
                      beta_binomial = F, beta_binomial_disp = 0.01
                      subclonal = FALSE, ccf_clonal = 0, k_subclonal = "1:1") {

  if (beta_binomial) {
    model <- cmdstanr::cmdstan_model("models/mixture_CNA_timing_betabinomial.stan")
    # model <- rstan::stan_model("models/mixture_CNA_timing_betabinomial.stan")
  } else {
    # model <- rstan::stan_model("models/mixture_CNA_timing_binomial.stan")
    model <- cmdstanr::cmdstan_model("models/mixture_CNA_timing_binomial.stan")
  }
  
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
    segment_name <- segment$segment_name

    segment_id <- paste(chr, segment$from, segment$to, sep = "_")

    # Get karyotype
    Major <- segment$Major
    minor <- segment$minor

    k <- paste(Major, minor, sep=':') 

    
    

    peaks <- get_clonal_peaks(k, purity, subclonal = FALSE, ccf_clonal, k_subclonal)
    
    #peaks <- get_subclonal_peaks(k, purity, clonal_weigth)
    
    

    if (k %in% possible_k) {
      # Get info for mutations
      segment_mutations <- mutations %>%
        filter(chr == segment$chr,from > segment$from, to < segment$to) %>%
        drop_na(DP)

      accepted_mutations <- data.frame()
      if (nrow(segment_mutations) > 0) {
        # Check if mutation is inside CI
        probs <- c(alpha/2, 1 - alpha/2)

        DP <- segment_mutations$DP
        NV <- segment_mutations$NV

        accepted_idx <- lapply(1:length(DP), function(i) {
          for (p in peaks) {
            if (beta_binomial) {
              quantiles <- TailRank::qbb(probs, DP[i], p * (1 - beta_binomial_disp) / beta_binomial_disp, (1 - p) * (1 - beta_binomial_disp) / beta_binomial_disp)
            } else {
              quantiles <- qbinom(probs, DP[i], p)
            }

            if ((NV[i] >= quantiles[1]) && (NV[i] <= quantiles[2])) {
              return(i)
            }
          }
        }) %>% unlist()
 
        # Get only good mutations
        accepted_mutations <- data.frame(DP = DP[accepted_idx], NV = NV[accepted_idx])
        
        
        
      }

      if (nrow(accepted_mutations) >= min_mutations_number) {
        accepted_mutations_all = dplyr::tibble(accepted_mutations=accepted_mutations, segment=segment_name)
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

        inference_results <- dplyr::bind_rows(inference_results, dplyr::tibble(tau = tau_posteriors, segment = segment_name, karyotype = k, chr = chr, segment_id = segment_id, segment_name = segment_name)) #substitute segment with segment_name instead of segment_idx for plotting

        summarized_results <- dplyr::bind_rows(summarized_results, dplyr::tibble(tau_low = tau_low, tau_mean = tau_mean, tau_high = tau_high, segment = segment_name, karyotype = k, chr = chr, segment_id = segment_id))
      }
    }
  }

  if (nrow(inference_results) == 0) {
    cli::cli_alert_danger("Inference concluded without errors but with no results.")
    return(NULL)
  }

  return(list(inference_results = inference_results, summarized_results = summarized_results, accepted_mutations=accepted_mutations_all))
}
