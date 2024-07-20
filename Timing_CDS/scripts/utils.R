get_clonal_peaks = function(k, purity) {
  multiplicities <- strsplit(k, ":") %>% unlist() %>% as.numeric()
  major <- multiplicities[1]
  n_tot <- sum(multiplicities)

  # get only Major and 1
  multiplicities <- c(1, major)

  peaks <- unique(multiplicities) * purity / (n_tot * purity + 2 * (1 - purity))
  return(sort(peaks))
}

get_tau_posteriors = function(fit, k) {
  # omega1 <- fit$draws("omega[1]", format = 'matrix') %>% as_tibble() %>% `colnames<-`('value') %>% dplyr::mutate(value = as.numeric(value)) %>% dplyr::pull(value)
  # omega2 <- fit$draws("omega[2]", format = 'matrix') %>% as_tibble() %>% `colnames<-`('value') %>% dplyr::mutate(value = as.numeric(value)) %>% dplyr::pull(value)
  omega1 <- fit$draws("omega[1]", format = "matrix") %>% as.numeric()   
  omega2 <- fit$draws("omega[2]", format = "matrix") %>% as.numeric()

  if (k == '2:1') {
    tau_posterior <- 3 * omega2 / (2*omega2 + omega1)
  } else {
    tau_posterior <- 2 * omega2 / (2*omega2 + omega1)
  }

  tau_posteriors <- dplyr::tibble(tau = tau_posterior)
  tau_posteriors
}

parse_gerstrung = function(vcf, cna) {
  # Mutations
  NV <- vcf@info$t_alt_count
  DP <- vcf@info$t_alt_count + vcf@info$t_ref_count
  from <- vcf@rowRanges@ranges@start
  to <- vcf@rowRanges@ranges@start + vcf@rowRanges@ranges@width - 1
  chr <- paste0("chr", vcf@rowRanges@seqnames)

  mutations <- dplyr::tibble(chr = chr, from = from, to = to, NV = NV, DP = DP)

  # Segments
  chr <- paste0("chr", cna@seqnames %>% as.numeric())
  from <- cna@ranges@start
  to <- cna@ranges@start + cna@ranges@width - 1

  Major <- cna$major_cn
  minor <- cna$minor_cn
  purity <- cna$clonal_frequency

  segments <- dplyr::tibble(chr = chr, from = from, to = to, Major = Major, minor = minor, purity = purity)

  return(list(mutations = mutations, segments = segments))
}
