#sample <- readRDS("~/Desktop/projects/CNA/50515723-b495-42a9-8750-e3da288bf6a3.rds")
setwd("~/Desktop/dottorato/CNAqc_project/2. Simple CN timing")
source("scripts/utils.R")
source("scripts/fit.R")
source("scripts/plot.R")
source("scripts/cluster_tau.R")
library(tidyverse)
library(MutationTimeR)

set.seed(27)
cn <- MutationTimeR:::refLengths[1:23]
t1 <- 0.7
purity <- 0.85
clusters <- data.frame(cluster=1:1, proportion=c(purity), n_ssms=c(100))
time2pi <- function(N,n,t1,t2){
  if(N==2 & n ==1)
    pi <-  c(3-2*t1, t1)
  else if(N==2 & n %in% c(0,2))
    pi <- c(2 -2*t1, t1)
  else pi <- 1
  pi <- pi/sum(pi)
}

cn$major_cn <- 2 #sample(1:2, 23, replace=TRUE)
cn$minor_cn <- sample(0:2, 23, replace=TRUE)
cn$clonal_frequency <- purity
cn$timing_param <- MutationTimeR:::defineMcnStates(cn, purity=purity, clusters=clusters, deltaFreq = 0)

for(i in seq_along(cn)){
  pi <- time2pi(cn$major_cn[i], cn$minor_cn[i], t1, t2)
  pi_sub <- clusters$n_ssms
  pi_sub[1] <- pi_sub[1] * (t1*2 + (1-t1)*(cn$major_cn[i]+ cn$minor_cn[i])) / 2
  pi_sub[-1] <- pi_sub[-1] * (cn$major_cn[i]+ cn$minor_cn[i]) / 2
  pi_sub <- pi_sub/sum(pi_sub)
  pi <- c(pi * pi_sub[1], pi_sub[-1])
  cn$timing_param[[i]][,"P.m.sX"] <- pi
  cn$timing_param[[i]][,"power.m.s"] <- rep(1, length(pi))
}
cn$n.snv_mnv <- width(MutationTimeR:::refLengths[1:23])/1e6 * 10

vcf <- MutationTimeR:::simulateMutations(cn, rho=0.01, n=40)

cn_timed <- cn[,c("major_cn","minor_cn","clonal_frequency")]


mt <- mutationTime(vcf, cn_timed, clusters=clusters, n.boot=10)
mcols(cn_timed) <- cbind(mcols(cn_timed),mt$T)

info(header(vcf)) <- rbind(info(header(vcf)),MutationTimeR:::mtHeader())
info(vcf) <- cbind(info(vcf), mt$V)

plotSample(vcf,cn_timed)

# Run our model
data <- parse_gerstrung(vcf = vcf, cna = cn)

fit <- fit_timing(segments = data$segments, mutations = data$mutations, purity = purity,
                  alpha = .05, min_mutations_number = 3,
                  beta_binomial = F, beta_binomial_disp = 0.01)
plot_inference(fit, data$segments)

# Check number of our timed segments versus theirs
gerstung_n_segs <- length(cn_timed)
cvg_n_segs <- nrow(fit$summarized_results)

cn_timed

s <- start(cn_timed)
e <- end(cn_timed)
chr <- seq(cn_timed)

timings <- mcols(cn_timed)
t_low <- timings$time.lo
t <- timings$time
t_high <- timings$time.up

comparison_results <- dplyr::tibble()
for (i in 1:gerstung_n_segs) {
  seg_id_i <- paste0("chr", chr[i], "_", s[i], "_", e[i])
  cvg_seg <- fit$summarized_results %>%
    filter(segment_id == seg_id_i)

  if (nrow(cvg_seg) == 1) {
    cvg_t <- cvg_seg %>% pull(tau_mean)
    cvg_l <- cvg_seg %>% pull(tau_low)
    cvg_h <- cvg_seg %>% pull(tau_high)
    spotted <- ((cvg_t >= t_low[i]) & (cvg_t <= t_high[i]))
    iou <- IOU(first = c(cvg_l, cvg_h), second = c(t_low[i], t_high[i]))
  } else {
    spotted <- NULL
    iou <- NULL
  }

  comparison_results <- dplyr::bind_rows(comparison_results, dplyr::tibble(seg_id = seg_id_i, spotted = spotted, iou = iou))
}

comparison_results

ggplot(data = comparison_results) +
  geom_histogram(mapping = aes(x=spotted), stat="count")

ggplot(data = comparison_results) +
  geom_histogram(mapping = aes(x=iou), bins = 10) +
  lims(x = c(0,1))

