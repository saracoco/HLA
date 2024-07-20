# PLOT
plot_inference = function(fit_results, segments, colour_by = "karyotype") {
  reference_genome <- CNAqc::chr_coordinates_GRCh38

  vfrom = reference_genome$from
  names(vfrom) = reference_genome$chr

  absoulte_segments <- segments %>%
    mutate(from = as.integer(from) + vfrom[chr],
           to = as.integer(to) + vfrom[chr])

  summarized_results <- fit_results$summarized_results %>%
    mutate(from = absoulte_segments[segment,]$from) %>%
    mutate(to = absoulte_segments[segment,]$to) %>%
    mutate(tau_mean = ifelse(tau_mean < 1, tau_mean, 1)) %>%
    mutate(tau_high = ifelse(tau_high < 1, tau_high, 1))

  p <- ggplot() +
    geom_rect(data = summarized_results, aes(xmin=from, xmax=to, ymin=tau_low, ymax=tau_high, fill = as.factor(.data[[colour_by]])), alpha = .5) +
    geom_segment(data = summarized_results, aes(y = tau_mean, yend = tau_mean, x = from, xend = to)) +
    scale_x_continuous(breaks = reference_genome$to, labels = gsub("chr", "", reference_genome$chr)) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 90)) +
    lims(y = c(0,1)) +
    labs(x = "chromosome", y = bquote(tau))
    # scale_fill_manual(values = c("forestgreen", "indianred3", "steelblue"), name = "")

  return(p)
}
