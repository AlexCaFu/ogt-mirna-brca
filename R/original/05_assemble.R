# =====================================================================
# R/05_assemble.R    merge layers into the master table and draw Figure 2
# =====================================================================
source("R/config.R")
suppressMessages({library(ggplot2); library(patchwork)})

l1 <- readr::read_csv("results/01_target_consensus.csv", show_col_types = FALSE)
l2 <- readr::read_csv("results/02_tcga_correlation.csv", show_col_types = FALSE)
s3 <- readr::read_csv("results/03_survival.csv",         show_col_types = FALSE)
geo <- tryCatch(readr::read_csv("results/04_geo_foldchange.csv", show_col_types = FALSE),
                error = function(e) NULL)

## per-microRNA survival (feature column holds the mature name for microRNAs)
mir_surv <- s3 %>% filter(feature %in% CANDIDATES$mirna) %>%
  transmute(mirna = feature, miR_HR = HR, miR_HR_p = HR_p)

master <- l2 %>%
  left_join(select(l1, mirna, n_pred_db, validated_in), by = "mirna") %>%
  left_join(mir_surv, by = "mirna")

## attach GEO fold changes if available (miRNA side and enzyme side)
if (!is.null(geo)) {
  mir_fc <- geo %>% filter(measures == "miRNA") %>% group_by(target) %>%
    summarise(miRNA_resistance = paste(sprintf("%s:%+.2f(adjP %.2g)", series, log2FC, adj_p), collapse = "; "), .groups = "drop") %>%
    rename(mirna = target)
  enz_fc <- geo %>% filter(measures == "mRNA") %>% group_by(target) %>%
    summarise(enzyme_resistance = paste(sprintf("%s:%+.2f(adjP %.2g)", series, log2FC, adj_p), collapse = "; "), .groups = "drop") %>%
    rename(target = target)
  master <- master %>% left_join(mir_fc, by = "mirna") %>% left_join(enz_fc, by = "target")
}

## tiering (same rule as the manuscript)
validated <- l1$validated_in[match(master$mirna, l1$mirna)] != ""
master$tier <- ifelse(master$inverse_sig & validated, "Tier 1",
               ifelse(master$inverse_sig & !validated, "Tier 2", "Tier 3"))
master <- master %>% arrange(factor(tier), desc(abs(rho)))
readr::write_csv(master, "results/TABLE2_integrated_master.csv")
message("Wrote results/TABLE2_integrated_master.csv")

## ---------- Figure 2A: correlation heatmap ----------
hm <- l2
hm$label <- factor(hm$mirna, levels = rev(sort(unique(hm$mirna))))
hm$target <- factor(hm$target, levels = c("OGT", "OGA", "EOGT"))
hm$star <- ifelse(hm$inverse_sig, "*", "")
pA <- ggplot(hm, aes(target, label, fill = rho)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f%s", rho, star)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Spearman rho") +
  labs(title = sprintf("A  miRNA-to-enzyme correlation (TCGA-BRCA, n=%d)", l2$n[1]),
       subtitle = "* FDR < 0.05 and rho < 0", x = NULL, y = NULL) +
  theme_minimal(base_size = 11) + theme(panel.grid = element_blank())

## ---------- Figure 2B: Cox HR forest ----------
fb <- s3
fb$kind <- ifelse(fb$feature %in% GENES | grepl("\\(", fb$feature), "enzyme", "miRNA")
fb <- fb %>% arrange(kind, HR)
fb$feature <- factor(fb$feature, levels = fb$feature)
fb$nominal <- fb$HR_p < 0.05
pB <- ggplot(fb, aes(HR, feature, color = nominal)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey50") +
  geom_point(size = 2) +
  scale_color_manual(values = c("FALSE" = "grey40", "TRUE" = "#B2182B"),
                     labels = c("ns", "nominal p < 0.05"), name = NULL) +
  scale_x_log10() +
  labs(title = "B  Overall-survival hazard ratio (Cox)", x = "Hazard ratio", y = NULL) +
  theme_minimal(base_size = 11)

fig <- pA / pB + plot_layout(heights = c(1, 1.2))
ggsave("figures/Figure2_combined.pdf", fig, width = 8.5, height = 11)
ggsave("figures/Figure2_combined.png", fig, width = 8.5, height = 11, dpi = 300)
ggsave("figures/Figure2A_correlation_heatmap.png", pA, width = 7, height = 5, dpi = 300)
ggsave("figures/Figure2B_HR_forest.png", pB, width = 7, height = 6, dpi = 300)

writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
message("Wrote figures/ and results/sessionInfo.txt")
message("Compare results/TABLE2_integrated_master.csv against the manuscript before publishing.")
