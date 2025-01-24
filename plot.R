library(qtl2)
library(tidyverse)
library(patchwork)

# lod <- read.csv("data/scan.csv", row.names = 1)
# plot_scan1(lod)

scan <- readRDS("data/scan.rds")

png("plots/lod.png", width = 8, height = 12, units = "in", res = 150)
par(mfrow = c(ncol(scan$out), 1), mar = c(1, 1, 1, 1))
for (pheno in 1:ncol(scan$out)) {
    # png(str_glue("plots/lod.{colnames(scan$out)[pheno]}.png"))
    plot(scan$out, scan$pmap, lodcolumn = pheno)
    # dev.off()
    # ggsave(str_glue("plots/lod.{colnames(scan$out)[pheno]}.png"))
}
dev.off()

png("plots/coef.retrofat.1.png", width = 6, height = 6, units = "in", res = 150)
par(mar = c(4, 4, 1, 1))
plot_coef(
    scan$coef[["1"]][["retrofat"]],
    scan$pmap["1"],
    scan1_output = scan$out[, "retrofat", drop = FALSE],
    bgcolor = "gray95",
    col = RColorBrewer::brewer.pal(9, "Set1"),
    legend = "bottomleft",
    ylim = c(-5, 5),
    xlim = c(250, 285),
)
dev.off()

png("plots/coef.bodylength_w_tail.7.png", width = 6, height = 6, units = "in", res = 150)
par(mar = c(4, 4, 1, 1))
plot_coef(
    scan$coef[["7"]][["bodylength_w_tail"]],
    scan$pmap["7"],
    scan1_output = scan$out[, "bodylength_w_tail", drop = FALSE],
    bgcolor = "gray95",
    col = RColorBrewer::brewer.pal(9, "Set1"),
    legend = "bottomleft",
    ylim = c(-10, 10),
    xlim = c(20, 50),
)
dev.off()

loci <- scan$pmap |>
    map_dfr(~enframe(.x, name = "locus", value = "pos"), .id = "chrom") |>
    mutate(chrom = fct_inorder(chrom))

df <- as.data.frame(scan$out) |>
    as_tibble(rownames = "locus") |>
    left_join(loci, by = "locus", relationship = "one-to-one") |>
    pivot_longer(-c(locus, chrom, pos), names_to = "phenotype", values_to = "LOD") |>
    mutate(phenotype = phenotype |>
               fct_relevel("bmi_bodylength_wo_tail", "bmi_bodylength_w_tail",
                           "bodylength_wo_tail", "bodylength_w_tail",
                           "bodyweight", "epifat", "parafat", "retrofat",
                           "fasting_glucose", "tail_length"))

df |>
    ggplot(aes(x = pos, y = LOD, color = phenotype)) +
    facet_grid(cols = vars(chrom), scales = "free_x", space = "free_x", switch = "x") +
    geom_point(size = 0.5) +
    scale_x_continuous(expand = c(0, 10)) +
    scale_color_manual(values = c("#7873B4", "#0373AF", "#E58525", "#FBDF8D", "#1F8650",
                                  "#BC3B30", "#6E9AAE", "#E94E98", "#777777", "#000000")) +
    xlab(NULL) +
    guides(color = guide_legend(override.aes = list(size = NULL))) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(0, "pt"),
          axis.text.x = element_blank())

ggsave("plots/porcupine.png", bg = "white", width = 10, height = 4)

# For comparison with paper figure, highlight top LOD per phenotype per chromosome

top <- df |>
    group_by(phenotype, chrom) |>
    slice_max(LOD, n = 1)

df |>
    ggplot(aes(x = pos, y = LOD)) +
    facet_grid(cols = vars(chrom), scales = "free_x", space = "free_x", switch = "x") +
    geom_point(size = 0.5, color = "#cccccc") +
    geom_point(aes(color = phenotype), data = top, size = 2) +
    scale_x_continuous(expand = c(0, 10)) +
    scale_color_manual(values = c("#7873B4", "#0373AF", "#E58525", "#FBDF8D", "#1F8650",
                                  "#BC3B30", "#6E9AAE", "#E94E98", "#777777", "#000000")) +
    xlab(NULL) +
    guides(color = guide_legend(override.aes = list(size = NULL))) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(0, "pt"),
          axis.text.x = element_blank())

ggsave("plots/porcupine_top.png", bg = "white", width = 10, height = 4)
