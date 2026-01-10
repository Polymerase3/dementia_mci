# load required packages
library(phiper)
library(rlang)
library(ggplot2)
library(Cairo)
library(openxlsx)
library(dplyr)
library(htmlwidgets)

set.seed(632961)

# create PHIP data object from input data with reproducible seed
withr::with_preserve_seed({
  ps <- phip_convert(
    data_long_path    = "data_original/other_data/dementia_full.parquet",
    sample_id         = "sampleID",
    peptide_id        = "peptideID",
    exist             = "exist",
    fold_change       = "fold_change",
    counts_input      = "input",
    counts_hit        = "count",
    backend           = "duckdb",
    peptide_library   = TRUE,
    materialise_table = TRUE,
    auto_expand       = TRUE,
    n_cores           = 10
  )
})

# create base results directory and save peptide library
out_dir <- file.path("results", "type_person")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
get_peptide_library(ps) %>%
  collect() %>%
  as.data.frame() %>%
  saveRDS(file = file.path(out_dir, "peptide_library.rds"))

# columns always to keep in saved data
base_cols <- c("sample_id", "peptide_id", "type_person", "fold_change", "exist")

# helper function to safely save an RDS file after creating directories
save_rds_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, path)
}

# helper function: select base columns (and any extras), drop NA rows, collect and save to RDS
make_and_save <- function(data, out_path, extra_vars = NULL, drop_na = TRUE) {
  cols_to_select <- unique(c(base_cols, extra_vars %||% character(0)))
  avail <- names(data)
  missing <- setdiff(cols_to_select, avail)
  if (length(missing)) {
    message("Skipping missing columns: ", paste(missing, collapse = ", "))
  }
  df <- data %>% dplyr::select(dplyr::any_of(cols_to_select))
  if (isTRUE(drop_na)) {
    df <- df %>% dplyr::filter(dplyr::if_all(dplyr::all_of(intersect(cols_to_select, names(.))), ~ !is.na(.)))
  }
  df <- df %>% dplyr::collect()
  save_rds_safe(df, out_path)
  invisible(df)
}

# filter to samples with type_person and create a group column
ps_cmp <- ps %>%
  filter(!is.na(type_person)) %>%
  mutate(
    group_char = dplyr::case_when(
      type_person == "MCI"      ~ "a_MCI",
      type_person == "Dementia" ~ "b_Dementia",
      TRUE                      ~ as.character(type_person)
    )
  )

# save the filtered dataset for this comparison
make_and_save(
  data     = ps_cmp,
  out_path = file.path(out_dir, "type_person_data.rds")
)

# ------------------ enrichment counts ------------------
CairoSVG(file.path(out_dir, "enrichment_counts.svg"), dpi = 300,
         height = 30, width = 40, unit = "cm", bg = "white")
p_enrich <- plot_enrichment_counts(ps_cmp, group_cols = "group_char",
                                   prevalence_threshold = 0.07) +
  theme(text = element_text(family = "DejaVu Sans"))
print(p_enrich)
dev.off()

# ------------------ alpha diversity --------------------
alpha_div <- compute_alpha_diversity(ps_cmp, group_cols = "group_char", carry_cols = c("sex", "age"))
dir.create(file.path(out_dir, "alpha_diversity"), recursive = TRUE, showWarnings = FALSE)
write.xlsx(alpha_div, file.path(out_dir, "alpha_diversity", "table.xlsx"))

CairoSVG(file.path(out_dir, "alpha_diversity", "richness.svg"), dpi = 300,
         height = 30, width = 40, unit = "cm", bg = "white")
p_alpha <- plot_alpha_diversity(alpha_div, metric = "richness", group_col = "group_char")
print(p_alpha)
dev.off()

CairoSVG(file.path(out_dir, "alpha_diversity", "shannon.svg"), dpi = 300,
         height = 30, width = 40, unit = "cm", bg = "white")
p_alpha <- plot_alpha_diversity(alpha_div, metric = "shannon", group_col = "group_char")
print(p_alpha)
dev.off()

# ------------------ beta diversity ---------------------
dist_bc <- phiper:::compute_distance(ps_cmp, value_col = "exist",
                                     method_normalization = "hellinger",
                                     distance = "bray", n_threads = 10)
dir.create(file.path(out_dir, "beta_diversity"), recursive = TRUE, showWarnings = FALSE)
dist_mat <- as.matrix(dist_bc)
openxlsx::write.xlsx(dist_mat, file = file.path(out_dir, "beta_diversity", "distance_matrix.xlsx"), rowNames = TRUE)

pcoa_res <- phiper:::compute_pcoa(dist_bc, neg_correction = "none", n_axes = 109, top_features = 30)
saveRDS(pcoa_res, file.path(out_dir, "beta_diversity", "pcoa_results.rds"))

cap_res <- phiper:::compute_capscale(dist_bc, ps = ps_cmp, formula = ~ group_char)
saveRDS(cap_res, file.path(out_dir, "beta_diversity", "capscale_results.rds"))

permanova_res <- phiper:::compute_permanova(dist_bc, ps = ps_cmp, group_col = "group_char", contrasts = "pairwise")
saveRDS(permanova_res, file.path(out_dir, "beta_diversity", "permanova_results.rds"))

# dispersion uses pairwise contrasts; keep the same output format
disp_res <- phiper:::compute_dispersion(dist_bc, ps = ps_cmp, group_col = "group_char", contrasts = "pairwise")
saveRDS(disp_res, file.path(out_dir, "beta_diversity", "dispersion_results.rds"))
print(disp_res)

# t-SNE 2d
tsne_res <- phiper:::compute_tsne(ps = ps_cmp, dist_obj = dist_bc, dims = 2L,
                                   perplexity = 15, meta_cols = c("group_char"))
openxlsx::write.xlsx(tsne_res, file = file.path(out_dir, "beta_diversity", "tsne2d_results.xlsx"), rowNames = TRUE)

CairoSVG(file.path(out_dir, "beta_diversity", "tsne2d_plot.svg"), dpi = 300,
         height = 30, width = 30, unit = "cm", bg = "white")
p_tsne2d <- phiper:::plot_tsne(tsne_res, view = "2d", colour = "group_char")
print(p_tsne2d)
dev.off()

# t-SNE 3d
tsne_res <- phiper:::compute_tsne(ps = ps_cmp, dist_obj = dist_bc, dims = 3L,
                                   perplexity = 20, meta_cols = c("group_char"))
openxlsx::write.xlsx(tsne_res, file = file.path(out_dir, "beta_diversity", "tsne3d_results.xlsx"), rowNames = TRUE)

p3d <- phiper:::plot_tsne(tsne_res, view = "3d", colour = "group_char")
htmlwidgets::saveWidget(p3d, file = file.path(out_dir, "beta_diversity", "tsne3d_plot.html"), selfcontained = TRUE)

# add group information to PCoA sample coordinates
pcoa_res$sample_coords <- pcoa_res$sample_coords %>%
  dplyr::left_join(
    ps_cmp$data_long %>% dplyr::select(sample_id, group_char) %>% dplyr::distinct(),
    by = "sample_id",
    copy = TRUE
  )

# PCoA plot with group centroids and ellipses
CairoSVG(file.path(out_dir, "beta_diversity", "pcoa_plot.svg"), dpi = 300,
         height = 30, width = 30, unit = "cm", bg = "white")
p_pcoa <- phiper:::plot_pcoa(pcoa_res, axes = c(1, 2), group_col = "group_char",
                             ellipse_by = "group", show_centroids = TRUE)
print(p_pcoa)
dev.off()

# scree plot for first 15 axes of PCoA
CairoSVG(file.path(out_dir, "beta_diversity", "scree_plot.svg"), dpi = 300,
         height = 30, width = 30, unit = "cm", bg = "white")
p_scree <- phiper:::plot_scree(pcoa_res, n_axes = 15, type = "line") + theme()
print(p_scree)
dev.off()
