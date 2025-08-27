##################################################
#
# Martínez et al. Evidence for a second population of M. polymorpha
# v1: Verbania, 21/08/2025
#
##################################################


### Libraries -----------
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(scatterpie)
library(patchwork)
library(readxl)
library(ggforce)
library(tibble)
library(scales)
library(adegenet)
library(hierfstat)
library(HardyWeinberg)



### Helper functions -----------

normalize_locality <- function(x){
  x <- trimws(gsub("\\s+"," ", x))
  dplyr::recode(x,
                "Túnel de la Atlántida" = "Tunel de la Atlantida",
                "Túnel de la Atlantida" = "Tunel de la Atlantida",
                "Tunel de la Atlántida" = "Tunel de la Atlantida",
                .default = x
  )
}

# combine two allele columns into one "a/b" diploid genotype; NA if either missing
ab_geno <- function(df, locus){
  a1 <- suppressWarnings(as.integer(df[[paste0(locus,".1")]]))
  a2 <- suppressWarnings(as.integer(df[[paste0(locus,".2")]]))
  ifelse(is.na(a1) | is.na(a2), NA_character_, paste(a1, a2, sep = "/"))
}


# matrix -> long with ordered factors
mat_to_long <- function(M, nm, levels){
  as.data.frame(M) |>
    rownames_to_column("Row") |>
    pivot_longer(-Row, names_to = "Col", values_to = all_of(nm)) |>
    mutate(Row = factor(Row, levels = levels),
           Col = factor(Col, levels = levels))
}


# canonicalize "a/b" (numeric alleles) -> "min/max"
canon_gt <- function(x, sep = "/") {
  x <- as.character(x)
  res <- rep(NA_character_, length(x))
  ok <- !is.na(x)
  parts <- strsplit(x[ok], split = sep, fixed = TRUE)
  res[ok] <- vapply(parts, function(p) {
    if (length(p) != 2) return(NA_character_)
    a <- suppressWarnings(as.numeric(trimws(p[1])))
    b <- suppressWarnings(as.numeric(trimws(p[2])))
    if (is.na(a) || is.na(b)) return(NA_character_)
    paste(sort(c(a,b)), collapse = sep)
  }, character(1))
  res
}

# build HardyWeinberg-style genotype counts vector for one locus
gt_to_HWcounts <- function(gt_char, sep = "/") {
  gt <- canon_gt(gt_char, sep = sep)
  gt <- gt[!is.na(gt)]
  n  <- length(gt)
  if (n < 5) return(NULL)                # too few to test
  
  alle <- suppressWarnings(as.numeric(unique(unlist(strsplit(gt, sep)))))
  alle <- sort(alle[!is.na(alle)])
  if (length(alle) < 2) return(NULL)     # monomorphic/undetermined
  
  # triangular genotype order
  labs <- character()
  for (i in seq_along(alle)) for (j in i:length(alle)) labs <- c(labs, paste(alle[i], alle[j], sep = sep))
  counts <- table(factor(gt, levels = labs))
  v <- as.numeric(counts)
  if (sum(v) < 5) return(NULL)
  v
}

# LD screen p-value (Monte-Carlo χ² on unordered genotypes)
ld_mc_p <- function(gt1_char, gt2_char, B = 5000) {
  g1 <- canon_gt(gt1_char); g2 <- canon_gt(gt2_char)
  ok <- !(is.na(g1) | is.na(g2)); g1 <- g1[ok]; g2 <- g2[ok]
  if (length(g1) < 5) return(NA_real_)
  tab <- table(g1, g2)
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = B)$p.value)
}


### Globals (palette, CRS, factor order) -----------
pal <- list(
  orange    = "#D55E00",
  grey_lo   = "#F2F2F2",
  grey_hi   = "#4D4D4D",
  grey_land = "grey97",
  grey_line = "grey80",
  grey_coast= "grey60"
)
crs_utm28  <- 32628

loc_levels   <- c("Cueva de los Lagos","Tunel de la Atlantida","Jameos del Agua","Charcos de Luis")
charcos_name <- "Charcos de Luis"

theme_set(theme_minimal(base_size = 10))


### Working directory -----------

setwd("/Users/amartinez/Dropbox/_Papers/_MunidopsisMorlockia protection/__MS_Microsatellites/")

######### Part 0: data preparation -----------

geno     <- read_xlsx("Data report/Table_Alelles_Supplementary_table.xlsx")
metadata <- read_xlsx("Data report/samples_Microsatellites_300517.xlsx")


# Clean genotype column names
colnames(geno) <- c("id",
                    "Mp-8.1","Mp-8.2","Mp-1.1","Mp-1.2","Mp-4.1","Mp-4.2",
                    "Mp-5.1","Mp-5.2","Mp-2.1","Mp-2.2","Mp-3.1","Mp-3.2",
                    "Mp-6.1","Mp-6.2","Mp-7.1","Mp-7.2"
)

geno <- merge(geno, metadata[c("id","Locality")])
geno$Locality <- normalize_locality(geno$Locality)

locus_cols <- grep("^Mp-", names(geno), value = TRUE)
loci <- sort(unique(sub("\\.\\d+$","", locus_cols)))

# long table of alleles (one row per allele call)
alle_long <- geno |>
  tidyr::pivot_longer(all_of(locus_cols), names_to = "Locus", values_to = "Allele") |>
  dplyr::mutate(
    Locus    = sub("\\.\\d+$","", Locus),
    Allele   = suppressWarnings(as.integer(Allele))
  ) |>
  dplyr::filter(!is.na(Allele), !is.na(Locality))


# Build "a/b" genotype columns and genind (needed for QC + later diversity)
geno_wide <- geno |>
  dplyr::transmute(id, Locality, !!!setNames(lapply(loci, \(L) ab_geno(geno, L)), loci))

gi <- adegenet::df2genind(geno_wide[, loci],
                          ploidy = 2,
                          ind.names = geno_wide$id,
                          pop = as.factor(geno_wide$Locality),
                          sep = "/"
)


######### Part 1: Genotyping quality control (QC) -----------


## 1.1 Call rate per site × locus ----------------
X_all <- adegenet::tab(gi, NA.method = "asis")
loc_names <- adegenet::locNames(gi)
site_fac  <- adegenet::pop(gi)

missing_geno_by_locus <- lapply(loc_names, function(L){
  cols <- grep(paste0("^", L, "\\."), colnames(X_all))
  if (length(cols) == 0) return(rep(TRUE, nrow(X_all)))
  apply(X_all[, cols, drop = FALSE], 1, function(r) all(is.na(r)))
})
names(missing_geno_by_locus) <- loc_names

qc_callrate <- lapply(loc_names, function(L){
  miss_vec <- missing_geno_by_locus[[L]]
  tibble(Locality = as.character(site_fac), Missing = miss_vec) |>
    dplyr::group_by(Locality) |>
    dplyr::summarise(
      N_attempted = dplyr::n(),
      N_genotyped = sum(!Missing),
      CallRate = round(100*N_genotyped/N_attempted, 1),
      .groups = "drop"
    ) |>
    dplyr::mutate(Locus = L, .before = 1)
}) |>
  dplyr::bind_rows() |>
  dplyr::mutate(Locality = factor(Locality, levels = loc_levels)) |>
  dplyr::arrange(Locality, Locus)


## 1.2 HWE exact tests per locus × site ----------
hw_min_n <- 5

qc_hwe <- lapply(split(geno_wide, geno_wide$Locality), function(dd){
  loci_here <- setdiff(names(dd), c("id","Locality"))
  out <- lapply(loci_here, function(L){
    cnt <- gt_to_HWcounts(dd[[L]])
    if (is.null(cnt)) {
      tibble::tibble(Locus = L, HWE_p = NA_real_, HWE_note = "n<5 or monomorphic/undetermined")
    } else {
      pval <- tryCatch({ HardyWeinberg::HWExact(cnt, verbose = FALSE)$pval }, error = function(e) NA_real_)
      if (is.na(pval)) {
        pval <- tryCatch({ HardyWeinberg::HWChisq(cnt, cc = 0, verbose = FALSE)$pval }, error = function(e) NA_real_)
      }
      tibble::tibble(Locus = L,
                     HWE_p = ifelse(is.finite(pval), round(pval, 4), NA_real_),
                     HWE_note = ifelse(is.na(pval), "test failed", "")
      )
    }
  }) |> dplyr::bind_rows()
  out$Locality <- unique(dd$Locality)
  out
}) |> dplyr::bind_rows() |>
  dplyr::left_join(qc_callrate |> dplyr::select(Locality, Locus, N_genotyped),
                   by = c("Locality","Locus")) |>
  dplyr::mutate(
    HWE_note = dplyr::case_when(
      is.na(HWE_p) & (is.na(N_genotyped) | N_genotyped < hw_min_n) ~ "n<5; HWE not assessed",
      TRUE ~ HWE_note
    )
  ) |>
  dplyr::arrange(factor(Locality, levels = loc_levels), Locus)

readr::write_csv(qc_hwe, "HWE_by_site.csv")

## 1.3 LD screen per site × locus pair ----------
ld_list <- lapply(split(geno_wide, geno_wide$Locality), function(dd){
  loci_here <- setdiff(names(dd), c("id","Locality"))
  if (length(loci_here) < 2) return(NULL)
  pairs <- t(combn(loci_here, 2))
  out <- apply(pairs, 1, function(p){
    c(Locus1 = p[1], Locus2 = p[2], p = ld_mc_p(dd[[p[1]]], dd[[p[2]]], B = 5000))
  })
  as.data.frame(t(out), stringsAsFactors = FALSE) |>
    dplyr::mutate(Locality = unique(dd$Locality)) |>
    dplyr::select(Locality, Locus1, Locus2, p)
}) |> dplyr::bind_rows()

if (nrow(ld_list) > 0) {
  qc_ld_pairs <- ld_list |>
    dplyr::mutate(p = suppressWarnings(as.numeric(p)),
                  p_adj = ifelse(is.na(p), NA_real_, p.adjust(p, method = "BH"))) |>
    dplyr::arrange(factor(Locality, levels = loc_levels), Locus1, Locus2)
  readr::write_csv(qc_ld_pairs, "Table_S3_LD_pairs_by_site.csv")
  
  ld_sig <- qc_ld_pairs |>
    dplyr::filter(!is.na(p_adj), p_adj < 0.05) |>
    dplyr::mutate(pair = paste(Locus1, "×", Locus2)) |>
    tidyr::pivot_longer(c(Locus1, Locus2), names_to = "which", values_to = "Locus") |>
    dplyr::group_by(Locality, Locus) |>
    dplyr::summarise(`LD pairs (BH<0.05)` = paste(pair, collapse = "; "), .groups = "drop")
} else {
  ld_sig <- tibble::tibble(Locality = character(), Locus = character(), `LD pairs (BH<0.05)` = character())
}

## 1.4 Build compact QC table (Table S1) --------
qc_table <- qc_callrate |>
  dplyr::left_join(qc_hwe |> dplyr::select(Locality, Locus, HWE_p, HWE_note),
                   by = c("Locality","Locus")) |>
  dplyr::left_join(ld_sig, by = c("Locality","Locus")) |>
  dplyr::mutate(
    `LD pairs (BH<0.05)` = dplyr::if_else(is.na(`LD pairs (BH<0.05)`) | `LD pairs (BH<0.05)` == "", "none", `LD pairs (BH<0.05)`),
    HWE_p = ifelse(is.na(HWE_p), NA, round(HWE_p, 4)),
    HWE_category = dplyr::case_when(
      !is.na(HWE_p) ~ "HWE tested",
      grepl("n<5", HWE_note, ignore.case = TRUE) ~ "HWE n<5",
      grepl("monomorphic", HWE_note, ignore.case = TRUE) ~ "HWE monomorphic",
      grepl("test failed", HWE_note, ignore.case = TRUE) ~ "HWE sparse table",
      TRUE ~ "HWE not assessed"
    ),
    LD_category = dplyr::if_else(`LD pairs (BH<0.05)` == "none", "no significant pairs", "≥1 significant pair")
  ) |>
  dplyr::select(Locus, Locality, N_attempted, N_genotyped, CallRate, HWE_p, HWE_category, HWE_note, `LD pairs (BH<0.05)`, LD_category) |>
  dplyr::arrange(factor(Locality, levels = loc_levels), Locus)

readr::write_csv(qc_table, "Table_S1_QC_by_site_locus.csv")



######## Part 2: Diversity indices & between-site stats -----------

hf <- genind2hierfstat(gi)
bs <- hierfstat::basic.stats(hf)

Ho_pop  <- colMeans(bs$Ho, na.rm = TRUE)
He_pop  <- colMeans(bs$Hs, na.rm = TRUE)
Fis_pop <- colMeans(bs$Fis, na.rm = TRUE)

ar   <- hierfstat::allelic.richness(hf)
AR_pop <- colMeans(ar$Ar, na.rm = TRUE)

NA_total_true <- alle_long |>
  dplyr::distinct(Locality, Locus, Allele) |>
  dplyr::count(Locality, name = "NA_total")

div_table <- data.frame(
  Locality = names(Ho_pop),
  N        = as.integer(table(adegenet::pop(gi))[names(Ho_pop)]),
  AR       = AR_pop[names(Ho_pop)],
  Ho       = Ho_pop,
  He       = He_pop,
  FIS      = Fis_pop,
  row.names = NULL
) |>
  dplyr::left_join(NA_total_true, by = "Locality") |>
  dplyr::mutate(Locality = factor(Locality, levels = loc_levels)) |>
  dplyr::arrange(Locality)

# Pairwise FST + 95% CI
pair_fst <- hierfstat::pairwise.WCfst(hf)
set.seed(99)
boot     <- hierfstat::boot.ppfst(dat = hf, nboot = 1000)
fst_mat  <- as.matrix(pair_fst)
ll_mat   <- boot$ll
ul_mat   <- boot$ul

# Charcos vs others numbers for prose (if present)
if (charcos_name %in% rownames(fst_mat)) {
  other_sites <- setdiff(rownames(fst_mat), charcos_name)
  fst_prose <- tibble::tibble(
    Contrast = paste(charcos_name, "vs", other_sites),
    FST      = round(as.numeric(fst_mat[charcos_name, other_sites]), 3),
    CI_low   = round(as.numeric(ll_mat[charcos_name, other_sites]), 3),
    CI_high  = round(as.numeric(ul_mat[charcos_name, other_sites]), 3)
  )
  readr::write_csv(fst_prose, "FST_Charcos_vs_others.csv")
}

# shared-allele counts (presence/absence per locus-allele)
M <- alle_long |>
  dplyr::distinct(Locality, Locus, Allele) |>
  tidyr::unite(locus_allele, Locus, Allele) |>
  dplyr::mutate(v = 1L) |>
  tidyr::pivot_wider(names_from = locus_allele, values_from = v, values_fill = 0L) |>
  dplyr::arrange(Locality) |>
  tibble::column_to_rownames("Locality") |>
  as.matrix()

shared_mat <- tcrossprod(M)



######## Part 3: Figures -----------


## 3.1 Map inputs --------------------------------
site_coords <- tibble::tribble(
  ~Locality,               ~lon,     ~lat,
  "Cueva de los Lagos",     -13.439,  29.157,
  "Tunel de la Atlantida",  -13.429,  29.155,
  "Jameos del Agua",        -13.433,  29.157,
  "Charcos de Luis",        -13.423,  29.203
)

# Prepare per-site summaries used in labels
n_by_site <- geno |>
  dplyr::distinct(id, Locality) |>
  dplyr::count(Locality, name = "N")

ar_by_site <- alle_long |>
  dplyr::group_by(Locality, Locus) |>
  dplyr::summarise(k = dplyr::n_distinct(Allele), .groups = "drop") |>
  dplyr::group_by(Locality) |>
  dplyr::summarise(mean_AR = mean(k, na.rm = TRUE), .groups = "drop")

priv_tbl <- alle_long |>
  dplyr::distinct(Locus, Allele, Locality) |>
  dplyr::add_count(Locus, Allele, name = "n_sites") |>
  dplyr::filter(n_sites == 1) |>
  dplyr::count(Locality, name = "n_private")

shared_tbl <- alle_long |>
  dplyr::distinct(Locality, Locus, Allele) |>
  dplyr::count(Locality, name = "n_alleles_total") |>
  dplyr::left_join(priv_tbl, by = "Locality") |>
  dplyr::mutate(n_private = tidyr::replace_na(n_private, 0L),
                n_shared  = pmax(n_alleles_total - n_private, 0L))

site_sum <- n_by_site |>
  dplyr::full_join(ar_by_site, by = "Locality") |>
  dplyr::full_join(shared_tbl, by = "Locality") |>
  dplyr::mutate(Locality = factor(Locality, levels = loc_levels)) |>
  dplyr::arrange(Locality) |>
  dplyr::left_join(site_coords, by = "Locality") |>
  dplyr::mutate(across(c(n_private,n_shared,n_alleles_total,N,mean_AR), as.numeric))

sites_sf <- site_sum |>
  sf::st_as_sf(coords = c("lon","lat"), crs = 4326, remove = FALSE) |>
  sf::st_transform(crs_utm28) |>
  dplyr::mutate(
    x   = sf::st_coordinates(geometry)[,1],
    y   = sf::st_coordinates(geometry)[,2],
    r_m = scales::rescale(n_alleles_total, to = c(500,1100))
  ) |>
  sf::st_drop_geometry()

fan <- tibble::tribble(
  ~Locality,               ~angle_deg, ~k_pie, ~k_lab, ~hjust,
  "Cueva de los Lagos",         205,       7.0,   9.2,    1,
  "Tunel de la Atlantida",      230,       7.8,   9.8,    1,
  "Jameos del Agua",            255,       7.0,   9.0,    1,
  "Charcos de Luis",            110,       4.6,   6.4,    0
)

sites_sf <- sites_sf |>
  dplyr::left_join(fan, by = "Locality") |>
  dplyr::mutate(
    angle_rad = angle_deg*pi/180,
    x_pie = x + k_pie*r_m*cos(angle_rad),
    y_pie = y + k_pie*r_m*sin(angle_rad),
    x_lab = x + k_lab*r_m*cos(angle_rad),
    y_lab = y + k_lab*r_m*sin(angle_rad),
    label_txt = sprintf("%s\nN=%d | AR=%.2f | Ho=%.2f | He=%.2f | FIS=%.2f",
                        Locality,
                        div_table$N[match(Locality, div_table$Locality)],
                        div_table$AR[match(Locality, div_table$Locality)],
                        div_table$Ho[match(Locality, div_table$Locality)],
                        div_table$He[match(Locality, div_table$Locality)],
                        div_table$FIS[match(Locality, div_table$Locality)])
  )

# Basemap
world <- rnaturalearth::ne_download(scale = 10, category = "cultural", type = "admin_0_countries", returnclass = "sf") |>
  sf::st_make_valid()
sites_ll <- site_sum |>
  sf::st_as_sf(coords = c("lon","lat"), crs = 4326, remove = FALSE)
bbox_ll <- sites_ll |>
  sf::st_union() |> sf::st_transform(crs_utm28) |> sf::st_buffer(15000) |>
  sf::st_transform(4326) |> sf::st_bbox() |> sf::st_as_sfc()

old_s2 <- sf::sf_use_s2(FALSE)
lanza  <- sf::st_crop(world, bbox_ll)
coast  <- rnaturalearth::ne_download(scale = 10, category = "physical", type = "coastline", returnclass = "sf") |>
  sf::st_make_valid() |>
  sf::st_filter(bbox_ll, .predicate = sf::st_intersects) |>
  sf::st_intersection(bbox_ll)
sf::sf_use_s2(old_s2)

lanza_utm <- sf::st_transform(lanza,  crs_utm28)
coast_utm <- sf::st_transform(coast,  crs_utm28)

## 3.2 Figure panels -----------------------------

# (A) map
p_map <- ggplot() +
  geom_sf(data = lanza_utm, fill = pal$grey_land, color = pal$grey_line, linewidth = 0.25) +
  geom_sf(data = coast_utm, color = pal$grey_coast, linewidth = 0.25) +
  geom_curve(data = sites_sf,
             aes(x = x, y = y, xend = x_pie, yend = y_pie),
             curvature = -0.15, color = "grey35", linewidth = 0.45,
             arrow = grid::arrow(length = unit(2.2,"mm"), type = "closed")) +
  ggforce::geom_circle(data = sites_sf,
                       aes(x0 = x_pie, y0 = y_pie, r = r_m*1.10),
                       inherit.aes = FALSE, fill = "white", color = NA, alpha = 0.85) +
  scatterpie::geom_scatterpie(
    data = sites_sf, aes(x = x_pie, y = y_pie, r = r_m),
    cols = c("n_private","n_shared"),
    color = "white", linewidth = 0.25, alpha = 0.95) +
  geom_curve(data = sites_sf,
             aes(x = x_pie, y = y_pie, xend = x_lab, yend = y_lab),
             curvature = 0.10, color = "grey35", linewidth = 0.35) +
  geom_label(data = sites_sf,
             aes(x = x_lab, y = y_lab, label = label_txt, hjust = hjust),
             size = 3, label.size = 0, label.r = unit(0.15,"lines"),
             fill = "white", alpha = 0.9) +
  geom_point(data = sites_sf, aes(x, y),
             shape = 21, fill = "white", color = "black", size = 1.7, stroke = 0.3) +
  scale_fill_manual(values = c(n_private = pal$orange, n_shared = "#999999"),
                    breaks = c("n_private","n_shared"),
                    labels = c("Private alleles","Shared alleles"), name = NULL) +
  coord_sf(crs = crs_utm28, expand = FALSE) +
  theme(legend.position = "bottom") +
  labs(title = "A. Map with private/shared allele pies and diversity labels")

# (B) Ho vs He dumbbell
div_long <- div_table |>
  dplyr::select(Locality, Ho, He, AR, N, FIS) |>
  dplyr::mutate(Locality = factor(Locality, levels = loc_levels)) |>
  tidyr::pivot_longer(c(Ho, He), names_to = "Metric", values_to = "Value")
seg_B <- tidyr::pivot_wider(div_long, names_from = Metric, values_from = Value)
x_max <- max(div_long$Value, na.rm = TRUE)
x_lab <- x_max + 0.10

p_B <- ggplot() +
  geom_segment(data = seg_B,
               aes(y = Locality, yend = Locality, x = He, xend = Ho),
               color = "grey55", linewidth = 0.6) +
  geom_point(data = div_long,
             aes(x = Value, y = Locality, color = Metric), size = 2.2) +
  scale_color_manual(values = c(He = pal$grey_hi, Ho = pal$orange), name = NULL) +
  scale_x_continuous(limits = c(0, x_lab),
                     expand = expansion(mult = c(0.02, 0.02)),
                     name = "Heterozygosity (He vs Ho)") +
  geom_text(data = div_table,
            aes(y = factor(Locality, levels = loc_levels),
                x = x_lab,
                label = sprintf("AR=%.2f | N=%d | FIS=%.2f", AR, N, FIS)),
            hjust = 1, size = 3, color = "grey20") +
  coord_cartesian(clip = "off") +
  theme(legend.position = "bottom",
        panel.grid.minor = element_blank(),
        plot.margin = margin(5.5, 30, 5.5, 5.5, "pt")) +
  labs(title = "B. Within-site diversity (Ho vs He); AR shown at right label", y = NULL)

# (C) FST lower triangle + shared alleles upper
fst_long   <- mat_to_long(pmax(fst_mat, 0), "FST",    loc_levels)
ll_long    <- mat_to_long(ll_mat,               "ll",  loc_levels)
ul_long    <- mat_to_long(ul_mat,               "ul",  loc_levels)
shared_long<- mat_to_long(shared_mat,       "Shared",  loc_levels)

grid_C <- fst_long |>
  dplyr::left_join(ll_long, by = c("Row","Col")) |>
  dplyr::left_join(ul_long, by = c("Row","Col")) |>
  dplyr::left_join(shared_long, by = c("Row","Col")) |>
  dplyr::mutate(i = as.integer(Row), j = as.integer(Col),
                in_lower = i > j, in_upper = i < j,
                FST_plot    = dplyr::if_else(in_lower, FST, NA_real_),
                Shared_plot = dplyr::if_else(in_upper, Shared, NA_integer_),
                sig_plot    = in_lower & !is.na(ll) & ll > 0)

p_C <- ggplot(grid_C, aes(Col, Row)) +
  geom_tile(aes(fill = FST_plot), colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(!is.na(FST_plot), sprintf("%.2f", FST_plot), "")),
            size = 3, colour = "black") +
  geom_point(data = dplyr::filter(grid_C, sig_plot),
             shape = 21, size = 2.2, stroke = 0.9,
             fill = NA, colour = pal$orange) +
  scale_fill_gradient(limits = c(0, 0.6), na.value = "white",
                      low = pal$grey_lo, high = pal$grey_hi, name = "FST (WC84)") +
  geom_tile(data = dplyr::filter(grid_C, in_upper),
            fill = "white", colour = "white", linewidth = 0.5) +
  geom_text(data = dplyr::filter(grid_C, in_upper),
            aes(label = Shared_plot), size = 3.2, fontface = "bold") +
  coord_equal() +
  labs(x = NULL, y = NULL,
       title = "C. Pairwise: lower = FST (ring if 95% CI > 0), upper = # shared alleles") +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))

# Compose
final_plot <- p_map | (p_B / p_C) + patchwork::plot_layout(widths = c(1.25, 1))
final_plot


######## Part 4: Main Table 1: Site-level genetic summary -----------

# Site-level QC summary from the already-built qc_table
site_qc <- qc_table %>%
  group_by(Locality) %>%
  summarise(
    N                 = suppressWarnings(max(N_attempted, na.rm = TRUE)),   # same across loci
    CallRate_mean     = round(mean(CallRate, na.rm = TRUE), 1),
    HWE_tested        = sum(HWE_category == "HWE tested", na.rm = TRUE),
    LD_sig_pairs_n    = sum(LD_category == "≥1 significant pair", na.rm = TRUE),
    .groups = "drop"
  )

# Private & total alleles per site (you already computed shared_tbl
alle_site <- shared_tbl %>%
  transmute(Locality, Private_alleles = n_private, Total_alleles = n_alleles_total)

# Charcos-vs-others FST with CI (string), "-" for Charcos row
fst_vs_charcos <- NULL
if (exists("fst_mat") && charcos_name %in% rownames(fst_mat)) {
  # Pull column for Charcos (or row; symmetric anyway)
  # We’ll format as "0.23 [0.08–0.38]"
  sites <- rownames(fst_mat)
  fst_col <- fst_mat[, charcos_name]
  ci_lo   <- ll_mat[, charcos_name]
  ci_hi   <- ul_mat[, charcos_name]
  fst_vs_charcos <- tibble::tibble(
    Locality = sites,
    `FST vs Charcos (95% CI)` = dplyr::if_else(
      sites == charcos_name,
      "-",
      sprintf("%.2f [%.2f–%.2f]", pmax(0, as.numeric(fst_col)), as.numeric(ci_lo), as.numeric(ci_hi))
    )
  )
}

# Assemble main table from div_table + QC + alleles + FST
tbl <- div_table %>%
  dplyr::select(Locality, N, AR, Ho, He, FIS) %>%
  dplyr::left_join(site_qc %>% dplyr::select(Locality, CallRate_mean, HWE_tested, LD_sig_pairs_n),
                   by = "Locality") %>%
  dplyr::left_join(alle_site %>% dplyr::select(Locality, Private_alleles),
                   by = "Locality")


# Add FST vs Charcos column safely (NA if not available)
if (exists("fst_vs_charcos") && !is.null(fst_vs_charcos)) {
  tbl <- tbl %>%
    dplyr::left_join(
      fst_vs_charcos %>%
        dplyr::rename(FST_vs_Charcos_95CI = `FST vs Charcos (95% CI)`),
      by = "Locality"
    )
} else {
  tbl <- dplyr::mutate(tbl, FST_vs_Charcos_95CI = NA_character_)
}

# Format, round, and order columns
Main_Table1 <- tbl %>%
  dplyr::mutate(
    `Ho/He`        = sprintf("%.2f/%.2f", Ho, He),
    AR             = round(AR, 2),
    FIS            = round(FIS, 2),
    CallRate_mean  = round(CallRate_mean, 1),
    Private_alleles = dplyr::coalesce(Private_alleles, 0L),
    Locality = factor(Locality, levels = loc_levels)
  ) %>%
  dplyr::select(
    Locality, N, AR, `Ho/He`, FIS,
    Private_alleles, CallRate_mean, HWE_tested, LD_sig_pairs_n,
    FST_vs_Charcos_95CI
  ) %>%
  dplyr::arrange(Locality) %>%
  dplyr::mutate(Locality = as.character(Locality)) %>%
  dplyr::rename(`FST vs Charcos (95% CI)` = FST_vs_Charcos_95CI)

# Write outputs
readr::write_csv(Main_Table1, "Table1_SiteSummary.csv")
