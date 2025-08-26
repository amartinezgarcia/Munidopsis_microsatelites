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

loc_levels <- c("Cueva de los Lagos","Tunel de la Atlantida","Jameos del Agua","Charcos de Luis")

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

locus_cols <- grep("^Mp-", names(geno), value = TRUE)

# long table of alleles (one row per allele call)
alle_long <- geno |>
  pivot_longer(all_of(locus_cols), names_to = "Locus", values_to = "Allele") |>
  mutate(
    Locus    = sub("\\.\\d+$","", Locus),
    Allele   = suppressWarnings(as.integer(Allele))
  ) |>
  filter(!is.na(Allele), !is.na(Locality))




######### Part 1: Summary data -----------


### 1.1. Summary of data by site --------------


## Number of individuals genotyped by site
n_by_site <- geno |> distinct(id, Locality) |> count(Locality, name = "N")


## Mean allellic richness by site
ar_by_site <- alle_long |>
  group_by(Locality, Locus) |>
  summarise(k = n_distinct(Allele), .groups = "drop") |>
  group_by(Locality) |>
  summarise(mean_AR = mean(k, na.rm = TRUE), .groups = "drop")

## Number of private (=unique) allelles by site
priv_tbl <- alle_long |>
  distinct(Locus, Allele, Locality) |>
  add_count(Locus, Allele, name = "n_sites") |>
  filter(n_sites == 1) |>
  count(Locality, name = "n_private")

## Number of shared allelles by site
shared_tbl <- alle_long |>
  distinct(Locality, Locus, Allele) |>
  count(Locality, name = "n_alleles_total") |>
  left_join(priv_tbl, by = "Locality") |>
  mutate(n_private = replace_na(n_private, 0L),
         n_shared  = pmax(n_alleles_total - n_private, 0L))


## Summary data per site
site_sum <- n_by_site |>
  full_join(ar_by_site, by = "Locality") |>
  full_join(shared_tbl, by = "Locality") |>
  mutate(Locality = factor(Locality, levels = loc_levels)) |>
  arrange(Locality)


### 1.2. Diversity indices (Ho, He, FIS, AR, FST) --------------

# Build "a/b" genotype columns for genind
loci <- sort(unique(sub("\\.\\d+$","", locus_cols)))

geno_wide <- geno |>
  transmute(id, Locality,
            !!!setNames(lapply(loci, \(L) ab_geno(geno, L)), loci))

gi <- df2genind(geno_wide[, loci],
                ploidy = 2,
                ind.names = geno_wide$id,
                pop = as.factor(geno_wide$Locality),
                sep = "/")

hf <- genind2hierfstat(gi)

bs <- hierfstat::basic.stats(hf)          # Ho, Hs(=He), Fis

Ho_pop  <- colMeans(bs$Ho, na.rm = TRUE)
He_pop  <- colMeans(bs$Hs, na.rm = TRUE)
Fis_pop <- colMeans(bs$Fis, na.rm = TRUE)

# Allelic richness (rarefied to min N across pops)
ar   <- hierfstat::allelic.richness(hf)   # locus × pop
AR_pop <- colMeans(ar$Ar, na.rm = TRUE)

# total distinct alleles per site 
NA_total_true <- alle_long |>
  distinct(Locality, Locus, Allele) |>
  count(Locality, name = "NA_total")

## Summary table
div_table <- data.frame(
  Locality = names(Ho_pop),
  N        = as.integer(table(pop(gi))[names(Ho_pop)]),
  AR       = AR_pop[names(Ho_pop)],
  Ho       = Ho_pop,
  He       = He_pop,
  FIS      = Fis_pop,
  row.names = NULL
) |>
  left_join(NA_total_true, by = "Locality") |>
  mutate(Locality = factor(Locality, levels = loc_levels)) |>
  arrange(Locality)


# Pairwise FST + 95% CI from bootstrap over loci
pair_fst <- hierfstat::pairwise.WCfst(hf)

boot     <- hierfstat::boot.ppfst(dat = hf, nboot = 1000)

fst_mat  <- as.matrix(pair_fst)
ll_mat   <- boot$ll
ul_mat   <- boot$ul


# shared-allele counts (presence/absence per locus-allele)
M <- alle_long |>
  distinct(Locality, Locus, Allele) |>
  unite(locus_allele, Locus, Allele) |>
  mutate(v = 1L) |>
  pivot_wider(names_from = locus_allele, values_from = v, values_fill = 0L) |>
  arrange(Locality) |>
  tibble::column_to_rownames("Locality") |>
  as.matrix()

shared_mat <- tcrossprod(M)




### 1.3. Coordinates for the 4 sites (WGS84) --------------

site_coords <- tribble(
  ~Locality,               ~lon,     ~lat,
  "Cueva de los Lagos",     -13.439,  29.157,
  "Tunel de la Atlantida",  -13.429,  29.155,
  "Jameos del Agua",        -13.433,  29.157,
  "Charcos de Luis",        -13.423,  29.203
)

site_sum <- site_sum |>
  left_join(site_coords, by = "Locality") |>
  mutate(across(c(n_private,n_shared,n_alleles_total,N,mean_AR), as.numeric))


# Site points (UTM) + pie/label offsets by angle
sites_sf <- site_sum |>
  st_as_sf(coords = c("lon","lat"), crs = 4326, remove = FALSE) |>
  st_transform(crs_utm28) |>
  mutate(
    x   = st_coordinates(geometry)[,1],
    y   = st_coordinates(geometry)[,2],
    r_m = rescale(n_alleles_total, to = c(500,1100))
  ) |>
  st_drop_geometry()


fan <- tribble(
  ~Locality,               ~angle_deg, ~k_pie, ~k_lab, ~hjust,
  "Cueva de los Lagos",         205,       7.0,   9.2,    1,
  "Tunel de la Atlantida",      230,       7.8,   9.8,    1,
  "Jameos del Agua",            255,       7.0,   9.0,    1,
  "Charcos de Luis",            110,       4.6,   6.4,    0
)


sites_sf <- sites_sf |>
  left_join(fan, by = "Locality") |>
  mutate(
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


### 1.4. Preparation of the basemap  --------------

# (lots of warning—no worry in principle)

world <- rnaturalearth::ne_download(scale = 10,
                                    category = "cultural", type = "admin_0_countries", returnclass = "sf") |>
  st_make_valid()
sites_ll <- site_sum |>
  st_as_sf(coords = c("lon","lat"), crs = 4326, remove = FALSE)
bbox_ll <- sites_ll |>
  st_union() |> st_transform(crs_utm28) |> st_buffer(15000) |>
  st_transform(4326) |> st_bbox() |> st_as_sfc()

old_s2 <- sf::sf_use_s2(FALSE)
lanza  <- st_crop(world, bbox_ll)
coast  <- rnaturalearth::ne_download(scale = 10, category = "physical",
                                     type = "coastline", returnclass = "sf") |>
  st_make_valid() |>
  st_filter(bbox_ll, .predicate = st_intersects) |>
  st_intersection(bbox_ll)
sf::sf_use_s2(old_s2)

lanza_utm <- st_transform(lanza,  crs_utm28)
coast_utm <- st_transform(coast,  crs_utm28)



######### Part 2: Preparation of the figure -----------

###  (Fig 1a) Map with pies + labels ---------------

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


###  (Fig 1b) Ho vs He dumbbell + AR inset ---------------

div_long <- div_table |>
  dplyr::select(Locality, Ho, He, AR, N, FIS) |>
  mutate(Locality = factor(Locality, levels = loc_levels)) |>
  pivot_longer(c(Ho, He), names_to = "Metric", values_to = "Value")

# dataset for the connecting segments (He -> Ho)
seg_B <- div_long |>
  tidyr::pivot_wider(names_from = Metric, values_from = Value)  # columns He, Ho

# place the text a little to the right of the max heterozygosity
x_max <- max(div_long$Value, na.rm = TRUE)
x_lab <- x_max + 0.10  # how far to the right the label sits

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


###  (Fig 1c) Pairwise shared-alleles heatmap ---------------


fst_long   <- mat_to_long(pmax(fst_mat, 0), "FST",    loc_levels)
ll_long    <- mat_to_long(ll_mat,               "ll",  loc_levels)
ul_long    <- mat_to_long(ul_mat,               "ul",  loc_levels)
shared_long<- mat_to_long(shared_mat,       "Shared",  loc_levels)

grid_C <- fst_long |>
  left_join(ll_long, by = c("Row","Col")) |>
  left_join(ul_long, by = c("Row","Col")) |>
  left_join(shared_long, by = c("Row","Col")) |>
  mutate(i = as.integer(Row), j = as.integer(Col),
         in_lower = i > j, in_upper = i < j,
         FST_plot    = if_else(in_lower, FST, NA_real_),
         Shared_plot = if_else(in_upper, Shared, NA_integer_),
         sig_plot    = in_lower & !is.na(ll) & ll > 0)

p_C <- ggplot(grid_C, aes(Col, Row)) +
  # lower triangle: FST (tile + value)
  geom_tile(aes(fill = FST_plot), colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(!is.na(FST_plot), sprintf("%.2f", FST_plot), "")),
            size = 3, colour = "black") +
  # mark significant pairs (95% CI excludes 0)
  geom_point(data = filter(grid_C, sig_plot),
             shape = 21, size = 2.2, stroke = 0.9,
             fill = NA, colour = pal$orange) +
  scale_fill_gradient(limits = c(0, 0.6), na.value = "white",
                      low = pal$grey_lo, high = pal$grey_hi, name = "FST (WC84)") +
  # upper triangle: # shared alleles on white
  geom_tile(data = filter(grid_C, in_upper),
            fill = "white", colour = "white", linewidth = 0.5) +
  geom_text(data = filter(grid_C, in_upper),
            aes(label = Shared_plot), size = 3.2, fontface = "bold") +
  coord_equal() +
  labs(x = NULL, y = NULL,
       title = "C. Pairwise: lower = FST (ring if 95% CI > 0), upper = # shared alleles") +
  theme(panel.grid = element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))



# ---- Compose & (optional) export ----
final_plot <- p_map | (p_B / p_C) + plot_layout(widths = c(1.25, 1))
final_plot