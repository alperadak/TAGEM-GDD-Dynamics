#======================================================================
# LONGITUDINAL FUNCTIONAL PCA OF DAILY GDD
#======================================================================

#======================================================================
# 1. SETUP
#======================================================================
rm(list = ls()); gc()
library(dplyr); library(tidyr); library(lubridate); library(refund)
library(ggplot2); library(ggrepel); library(purrr); library(factoextra)
library(cluster); library(viridis)

#======================================================================
# 2. LOAD AND PREPARE DAILY CLIMATE DATA
#======================================================================
env.data <- readRDS("E:\\ALPER_BURAK_ENV\\env_data_1981_2025.rds")

env.data2 <- env.data %>%
  mutate(YYYYMMDD = as.Date(YYYYMMDD), Year = year(YYYYMMDD), Month = month(YYYYMMDD), Day = day(YYYYMMDD),
         env = as.character(env), GDD_10 = ((T2M_MAX + T2M_MIN) / 2) - 10, GDD_10 = pmax(GDD_10, 0), GDD = GDD_10) %>%
  filter(between(Year, 1981, 2025), !(Month == 2 & Day == 29)) %>%
  group_by(env, Year) %>%
  arrange(YYYYMMDD, .by_group = TRUE) %>%
  mutate(DOY = row_number(), CumGDD_10 = cumsum(GDD_10)) %>%
  ungroup() %>%
  mutate(ENVwithYear = paste(env, Year, sep = "_")) %>%
  arrange(env, Year, DOY)

# Analysis dataset
gdd_dat <- env.data2 %>%
  transmute(env = as.character(env), Year = as.integer(Year), DOY = as.integer(DOY),
            GDD = as.numeric(GDD), CumGDD_10 = as.numeric(CumGDD_10)) %>%
  filter(!is.na(env), !is.na(Year), !is.na(DOY), !is.na(GDD), !is.na(CumGDD_10))

#======================================================================
# 3. CREATE ENVIRONMENT × YEAR FUNCTIONAL DATA STRUCTURE
#======================================================================
# Rows = environment-year curves; columns = DOY-specific GDD

gdd_wide <- gdd_dat %>%
  select(env, Year, DOY, GDD) %>%
  pivot_wider(names_from = DOY, values_from = GDD, names_prefix = "DOY_") %>%
  arrange(env, Year)

# Repeated visit = year within environment
gdd_wide <- gdd_wide %>%
  group_by(env) %>%
  arrange(Year, .by_group = TRUE) %>%
  mutate(visit_order = row_number()) %>%
  ungroup()

#======================================================================
# 4. VISUALIZE LONG-TERM CUMULATIVE GDD TRAJECTORIES
#======================================================================
month_df <- data.frame(
  Month = month.abb,
  start = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
  end   = c(31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)
) %>%
  mutate(mid = (start + end) / 2, shade = rep(c(TRUE, FALSE), 6))

p_gdd <- ggplot() +
  geom_rect(data = filter(month_df, shade), aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey90", inherit.aes = FALSE) +
  geom_line(data = gdd_dat, aes(x = DOY, y = CumGDD_10, group = Year, color = Year),
            linewidth = 0.5, alpha = 0.75) +
  facet_wrap(~env, ncol = 5) +
  scale_color_viridis_c(option = "plasma", name = "Year", breaks = c(seq(1981, 2021, by = 5), 2025)) +
  scale_x_continuous(breaks = month_df$mid, labels = month_df$Month, limits = c(1, 365), expand = c(0, 0)) +
  labs(x = "Month", y = "Cumulative GDD") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 9),
        axis.text.x = element_text(angle = 90, vjust = 0.5))

p_gdd

ggsave("E:/ALPER_BURAK_ENV/Long_term_cumulative_GDD_trajectories.png",
       plot = p_gdd, width = 9.5, height = 8, units = "in", dpi = 600)

#======================================================================
# 5. CREATE FUNCTIONAL RESPONSE MATRIX
#======================================================================
# Y: rows = environment-year curves; columns = DOY-specific daily GDD

Y <- gdd_wide %>% select(starts_with("DOY_")) %>% as.matrix()
storage.mode(Y) <- "numeric"

cat("Dimension of Y:", dim(Y), "\n")
cat("Missing values:", sum(is.na(Y)), "\n")
cat("Infinite values:", sum(!is.finite(Y)), "\n")

#======================================================================
# 6. HANDLE MISSING AND INFINITE VALUES
#======================================================================
Y[!is.finite(Y)] <- NA
Y_imp <- Y

for(j in seq_len(ncol(Y_imp))) {
  if(anyNA(Y_imp[, j])) {
    col_mean <- mean(Y_imp[, j], na.rm = TRUE)
    if(!is.finite(col_mean)) stop(paste("DOY column", j, "contains no usable GDD observations."))
    Y_imp[is.na(Y_imp[, j]), j] <- col_mean
  }
}

cat("Missing values after processing:", sum(is.na(Y_imp)), "\n")

#======================================================================
# 7. DEFINE LFPCA DATA LEVELS AND FUNCTIONAL DOMAINS
#======================================================================
# Subject = research location
subject_id <- gdd_wide$env

# Visit = repeated year within research location
visit_id <- gdd_wide$visit_order

# Longitudinal time = calendar year scaled to 0-1
year_min <- min(gdd_wide$Year)
year_max <- max(gdd_wide$Year)
obsT <- (gdd_wide$Year - year_min) / (year_max - year_min)

# Functional domain = DOY scaled to 0-1
funcArg_raw <- as.numeric(sub("DOY_", "", colnames(Y_imp)))
funcArg <- (funcArg_raw - min(funcArg_raw)) / (max(funcArg_raw) - min(funcArg_raw))

#======================================================================
# 8. BASIC DATA SANITY CHECK
#======================================================================
cat("Number of environments:", length(unique(subject_id)), "\n")
cat("Number of observed env-year curves:", nrow(Y_imp), "\n")
cat("Number of DOY points:", ncol(Y_imp), "\n")
cat("Years:", min(gdd_wide$Year), "-", max(gdd_wide$Year), "\n")

#======================================================================
# 9. FIT LONGITUDINAL FPCA
#======================================================================
set.seed(123)

lfpca_gdd <- fpca.lfda(
  Y = Y_imp, subject.index = subject_id, visit.index = visit_id,
  obsT = obsT, funcArg = funcArg,
  numTEvalPoints = length(sort(unique(gdd_wide$Year))),
  LongiModel.method = "fpca.sc",
  mFPCA.pve = 0.95, sFPCA.pve = 0.95,
  mFPCA.knots = 35, sFPCA.nbasis = 10
)

#======================================================================
# 10. INSPECT LFPCA MODEL
#======================================================================
print(lfpca_gdd$mFPCA.npc)
print(lfpca_gdd$sFPCA.npc)
summary(lfpca_gdd)

cat("\nNumber of marginal FPCA dimensions retained:", lfpca_gdd$mFPCA.npc, "\n")
cat("\nNumber of longitudinal PCs retained within each marginal component:\n")
print(lfpca_gdd$sFPCA.npc)

longitudinal_pc_df <- data.frame(
  mFPCA = seq_along(lfpca_gdd$sFPCA.npc),
  n_longitudinal_PC = lfpca_gdd$sFPCA.npc
)

longitudinal_pc_df

#======================================================================
# 11. MARGINAL FPCA VARIANCE EXPLAINED
#======================================================================
marginal_eigenvalues <- lfpca_gdd$mFPCA.evalues
marginal_pve <- (marginal_eigenvalues / sum(marginal_eigenvalues)) * 100
marginal_cum_pve <- cumsum(marginal_pve)

marginal_variance <- data.frame(
  mFPCA = seq_along(marginal_eigenvalues),
  Eigenvalue = marginal_eigenvalues,
  PVE = marginal_pve,
  Cumulative_PVE = marginal_cum_pve
)

print(marginal_variance)

#======================================================================
# 12. FULL SCREE / CUMULATIVE PVE
#======================================================================
scree_eval <- lfpca_gdd$mFPCA.scree.eval

scree_df <- data.frame(
  mFPCA = seq_along(scree_eval),
  Eigenvalue = scree_eval,
  PVE = scree_eval / sum(scree_eval) * 100
) %>%
  mutate(Cumulative_PVE = cumsum(PVE))

p_scree <- ggplot(scree_df, aes(x = mFPCA, y = Cumulative_PVE)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 95, linetype = "dashed") +
  geom_vline(xintercept = lfpca_gdd$mFPCA.npc, linetype = "dashed") +
  scale_x_continuous(breaks = seq(1, length(scree_eval), by = 2)) +
  theme_bw(base_size = 14) +
  labs(x = "Marginal functional principal component",
       y = "Cumulative variance explained (%)") +
  theme(panel.grid = element_blank())

p_scree

ggsave(
  "E:/ALPER_BURAK_ENV/LFPCA_Marginal_FPCA_Cumulative_Variance_Explained.png",
  plot = p_scree,
  width = 6,
  height = 4,
  units = "in",
  dpi = 600
)

length(lfpca_gdd$mFPCA.evalues)
length(lfpca_gdd$mFPCA.scree.eval)
sum(lfpca_gdd$mFPCA.evalues)
sum(lfpca_gdd$mFPCA.scree.eval)

#======================================================================
# 13. BIVARIATE SMOOTH MEAN GDD SURFACE
#======================================================================
# Rows = longitudinal time / year (45); columns = functional time / DOY (365)

mean_surface <- lfpca_gdd$bivariateSmoothMeanFunc
dim(mean_surface)

mean_surface_df <- as.data.frame(mean_surface)
colnames(mean_surface_df) <- 1:365
mean_surface_df$Year <- 1981:2025

mean_surface_long <- mean_surface_df %>%
  pivot_longer(cols = -Year, names_to = "DOY", values_to = "Mean_GDD") %>%
  mutate(DOY = as.integer(DOY))

#======================================================================
# 14. DAILY BIVARIATE MEAN GDD SURFACE
#======================================================================
p_mean_surface <- ggplot(mean_surface_long, aes(x = DOY, y = Year, fill = Mean_GDD)) +
  geom_raster() +
  scale_fill_viridis_c(option = "plasma", name = "Mean GDD") +
  scale_x_continuous(breaks = month_df$mid, labels = month_df$Month,
                     limits = c(1, 365), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(1981, 1985, 1990, 1995, 2000, 2005, 2010, 2015, 2020, 2025),
                     expand = c(0, 0)) +
  labs(x = "Month", y = "Year", title = "Long-term mean seasonal GDD surface") +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

p_mean_surface

#======================================================================
# 15. MONTHLY MEAN GDD FROM BIVARIATE MEAN SURFACE
#======================================================================
month_limits <- c(0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)

monthly_mean_surface <- mean_surface_long %>%
  mutate(Month = cut(DOY, breaks = month_limits, labels = month.abb, include.lowest = TRUE)) %>%
  group_by(Year, Month) %>%
  summarise(Mean_GDD = mean(Mean_GDD, na.rm = TRUE), .groups = "drop") %>%
  mutate(Month = factor(Month, levels = month.abb),
         Year = factor(Year, levels = 1981:2025))

#======================================================================
# 16. MONTHLY MEAN GDD HEATMAP
#======================================================================
p_monthly_surface <- ggplot(monthly_mean_surface, aes(x = Month, y = Year, fill = Mean_GDD)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.1f", Mean_GDD)), size = 1.5, color = "white") +
  scale_fill_viridis_c( name = "Mean GDD") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Month", y = "Year") +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        axis.text.y = element_text(size = 6),
        strip.text = element_text(face = "bold"))

p_monthly_surface

ggsave("E:/ALPER_BURAK_ENV/Monthly_mean_GDD_surface.png",
       plot = p_monthly_surface, width = 6, height = 4, units = "in", dpi = 600)

#======================================================================
# 17. EXTRACT AND VISUALIZE MARGINAL EIGENFUNCTIONS
#======================================================================
n_mfpca <- lfpca_gdd$mFPCA.npc

mfpca_functions <- as.data.frame(
  lfpca_gdd$mFPCA.efunctions[, seq_len(n_mfpca), drop = FALSE]
)

colnames(mfpca_functions) <- paste0("mFPCA", seq_len(n_mfpca))
mfpca_functions$DOY <- funcArg_raw

mfpca_functions_long <- mfpca_functions %>%
  pivot_longer(cols = starts_with("mFPCA"), names_to = "Component", values_to = "Eigenfunction")

p_mfpca_functions <- ggplot(mfpca_functions_long, aes(x = DOY, y = Eigenfunction)) +
  geom_line(linewidth = 1) +
  facet_wrap(~Component) +
  theme_bw(base_size = 14) +
  labs(x = "Day of year", y = "Eigenfunction",
       title = "Marginal FPCA eigenfunctions for daily GDD") +
  theme(panel.grid = element_blank())

p_mfpca_functions

#======================================================================
# 18. EXTRACT LOCATION-SPECIFIC LONGITUDINAL SCORE FUNCTIONS
#======================================================================
env_names <- unique(gdd_wide$env)

if(length(lfpca_gdd$sFPCA.xiHat.bySubj) != length(env_names)) {
  stop("Number of subject-level score matrices does not match number of environments.")
}

names(lfpca_gdd$sFPCA.xiHat.bySubj) <- env_names

score_dimensions <- lapply(lfpca_gdd$sFPCA.xiHat.bySubj, dim)
print(score_dimensions[1:3])

#======================================================================
# 19. CONVERT LONGITUDINAL SCORE FUNCTIONS TO LONG FORMAT
#======================================================================
n_time_eval <- nrow(lfpca_gdd$sFPCA.xiHat.bySubj[[1]])
year_grid <- seq(from = year_min, to = year_max, length.out = n_time_eval)

lfpca_scores_long <- imap_dfr(
  lfpca_gdd$sFPCA.xiHat.bySubj,
  function(mat, env_name) {
    mat <- mat[, seq_len(n_mfpca), drop = FALSE]
    temp <- as.data.frame(mat)
    colnames(temp) <- paste0("mFPCA", seq_len(ncol(temp)))
    temp$env <- env_name
    temp$Year <- year_grid
    temp %>% pivot_longer(cols = starts_with("mFPCA"),
                          names_to = "Component",
                          values_to = "Score")
  }
)

head(lfpca_scores_long)

#======================================================================
# 20. VISUALIZE LONG-TERM LFPCA SCORE TRAJECTORIES
#======================================================================
p_lfpca_scores <- ggplot(lfpca_scores_long, aes(x = Year, y = Score, group = env)) +
  geom_line(alpha = 0.65, linewidth = 0.7) +
  facet_wrap(~Component) +
  theme_bw(base_size = 14) +
  labs(x = "Year", y = "Longitudinal FPCA score",
       title = "Long-term thermal trajectories of environments") +
  theme(panel.grid = element_blank())

p_lfpca_scores

#======================================================================
# 21. PREPARE LOCATION-LEVEL LONGITUDINAL FEATURES FOR CLUSTERING
#======================================================================
# Rows = research locations
# Columns = complete longitudinal score trajectories across retained mFPCs

X_longitudinal <- do.call(
  rbind,
  lapply(
    lfpca_gdd$sFPCA.xiHat.bySubj,
    function(mat) {
      mat_use <- mat[, seq_len(n_mfpca), drop = FALSE]
      as.vector(mat_use)
    }
  )
)

#======================================================================
# 22. STANDARDIZE LOCATION-LEVEL LONGITUDINAL FEATURES
#======================================================================
# Rows = research locations
# Columns = longitudinal mFPCA score features

rownames(X_longitudinal) <- env_names
X_scaled <- scale(X_longitudinal)

cat("Dimension of X_longitudinal:", dim(X_longitudinal), "\n")
cat("Dimension of X_scaled:", dim(X_scaled), "\n")
cat("Missing values:", sum(is.na(X_scaled)), "\n")
cat("Infinite values:", sum(!is.finite(X_scaled)), "\n")

#======================================================================
# 23. SECONDARY PCA OF LONGITUDINAL LFPCA FEATURES
#======================================================================
pca_cluster <- prcomp(X_scaled, center = FALSE, scale. = FALSE)

pca_eigenvalues <- pca_cluster$sdev^2
pca_pve <- pca_eigenvalues / sum(pca_eigenvalues)
pca_cum_pve <- cumsum(pca_pve)

pca_variance <- data.frame(
  PC = seq_along(pca_pve),
  PVE = pca_pve * 100,
  Cumulative_PVE = pca_cum_pve * 100
)

print(pca_variance)

#======================================================================
# 24. SELECT PCs EXPLAINING AT LEAST 95% OF VARIANCE
#======================================================================
n_pc_keep <- which(pca_cum_pve >= 0.95)[1]

cat("Number of secondary PCs retained:", n_pc_keep, "\n")
cat("Cumulative variance explained:", round(pca_cum_pve[n_pc_keep] * 100, 2), "%\n")

X_C <- pca_cluster$x[, seq_len(n_pc_keep), drop = FALSE]

cat("Dimension of clustering matrix X_C:", dim(X_C), "\n")

#======================================================================
# 25. VISUALIZE SECONDARY PCA VARIANCE EXPLAINED
#======================================================================
p_pca_variance <- ggplot(
  pca_variance,
  aes(x = PC, y = Cumulative_PVE)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 95, linetype = "dashed") +
  geom_vline(xintercept = n_pc_keep, linetype = "dashed") +
  
  scale_x_continuous(
    breaks = sort(unique(c(
      seq(1, nrow(pca_variance), by = 2),
      6
    )))
  ) +
  
  scale_y_continuous(
    breaks = seq(30, 100, by = 5),
    limits = c(30, 100)
  ) +
  
  theme_bw(base_size = 14) +
  labs(
    x = "Principal component",
    y = "Cumulative variance explained (%)"
  ) +
  theme(panel.grid = element_blank())

p_pca_variance


ggsave(
  filename = "E:/ALPER_BURAK_ENV/Secondary_PCA_Cumulative_Variance_Explained.png",
  plot = p_pca_variance,
  width = 6,
  height = 4,
  units = "in",
  dpi = 600
)

#======================================================================
# 26. DETERMINE OPTIMAL NUMBER OF K-MEANS CLUSTERS
#======================================================================
set.seed(123)

k_range <- 2:20

sil_kmeans <- data.frame(
  K = k_range,
  Mean_Silhouette = sapply(k_range, function(k) {
    km_temp <- kmeans(X_C, centers = k, nstart = 100, iter.max = 1000)
    mean(cluster::silhouette(km_temp$cluster, dist(X_C))[, "sil_width"])
  })
)

print(sil_kmeans)

best_k <- sil_kmeans$K[which.max(sil_kmeans$Mean_Silhouette)]
best_silhouette <- max(sil_kmeans$Mean_Silhouette)

cat("Optimal number of k-means clusters:", best_k, "\n")
cat("Maximum mean silhouette width:", round(best_silhouette, 3), "\n")

#======================================================================
# 27. VISUALIZE K-MEANS CLUSTER NUMBER SELECTION
#======================================================================
p_k_selection <- ggplot(sil_kmeans, aes(x = K, y = Mean_Silhouette)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_point(data = dplyr::filter(sil_kmeans, K == best_k), size = 5, shape = 21, stroke = 1.2) +
  geom_vline(xintercept = best_k, linetype = "dashed") +
  scale_x_continuous(breaks = k_range) +
  theme_bw(base_size = 14) +
  labs(x = "Number of clusters (K)",
       y = "Mean silhouette width") +
  theme(panel.grid = element_blank())

p_k_selection

ggsave(
  filename = "E:/ALPER_BURAK_ENV/LFPCA_PCA_kmeans_cluster_selection.png",
  plot = p_k_selection,
  width = 6,
  height = 5,
  units = "in",
  dpi = 600
)
#======================================================================
# 28. FIT FINAL K-MEANS CLUSTERING
#======================================================================
set.seed(123)

km_lfpca <- kmeans(
  X_C,
  centers = best_k,
  nstart = 500,
  iter.max = 1000
)

table(km_lfpca$cluster)

#======================================================================
# 29. CREATE FINAL LOCATION CLUSTER TABLE
#======================================================================
env_cluster_lfpca <- data.frame(
  env = rownames(X_C),
  Cluster = factor(km_lfpca$cluster)
) %>%
  arrange(Cluster, env)

env_cluster_lfpca

cluster_size_lfpca <- env_cluster_lfpca %>%
  dplyr::count(Cluster, name = "N_locations")

cluster_size_lfpca

#======================================================================
# 30. DEFINE CLUSTER COLOR PALETTE
#======================================================================
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "gold", "#0072B2",
  "#D55E00", "#CC79A7", "#000000", "#7F7F7F", "#8C510A"
)

cluster_colors <- okabe_ito[seq_len(best_k)]

#======================================================================
# 31. PCA SCATTER PLOT OF K-MEANS CLASSIFICATION
#======================================================================
cluster_plot_df <- data.frame(
  env = rownames(pca_cluster$x),
  PC1 = pca_cluster$x[, 1],
  PC2 = pca_cluster$x[, 2]
) %>%
  left_join(env_cluster_lfpca, by = "env") %>%
  mutate(Label = paste0(env, " (", Cluster, ")"))

p_cluster_pca <- ggplot(cluster_plot_df,
                        aes(x = PC1, y = PC2, color = Cluster, label = Label)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5, max.overlaps = Inf) +
  scale_color_manual(values = cluster_colors, name = "LFPCA-informed\nk-means cluster") +
  theme_bw(base_size = 14) +
  labs(x = paste0("PC1 (", round(100 * pca_pve[1], 1), "%)"),
       y = paste0("PC2 (", round(100 * pca_pve[2], 1), "%)"),
       #title = "Longitudinal GDD-based classification of research locations",
       #subtitle = paste0("LFPCA + secondary PCA + k-means; K = ", best_k)
       ) +
  theme(panel.grid = element_blank())

p_cluster_pca

ggsave(
  filename = "E:/ALPER_BURAK_ENV/LFPCA_PCA_kmeans_classification.png",
  plot = p_cluster_pca,
  width = 9,
  height = 6,
  units = "in",
  dpi = 600
)

#======================================================================
# 32. FINAL CLUSTER SILHOUETTE PLOT
#======================================================================
sil_final <- cluster::silhouette(km_lfpca$cluster, dist(X_C))

p_final_silhouette <- factoextra::fviz_silhouette(sil_final) +
  theme_bw(base_size = 14) +
  labs(title = "Silhouette structure of final k-means clusters",
       subtitle = paste0("K = ", best_k,
                         "; mean silhouette = ", round(mean(sil_final[, "sil_width"]), 3))) +
  theme(panel.grid = element_blank())

p_final_silhouette

#======================================================================
# 33. JOIN K-MEANS CLUSTERS BACK TO ORIGINAL GDD DATA
#======================================================================
gdd_clustered <- gdd_dat %>%
  left_join(env_cluster_lfpca, by = "env")

#======================================================================
# 34. LONG-TERM DAILY GDD TRAJECTORIES OF K-MEANS CLUSTERS
#======================================================================
cluster_daily_gdd <- gdd_clustered %>%
  group_by(Cluster, DOY) %>%
  summarise(
    Mean_GDD = mean(GDD, na.rm = TRUE),
    SD_GDD = sd(GDD, na.rm = TRUE),
    .groups = "drop"
  )

p_cluster_gdd <- ggplot() +
  geom_rect(data = dplyr::filter(month_df, shade),
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "grey90", alpha = 0.5, inherit.aes = FALSE) +
  geom_line(data = cluster_daily_gdd,
            aes(x = DOY, y = Mean_GDD, color = Cluster),
            linewidth = 1.2) +
  scale_color_manual(values = cluster_colors, name = "k-means cluster") +
  scale_x_continuous(breaks = month_df$mid, labels = month_df$Month,
                     limits = c(1, 365), expand = c(0, 0)) +
  theme_bw(base_size = 14) +
  labs(x = "Month",
       y = "Mean daily GDD",
       title = "Long-term daily GDD trajectories of research location clusters") +
  theme(panel.grid = element_blank())

p_cluster_gdd

#======================================================================
# 35. LOAD TURKEY MAP AND RESEARCH LOCATION COORDINATES
#======================================================================
library(sf); library(readxl)

tur_prov <- st_read("C:/Users/a-ada/Downloads/gadm41_TUR_shp/gadm41_TUR_1.shp", quiet = TRUE)
env_locations <- read_excel("E:/ALPER_BURAK_ENV/Institues_with_coords.xlsx")

head(env_locations)
names(env_locations)

#======================================================================
# 36. JOIN K-MEANS CLUSTERS WITH LOCATION COORDINATES
#======================================================================
# env_cluster_lfpca = env, Cluster
# env_locations = Env, Latitude, Longitude

env_locations_cluster <- env_locations %>%
  left_join(env_cluster_lfpca, by = c("Env" = "env"))

print(env_locations_cluster)

#======================================================================
# 37. CHECK UNMATCHED RESEARCH LOCATIONS
#======================================================================
env_locations_cluster %>% filter(is.na(Cluster))

#======================================================================
# 38. CREATE MAP LABELS
#======================================================================
# Example: BATEM (2)

env_locations_cluster <- env_locations_cluster %>%
  mutate(Label = paste0(Env, "\n(", Cluster, ")"))

#======================================================================
# 39. DEFINE MAP CLUSTER COLORS
#======================================================================
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "gold", "#0072B2",
               "#D55E00", "#CC79A7", "#000000", "#7F7F7F", "#8C510A")

n_clusters <- nlevels(factor(env_locations_cluster$Cluster))
cluster_colors <- okabe_ito[seq_len(n_clusters)]

#======================================================================
# 40. MAP LFPCA-INFORMED K-MEANS CLUSTERS ACROSS TÜRKİYE
#======================================================================
p_map <- ggplot() +
  geom_sf(data = tur_prov, fill = "gray95", color = "gray80", linewidth = 0.25) +
  geom_point(data = env_locations_cluster,
             aes(x = Longitude, y = Latitude, color = Cluster), size = 3) +
  geom_text_repel(data = env_locations_cluster,
                  aes(x = Longitude, y = Latitude, label = Label, color = Cluster),
                  size = 2.6, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = cluster_colors, name = "LFPCA-informed\nk-means cluster") +
  coord_sf(xlim = c(25, 45), ylim = c(35, 43), expand = FALSE) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = "right", panel.grid = element_blank())

p_map

library(sf); library(ggplot2); library(ggrepel); library(rnaturalearth); library(rnaturalearthdata)

world <- ne_countries(scale = "medium", returnclass = "sf")

p_map <- ggplot() +
  geom_sf(data = world, fill = "gray97", color = "gray75", linewidth = 0.20) +
  geom_sf(data = tur_prov, fill = "gray88", color = "gray40", linewidth = 0.35) +
  geom_point(data = env_locations_cluster, aes(x = Longitude, y = Latitude, color = Cluster), size = 3) +
  geom_text_repel(data = env_locations_cluster, aes(x = Longitude, y = Latitude, label = Label, color = Cluster), size = 2.5, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.3, min.segment.length = 0, show.legend = FALSE) +
  scale_color_manual(values = cluster_colors, name = "LFPCA-informed\nk-means cluster") +
  coord_sf(xlim = c(20, 47), ylim = c(34, 45), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right", panel.grid = element_blank(), axis.title = element_text(size = 11), axis.text = element_text(size = 9))

p_map


install.packages("ggspatial")

library(sf); library(readxl); library(dplyr); library(ggplot2); library(ggrepel); library(ggspatial)

#======================================================================
# LOAD TÜRKİYE AND RESEARCH LOCATIONS
#======================================================================

tur_prov <- st_read("C:/Users/a-ada/Downloads/gadm41_TUR_shp/gadm41_TUR_1.shp", quiet = TRUE)
env_locations <- read_excel("E:/ALPER_BURAK_ENV/Institues_with_coords.xlsx")

env_locations_cluster <- env_locations %>% left_join(env_cluster_lfpca, by = c("Env" = "env"))
env_locations_cluster <- env_locations_cluster %>% mutate(Label = paste0(Env, "\n(", Cluster, ")"))

okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "gold", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#7F7F7F", "#8C510A")
n_clusters <- nlevels(factor(env_locations_cluster$Cluster))
cluster_colors <- okabe_ito[seq_len(n_clusters)]

#======================================================================
# MAP WITH AUTOMATIC GEOGRAPHIC LABELS
#======================================================================

p_map <- ggplot() +
  annotation_map_tile(type = "osm", zoomin = 0) +
  geom_sf(data = tur_prov, fill = NA, color = "gray40", linewidth = 0.35) +
  geom_point(data = env_locations_cluster, aes(x = Longitude, y = Latitude, color = Cluster), size = 3) +
  geom_text_repel(data = env_locations_cluster, aes(x = Longitude, y = Latitude, label = Label, color = Cluster),
                  size = 2.6, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = cluster_colors, name = "LFPCA-informed\nk-means cluster") +
  coord_sf(xlim = c(20, 47), ylim = c(34, 45), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(legend.position = "right", panel.grid = element_blank())

p_map

ggsave("E:/ALPER_BURAK_ENV/Turkey_Institute_Map_LFPCA_Cluster_Basemap.png",
       plot = p_map, width = 11, height = 6, units = "in", dpi = 600)
#======================================================================
# 41. SAVE TURKEY CLUSTER MAP
#======================================================================
ggsave("E:/ALPER_BURAK_ENV/Turkey_Institute_Map_LFPCA_Cluster2.png",
       plot = p_map, width = 9, height = 5, units = "in", dpi = 600)

#======================================================================
# 42. PREPARE DAILY GDD DATA FOR TEMPORAL TREND ANALYSIS
#======================================================================
gdd_daily_trend <- gdd_dat %>%
  mutate(Year_centered = Year - 1981) %>%
  arrange(env, DOY, Year)

#======================================================================
# 43. FIT DAILY LINEAR GDD TREND FOR EACH LOCATION × DOY
#======================================================================
# Model: GDD = Intercept + Slope × (Year - 1981)
# Intercept = estimated GDD in 1981
# Slope = annual change in daily GDD

daily_models <- gdd_daily_trend %>%
  group_by(env, DOY) %>%
  nest() %>%
  mutate(model = purrr::map(data, ~lm(GDD ~ Year_centered, data = .x)))

daily_models

#======================================================================
# 44. EXTRACT DAILY LINEAR TREND STATISTICS
#======================================================================
daily_trend_results <- daily_models %>%
  mutate(
    coef_info = purrr::map(model, broom::tidy),
    model_info = purrr::map(model, broom::glance),
    Intercept = purrr::map_dbl(model, ~coef(.x)[1]),
    Slope = purrr::map_dbl(model, ~coef(.x)[2]),
    R2 = purrr::map_dbl(model, ~summary(.x)$r.squared),
    Adjusted_R2 = purrr::map_dbl(model, ~summary(.x)$adj.r.squared),
    RMSE = purrr::map_dbl(model, ~sqrt(mean(residuals(.x)^2))),
    P_value = purrr::map_dbl(model, ~coef(summary(.x))[2, 4]),
    N_years = purrr::map_int(data, nrow)
  ) %>%
  select(env, DOY, N_years, Intercept, Slope, R2, Adjusted_R2, RMSE, P_value) %>%
  mutate(Significance = case_when(P_value < 0.001 ~ "***", P_value < 0.01 ~ "**",
                                  P_value < 0.05 ~ "*", TRUE ~ "ns"))

print(daily_trend_results, n = 50)

#======================================================================
# 45. SAVE DAILY GDD LINEAR TREND STATISTICS
#======================================================================
write.csv(daily_trend_results,
          "E:/ALPER_BURAK_ENV/Daily_GDD_Linear_Trend_Results.csv",
          row.names = FALSE)

#======================================================================
# 46. EXAMPLE LOCATION × DOY LINEAR TREND
#======================================================================
selected_env <- "FAEM"
selected_DOY <- 216

plot_dat <- gdd_daily_trend %>%
  filter(env == selected_env, DOY == selected_DOY)

fit_lm <- lm(GDD ~ Year, data = plot_dat)
lm_sum <- summary(fit_lm)

intercept <- coef(fit_lm)[1]
slope <- coef(fit_lm)[2]
r2 <- lm_sum$r.squared
pval <- coef(lm_sum)[2, 4]

subtitle_text <- paste0(
  "Linear model: GDD = ", round(intercept, 2),
  ifelse(slope >= 0, " + ", " - "), round(abs(slope), 2), " × Year",
  "   |   Slope = ", round(slope, 2),
  "   |   R² = ", round(r2, 2),
  "   |   P = ", ifelse(pval < 0.001, "< 0.001", format.pval(pval, digits = 3))
)

p_example_trend <- ggplot(plot_dat, aes(x = Year, y = GDD)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 1) +
  labs(x = "Year", y = "Daily GDD",
       title = paste0(selected_env, " — DOY ", selected_DOY),
       subtitle = subtitle_text) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 11))

p_example_trend

#======================================================================
# 47. DEFINE MONTH INFORMATION FOR DAILY TREND VISUALIZATION
#======================================================================
month_trend_df <- data.frame(
  Month = month.abb,
  xmin = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
  xmax = c(31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365)
) %>%
  mutate(xmid = (xmin + xmax) / 2, Shade = rep(c("A", "B"), 6))

#======================================================================
# 48. DAILY GDD TEMPORAL SLOPE PROFILES WITH MONTHLY SHADING
#======================================================================
p_daily_slope <- ggplot() +
  geom_rect(data = month_trend_df,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Shade),
            inherit.aes = FALSE, alpha = 0.40) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_line(data = daily_trend_results, aes(x = DOY, y = Slope), linewidth = 0.7) +
  facet_wrap(~env, ncol = 5) +
  scale_fill_manual(values = c("A" = "grey70", "B" = "grey90"), guide = "none") +
  scale_x_continuous(breaks = month_trend_df$xmid, labels = month_trend_df$Month,
                     limits = c(1, 365), expand = c(0, 0)) +
  theme_bw(base_size = 11) +
  labs(x = "Month",
       y = expression("Linear trend in daily GDD (GDD "*year^{-1}*")")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))

p_daily_slope

ggsave("E:/ALPER_BURAK_ENV/Daily_GDD_Linear_Slope_All_Locations_month.png",
       plot = p_daily_slope, width = 10, height = 8, units = "in", dpi = 600)

p_daily_slope <- ggplot() +
  geom_rect(data = month_trend_df,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Shade),
            inherit.aes = FALSE, alpha = 0.40) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_point(data = daily_trend_results,
             aes(x = DOY, y = Slope, color = R2),
             size = 0.8, alpha = 0.85) +
  facet_wrap(~env, ncol = 5) +
  scale_fill_manual(values = c("A" = "grey70", "B" = "grey90"), guide = "none") +
  scale_color_viridis_c(option = "plasma", name = expression(R^2)) +
  scale_x_continuous(breaks = month_trend_df$xmid,
                     labels = month_trend_df$Month,
                     limits = c(1, 365), expand = c(0, 0)) +
  theme_bw(base_size = 11) +
  labs(x = "Month",
       y = expression("Daily GDD trend slope (GDD "*year^{-1}*")")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))

p_daily_slope

ggsave("E:/ALPER_BURAK_ENV/Daily_GDD_Linear_Slope_R2_All_Locations.png",
       plot = p_daily_slope, width = 10, height = 8, units = "in", dpi = 600)

#======================================================================
# 49. ASSIGN DAILY GDD TREND SLOPES TO MONTHS
#======================================================================
# 2001 is used as a non-leap reference year

daily_trend_auc <- daily_trend_results %>%
  mutate(
    Date_ref = as.Date(DOY - 1, origin = "2001-01-01"),
    Month_num = as.integer(format(Date_ref, "%m")),
    Month = factor(format(Date_ref, "%b"), levels = month.abb)
  )

#======================================================================
# 50. DEFINE TRAPEZOIDAL INTEGRATION FUNCTION
#======================================================================
trapz_auc <- function(x, y) {
  ord <- order(x); x <- x[ord]; y <- y[ord]
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2, na.rm = TRUE)
}

#======================================================================
# 51. CALCULATE MONTHLY INTEGRATED DAILY GDD TREND AUC
#======================================================================
monthly_slope_auc <- daily_trend_auc %>%
  group_by(env, Month_num, Month) %>%
  arrange(DOY, .by_group = TRUE) %>%
  summarise(
    Net_AUC = trapz_auc(DOY, Slope),
    Positive_AUC = trapz_auc(DOY, pmax(Slope, 0)),
    Negative_AUC = trapz_auc(DOY, pmin(Slope, 0)),
    Absolute_AUC = trapz_auc(DOY, abs(Slope)),
    Mean_Slope = mean(Slope, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(env, Month_num)

monthly_slope_auc

#======================================================================
# 52. RANK RESEARCH LOCATIONS WITHIN EACH MONTH BY NET AUC
#======================================================================
monthly_slope_auc_ranked <- monthly_slope_auc %>%
  mutate(Month = factor(Month, levels = month.abb)) %>%
  group_by(Month) %>%
  arrange(desc(Net_AUC), .by_group = TRUE) %>%
  mutate(Rank = row_number(),
         Rank_label = ifelse(Rank <= 5, as.character(Rank), "")) %>%
  ungroup()

#======================================================================
# 53. MONTHLY NET SLOPE AUC HEATMAP WITH TOP-FIVE LOCATIONS
#======================================================================
p_auc_heatmap <- ggplot(monthly_slope_auc_ranked,
                        aes(x = Month, y = env, fill = Net_AUC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = Rank_label), fontface = "bold", size = 3.5) +
  scale_fill_gradient(low = "cyan", high = "red", name = "Net slope AUC") +
  labs(x = "Month", y = "TAGEM research location"
       #title = "Monthly integrated long-term trends in daily GDD",
       #subtitle = "Net AUC was calculated from daily GDD trend slopes\nNumbers indicate the five highest-ranked locations within each month"
       ) +
  theme_bw(base_size = 13) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(face = "bold"))

p_auc_heatmap

ggsave("E:/ALPER_BURAK_ENV/Monthly_Net_Slope_AUC_Heatmap.png",
       plot = p_auc_heatmap, width = 9, height = 7, units = "in", dpi = 600)
#======================================================================
# 54. ASSIGN DAILY SLOPES TO MONTHS FOR DISTRIBUTION ANALYSIS
#======================================================================
daily_trend_results2 <- daily_trend_results %>%
  mutate(
    Month = format(as.Date(DOY - 1, origin = "2001-01-01"), "%b"),
    Month = factor(Month, levels = month.abb)
  )

#======================================================================
# 55. CALCULATE MONTHLY MEAN SLOPE BY RESEARCH LOCATION
#======================================================================
plot_data <- daily_trend_results2 %>%
  group_by(env, Month) %>%
  mutate(mean_slope = mean(Slope, na.rm = TRUE)) %>%
  ungroup()

rank_data <- daily_trend_results2 %>%
  group_by(Month, env) %>%
  summarise(mean_slope = mean(Slope, na.rm = TRUE),
            ymax = max(Slope, na.rm = TRUE),
            .groups = "drop") %>%
  group_by(Month) %>%
  arrange(desc(mean_slope), .by_group = TRUE) %>%
  mutate(Rank = row_number()) %>%
  filter(Rank <= 5)

#======================================================================
# 56. MONTHLY DISTRIBUTION OF DAILY GDD TREND SLOPES
#======================================================================
p_monthly_slope_box <- ggplot(plot_data, aes(x = env, y = Slope, fill = mean_slope)) +
  geom_boxplot(outlier.size = 0.5, width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_text(data = rank_data,
            aes(x = env, y = ymax + 0.02 * diff(range(plot_data$Slope, na.rm = TRUE)), label = Rank),
            inherit.aes = FALSE, size = 4, fontface = "bold") +
  facet_wrap(~Month, ncol = 3) +
  scale_fill_viridis_c(option = "plasma", name = "Mean\nSlope") +
  theme_bw(base_size = 12) +
  labs(x = "TAGEM research location",
       y = expression("Daily GDD trend slope (GDD "*year^{-1}*")")) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5))

p_monthly_slope_box




#======================================================================
# 57. FPCA OF DAILY GDD TREND SLOPE CURVES
#======================================================================

library(dplyr)
library(tidyr)
library(fdapace)
library(ggplot2)
library(ggrepel)

#----------------------------------------------------------------------
# 57.1 Wide matrix: rows = locations, columns = DOY, values = Slope
#----------------------------------------------------------------------

slope_wide <- daily_trend_auc %>%
  dplyr::select(env, DOY, Slope) %>%
  pivot_wider(names_from = DOY, values_from = Slope, names_prefix = "DOY_") %>%
  arrange(env)

ids <- slope_wide$env

Ymat_slope <- slope_wide %>%
  dplyr::select(-env) %>%
  as.matrix()

storage.mode(Ymat_slope) <- "numeric"

t_points_slope <- as.numeric(sub("DOY_", "", colnames(Ymat_slope)))

dim(Ymat_slope)   # Expected: 30 x 365


#----------------------------------------------------------------------
# 57.2 Convert to FPCA input
#----------------------------------------------------------------------

Ly_slope <- lapply(seq_len(nrow(Ymat_slope)), function(i) Ymat_slope[i, ])
Lt_slope <- replicate(nrow(Ymat_slope), t_points_slope, simplify = FALSE)

names(Ly_slope) <- ids
names(Lt_slope) <- ids


#----------------------------------------------------------------------
# 57.3 Run dense FPCA
#----------------------------------------------------------------------

FPCA_slope <- FPCA(
  Ly_slope,
  Lt_slope,
  list(
    dataType = "Dense",
    error = FALSE,
    methodSelectK = "FVE",
    FVEthreshold = 0.95,
    maxK = 5,
    plot = TRUE
  )
)

plot(FPCA_slope)
print(FPCA_slope)
FPCA_slope$selectK
FPCA_slope$lambda
head(FPCA_slope$xiEst)


#----------------------------------------------------------------------
# 57.4 Variance explained
#----------------------------------------------------------------------

eigenvalues_slope <- FPCA_slope$lambda
variance_explained_slope <- eigenvalues_slope / sum(eigenvalues_slope) * 100
cumulative_variance_slope <- cumsum(variance_explained_slope)

fpca_variance_slope <- data.frame(
  Component = paste0("FPCA", seq_along(eigenvalues_slope)),
  Eigenvalue = eigenvalues_slope,
  Variance = variance_explained_slope,
  Cumulative = cumulative_variance_slope
)

print(fpca_variance_slope)

#----------------------------------------------------------------------
# 57.4B Cumulative variance explained by slope-FPCA
#----------------------------------------------------------------------

fpca_variance_slope <- fpca_variance_slope %>%
  mutate(
    FPCA_num = seq_len(n()),
    Cumulative = cumsum(Variance)
  )

selected_K_slope <- FPCA_slope$selectK

p_fpca_cumulative_slope <- ggplot(
  fpca_variance_slope,
  aes(x = FPCA_num, y = Cumulative)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  
  geom_hline(
    yintercept = 95,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  geom_vline(
    xintercept = selected_K_slope,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  scale_x_continuous(
    breaks = seq_len(nrow(fpca_variance_slope))
  ) +
  
  scale_y_continuous(
    breaks = seq(40, 100, by = 10),
    limits = c(40, 100)
  ) +
  
  theme_bw(base_size = 14) +
  
  labs(
    x = "Functional principal component",
    y = "Cumulative variance explained (%)"
  ) +
  
  theme(
    panel.grid = element_blank()
  )

p_fpca_cumulative_slope


ggsave(
  "E:/ALPER_BURAK_ENV/FPCA_Slope_Cumulative_Variance_Explained.png",
  plot = p_fpca_cumulative_slope,
  width = 5,
  height = 4,
  units = "in",
  dpi = 600
)

#----------------------------------------------------------------------
# 57.5 Extract location FPCA scores
#----------------------------------------------------------------------

n_pc_slope <- min(5, ncol(FPCA_slope$xiEst))

fpca_slope_scores <- as.data.frame(FPCA_slope$xiEst[, seq_len(n_pc_slope), drop = FALSE])
colnames(fpca_slope_scores) <- paste0("FPCA", seq_len(n_pc_slope))
fpca_slope_scores$env <- ids

fpca_slope_scores <- fpca_slope_scores %>%
  dplyr::select(env, everything())

print(fpca_slope_scores)


#----------------------------------------------------------------------
# 57.6 FPCA1 vs FPCA2
#----------------------------------------------------------------------

pc1_var <- round(variance_explained_slope[1], 1)
pc2_var <- round(variance_explained_slope[2], 1)

p_fpca_slope <- ggplot(fpca_slope_scores, aes(x = FPCA1, y = FPCA2, label = env)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3.5, max.overlaps = Inf) +
  theme_bw(base_size = 14) +
  labs(
    x = paste0("FPCA1 (", pc1_var, "%)"),
    y = paste0("FPCA2 (", pc2_var, "%)"),
    title = "FPCA of daily long-term GDD trend slopes",
    subtitle = "Research locations characterized by their 365-day slope trajectories"
  )

p_fpca_slope


#----------------------------------------------------------------------
# 57.7 Scree plot
#----------------------------------------------------------------------

p_scree_slope <- ggplot(fpca_variance_slope, aes(x = Component, y = Variance)) +
  geom_col() +
  geom_text(aes(label = paste0(round(Variance, 1), "%")), vjust = -0.4) +
  theme_bw(base_size = 14) +
  labs(
    x = "Functional principal component",
    y = "Variance explained (%)",
    title = "Variance explained by FPCA components"
  )

p_scree_slope


#----------------------------------------------------------------------
# 57.8 FPCA eigenfunctions across DOY
#----------------------------------------------------------------------

phi_slope <- as.data.frame(FPCA_slope$phi)
phi_slope$DOY <- FPCA_slope$workGrid

phi_long_slope <- phi_slope %>%
  pivot_longer(cols = -DOY, names_to = "Component", values_to = "Eigenfunction")

p_eigen_slope <- ggplot(phi_long_slope, aes(x = DOY, y = Eigenfunction, color = Component)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  theme_bw(base_size = 14) +
  labs(
    x = "Day of year",
    y = "Eigenfunction",
    title = "Functional principal components of daily GDD trend slopes",
    color = "Component"
  )

p_eigen_slope


#======================================================================
# 58. K-MEANS CLUSTERING USING FPCA SCORES
#======================================================================

library(dplyr)
library(ggplot2)
library(ggrepel)
library(cluster)

#----------------------------------------------------------------------
# 58.1 Prepare FPCA score matrix
#----------------------------------------------------------------------

K_fpca <- FPCA_slope$selectK

X_kmeans <- fpca_slope_scores %>%
  dplyr::select(all_of(paste0("FPCA", 1:K_fpca))) %>%
  as.matrix()

rownames(X_kmeans) <- fpca_slope_scores$env

dim(X_kmeans)
head(X_kmeans)


#----------------------------------------------------------------------
# 58.2 Evaluate K = 2:10
#----------------------------------------------------------------------

set.seed(123)

k_values <- 2:20

wss <- sapply(k_values, function(k) kmeans(X_kmeans, centers = k, nstart = 100)$tot.withinss)

silhouette_mean <- sapply(k_values, function(k) {
  km <- kmeans(X_kmeans, centers = k, nstart = 100)
  mean(silhouette(km$cluster, dist(X_kmeans))[, 3])
})

cluster_eval <- data.frame(K = k_values, WSS = wss, Silhouette = silhouette_mean)

print(cluster_eval)


#----------------------------------------------------------------------
# 58.3 Elbow plot
#----------------------------------------------------------------------

p_elbow <- ggplot(cluster_eval, aes(K, WSS)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_x_continuous(breaks = k_values) +
  theme_bw(base_size = 14) +
  labs(x = "Number of clusters (K)", y = "Total within-cluster sum of squares",
       title = "Elbow method for FPCA-based clustering")

p_elbow


#----------------------------------------------------------------------
# 58.4 Silhouette plot
#----------------------------------------------------------------------
# Find optimal K
best_k <- cluster_eval$K[which.max(cluster_eval$Silhouette)]

# Silhouette plot
p_silhouette <- ggplot(cluster_eval, aes(K, Silhouette)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_vline(
    xintercept = best_k,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  scale_x_continuous(breaks = k_values) +
  theme_bw(base_size = 14) +
  labs(
    x = "Number of clusters (K)",
    y = "Mean silhouette width",
    #title = "Silhouette analysis for FPCA-based clustering"
  )

p_silhouette

# Check selected K
best_k

ggsave(
  "E:/ALPER_BURAK_ENV/FPCA_Slope_Silhouette.png",
  plot = p_silhouette,
  width = 6,
  height = 4,
  units = "in",
  dpi = 600
)


#----------------------------------------------------------------------
# 58.5 Best K and final k-means
#----------------------------------------------------------------------

best_k <- cluster_eval$K[which.max(cluster_eval$Silhouette)]
best_k

set.seed(123)
km_fpca_slope <- kmeans(X_kmeans, centers = best_k, nstart = 100)

fpca_slope_scores$Cluster <- factor(km_fpca_slope$cluster)

env_cluster_fpca_slope <- data.frame(
  env = rownames(X_kmeans),
  Cluster = factor(km_fpca_slope$cluster)
) %>% arrange(Cluster, env)

print(env_cluster_fpca_slope)


#----------------------------------------------------------------------
# 58.6 FPCA1 vs FPCA2 with clusters
#----------------------------------------------------------------------

pc1_var <- round(variance_explained_slope[1], 1)
pc2_var <- round(variance_explained_slope[2], 1)

okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "gold", "#0072B2",
  "#D55E00", "#CC79A7", "#000000", "#7F7F7F", "#8C510A"
)

p_fpca_cluster <- ggplot(
  fpca_slope_scores,
  aes(FPCA1, FPCA2, color = Cluster, label = env)
) +
  geom_point(size = 4) +
  geom_text_repel(size = 3.5, max.overlaps = Inf) +
  
  scale_color_manual(values = okabe_ito) +
  
  theme_bw(base_size = 14) +
  labs(
    x = paste0("FPCA1 (", pc1_var, "%)"),
    y = paste0("FPCA2 (", pc2_var, "%)"),
    color = "Cluster"
  )

p_fpca_cluster

ggsave(
     "E:/ALPER_BURAK_ENV/FPCA_Slope_based_cluster.png",
     plot = p_fpca_cluster,
     width = 8.5,
     height = 4.5,
     units = "in",
     dpi = 600
   )

#----------------------------------------------------------------------
# 58.7 Join clusters to original daily slope data
#----------------------------------------------------------------------

daily_trend_auc_clustered <- daily_trend_auc %>%
  left_join(env_cluster_fpca_slope, by = "env")


#----------------------------------------------------------------------
# 58.8 Mean daily slope trajectory by cluster
#----------------------------------------------------------------------

cluster_slope_profile <- daily_trend_auc_clustered %>%
  group_by(Cluster, DOY) %>%
  summarise(Mean_Slope = mean(Slope, na.rm = TRUE),
            SD_Slope = sd(Slope, na.rm = TRUE), .groups = "drop")

p_cluster_slope <- ggplot(cluster_slope_profile, aes(DOY, Mean_Slope, color = Cluster, fill = Cluster)) +
  geom_ribbon(aes(ymin = Mean_Slope - SD_Slope, ymax = Mean_Slope + SD_Slope),
              alpha = 0.15, color = NA) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  theme_bw(base_size = 14) +
  labs(x = "Day of year", y = "Mean GDD trend slope", color = "Cluster", fill = "Cluster",
       title = "Seasonal GDD trend profiles of FPCA-based clusters")

p_cluster_slope


#======================================================================
# 59. MAP FPCA-SLOPE K-MEANS CLUSTERS ACROSS TÜRKİYE
#======================================================================

library(sf); library(readxl); library(dplyr); library(ggplot2); library(ggrepel)

#----------------------------------------------------------------------
# 59.1 Load Türkiye map and research location coordinates
#----------------------------------------------------------------------

tur_prov <- st_read("C:/Users/a-ada/Downloads/gadm41_TUR_shp/gadm41_TUR_1.shp", quiet = TRUE)
env_locations <- read_excel("E:/ALPER_BURAK_ENV/Institues_with_coords.xlsx")

head(env_locations)
names(env_locations)

#----------------------------------------------------------------------
# 59.2 Join FPCA-slope k-means clusters with location coordinates
#----------------------------------------------------------------------

env_locations_cluster_slope <- env_locations %>%
  left_join(env_cluster_fpca_slope, by = c("Env" = "env"))

print(env_locations_cluster_slope)

#----------------------------------------------------------------------
# 59.3 Check unmatched research locations
#----------------------------------------------------------------------

env_locations_cluster_slope %>% filter(is.na(Cluster))

#----------------------------------------------------------------------
# 59.4 Create map labels
# Example: BATEM (2)
#----------------------------------------------------------------------

env_locations_cluster_slope <- env_locations_cluster_slope %>%
  mutate(Label = paste0(Env, "\n(", Cluster, ")"))

#----------------------------------------------------------------------
# 59.5 Define cluster colors
#----------------------------------------------------------------------

okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "gold", "#0072B2",
               "#D55E00", "#CC79A7", "#000000", "#7F7F7F", "#8C510A")

n_clusters_slope <- nlevels(factor(env_locations_cluster_slope$Cluster))
cluster_colors_slope <- okabe_ito[seq_len(n_clusters_slope)]

#----------------------------------------------------------------------
# 59.6 Map FPCA-informed GDD-change k-means clusters
#----------------------------------------------------------------------

p_map_fpca_slope <- ggplot() +
  geom_sf(data = tur_prov, fill = "gray95", color = "gray80", linewidth = 0.25) +
  geom_point(data = env_locations_cluster_slope,
             aes(x = Longitude, y = Latitude, color = Cluster), size = 3) +
  geom_text_repel(data = env_locations_cluster_slope,
                  aes(x = Longitude, y = Latitude, label = Label, color = Cluster),
                  size = 2.6, max.overlaps = Inf, show.legend = FALSE) +
  scale_color_manual(values = cluster_colors_slope,
                     name = "GDD slope-informed\nk-means cluster") +
  coord_sf(xlim = c(25, 45), ylim = c(35, 43), expand = FALSE) +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude") +
  theme(legend.position = "right", panel.grid = element_blank())

p_map_fpca_slope

#----------------------------------------------------------------------
# 59.7 Save map
#----------------------------------------------------------------------

ggsave("E:/ALPER_BURAK_ENV/Turkey_Institute_Map_FPCA_Slope_Cluster.png",
       plot = p_map_fpca_slope, width = 9, height = 5, units = "in", dpi = 600)



