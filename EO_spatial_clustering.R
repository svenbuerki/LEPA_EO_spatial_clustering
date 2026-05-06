# ============================================================
# EO Spatial Connectivity & Clustering — Stage 1 C3 Hypothesis
# Goal: Identify geographic groups of locations based on 500m
#       pollinator-dispersal connectivity to proxy shared vs.
#       independent demographic history (ancestral bottleneck)
# Input:  Peggy_EOs_Germplasm_w_lat_long_from_Events_30Apr2026.csv
# Output: EO_location_groups.csv, EO_connectivity_summary.csv,
#         EO_connectivity_map.pdf/.png, EO_clustering_dendrogram.pdf
# ============================================================

library(sf)
library(igraph)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(scales)
library(maps)
library(cluster)   # silhouette
library(ggnewscale) # dual fill scales in ggplot2

set.seed(42)

out_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) ""
)
if (!nzchar(out_dir) || out_dir == ".") {
  out_dir <- "/Users/sven/Documents/Current_projects/NSF24-543_Self-incompatibility/Brainstorm_NSF/Data"
}

# ============================================================
# 1. LOAD DATA
# ============================================================
dat <- read.csv(
  file.path(out_dir, "Peggy_EOs_Germplasm_w_lat_long_from_Events_30Apr2026.csv"),
  stringsAsFactors = FALSE,
  fileEncoding = "UTF-8-BOM"
)

# Remove rows with missing coordinates
dat <- dat[!is.na(dat$eventDecimalLatitude) & !is.na(dat$eventDecimalLongitude), ]
dat$lat <- as.numeric(dat$eventDecimalLatitude)
dat$lon <- as.numeric(dat$eventDecimalLongitude)
dat <- dat[!is.na(dat$lat) & !is.na(dat$lon), ]

cat(sprintf("Records after removing missing coords: %d\n", nrow(dat)))
cat(sprintf("Unique EOs (EOCode):     %d\n", length(unique(dat$EOCode))))
cat(sprintf("Unique locations (locationID): %d\n\n", length(unique(dat$locationID))))

# ============================================================
# 2. BUILD SF OBJECTS — event points, then convex hulls per location
# ============================================================

# Event points as sf (WGS84)
pts_sf <- st_as_sf(dat, coords = c("lon", "lat"), crs = 4326)

# For each locationID: collect all event points → convex hull
# st_convex_hull on a MULTIPOINT gives a POLYGON (or POINT/LINESTRING for 1-2 pts)
loc_hulls_wgs <- lapply(split(pts_sf, pts_sf$locationID), function(sub) {
  eo   <- unique(sub$EOCode)[1]
  loc  <- unique(sub$locationID)[1]
  n    <- nrow(sub)
  geom <- st_union(sub)              # MULTIPOINT
  hull <- st_convex_hull(geom)       # POLYGON / LINESTRING / POINT
  st_sf(
    locationID = loc,
    EOCode     = eo,
    n_events   = n,
    geometry   = hull
  )
})
loc_hulls_wgs <- do.call(rbind, loc_hulls_wgs)
rownames(loc_hulls_wgs) <- NULL

cat("=== Location hull geometry types ===\n")
print(table(st_geometry_type(loc_hulls_wgs)))
cat("\n")

# Project to UTM Zone 11N (EPSG:32611) for meter-accurate distances
loc_hulls <- st_transform(loc_hulls_wgs, crs = 32611)

# Also keep event centroids per location for labelling
loc_centroids <- aggregate(
  cbind(lat = dat$lat, lon = dat$lon) ~ locationID + EOCode,
  data = dat,
  FUN  = mean
)
event_counts <- table(dat$locationID)
loc_centroids$n_events <- as.integer(
  event_counts[as.character(loc_centroids$locationID)]
)
loc_centroids$n_events[is.na(loc_centroids$n_events)] <- 1L

# ============================================================
# 3. PAIRWISE MINIMUM DISTANCES BETWEEN LOCATION HULLS (metres)
# ============================================================
n_loc    <- nrow(loc_hulls)
loc_ids  <- loc_hulls$locationID

dist_mat <- st_distance(loc_hulls)           # units: metres (UTM)
dist_mat <- matrix(as.numeric(dist_mat),
                   nrow = n_loc, ncol = n_loc,
                   dimnames = list(loc_ids, loc_ids))

cat("=== Pairwise hull-to-hull distance range (m) ===\n")
off_diag <- dist_mat[dist_mat > 0]
cat(sprintf("Min (off-diagonal): %.0f m\n", min(off_diag)))
cat(sprintf("Max:                %.0f m\n\n", max(dist_mat)))

# ============================================================
# 4. 500m CONNECTIVITY → ADJACENCY MATRIX → CONNECTED COMPONENTS
# ============================================================
THRESHOLD_M <- 500   # pollinator dispersal limit

# Exclude diagonal (self-loops) but include dist=0 pairs (overlapping hulls = connected)
adj <- dist_mat <= THRESHOLD_M
diag(adj) <- FALSE

g      <- graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
V(g)$name <- loc_ids
comps  <- components(g)

loc_hulls$group <- factor(comps$membership[match(loc_hulls$locationID,
                                                  names(comps$membership))])
loc_centroids$group <- factor(
  comps$membership[match(loc_centroids$locationID, names(comps$membership))]
)

cat(sprintf("=== Connectivity (threshold = %d m) ===\n", THRESHOLD_M))
cat(sprintf("Number of geographic groups (connected components): %d\n", comps$no))
cat(sprintf("Largest group size: %d locations\n\n", max(comps$csize)))

# Summary per group
grp_summary <- data.frame(
  group      = seq_len(comps$no),
  n_locations = comps$csize,
  EOs        = sapply(seq_len(comps$no), function(g_id) {
    locs <- names(comps$membership)[comps$membership == g_id]
    paste(sort(unique(loc_hulls$EOCode[loc_hulls$locationID %in% locs])),
          collapse = ", ")
  }),
  stringsAsFactors = FALSE
)
cat("=== Group summary ===\n")
print(grp_summary)
cat("\n")

# ============================================================
# 4b. HULL AREA & GENETIC DRIFT PROXY PER GROUP
#     Union of member location hulls (UTM) gives total connected
#     habitat area. Smaller area → smaller Ne → stronger drift.
#     Drift index: 0 = weakest drift (largest group), 1 = strongest.
# ============================================================
loc_hulls$hull_area_m2 <- as.numeric(st_area(loc_hulls))

grp_area <- do.call(rbind, lapply(
  sort(unique(as.integer(as.character(loc_hulls$group)))), function(gid) {
    sub_h      <- loc_hulls[as.integer(as.character(loc_hulls$group)) == gid, ]
    union_area <- as.numeric(st_area(st_union(sub_h)))
    data.frame(
      group       = gid,
      n_locations = nrow(sub_h),
      EOs         = paste(sort(unique(sub_h$EOCode)), collapse = ", "),
      area_m2     = round(union_area, 1),
      area_ha     = round(union_area / 10000, 4),
      stringsAsFactors = FALSE
    )
  }
))

a_min <- min(grp_area$area_ha)
a_max <- max(grp_area$area_ha)
grp_area$drift_index <- round(
  if (a_max > a_min) 1 - (grp_area$area_ha - a_min) / (a_max - a_min)
  else rep(0.5, nrow(grp_area)),
  3
)

cat("=== Group areas and drift index (highest drift first) ===\n")
print(grp_area[order(grp_area$drift_index, decreasing = TRUE), ], row.names = FALSE)
cat("\n")

write.csv(grp_area, file.path(out_dir, "EO_group_areas.csv"), row.names = FALSE)

# Attach area and drift index to loc_centroids
loc_centroids <- merge(
  loc_centroids,
  grp_area[, c("group", "area_ha", "drift_index")],
  by.x = "group", by.y = "group", all.x = TRUE
)

# ============================================================
# 4c. HIERARCHICAL CLUSTERING — computed early for shared colours
#     Macro-cluster (BL) assignments used in both the network figure
#     and dendrogram so both figures share the same colour palette.
# ============================================================
if (comps$no >= 3) {
  grp_coords <- do.call(rbind, lapply(seq_len(comps$no), function(g_id) {
    locs <- loc_centroids[loc_centroids$group == g_id, ]
    data.frame(group = g_id, lat = mean(locs$lat), lon = mean(locs$lon))
  }))
  grp_sf     <- st_as_sf(grp_coords, coords = c("lon", "lat"), crs = 4326)
  grp_sf_utm <- st_transform(grp_sf, crs = 32611)
  grp_dist   <- matrix(as.numeric(st_distance(grp_sf_utm)),
                       nrow = comps$no,
                       dimnames = list(paste0("G", seq_len(comps$no)),
                                       paste0("G", seq_len(comps$no)))) / 1000
  hc    <- hclust(as.dist(grp_dist), method = "ward.D2")
  k_max <- min(6, comps$no - 1)
  if (k_max >= 2) {
    sil_scores <- sapply(2:k_max, function(k) {
      cl  <- cutree(hc, k = k)
      sil <- silhouette(cl, as.dist(grp_dist))
      mean(sil[, "sil_width"])
    })
    names(sil_scores) <- paste0("k=", 2:k_max)
    best_k <- which.max(sil_scores) + 1
    grp_coords$macro_cluster <- factor(cutree(hc, k = best_k))
    loc_centroids$macro_cluster <- grp_coords$macro_cluster[
      match(as.integer(as.character(loc_centroids$group)), grp_coords$group)
    ]
  } else {
    best_k <- 1
    grp_coords$macro_cluster <- factor(1)
    loc_centroids$macro_cluster <- factor(1)
  }
} else {
  grp_coords <- do.call(rbind, lapply(seq_len(comps$no), function(g_id) {
    locs <- loc_centroids[loc_centroids$group == g_id, ]
    data.frame(group = g_id, lat = mean(locs$lat), lon = mean(locs$lon),
               macro_cluster = factor(1))
  }))
  loc_centroids$macro_cluster <- factor(1)
  best_k <- 1
  k_max  <- 1
  hc     <- NULL
  grp_dist   <- NULL
  sil_scores <- NULL
}
clust_cols <- setNames(
  RColorBrewer::brewer.pal(max(best_k, 3), "Set1")[seq_len(best_k)],
  as.character(seq_len(best_k))
)

# ============================================================
# 5. CONNECTIVITY EDGE LIST (for mapping)
# ============================================================
edge_list <- as.data.frame(as_edgelist(g), stringsAsFactors = FALSE)
colnames(edge_list) <- c("from", "to")

# Attach coordinates for each endpoint
edge_list <- merge(edge_list,
                   loc_centroids[, c("locationID", "lon", "lat")],
                   by.x = "from", by.y = "locationID")
edge_list <- merge(edge_list,
                   loc_centroids[, c("locationID", "lon", "lat")],
                   by.x = "to", by.y = "locationID",
                   suffixes = c("_from", "_to"))

# ============================================================
# 6. OUTPUT TABLES
# ============================================================
loc_out <- loc_centroids[order(as.integer(as.character(loc_centroids$group)),
                               loc_centroids$EOCode, loc_centroids$locationID), ]
loc_out <- merge(loc_out, grp_summary[, c("group", "EOs")],
                 by = "group", suffixes = c("", "_all_EOs"))
write.csv(loc_out, file.path(out_dir, "EO_location_groups.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# 6a. FULL PAIRWISE LOCATION CONNECTIVITY TABLE
#     All n*(n-1)/2 location pairs with distance and connectivity status
# ------------------------------------------------------------------
loc_ids_vec <- loc_hulls$locationID
pair_rows <- do.call(rbind, lapply(seq_len(n_loc - 1), function(i) {
  do.call(rbind, lapply(seq(i + 1, n_loc), function(j) {
    li <- loc_ids_vec[i]
    lj <- loc_ids_vec[j]
    gi <- as.integer(as.character(loc_hulls$group[i]))
    gj <- as.integer(as.character(loc_hulls$group[j]))
    eoi <- loc_hulls$EOCode[i]
    eoj <- loc_hulls$EOCode[j]
    d   <- dist_mat[i, j]
    data.frame(
      locationID_A   = li,
      EOCode_A       = eoi,
      group_A        = gi,
      locationID_B   = lj,
      EOCode_B       = eoj,
      group_B        = gj,
      distance_m     = round(d, 1),
      connected_500m = d <= THRESHOLD_M,
      same_EO        = eoi == eoj,
      same_group     = gi == gj,
      link_type      = ifelse(d > THRESHOLD_M, "not connected",
                       ifelse(eoi == eoj, "within-EO", "between-EO")),
      stringsAsFactors = FALSE
    )
  }))
}))
pair_rows <- pair_rows[order(pair_rows$distance_m), ]

write.csv(pair_rows, file.path(out_dir, "EO_pairwise_connectivity.csv"),
          row.names = FALSE)

cat("=== Pairwise connectivity table (connected pairs only) ===\n")
print(pair_rows[pair_rows$connected_500m, ], row.names = FALSE)
cat("\n")

# ------------------------------------------------------------------
# 6b. GROUP-TO-GROUP MINIMUM DISTANCE TABLE
#     Minimum hull-to-hull distance between every pair of groups
# ------------------------------------------------------------------
grp_ids <- sort(unique(as.integer(as.character(loc_hulls$group))))
grp_pair_rows <- do.call(rbind, lapply(seq_len(length(grp_ids) - 1), function(gi) {
  do.call(rbind, lapply(seq(gi + 1, length(grp_ids)), function(gj) {
    g1 <- grp_ids[gi]; g2 <- grp_ids[gj]
    locs1 <- which(as.integer(as.character(loc_hulls$group)) == g1)
    locs2 <- which(as.integer(as.character(loc_hulls$group)) == g2)
    min_d  <- min(dist_mat[locs1, locs2, drop = FALSE])
    eos1   <- paste(sort(unique(loc_hulls$EOCode[locs1])), collapse = ", ")
    eos2   <- paste(sort(unique(loc_hulls$EOCode[locs2])), collapse = ", ")
    data.frame(
      group_A      = g1,
      EOs_A        = eos1,
      n_locs_A     = length(locs1),
      group_B      = g2,
      EOs_B        = eos2,
      n_locs_B     = length(locs2),
      min_dist_m   = round(min_d, 1),
      connected    = min_d <= THRESHOLD_M,
      stringsAsFactors = FALSE
    )
  }))
}))
grp_pair_rows <- grp_pair_rows[order(grp_pair_rows$min_dist_m), ]

write.csv(grp_pair_rows, file.path(out_dir, "EO_group_distances.csv"),
          row.names = FALSE)

cat("=== Closest group pairs (top 15) ===\n")
print(head(grp_pair_rows, 15), row.names = FALSE)
cat("\n")

# ------------------------------------------------------------------
# 6c. CONNECTIVITY SUMMARY
# ------------------------------------------------------------------
conn_summary <- data.frame(
  metric = c(
    "Total location pairs evaluated",
    "Connected pairs (<=500m)",
    "  Within-EO connections",
    "  Between-EO connections",
    "Isolated locations (no neighbour within 500m)",
    "Number of geographic groups (connected components)",
    "Groups with >1 location",
    "Closest unconnected group pair (m)",
    "Farthest connected pair (m)"
  ),
  value = c(
    nrow(pair_rows),
    sum(pair_rows$connected_500m),
    sum(pair_rows$connected_500m & pair_rows$same_EO),
    sum(pair_rows$connected_500m & !pair_rows$same_EO),
    sum(degree(g) == 0),
    comps$no,
    sum(comps$csize > 1),
    round(min(grp_pair_rows$min_dist_m[!grp_pair_rows$connected]), 1),
    round(max(pair_rows$distance_m[pair_rows$connected_500m]), 1)
  )
)
write.csv(conn_summary, file.path(out_dir, "EO_connectivity_summary.csv"),
          row.names = FALSE)
cat("=== Connectivity summary ===\n")
print(conn_summary, row.names = FALSE)
cat("\n")

NEAR_MISS_M <- 2000   # near-miss threshold used in BL panel figure (Section 11)

# ============================================================
# 9. HIERARCHICAL CLUSTERING OF GROUPS (group centroid distances)
#    Only if ≥ 3 groups exist
# ============================================================
if (comps$no >= 3) {

  # Clustering pre-computed in Section 4c; print silhouette summary here
  if (k_max >= 2) {
    cat("=== Silhouette scores (group-level clustering) ===\n")
    print(round(sil_scores, 3))
    cat(sprintf("\nOptimal k (group clusters): %d\n\n", best_k))
  } else {
    cat("Too few groups for silhouette analysis; single macro-cluster assigned.\n\n")
  }

  # Dendrogram — relabel leaves: Group | EO(s) | location IDs
  hc_plot <- hc
  hc_plot$labels <- sapply(seq_len(comps$no), function(i) {
    eos  <- grp_summary$EOs[grp_summary$group == i]
    locs <- sort(loc_centroids$locationID[
      as.integer(as.character(loc_centroids$group)) == i
    ])
    paste0("G", i, " | ", eos, " | pop.", paste(locs, collapse = ","))
  })

  # Height at which to draw the cut line for best_k clusters
  cut_h <- if (k_max >= 2 && best_k <= comps$no - 1) {
    mean(hc$height[c(comps$no - best_k, comps$no - best_k + 1)])
  } else NA

  clust_cols <- RColorBrewer::brewer.pal(max(best_k, 3), "Set1")[seq_len(best_k)]

  # dend_plot_annotated() is defined and called in Section 10, after
  # grp_area$pop_size is available.  A minimal plain version is kept here
  # so the variable hc_plot / cut_h / cl_assign objects remain in scope.
  cl_assign <- cutree(hc_plot, k = best_k)
  cl_order  <- hc_plot$order
  first_pos <- tapply(match(seq_len(comps$no), cl_order), cl_assign, min)
  lr_rank   <- rank(first_pos, ties.method = "first")


} else {
  cat("Fewer than 3 groups — skipping hierarchical clustering of groups.\n\n")
}

# ============================================================
# 9. FINAL CONSOLE SUMMARY
# ============================================================
cat("=== Location → Group assignments ===\n")
print(loc_out[, c("EOCode", "locationID", "lat", "lon",
                  "n_events", "group")])

message("\nDone. Outputs written to: ", out_dir)
message("  EO_location_groups.csv")
message("  EO_group_areas.csv")
message("  EO_connectivity_summary.csv")
if (comps$no >= 3) message("  EO_clustering_dendrogram.pdf / .png")

# ============================================================
# 10. POPULATION SIZE DATA, BL INTEGRATION & STAGE 1 BL COMPARISON
#     Goal: load plant counts per location, sum to group level,
#     attach BL assignments, produce BL-level summary table.
# ============================================================

pop_dat <- read.csv(
  file.path(out_dir, "Data_aggregation_by_location_fecundity.csv"),
  stringsAsFactors = FALSE
)
pop_dat$pop_size <- pop_dat$total_OrganismQuantityFertile +
                    pop_dat$total_OrganismQuantityVegetative

# Merge pop_size into loc_centroids by locationID
loc_centroids <- merge(
  loc_centroids,
  pop_dat[, c("locationID", "pop_size")],
  by = "locationID", all.x = TRUE
)
loc_centroids$pop_size[is.na(loc_centroids$pop_size)] <- 0L

# Derive BL assignments directly from clustering (lr_rank gives left-to-right
# dendrogram order so BL1 = leftmost cluster, BL5 = rightmost).
# This replaces reading EO_group_BL_summary.csv, which breaks after re-runs
# that change group numbering.
bl_name_map   <- setNames(paste0("BL", lr_rank), as.character(seq_len(best_k)))
clust_to_bl   <- bl_name_map                          # cluster index → BL name
bl_strip_cols <- setNames(clust_cols[seq_len(best_k)], bl_name_map)

# Attach BL to loc_centroids via group → cutree cluster → BL name
loc_centroids$BL <- bl_name_map[
  as.character(cl_assign[as.integer(as.character(loc_centroids$group))])
]

# Group-level pop size: sum of all location pop_sizes within the group
grp_pop <- aggregate(pop_size ~ group, data = loc_centroids, FUN = sum, na.rm = TRUE)
grp_area <- merge(grp_area, grp_pop, by = "group", all.x = TRUE)
grp_area$pop_size[is.na(grp_area$pop_size)] <- 0L
grp_area$BL <- bl_name_map[as.character(cl_assign[grp_area$group])]

# BL-level summary table
bl_summary <- do.call(rbind, lapply(sort(unique(bl_name_map)), function(bl) {
  g_sub <- grp_area[!is.na(grp_area$BL) & grp_area$BL == bl, ]
  l_sub <- loc_centroids[!is.na(loc_centroids$BL) & loc_centroids$BL == bl, ]
  data.frame(
    BL                 = bl,
    n_groups           = nrow(g_sub),
    n_locations        = nrow(l_sub),
    n_EOs              = length(unique(l_sub$EOCode)),
    total_area_ha      = round(sum(g_sub$area_ha,    na.rm = TRUE), 2),
    mean_area_ha       = round(mean(g_sub$area_ha,   na.rm = TRUE), 3),
    total_pop_size     = sum(g_sub$pop_size,          na.rm = TRUE),
    mean_pop_per_group = round(mean(g_sub$pop_size,   na.rm = TRUE)),
    mean_drift_index   = round(mean(g_sub$drift_index, na.rm = TRUE), 3),
    stringsAsFactors   = FALSE
  )
}))

write.csv(bl_summary, file.path(out_dir, "EO_BL_summary.csv"), row.names = FALSE)
cat("=== BL-level summary — Stage 1 (are BLs equal in size and population?) ===\n")
print(bl_summary, row.names = FALSE)
cat("\n")

# ============================================================
# 10b. ANNOTATED DENDROGRAM — replaces plain version
#      Two-panel layout (layout()):
#        Panel 1 (top, ~70%): dendrogram, no EO circles
#        Panel 2 (bot, ~30%): per-leaf DI color strip +
#                             census-population bars
#      Both panels share the same x range → leaf-perfect alignment.
# ============================================================
if (comps$no >= 3 && !is.null(hc)) {

  n_lv      <- comps$no
  leaf_grps <- hc_plot$order   # group index at each leaf position (L→R)
  leaf_di   <- grp_area$drift_index[match(leaf_grps, grp_area$group)]
  leaf_pop  <- grp_area$pop_size[match(leaf_grps,  grp_area$group)]
  leaf_pop[is.na(leaf_pop)] <- 0L
  max_pop   <- max(leaf_pop, na.rm = TRUE)

  # DI color ramp: blue→white→red (matches panel figure and network)
  di_pal  <- colorRampPalette(c("#2166AC", "#F7F7F7", "#D6604D"))(101)
  di_cols <- di_pal[pmax(1L, round(leaf_di * 100) + 1L)]

  # Soft BL background colors (very low alpha of clust_cols Set1 palette)
  bl_bg <- adjustcolor(clust_cols, alpha.f = 0.12)

  dend_plot_annotated <- function() {

    layout(matrix(c(1, 2), nrow = 2, ncol = 1), heights = c(3.6, 1.2))

    # ── Panel 1: Dendrogram ───────────────────────────────────────────
    par(mar = c(6.5, 4.5, 5, 4.5))   # right margin matches Panel 2 for alignment
    plot(hc_plot,
         main = "", xlab = "", ylab = "Ward's D2 linkage criterion",
         sub  = "", cex = 0.68, hang = -1, axes = TRUE)
    # BL background shading — drawn after plot() to get usr coordinates;
    # same low-alpha colors as annotation panel below for visual continuity.
    usr <- par("usr")
    for (k_i in seq_len(best_k)) {
      lk  <- which(cl_assign == k_i)
      pos <- which(cl_order  %in% lk)
      rect(min(pos) - 0.5, usr[3], max(pos) + 0.5, cut_h,
           col = bl_bg[k_i], border = NA)
    }
    # Redraw BL borders on top of shading
    if (k_max >= 2) rect.hclust(hc_plot, k = best_k, border = "grey30")
    if (!is.na(cut_h)) {
      abline(h = cut_h, lty = 2, col = "grey40", lwd = 1.2)
      text(x = 0.5, y = cut_h * 1.04,
           labels = sprintf("k = %d independent bottleneck lineages", best_k),
           adj = c(0, 0), cex = 0.78, col = "grey20", font = 3)
    }
    # BL labels placed in bottom margin (line 5), below leaf labels — no overlap
    for (k_i in seq_len(best_k)) {
      lk  <- which(cl_assign == k_i)
      pos <- which(cl_order  %in% lk)
      mtext(clust_to_bl[as.character(k_i)], side = 1, line = 5,
            at = mean(range(pos)), col = "black", cex = 0.82, font = 2)
    }
    # Save usr after dendrogram so annotation panel can match x coordinates exactly
    dend_usr <<- par("usr")
    mtext(
      "Independent ancestral bottleneck lineages of Lepidium papilliferum",
      side = 3, line = 3.2, cex = 1.00, font = 2, adj = 0.5
    )
    mtext(
      sprintf(
        "Ward's D2 hierarchical clustering of %d geographic groups  |  Silhouette-optimal k = %d  |  C3 Hypothesis Stages 1 & 4",
        comps$no, best_k
      ),
      side = 3, line = 1.9, cex = 0.75, col = "grey30", adj = 0.5
    )

    # ── Panel 2: Annotation strips ────────────────────────────────────
    # Use dend_usr[1:2] so leaf x-positions match exactly across panels
    par(mar = c(1, 4.5, 0.3, 4.5))   # wide right margin for DI legend bar
    plot(NA,
         xlim = dend_usr[1:2], ylim = c(0, 1),
         axes = FALSE, xlab = "", ylab = "",
         xaxs = "i", yaxs = "i")

    # Soft BL background shading — same colors as dendrogram panel above;
    # no BL text labels here (they are already in the dendrogram).
    for (k_i in seq_len(best_k)) {
      lk  <- which(cl_assign == k_i)
      pos <- which(cl_order  %in% lk)
      rect(min(pos) - 0.49, 0, max(pos) + 0.49, 1,
           col = bl_bg[k_i], border = NA)
    }

    # Connecting tick marks: vertical lines from top of annotation panel
    # upward into the margin, aligned with dendrogram leaf positions
    par(xpd = TRUE)
    for (i in seq_len(n_lv)) {
      segments(i, 1, i, 1.25, col = "grey55", lwd = 0.5)
    }
    par(xpd = FALSE)

    # ─ Row 1 (y 0.54–0.98): Drift index color strip ─────────────────
    # Expanded upward now that BL text labels are removed.
    for (i in seq_len(n_lv)) {
      rect(i - 0.44, 0.54, i + 0.44, 0.98,
           col = di_cols[i], border = "white", lwd = 0.3)
    }
    # Label only groups with DI <= 0.75 (the low-drift exceptions)
    for (i in seq_len(n_lv)) {
      if (leaf_di[i] <= 0.75) {
        text(i, 0.76, sprintf("%.2f", leaf_di[i]),
             cex = 0.50, col = "grey15", srt = 90, font = 2)
      }
    }
    axis(2, at = 0.76, labels = "Drift\nindex", tick = FALSE,
         las = 1, cex.axis = 0.62, col.axis = "grey25", line = 0.2)

    # ─ Row 2 (y 0.04–0.46): Census population bars ───────────────────
    for (i in seq_len(n_lv)) {
      bh  <- (leaf_pop[i] / max_pop) * 0.40
      # Darker grey for larger populations (dark at max, light grey at min)
      pg  <- 0.55 - 0.30 * (leaf_pop[i] / max_pop)
      rect(i - 0.38, 0.04, i + 0.38, 0.04 + bh,
           col = grey(pg), border = "white", lwd = 0.3)
      # Label bars with N >= 300
      if (leaf_pop[i] >= 300) {
        lbl <- if (leaf_pop[i] >= 1000)
          sprintf("%.1fk", leaf_pop[i] / 1000) else as.character(leaf_pop[i])
        text(i, 0.04 + bh + 0.01, lbl,
             cex = 0.42, col = "grey25", adj = c(0.5, 0), srt = 90)
      }
    }
    # Max-N reference line + label
    abline(h = 0.44, lty = 3, col = "grey65", lwd = 0.8)
    text(n_lv + 0.46, 0.46,
         sprintf("max N = %d", max_pop),
         cex = 0.48, col = "grey45", adj = c(1, 0))
    axis(2, at = 0.25, labels = "Census\nsize N", tick = FALSE,
         las = 1, cex.axis = 0.62, col.axis = "grey25", line = 0.2)

    # Dividing line between the two strips (white gap)
    abline(h = 0.52, col = "white", lwd = 3)

    # ─ DI color scale legend (vertical gradient bar in right margin) ──────
    par(xpd = TRUE)
    plt <- par("plt")   # plot region in NFC: c(x1, x2, y1, y2)
    # Place bar 15–50% of the way into the right margin
    bar_x0 <- grconvertX(plt[2] + (1 - plt[2]) * 0.15, from = "nfc", to = "user")
    bar_x1 <- grconvertX(plt[2] + (1 - plt[2]) * 0.50, from = "nfc", to = "user")
    lbl_x  <- grconvertX(plt[2] + (1 - plt[2]) * 0.58, from = "nfc", to = "user")
    n_grad <- 40
    gy     <- seq(0.54, 0.98, length.out = n_grad + 1)
    for (gi in seq_len(n_grad)) {
      di_val <- (gi - 0.5) / n_grad
      rect(bar_x0, gy[gi], bar_x1, gy[gi + 1],
           col = di_pal[round(di_val * 100) + 1L], border = NA)
    }
    rect(bar_x0, 0.54, bar_x1, 0.98, col = NA, border = "grey30", lwd = 0.6)
    text(lbl_x, 0.98,  "1",   adj = c(0, 1),   cex = 0.55, col = "grey20")
    text(lbl_x, 0.76,  "0.5", adj = c(0, 0.5), cex = 0.55, col = "grey20")
    text(lbl_x, 0.54,  "0",   adj = c(0, 0),   cex = 0.55, col = "grey20")
    text((bar_x0 + bar_x1) / 2, 0.52, "DI",
         adj = c(0.5, 1), cex = 0.60, col = "grey30", font = 3)
    par(xpd = FALSE)

    # Caption
    mtext(
      paste0(
        "Annotation strips aligned leaf-by-leaf with the dendrogram above.  ",
        "DI strip: blue = large habitat (DI = 0, weak drift),  ",
        "red = small habitat (DI = 1, strong drift);  labeled where DI \u2264 0.75.  ",
        sprintf("Census bar: fertile + vegetative plant count per group (max N = %d).", max_pop)
      ),
      side = 1, line = 0.2, cex = 0.56, col = "grey35", adj = 0.5
    )
  }

  pdf(file.path(out_dir, "EO_clustering_dendrogram.pdf"),
      width = 16, height = 9, onefile = FALSE)
  dend_plot_annotated()
  dev.off()

  png(file.path(out_dir, "EO_clustering_dendrogram.png"),
      width = 16, height = 9, units = "in", res = 300)
  dend_plot_annotated()
  dev.off()

  message("  EO_clustering_dendrogram.pdf / .png (annotated: DI + census N strips)")
}

# ============================================================
# 11. BL-PANEL NETWORK FIGURE — Stage 4: Genetic Drift Predictions
#     Five panels (one per BL); nodes = locations (not groups);
#     shape = isolated (circle) vs. connected (diamond);
#     fill = group drift index; size = location population size;
#     near-miss edges drawn behind nodes;
#     connected edges drawn ON TOP of nodes.
# ============================================================

# Tag pair_rows with BL for within-BL filtering
pair_rows$BL_A <- loc_centroids$BL[
  match(pair_rows$locationID_A, loc_centroids$locationID)
]
pair_rows$BL_B <- loc_centroids$BL[
  match(pair_rows$locationID_B, loc_centroids$locationID)
]

# Per-BL: FR layout + node/edge data frames
bl_node_list2 <- list()
bl_conn_list2 <- list()
bl_near_list2 <- list()
bl_far_list2  <- list()

for (bl in sort(unique(loc_centroids$BL[!is.na(loc_centroids$BL)]))) {

  bl_locs <- loc_centroids$locationID[
    !is.na(loc_centroids$BL) & loc_centroids$BL == bl
  ]

  # Induced subgraph (connected edges only, <=500 m)
  sub_g <- induced_subgraph(g, V(g)[V(g)$name %in% as.character(bl_locs)])

  # Connected and near-miss location pairs within this BL
  in_bl <- !is.na(pair_rows$BL_A) & !is.na(pair_rows$BL_B) &
           pair_rows$BL_A == bl & pair_rows$BL_B == bl
  bl_conn <- pair_rows[in_bl & pair_rows$connected_500m, ]
  bl_near <- pair_rows[in_bl & !pair_rows$connected_500m &
                       pair_rows$distance_m <= NEAR_MISS_M, ]

  # Far pairs (>2 km, within BL): keep only the closest location pair per
  # disconnected group-pair so labels stay readable
  bl_far_raw <- pair_rows[in_bl & !pair_rows$connected_500m &
                          pair_rows$distance_m > NEAR_MISS_M, ]
  if (nrow(bl_far_raw) > 0) {
    bl_far_raw$grp_pair <- paste0(
      pmin(bl_far_raw$group_A, bl_far_raw$group_B), "_",
      pmax(bl_far_raw$group_A, bl_far_raw$group_B)
    )
    bl_far_raw <- bl_far_raw[order(bl_far_raw$distance_m), ]
    bl_far <- bl_far_raw[!duplicated(bl_far_raw$grp_pair), ]
  } else {
    bl_far <- bl_far_raw
  }

  # Build layout graph: connected (weight 2) + near-miss (weight 0.4) +
  # far (weight 0.1 — pulls disconnected groups loosely together in layout)
  g_bl <- sub_g
  extra_edges <- rbind(
    if (nrow(bl_near) > 0) bl_near[, c("locationID_A", "locationID_B")] else NULL,
    if (nrow(bl_far)  > 0) bl_far[,  c("locationID_A", "locationID_B")] else NULL
  )
  if (!is.null(extra_edges) && nrow(extra_edges) > 0) {
    g_bl <- add_edges(
      g_bl,
      as.vector(rbind(as.character(extra_edges$locationID_A),
                      as.character(extra_edges$locationID_B)))
    )
  }
  E(g_bl)$weight <- c(rep(2.0, ecount(sub_g)),
                       rep(0.4, nrow(bl_near)),
                       rep(0.1, nrow(bl_far)))

  # FR layout, normalised to [0, 1]
  set.seed(42)
  if (vcount(g_bl) > 1) {
    lay_bl <- layout_with_fr(g_bl, weights = E(g_bl)$weight)
  } else {
    lay_bl <- matrix(c(0.5, 0.5), nrow = 1)
  }
  rownames(lay_bl) <- V(g_bl)$name
  xr <- range(lay_bl[, 1]); yr <- range(lay_bl[, 2])
  lay_bl[, 1] <- (lay_bl[, 1] - xr[1]) / max(diff(xr), 1e-3)
  lay_bl[, 2] <- (lay_bl[, 2] - yr[1]) / max(diff(yr), 1e-3)

  # Node data frame: one row per location
  nd <- data.frame(
    locationID = as.integer(rownames(lay_bl)),
    x          = lay_bl[, 1],
    y          = lay_bl[, 2],
    stringsAsFactors = FALSE
  )
  nd <- merge(nd,
              loc_centroids[, c("locationID", "EOCode", "group",
                                "n_events", "pop_size", "BL")],
              by = "locationID")
  nd <- merge(nd,
              grp_area[, c("group", "area_ha", "drift_index")],
              by = "group")
  nd$isolated <- degree(sub_g)[as.character(nd$locationID)] == 0
  nd$label    <- paste0(nd$EOCode, "\n(", nd$locationID, ")\n",
                        nd$pop_size, " ind.")
  bl_node_list2[[bl]] <- nd

  # Helper: attach layout coords to edge table
  add_lay2 <- function(pairs) {
    if (nrow(pairs) == 0) return(NULL)
    pairs$x_A   <- lay_bl[as.character(pairs$locationID_A), 1]
    pairs$y_A   <- lay_bl[as.character(pairs$locationID_A), 2]
    pairs$x_B   <- lay_bl[as.character(pairs$locationID_B), 1]
    pairs$y_B   <- lay_bl[as.character(pairs$locationID_B), 2]
    pairs$mid_x <- (pairs$x_A + pairs$x_B) / 2
    pairs$mid_y <- (pairs$y_A + pairs$y_B) / 2
    pairs$BL    <- bl
    pairs
  }
  bl_conn_list2[[bl]] <- add_lay2(bl_conn)
  bl_near_list2[[bl]] <- add_lay2(bl_near)
  bl_far_list2[[bl]]  <- add_lay2(bl_far)
}

all_nodes2 <- do.call(rbind, bl_node_list2)
all_conn2  <- do.call(rbind, Filter(Negate(is.null), bl_conn_list2))
if (is.null(all_conn2)) all_conn2 <- data.frame()
all_near2  <- do.call(rbind, Filter(Negate(is.null), bl_near_list2))
if (is.null(all_near2)) all_near2 <- data.frame()
all_far2   <- do.call(rbind, Filter(Negate(is.null), bl_far_list2))
if (is.null(all_far2))  all_far2  <- data.frame()

# Facet labels: BL name + summary statistics
bl_labels2 <- setNames(
  sprintf(
    "%s  |  %d groups \u00b7 %d pop. \u00b7 %d EOs\n%.2f ha total \u00b7 %d ind. total",
    bl_summary$BL, bl_summary$n_groups, bl_summary$n_locations,
    bl_summary$n_EOs, bl_summary$total_area_ha, bl_summary$total_pop_size
  ),
  bl_summary$BL
)
all_nodes2$BL_label <- factor(bl_labels2[all_nodes2$BL], levels = bl_labels2)
if (nrow(all_conn2) > 0)
  all_conn2$BL_label <- factor(bl_labels2[all_conn2$BL], levels = bl_labels2)
if (nrow(all_near2) > 0)
  all_near2$BL_label <- factor(bl_labels2[all_near2$BL], levels = bl_labels2)
if (nrow(all_far2) > 0)
  all_far2$BL_label  <- factor(bl_labels2[all_far2$BL],  levels = bl_labels2)

p_bl_net <- ggplot(all_nodes2, aes(x = x, y = y)) +

  # 0. Far edges (>2 km within BL) — only the single shortest edge per BL,
  #    drawn as a dashed grey line with a km label
  {if (nrow(all_far2) > 0) {
    all_far2_min <- do.call(rbind, lapply(split(all_far2, all_far2$BL_label), function(x)
      x[which.min(x$distance_m), ]
    ))
    list(
      geom_segment(data = all_far2_min,
                   aes(x = x_A, y = y_A, xend = x_B, yend = y_B),
                   color = "grey55", linewidth = 0.6,
                   linetype = "dashed", alpha = 0.80,
                   inherit.aes = FALSE),
      geom_label(data = all_far2_min,
                 aes(x = mid_x, y = mid_y,
                     label = paste0(round(distance_m / 1000, 1), " km")),
                 size = 3.2, fill = "white", color = "grey20",
                 linewidth = 0.3, alpha = 0.95, inherit.aes = FALSE,
                 label.padding = unit(0.12, "lines"), fontface = "bold")
    )
  }} +

  # 1. Near-miss edges — solid grey, behind nodes
  {if (nrow(all_near2) > 0)
    geom_segment(data = all_near2,
                 aes(x = x_A, y = y_A, xend = x_B, yend = y_B),
                 color       = "grey55", linewidth = 0.6,
                 linetype    = "solid", alpha = 0.70,
                 inherit.aes = FALSE)
  } +


  # 2. Nodes: shape = isolation status; fill = drift index; size = pop size
  geom_point(aes(fill = drift_index, size = pop_size, shape = isolated),
             color = "grey20", stroke = 0.65, alpha = 0.93) +

  scale_fill_gradientn(
    colours = c("#2166AC", "#F7F7F7", "#D6604D"),
    name    = "Drift index\n(area proxy)",
    limits  = c(0, 1),
    breaks  = c(0, 0.5, 1),
    labels  = c("0\n(large\nhabitat)", "0.5", "1\n(small\nhabitat)")
  ) +
  scale_size_area(
    name     = "Population\nsize (N)",
    max_size = 15,
    breaks   = c(50, 200, 500, 1000)
  ) +
  scale_shape_manual(
    values = c("TRUE" = 21, "FALSE" = 23),
    name   = "Status",
    labels = c("TRUE" = "Isolated", "FALSE" = "Connected")
  ) +

  # 3. Connected edges — ON TOP of nodes
  {if (nrow(all_conn2) > 0)
    geom_segment(data = all_conn2,
                 aes(x = x_A, y = y_A, xend = x_B, yend = y_B,
                     color = link_type),
                 linewidth = 1.6, alpha = 0.90, inherit.aes = FALSE)
  } +

  scale_color_manual(
    values = c("within-EO"  = "#1A5276",
               "between-EO" = "#B03A2E"),
    name   = "Connection",
    labels = c(
      "within-EO"  = paste0("Within-EO (\u2264", THRESHOLD_M, " m)"),
      "between-EO" = paste0("Between-EO (\u2264", THRESHOLD_M, " m)")
    ),
    drop = FALSE
  ) +

  # 4. Distance labels on connected edges
  {if (nrow(all_conn2) > 0)
    geom_label(data = all_conn2,
               aes(x = mid_x, y = mid_y,
                   label = paste0(round(distance_m), " m")),
               size = 3.0, fill = "white", color = "grey10",
               linewidth = 0.25, alpha = 0.92, inherit.aes = FALSE,
               label.padding = unit(0.10, "lines"))
  } +

  # 5. Node labels (white background, repelled)
  geom_label_repel(
    aes(label = label),
    size          = 2.1,
    fill          = "white",
    color         = "grey15",
    box.padding   = 0.50,
    point.padding = 0.35,
    label.size    = NA,
    alpha         = 0.88,
    max.overlaps  = 60,
    lineheight    = 0.80,
    min.segment.length = 0.2
  ) +

  guides(
    fill  = guide_colorbar(order = 1),
    size  = guide_legend(order = 2,
                         override.aes = list(shape = 21, fill = "grey60")),
    shape = guide_legend(order = 3,
                         override.aes = list(fill = "grey60", size = 4)),
    color = guide_legend(order = 4,
                         override.aes = list(linewidth = 1.5))
  ) +

  facet_wrap(~ BL_label, nrow = 1, scales = "free") +

  labs(
    title    = paste0(
      "Predicted genetic drift intensity across five independent ",
      "bottleneck lineages of Lepidium papilliferum"
    ),
    subtitle = paste0(
      "C3 Hypothesis \u2014 Stage 4  |  39 locations in 5 BLs  |  ",
      "Fill = drift index (red = strong drift; blue = weak drift)  |  ",
      "Node size = population size  |  \u25c6 connected  \u25cf isolated\n",
      "Colored solid: \u2264500 m (connected, pollinator range) \u2014 ",
      "Grey solid: >500 m to \u22642 km (near-miss, no gene flow) \u2014 ",
      "Dashed + label: shortest distance >2 km between disconnected groups within BL"
    ),
    caption = paste0(
      "No between-BL connections exist (all BLs separated by >10 km).  ",
      "Drift index: DI = 1 \u2212 (area \u2212 min) / (max \u2212 min); ",
      "0 = largest habitat (weakest drift), 1 = smallest habitat (strongest drift).  ",
      "Population size = fertile + vegetative plant counts per location.  ",
      "C3 Hypothesis \u2014 Stage 4."
    ),
    x = NULL, y = NULL
  ) +
  theme_void(base_size = 10) +
  theme(
    strip.text       = element_text(face = "bold", size = 8, lineheight = 1.15,
                                    margin = margin(b = 5, t = 5)),
    strip.background = element_rect(fill = "white", color = "grey55",
                                    linewidth = 0.6),
    panel.border     = element_rect(color = "grey55", fill = NA,
                                    linewidth = 0.9),
    panel.spacing    = unit(0.9, "cm"),
    plot.title       = element_text(face = "bold", size = 11, hjust = 0.5,
                                    margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 8.5, color = "grey35", hjust = 0.5,
                                    margin = margin(b = 6)),
    plot.caption     = element_text(size = 7, color = "grey50", hjust = 0.5,
                                    margin = margin(t = 6)),
    legend.position  = "right",
    plot.background  = element_rect(fill = "white", color = NA),
    plot.margin      = margin(12, 12, 12, 12)
  )

# Color each BL strip background to match the dendrogram cluster colors.
# Strip grob tree: strip gtable → grobs[[1]] (gTree) → $children → background rect
# (named "strip.background.x..rect.*" in ggplot2 >= 3.4).
gt_bl <- ggplotGrob(p_bl_net)
strip_idx <- sort(which(grepl("^strip-t", gt_bl$layout$name)))
bl_facet_order <- levels(all_nodes2$BL_label)
bl_name_order  <- sub("\\s.*", "", bl_facet_order)   # extract "BL1", "BL2", …
for (i in seq_along(strip_idx)) {
  j     <- strip_idx[i]
  col_i <- adjustcolor(bl_strip_cols[bl_name_order[i]], alpha.f = 0.45)
  bdr_i <- adjustcolor(bl_strip_cols[bl_name_order[i]], alpha.f = 0.85)
  inner <- gt_bl$grobs[[j]]$grobs[[1]]          # the gTree inside the strip gtable
  bg_nm <- grep("^strip\\.background", names(inner$children), value = TRUE)[1]
  if (!is.na(bg_nm)) {
    inner$children[[bg_nm]]$gp$fill <- col_i
    inner$children[[bg_nm]]$gp$col  <- bdr_i
    gt_bl$grobs[[j]]$grobs[[1]] <- inner
  }
}

cairo_pdf(file.path(out_dir, "EO_BL_drift_panel.pdf"), width = 24, height = 8)
grid::grid.draw(gt_bl)
dev.off()
png(file.path(out_dir, "EO_BL_drift_panel.png"),
    width = 24, height = 8, units = "in", res = 300, type = "cairo")
grid::grid.draw(gt_bl)
dev.off()

message("  EO_BL_summary.csv")
message("  EO_BL_drift_panel.pdf / .png")
