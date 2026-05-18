library(data.table)
library(ggplot2)
library(tidyverse)
library(rpart)

sample_id <- "GM19462"
base_dir <- "~/cpgskew/female_samples_ready"

cat("Loading CpG BED...\n")

cpg_df <- fread("CGIs_CHM13v2.0_chrX.bed", header = FALSE) %>%
  as_tibble() %>%
  rename(
    chr = V1,
    start = V2,
    end = V3
  ) %>%
  filter(chr == "chrX") %>%
  mutate(
    CGI_id = paste0(chr, ":", start, "-", end),
    cpg_id = row_number(),
    cpg_start = start,
    cpg_end = end
  ) %>%
  arrange(start)

cat("Total CpGs:", nrow(cpg_df), "\n")

cat("Loading ps block info from whatshap...\n")

blocks_file <- file.path(base_dir, paste0(sample_id, "_blocks.tsv"))

ps_blocks <- fread(blocks_file) %>%
  as_tibble() %>%
  filter(chromosome == "chrX") %>%
  rename(PS = phase_set, ps_start = from, ps_end = to) %>%
  mutate(
    PS       = as.numeric(PS),
    ps_start = as.numeric(ps_start),
    ps_end   = as.numeric(ps_end)
  ) %>%
  dplyr::select(PS, ps_start, ps_end) %>%
  arrange(ps_start)

cat("Loaded PS blocks:", nrow(ps_blocks), "\n")

head(ps_blocks)

cat("Mapping CpG ↔ PS...\n")

cpg_ps_map <- cpg_df %>%
  dplyr::select(cpg_id, CGI_id, cpg_start, cpg_end) %>%
  left_join(ps_blocks, by = join_by(overlaps(cpg_start, cpg_end, ps_start, ps_end)))

# ADD BACKBONE FLAG HERE

cpg_ps_map <- cpg_ps_map %>%
  group_by(CGI_id) %>%
  mutate(
    ps_relation = case_when(
      is.na(PS) ~ "no_PS",
      n() > 1 ~ "multi_PS",
      TRUE ~ "single_PS"
    )
  ) %>%
  ungroup()

################################################################################
#### sanity check for backbone
##############

cat("\n--- CpG ↔ PS Mapping Summary ---\n")

total_cpg  <- nrow(cpg_df)
mapped_cpg <- cpg_ps_map %>% filter(ps_relation != "no_PS") %>% distinct(CGI_id) %>% nrow()

cat("Total CpGs:", total_cpg, "\n")
cat("Mapped CpGs:", mapped_cpg, "\n")
cat("Unmapped CpGs:", total_cpg - mapped_cpg, "\n")
cat("Multi-PS CpGs:", cpg_ps_map %>% filter(ps_relation == "multi_PS") %>% distinct(CGI_id) %>% nrow(), "\n")

cat("Unmapped CpG IDs:\n")
cpg_ps_map %>% filter(ps_relation == "no_PS") %>% distinct(CGI_id) %>% head(15) %>% print()

ps_gaps <- ps_blocks %>%
  arrange(ps_start) %>%
  mutate(next_start = lead(ps_start), gap = next_start - ps_end) %>%
  filter(gap > 0)

cat("Number of PS gaps:", nrow(ps_gaps), "\n")
ps_gaps %>% head(10)

################################################################################

cat("Loading reads...\n")

hp_file <- file.path(base_dir, paste0(sample_id, "_chrX_Y_18_21.hp_CGIX_hpreads.tsv.gz"))
hp_df <- fread(hp_file) %>%
  as_tibble() %>%
  rename(
    read_id = V1,
    HP = V2,
    PS = V3
  )

meth_file <- file.path(base_dir, paste0(sample_id, "_chrX_CGIX_clustered_reads.tsv.gz"))
meth_df <- fread(meth_file) %>%
  as_tibble() %>%
  rename(
    read_id = read_name,
    cluster = cluster_id,
    meth_start = start,
    meth_end = end
  ) %>%
  dplyr::select(read_id, cluster, meth_start, meth_end, CGI_id, assigned_X)

reads_df <- meth_df %>%
  left_join(
    cpg_ps_map %>% distinct(CGI_id, PS),
    by = "CGI_id"
  ) %>%
  left_join(hp_df, by = c("read_id", "PS")) %>%
  filter(!is.na(HP))

################################################################################

reads_df <- reads_df %>%
  group_by(CGI_id) %>%
  mutate(
    n_clusters    = n_distinct(cluster[!is.na(PS)]),
    cluster_valid = (n_clusters == 2)
  ) %>%
  ungroup()

reads_df <- reads_df %>%
  mutate(
    usable_for_skew = case_when(
      !is.na(PS) &
        cluster_valid &
        assigned_X %in% c("Xa", "Xi") ~ TRUE,
      TRUE ~ FALSE
    )
  )

# Sanity check
cat("CGIs usable for skew:", reads_df %>% filter(usable_for_skew) %>% distinct(CGI_id) %>% nrow(), "\n")
cat("Cluster distribution:\n")
reads_df %>%
  filter(!is.na(PS)) %>%
  distinct(CGI_id, n_clusters) %>%
  dplyr::count(n_clusters) %>%
  print()
##############

################################################################################

cat("Computing skew...\n")

counts_cpg <- reads_df %>%
  filter(usable_for_skew) %>%
  group_by(CGI_id, assigned_X, HP) %>%
  summarise(count = n(), .groups = "drop")

combi_map <- c(
  "Xa_1" = "H1_Xa",
  "Xa_2" = "H2_Xa",
  "Xi_1" = "H1_Xi",
  "Xi_2" = "H2_Xi"
)

skew_cpg <- counts_cpg %>%
  unite(combi, assigned_X, HP) %>%
  mutate(combi = combi_map[combi]) %>%
  pivot_wider(
    id_cols = CGI_id,
    names_from = combi,
    values_from = count,
    values_fill = 0
  )

for(col in c("H1_Xa","H1_Xi","H2_Xa","H2_Xi")){
  if(!col %in% colnames(skew_cpg)){
    skew_cpg[[col]] <- 0
  }
}

skew_cpg <- skew_cpg %>%
  mutate(
    total = H1_Xa + H1_Xi + H2_Xa + H2_Xi,
    skew = (H1_Xa + H2_Xi) / total
  )

######

final_df <- cpg_df %>%
  dplyr::select(CGI_id, cpg_start, cpg_end) %>%
  left_join(skew_cpg, by = "CGI_id") %>%
  left_join(                                      
    cpg_ps_map %>%                                
      filter(ps_relation == "single_PS") %>%                      
      distinct(CGI_id, PS, ps_start, ps_end),     
    by = "CGI_id"                                 
  ) %>%                                           
  mutate(skew = replace_na(skew, -0.1))


################################################################################

######### splitting ps blocks using R part regression tree
##############################################################

block_rpart <- final_df %>%
  filter(!is.na(PS), skew >= 0) %>%
  group_by(PS, ps_start, ps_end) %>%
  filter(dplyr::n() >= 10) %>%
  group_modify(~ {
    
    df <- .x %>%
      arrange(cpg_start)
    
    fit <- rpart(
      skew ~ cpg_start,
      data = df,
      method = "anova",
      control = rpart.control(
        cp = 0.02,
        minsplit = 8,
        maxdepth = 1
      )
    )
    
    # no split happened
    if(nrow(fit$frame) == 1){
      return(tibble(
        needs_split = FALSE,
        boundary_pos = NA_real_,
        left_mean = NA_real_,
        right_mean = NA_real_
      ))
    }
    
    split_val <- fit$splits[1, "index"]
    
    left_mean  <- mean(df$skew[df$cpg_start <= split_val], na.rm = TRUE)
    right_mean <- mean(df$skew[df$cpg_start >  split_val], na.rm = TRUE)
    
    # require real difference
    if (!((left_mean > 0.7 & right_mean < 0.3) | 
          (left_mean < 0.3 & right_mean > 0.7))){
      return(tibble(
        needs_split = FALSE,
        boundary_pos = NA_real_,
        left_mean = left_mean,
        right_mean = right_mean
      ))
    }
    
    tibble(
      needs_split = TRUE,
      boundary_pos = as.numeric(split_val),
      left_mean = left_mean,
      right_mean = right_mean
    )
    
  }) %>%
  ungroup()

# inspect splits
block_rpart %>% filter(needs_split) %>% arrange(ps_start)

################################################################################

# -------------------------
# 1. Prepare CpG plot data
# -------------------------
plot_df <- final_df %>%
  mutate(
    midpoint = (cpg_start + cpg_end) / 2,
    pos_mb = midpoint / 1e6,
    skew_type = case_when(
      skew == -0.1 ~ "invalid",
      TRUE ~ "valid"
    )
  )

# -------------------------
# 2. Prepare PS blocks
# -------------------------
ps_plot <- ps_blocks %>%
  mutate(
    ps_start_mb = ps_start / 1e6,
    ps_end_mb   = ps_end   / 1e6
  )

# Get PS blocks that have at least one valid skew point within them
ps_with_data <- ps_plot %>%
  inner_join(
    plot_df %>% 
      filter(skew_type == "valid") %>%
      distinct(pos_mb),
    by = join_by(ps_start_mb <= pos_mb, ps_end_mb >= pos_mb)
  ) %>%
  distinct(PS, ps_start_mb, ps_end_mb)

# -------------------------
# 3. CpG index mapping (for secondary axis)
# -------------------------
cpg_axis_map <- plot_df %>%
  arrange(pos_mb)%>%
  mutate(cpg_id = row_number())


split_boundaries <- block_rpart %>%
  filter(needs_split) %>%
  mutate(boundary_mb = boundary_pos / 1e6)


# -------------------------
# 4. Plot
# -------------------------
idx <- unique(round(seq(1, nrow(cpg_axis_map), length.out = 10)))
max_x <- max(
  cpg_df$cpg_end,
  ps_blocks$ps_end
) / 1e6


p <- ggplot() +
  
  # --- PS BLOCK BOUNDARIES (only blocks with valid data) ---
  geom_vline(
    data = ps_with_data %>%
      pivot_longer(
        c(ps_start_mb, ps_end_mb),
        names_to = "boundary_type",
        values_to = "boundary"
      ),
    aes(xintercept = boundary, color = boundary_type),
    linewidth = 0.3,
    alpha = 0.6,
    linetype = "dashed"
  ) +
  scale_color_manual(
    values = c(
      "ps_start_mb" = "steelblue",
      "ps_end_mb"   = "tomato"
    ),
    labels = c(
      "ps_start_mb" = "PS Start",
      "ps_end_mb"   = "PS End"
    ),
    name = "PS Boundary"
  ) +
  
  # --- invalid points ---
  geom_point(
    data = plot_df %>% filter(skew_type == "invalid"),
    aes(x = pos_mb, y = skew),
    color = "grey80",
    size = 0.5,
    alpha = 0.4
  ) +
  
  # --- valid points ---
  geom_point(
    data = plot_df %>% filter(skew_type == "valid"),
    aes(x = pos_mb, y = skew),
    color = "blue",
    size = 1.2,
    alpha = 0.9
  ) +
  
  labs(
    x = "Genomic Position (Mb)",
    y = "Skew",
    title = "CpG-level Skew across chrX with PS Gap Regions and split using Rpart"
  ) +
  
  scale_y_continuous(
    limits = c(-0.2, 1.1),
    breaks = c(-0.1, 0, 0.5, 1)
  ) +
  
  scale_x_continuous(
    name = "Genomic Position (Mb)",
    limits = c(0, max_x),
    breaks = unique(c(
      pretty(c(0, max_x), n = 10),
      max_x
    )),
    expand = c(0, 0),
    sec.axis = sec_axis(
      transform = ~ .,
      breaks = cpg_axis_map$pos_mb[idx],
      labels = cpg_axis_map$cpg_id[idx],
      name = "CpG Index"
    )
  ) +
  
  theme_minimal() +
  theme(
    axis.text.x.top = element_text(angle = 45, hjust = 0, size = 8)
  )+
  geom_vline(                  # add this as a new layer
    data = split_boundaries,
    aes(xintercept = boundary_mb),
    color     = "darkgreen",
    linewidth = 0.4,
    linetype  = "solid",
    alpha     = 0.8
  )

ggsave("skew_plot.png", plot = p, width = 20, height = 8, dpi = 300)

# -------------------------
# 5. Show plot
# -------------------------
print(p)


################################################################################
######### Build updated PS blocks
################################################################################

cat("Building updated PS blocks...\n")

# blocks that don't need splitting — keep as is
ps_unsplit <- ps_blocks %>%
  filter(!PS %in% (block_rpart %>% filter(needs_split) %>% pull(PS)))

# find which CGI the boundary falls in or between
boundary_cgi <- final_df %>%
  filter(!is.na(PS)) %>%
  inner_join(
    block_rpart %>% filter(needs_split) %>% dplyr::select(PS, boundary_pos),
    by = "PS"
  ) %>%
  group_by(PS) %>%
  arrange(cpg_start) %>%
  mutate(
    boundary_inside_cgi   = boundary_pos >= cpg_start & boundary_pos <= cpg_end,
    last_left_cgi_end     = cpg_end[max(which(cpg_start <= boundary_pos))],
    first_right_cgi_start = cpg_start[min(which(cpg_start > boundary_pos))]
  ) %>%
  distinct(PS, boundary_pos, boundary_inside_cgi, last_left_cgi_end, first_right_cgi_start) %>%
  ungroup()
  
# blocks that need splitting — create two rows each
ps_split_left <- boundary_cgi %>%
  left_join(block_rpart %>% filter(needs_split) %>% dplyr::select(PS, ps_start, ps_end), by = "PS") %>%
  transmute(
    PS       = PS,
    ps_start = ps_start,
    ps_end   = case_when(
      boundary_inside_cgi ~ as.numeric(boundary_pos),
      TRUE                ~ last_left_cgi_end
    ),
    split = "left",
    boundary_inside_cgi = boundary_inside_cgi
  )

ps_split_right <- boundary_cgi %>%
  left_join(block_rpart %>% filter(needs_split) %>% dplyr::select(PS, ps_start, ps_end), by = "PS") %>%
  transmute(
    PS       = case_when(
      boundary_inside_cgi ~ as.numeric(boundary_pos) + 1,
      TRUE                ~ first_right_cgi_start
    ),
    ps_start = case_when(
      boundary_inside_cgi ~ as.numeric(boundary_pos) + 1,
      TRUE                ~ first_right_cgi_start
    ),
    ps_end   = ps_end,
    split    = "right",
    boundary_inside_cgi = boundary_inside_cgi
  )

ps_blocks_updated <- bind_rows(
  ps_unsplit %>% mutate(split = "none"),
  ps_split_left,
  ps_split_right
) %>%
  arrange(ps_start)

cat("Original PS blocks:", nrow(ps_blocks), "\n")
cat("Updated PS blocks:", nrow(ps_blocks_updated), "\n")

################################################################################
######### Reassign reads to updated PS blocks via CGI midpoint
################################################################################

cat("Reassigning reads to updated PS blocks...\n")

reads_for_ps <- reads_df %>%
  dplyr::select(read_id, CGI_id, meth_start, meth_end) %>%
  distinct()

reads_ps_updated <- reads_for_ps %>%
  left_join(
    ps_blocks_updated,
    by = join_by(overlaps(meth_start, meth_end, ps_start, ps_end))
  )

cat("Reads reassigned:", nrow(reads_ps_updated), "\n")

################################################################################
######### Update hp_df with new PS tags
################################################################################

cat("Updating HP reads with new PS tags...\n")

read_ps_lookup <- reads_ps_updated %>%
  filter(!is.na(PS)) %>%
  distinct(read_id, PS) %>%
  rename(PS_updated = PS)

hp_df_updated <- hp_df %>%
  left_join(read_ps_lookup, by = "read_id") %>%
  mutate(
    PS = case_when(
      !is.na(PS_updated) ~ PS_updated,
      TRUE ~ PS
    )
  ) %>%
  dplyr::select(read_id, HP, PS)

cat("HP reads updated:", nrow(hp_df_updated), "\n")

################################################################################
######### Run calculate_skew_by_block with updated HP reads
################################################################################

# define function
calculate_skew_by_block <- function(clustered_reads, haplotyped_reads){
  clustered_reads <- clustered_reads %>% filter(assigned_X %in% c("Xa","Xi"))
  clustered_reads <- clustered_reads %>% distinct(read_name, .keep_all = TRUE)
  df2 <- left_join(clustered_reads, haplotyped_reads)
  df2 <- df2 %>% filter(!is.na(HP))
  counts_by_block <- df2 %>% group_by(PS, assigned_X, HP) %>% summarise(counts = n(), .groups = "drop")
  skew_by_block <- counts_by_block %>%
    unite(combi, assigned_X, HP) %>%
    mutate(combi = recode(combi,
                          "Xa_1" = "H1_Xa", "Xa_2" = "H2_Xa",
                          "Xi_1" = "H1_Xi", "Xi_2" = "H2_Xi")) %>%
    pivot_wider(id_cols = PS, names_from = combi, values_from = counts, values_fill = 0)
  if(!"H1_Xa" %in% colnames(skew_by_block)) skew_by_block$H1_Xa <- 0
  if(!"H1_Xi" %in% colnames(skew_by_block)) skew_by_block$H1_Xi <- 0
  if(!"H2_Xa" %in% colnames(skew_by_block)) skew_by_block$H2_Xa <- 0
  if(!"H2_Xi" %in% colnames(skew_by_block)) skew_by_block$H2_Xi <- 0
  skew_by_block <- skew_by_block %>%
    mutate(H1_Xa_skew = (H1_Xa + H2_Xi) / (H1_Xa + H1_Xi + H2_Xa + H2_Xi))
  return(skew_by_block)
}

##########

cat("Computing skew per updated PS block...\n")

clustered_reads <- fread(meth_file) %>%
  as_tibble() %>%
  rename(
    cluster    = cluster_id,
    meth_start = start,
    meth_end   = end
  )

skew_by_block_updated <- calculate_skew_by_block(
  clustered_reads,
  hp_df_updated %>% rename(read_name = read_id)
)

# add positional info for downstream use
skew_by_block_updated <- skew_by_block_updated %>%
  left_join(
    ps_blocks_updated %>% dplyr::select(PS, ps_start, ps_end, split, boundary_inside_cgi),
    by = "PS"
  )

cat("PS blocks with skew:", nrow(skew_by_block_updated), "\n")

# -------------------------
# Skew per PS block bar plot
# -------------------------
ps_plot_df <- skew_by_block_updated %>%
  mutate(
    ps_mid_mb  = (ps_start + ps_end) / 2e6,
    ps_width   = (ps_end - ps_start) / 1e6,
    block_type = case_when(
      split == "left"  ~ "split",
      split == "right" ~ "split",
      TRUE             ~ "unsplit"
    )
  ) %>%
  filter(!is.na(H1_Xa_skew))

p2 <- ggplot(ps_plot_df) +
  geom_col(
    aes(x = ps_mid_mb, y = H1_Xa_skew, width = ps_width, fill = block_type),
    alpha = 0.8
  ) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.4) +
  scale_fill_manual(
    values = c("unsplit" = "steelblue", "split" = "tomato"),
    name   = "Block type"
  ) +
  scale_x_continuous(
    name   = "Genomic Position (Mb)",
    limits = c(0, max_x),
    breaks = pretty(c(0, max_x), n = 10),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    name   = "Skew (H1_Xa_skew)",
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  labs(title = "Skew per PS block across chrX") +
  theme_minimal()

ggsave("skew_by_block_plot.png", plot = p2, width = 20, height = 6, dpi = 300)
print(p2)

################################################################################
######### Write outputs
################################################################################

write_tsv(
  skew_cpg,
  file.path(paste0(sample_id, "_skewPerCGI.tsv.gz"))
)

write_tsv(
  skew_by_block_updated %>%
    dplyr::select(PS, ps_start, ps_end, split, boundary_inside_cgi,
                  H1_Xa, H2_Xa, H1_Xi, H2_Xi, H1_Xa_skew),
  file.path(paste0(sample_id, "_skewPerPS.tsv.gz"))
)

cat("Outputs written.\n")