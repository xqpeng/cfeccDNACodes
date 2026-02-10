# =============================================================================
# Triplex Position Distribution Analysis (0-based Physical Distance Mapping)
# =============================================================================

library(tidyverse)
library(readr)
library(ggplot2)
library(Biostrings) 
library(parallel) 

# --- 1. Path Configuration ---
POS_RES_DIR <- "D:\\Triplex\\triplex(R)\\batch_results_Final_DocCorrect2"
POS_FASTA_DIR <- "D:\\Triplex\\Target_sequences\\output_fasta_5bp_filtered"
GROUP_FILE <- "D:\\Triplex\\Target_sequences\\5093_README.txt"

# Version 3 output directory for zero-based distance analysis
OUTPUT_DIR <- file.path(POS_RES_DIR, "Combined_Distance_Analysis_Strict_18Threads_V3_ZeroBased")
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR)

# Directory for caching intermediate batch processing results
TEMP_DIR <- file.path(OUTPUT_DIR, "temp_chunks_corrected")
if (!dir.exists(TEMP_DIR)) dir.create(TEMP_DIR)

# --- 2. Aesthetic and Font Settings ---
my_colors <- c(
  "Control" = "#525252",
  "HCC"     = "#D93636",
  "HNSCC"   = "#3478C4",
  "NPC"     = "#54A356",
  "CRC"     = "#9655AA",
  "LC"      = "#D6A629"
)

if (.Platform$OS.type == "windows") {
  windowsFonts(Arial = windowsFont("Arial"))
  font_family <- "Arial"
} else {
  font_family <- "sans"
}

# --- 3. Plot Theme Definition (Optimized for Publication) ---
custom_theme <- theme_bw() + theme(
  text = element_text(family = font_family),
  plot.title = element_blank(), 
  
  # Axis settings
  axis.title = element_text(size = 24, family = font_family), 
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 22, color = "black", family = font_family),
  axis.text.y = element_text(size = 22, color = "black", family = font_family),
  
  # Legend appearance and positioning
  legend.position = c(0.99, 0.99), 
  legend.justification = c(1, 1),
  legend.background = element_rect(fill = "white", color = "black", size = 0.5),
  legend.key = element_rect(fill = "transparent", color = NA),
  legend.title = element_blank(),
  legend.text = element_text(size = 14, family = font_family),
  
  panel.grid.major = element_line(size = 0.2), 
  panel.grid.minor = element_blank(),
  
  plot.margin = margin(t = 10, r = 60, b = 10, l = 10, unit = "pt")
)

# Legend aesthetics override for density plots
legend_fix <- guides(
  fill = "none", 
  color = guide_legend(override.aes = list(size = 1.5, linetype = 1, alpha = 1))
)

# =============================================================================
# Core Calculation: 0-based Physical Distance Algorithm
# =============================================================================
process_single_folder <- function(folder_name, res_dir, fasta_dir) {
  library(Biostrings); library(dplyr)
  sample_name_core <- gsub("^analysis_", "", folder_name)
  
  fasta_filename <- paste0(sample_name_core, ".fasta")
  fasta_path <- file.path(fasta_dir, fasta_filename)
  csv_path <- file.path(res_dir, folder_name, "best_triplex_sequences_detailed.csv")
   
  result_df <- NULL
  if (file.exists(fasta_path) && file.exists(csv_path)) {
    fasta_seqs <- tryCatch(readDNAStringSet(fasta_path), error = function(e) NULL)
    csv_data <- tryCatch(read.csv(csv_path, stringsAsFactors = FALSE), error = function(e) NULL)
    
    if (!is.null(fasta_seqs) && !is.null(csv_data) && nrow(csv_data) > 0) {
      ref_seqs_str <- as.character(fasta_seqs)
      names(ref_seqs_str) <- sapply(strsplit(names(ref_seqs_str), "\\|"), `[`, 1)
      
      n_rows <- nrow(csv_data)
      min_dists <- numeric(n_rows) 
      valid_rows <- logical(n_rows)
      
      for (k in 1:n_rows) {
        raw_id <- csv_data$SequenceID[k]
        short_csv_id <- strsplit(raw_id, "\\|")[[1]][1]
        
        if (short_csv_id %in% names(ref_seqs_str)) {
          original_seq <- ref_seqs_str[[short_csv_id]]
          triplex_seq <- csv_data$Sequence[k]
          
          # Locate the triplex sequence within the original read
          i <- as.integer(regexpr(triplex_seq, original_seq, fixed = TRUE))
          
          if (i < 0) {
            rc_seq <- as.character(reverseComplement(DNAString(triplex_seq)))
            i <- as.integer(regexpr(rc_seq, original_seq, fixed = TRUE))
          }
          
          if (i > 0) {
            length_seq <- nchar(original_seq)
            triplex_len <- nchar(triplex_seq)
            triplex_end_pos <- i + triplex_len - 1
            
            # Distance Calculation (0-based): Number of bases from the nearest end
            dist_to_5 <- i - 1
            dist_to_3 <- length_seq - triplex_end_pos
            
            # Select the minimum distance to either flanking end
            min_dist <- min(dist_to_5, dist_to_3)
            
            min_dists[k] <- min_dist
            valid_rows[k] <- TRUE
          }
        }
      }
      
      if (any(valid_rows)) {
        clean_sample <- gsub("^TS_", "", sample_name_core)
        result_df <- data.frame(
          Sample = clean_sample,
          MinDist = min_dists[valid_rows],
          stringsAsFactors = FALSE
        )
      }
    }
  }
  return(result_df)
}

# =============================================================================
# Parallel Processing Scheduler (18-Thread Support)
# =============================================================================
run_batch_processing <- function(res_dir, fasta_dir, n_cores, batch_size) {
  res_folders <- list.dirs(res_dir, full.names = FALSE, recursive = FALSE)
  res_folders <- res_folders[grep("^analysis_", res_folders)]
  
  total_folders <- length(res_folders)
  if(total_folders == 0) return(NULL)
  
  chunks <- split(res_folders, ceiling(seq_along(res_folders) / batch_size))
  cat(sprintf("\n>>> Initializing: %d folders in %d batches\n", total_folders, length(chunks)))
  
  all_results <- list()
  
  for (i in seq_along(chunks)) {
    chunk_folders <- chunks[[i]]
    cat(sprintf("   Processing Batch %d/%d (18 Threads)... ", i, length(chunks)))
    
    cl <- makeCluster(n_cores)
    clusterEvalQ(cl, {library(Biostrings); library(dplyr)})
    
    batch_res_list <- parLapply(cl, chunk_folders, process_single_folder, 
                                res_dir = res_dir, fasta_dir = fasta_dir)
    stopCluster(cl)
    
    batch_df <- do.call(rbind, batch_res_list)
    
    if (!is.null(batch_df) && nrow(batch_df) > 0) {
      temp_file <- file.path(TEMP_DIR, paste0("batch_", i, ".rds"))
      saveRDS(batch_df, temp_file)
      cat(sprintf("OK! (Saved %d rows)\n", nrow(batch_df)))
      all_results[[i]] <- batch_df
    } else {
      cat("No valid data.\n")
    }
    gc() 
  }
  return(do.call(rbind, all_results))
}

# =============================================================================
# Main Execution Pipeline
# =============================================================================

n_cores <- 18 
BATCH_SIZE <- 20 

cat("\n>>> [1/3] Loading Group Metadata...\n")
readme_lines <- readLines(GROUP_FILE)
header_idx <- which(grepl("Sample", readme_lines) & grepl("Group", readme_lines))[1]
group_data <- read.table(text = readme_lines[header_idx:length(readme_lines)], header = TRUE, sep = "\t", quote = "", fill = TRUE)

group_data <- group_data %>% 
  filter(library == "NonBS-seq") %>% 
  select(Sample, Group) %>% 
  mutate(Sample = trimws(gsub("\\.BS$", "", Sample)), Group = trimws(Group)) %>%
  filter(Group %in% names(my_colors)) 

cat("\n>>> [2/3] Executing Parallel Back-tracking Calculation...\n")
start_t <- Sys.time()
pos_result <- run_batch_processing(POS_RES_DIR, POS_FASTA_DIR, n_cores, batch_size = BATCH_SIZE)
print(difftime(Sys.time(), start_t))

if (!is.null(pos_result)) {
  all_data <- pos_result %>% inner_join(group_data, by = "Sample")
  all_data$Group <- factor(all_data$Group, levels = c(names(my_colors)))
  
  cat("\n>>> [3/3] Generating Visualization...\n")
  
  max_x <- 80 
  
  # Density plot of triplex distribution relative to sequence ends
  p <- ggplot(all_data, aes(x = MinDist, color = Group)) +
    geom_density(size = 1, alpha = 0, key_glyph = "path") + 
    
    scale_x_continuous(limits = c(0, max_x), expand = c(0, 0), breaks = seq(0, 80, 10)) +
    scale_color_manual(values = my_colors) +
    
    custom_theme + 
    legend_fix + 
    
    labs(
      title = NULL,
      x = "Triplex distance to nearest end (bp)",
      y = "Density"
    )
  
  out_file <- file.path(OUTPUT_DIR, "04_Distance_Distribution_ZeroBased.png")
  ggsave(out_file, p, width = 23.44, height = 18.36, units = "cm", dpi = 300)
  
  cat("\n✅ Visualization generated successfully.\n")
  cat("✅ Metric: 0-based physical distance (distance=0 at sequence boundaries).\n")
  cat("✅ Output File: ", out_file, "\n")
  
} else {
  cat("❌ Error: Calculation returned empty dataset.\n")
}