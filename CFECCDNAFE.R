#!/usr/bin/env Rscript
library(EccDNAFeature)

# parameters
args <- commandArgs(trailingOnly = TRUE)

# check the number of parameters
if (length(args) < 2) {
  cat("Usage: Rscript EccDNAFE.R <feature1,feature2,...> <dir1> <dir2> ...\n")
  cat("Example: Rscript EccDNAFE.R JNM,SBM,OLR /path/to/dir1 /path/to/dir2 /path/to/dir3\n")
  cat("Available features: BPM, EDM, JNM, SBM, OJM, OLR, CNV_onco, CNV_im\n")
  quit(status = 1)
}

# split feature parameters
feature_names <- unlist(strsplit(args[1], ","))
input_dirs <- args[-1]

# mapping feature function 
feature_mapping <- list(
  JNM = EccDNAFeature::runJNM,
  SBM = EccDNAFeature::runSBM,
  OJM = EccDNAFeature::runOJM,
  OLR = EccDNAFeature::runOLR
  BPM = EccDNAFeature::runBPM,
  EDM = EccDNAFeature::runEDM,
  CNV_onco = EccDNAFeature::runCNV_onco,
  CNV_im = EccDNAFeature::runCNV_immu
)

# check the input feature names
valid_features <- names(feature_mapping)
invalid_features <- setdiff(feature_names, valid_features)

if (length(invalid_features) > 0) {
  cat("Error: Invalid feature name(s):", paste(invalid_features, collapse = ", "), "\n")
  cat("Available features:", paste(valid_features, collapse = ", "), "\n")
  quit(status = 1)
}

# check the input 
missing_dirs <- input_dirs[!dir.exists(input_dirs)]
if (length(missing_dirs) > 0) {
  cat("Error: Directory not found:", paste(missing_dirs, collapse = ", "), "\n")
  quit(status = 1)
}

# select feature functions
feature_funcs <- feature_mapping[feature_names]

cat("=== EccDNA Feature Extraction ===\n")
cat("Features:", paste(feature_names, collapse = ", "), "\n")
cat("Directories:", paste(input_dirs, collapse = ", "), "\n")
cat("==================================\n\n")

# collect the bam files in each directory
for (input_dir in input_dirs) {
  # select all the bam files 
  bam_files <- list.files(
    path = input_dir,
    pattern = "\\.bam$",
    full.names = TRUE
  )
  
  # exclude possible index files of bam 
  bam_files <- bam_files[!grepl("\\.bai$", bam_files)]
  
  if (length(bam_files) == 0) {
    cat("No BAM files found in:", input_dir, "\n")
    next
  }
  
  # create output directory for each input directory 
  main_output_dir <- file.path(input_dir, "eccdna_results")
  dir.create(main_output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # create sub-directory for each feature
  for (fname in feature_names) {
    dir.create(file.path(main_output_dir, fname), showWarnings = FALSE)
  }
  
  cat("Processing directory:", input_dir, "\n")
  cat("Found", length(bam_files), "BAM files\n")
  
  # process each bam file 
  for (bam in bam_files) {
    basename_full <- basename(bam)
    prefix <- sub("\\.bam$", "", basename_full)
    
    cat("  Processing:", prefix, "\n")
    
    # create output file to store each feature profile  
    for (fname in feature_names) {
      out_subdir <- file.path(main_output_dir, fname)
      out_file <- file.path(out_subdir, paste0(prefix, "_", fname, ".txt"))
      
      result <- tryCatch({
        cat("    Calculating", fname, "... ")
        feature_funcs[[fname]](bam)
      }, error = function(e) {
        cat("ERROR:", conditionMessage(e), "\n")
        NULL
      })
      
      # output the feature profile 
      
      if (!is.null(result)) {
        if (is.data.frame(result)) {
          write.table(result, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
          cat("SUCCESS (", nrow(result), "rows )\n")
        } else if (is.vector(result) && length(result) > 0) {
          writeLines(as.character(result), out_file)
          cat("SUCCESS (", length(result), "elements )\n")
        } else {
          capture.output(print(result), file = out_file)
          cat("SUCCESS (printed output)\n")
        }
      } else {
        writeLines("# Calculation failed", out_file)
        cat("FAILED\n")
      }
    }
  }
  
  cat("Completed:", input_dir, "→", main_output_dir, "\n\n")
}

cat("=== All tasks completed! ===\n")