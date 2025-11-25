feature_CNV <- function(bamfile){

  # load hg38 genome object
  genome <- BSgenome.Hsapiens.UCSC.hg38

  # get chromosome lengths 
  chrom_lengths <- seqlengths(genome)[paste0("chr", c(1:22))]

  # load oncogene list ----
  data("Oncogene", package = "EccDNAFeature")

reads_gr <- granges(readGAlignments(bamfile))
  seqlevelsStyle(reads_gr) <- "UCSC" 
  
  results_df <- data.frame(
    gene_name = character(),
    log2_corrected_ratio = numeric(),
    stringsAsFactors = FALSE
  )
  
  
  
  cat("\nprocessing ", nrow(Oncogene), " oncogenes...\n")
  
  # process each oncogene ----
  for (i in 1:nrow(Oncogene)) {
    
    current_gene_info <- Oncogene[i, ]
    gene_name <- current_gene_info$geneName
    
    gene_gr <- GRanges(
      seqnames = current_gene_info$chr,
      ranges = IRanges(start = current_gene_info$start, end = current_gene_info$end)
    )
    seqlevelsStyle(gene_gr) <- "UCSC"
    gene_length <- width(gene_gr)
    
    if (is.na(gene_length) || gene_length <= 0) {
      warning(paste("Oncogene '", gene_name, "' has invalid length, skip.", sep = ""))
      next
    }
    
    # C2 (raw): reads counts in the oncogene region  
    C2_raw <- sum(countOverlaps(gene_gr, reads_gr))
    
    # split the genome into windows with length=the length of the oncogene 
    windows_gr <- tileGenome(chrom_lengths, tilewidth = gene_length, cut.last.tile.in.chrom = TRUE)
    
    # C1 (raw): reads counts in windows 
    C1_raw_counts <- countOverlaps(windows_gr, reads_gr)
    C1_raw_mean <- mean(C1_raw_counts)
    
    # compute the CG contents of windows
    window_seqs <- getSeq(genome, windows_gr)
    window_freqs <- alphabetFrequency(window_seqs, baseOnly = TRUE, as.prob = TRUE)
    window_gc_content <- window_freqs[, "C"] + window_freqs[, "G"]
    
    C2_corrected <- NA # 初始化
    
    # GC model for recalculating read counts under GC correction
    model_data <- data.frame(counts = C1_raw_counts, gc = window_gc_content)
    model_data <- model_data[model_data$counts > 0, ]
    
    if (nrow(model_data) > 100) {
      loess_fit <- loess(counts ~ gc, data = model_data, span = 0.3)
      
      gene_seq <- getSeq(genome, gene_gr)
      gene_freqs <- alphabetFrequency(gene_seq, baseOnly = TRUE, as.prob = TRUE)
      gene_gc_content <- gene_freqs[, "C"] + gene_freqs[, "G"]
      
      predicted_gene_count <- predict(loess_fit, newdata = data.frame(gc = gene_gc_content))
      
      correction_factor <- C1_raw_mean / predicted_gene_count
      
      if (!is.na(correction_factor) && !is.infinite(correction_factor)) {
        C2_corrected <- C2_raw * correction_factor
      }
    } else {
      warning("Oncogene '", gene_name, "' (in BAM '", basename(bamfile), "') does not have sufficient data, and it cann't build reliable GC model. Skip GC correction!")
      gene_gc_content <- NA 
    }
    
    # calculate log2 ratio
    log2_raw_ratio <- ifelse(C1_raw_mean == 0, 0, log2(C2_raw / C1_raw_mean))
    log2_corrected_ratio <- ifelse(is.na(C2_corrected) || C1_raw_mean == 0, NA, log2(C2_corrected / C1_raw_mean))
    
    # store result
    results_df[i, ] <- list(
      gene_name,
      log2_corrected_ratio
    )
    
    
  } 
  
  cat("---------------------------------------------------\n")
  cat("CNV_onco profile of  '", basename(bamfile), "'  is finished！\n\n")
  
   return(results_df)
  }