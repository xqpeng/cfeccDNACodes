runBPM <- function(bamfile){
  param <- ScanBamParam(flag = scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = NA,
                                           hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
                                           isFirstMateRead = NA, isSecondMateRead = NA,
                                           isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
                                           isDuplicate = NA, isSupplementaryAlignment = NA),
                        what = "seq")
  primary <- readGAlignments(bamfile, param = param,use.names = T)


  df1 <- data.frame(
  	 p.chr = seqnames(primary),
 	 p.start = start(primary),
 	 p.end = end(primary),
 	 p.strand = strand(primary),
 	 p.cigar = cigar(primary),
  	 id = names(primary)
  )
  df1$sequence <- as.character(mcols(primary)$seq)
  df1$p.signal <- ifelse(grepl("^[0-9]+M|[0-9]+M$", df1$p.cigar), "A", "B")

  hg.38 <- BSgenome.Hsapiens.UCSC.hg38
  non_grch38_chr <- as.character(setdiff(seqlevels(hg.38), c(paste0("chr", 1:22), "chrX", "chrY")))


  df1$p.chr <- as.character(df1$p.chr)
  df1 <- df1[!(df1$p.chr %in% non_grch38_chr), ]
  df1 <- df1[(as.character(df1$p.signal) == 'A'), ]
   
  bpm <- character(nrow(df1)) 
  positive_idx <- as.character(df1$p.strand) == '+'
  if(any(positive_idx)) {
    positive_df <- df1[positive_idx, ]
    upstream_seq <- as.character(getSeq(
      hg.38, 
      GRanges(
        seqnames = positive_df$p.chr,
        ranges = IRanges(start = positive_df$p.start - 2, end = positive_df$p.start - 1)
      )
    ))
    read_bpm <- substr(positive_df$sequence, 1, 2)
    bpm[positive_idx] <- paste0(upstream_seq, read_bpm)
  }
  
  # 负链处理
  negative_idx <- as.character(df1$p.strand) == '-'
  if(any(negative_idx)) {
    negative_df <- df1[negative_idx, ]
    downstream_seq <- as.character(getSeq(
      hg.38,
      GRanges(
        seqnames = negative_df$p.chr,
        ranges = IRanges(start = negative_df$p.end + 1, end = negative_df$p.end + 2)
      )
    ))
    read_bpm <- substr(negative_df$sequence, 
                      nchar(negative_df$sequence) - 1, 
                      nchar(negative_df$sequence))
    combined_bpm <- paste0(read_bpm, downstream_seq)
    bpm[negative_idx] <- chartr("ATCG", "TAGC", combined_bpm)
  }
  
  pa.df <- data.frame(
    Qname = df1$id,
    Chr = df1$p.chr,
    Motif = bpm
  )

  pa.df <- pa.df[!grepl("N", pa.df$Motif), ]

  sbpm<- as.data.frame(prop.table(table(pa.df$Motif)))
  names(sbpm) <- c('motif','frequency')
  number_motif<-4
  base_pairs <- c("A", "T", "C", "G")
  all_combinations <- expand.grid(rep(list(base_pairs), number_motif))
  all_combinations$motif <- apply(all_combinations, 1, paste, collapse = "")
  empty.df <- as.data.frame(all_combinations[,-c(1:number_motif)])
  names(empty.df)<-'motif'
  BPM<-merge(empty.df,sbpm,by='motif',all.x=TRUE)
  BPM[is.na(BPM)] <- 0

  return(BPM)
}