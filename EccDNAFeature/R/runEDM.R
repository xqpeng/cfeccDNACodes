runEDM <- function(bamfile){
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
   
  edm <- character(nrow(df1)) 
  #positive strand 
  positive_idx <- as.character(df1$p.strand) == '+'
  if(any(positive_idx)) {
    positive_df <- df1[positive_idx, ]
    read_edm <- substr(positive_df$sequence, 1, 4)
    edm[positive_idx] <- read_edm
  }
   #negative strand 
  negative_idx <- as.character(df1$p.strand) == '-'
  if(any(negative_idx)) {
    negative_df <- df1[negative_idx, ]
    read_edm <- substr(negative_df$sequence, 
                      nchar(negative_df$sequence) - 3, 
                      nchar(negative_df$sequence))
    read_edm <- stringi::stri_reverse(read_edm) 
    edm[negative_idx] <- chartr("ATCG", "TAGC", read_edm)
  }
  
  pa.df <- data.frame(
    Qname = df1$id,
    Chr = df1$p.chr,
    Motif = edm
  )

  pa.df <- pa.df[!grepl("N", pa.df$Motif), ]

  sedm<- as.data.frame(prop.table(table(pa.df$Motif)))
  names(sedm) <- c('motif','frequency')
  number_motif<-4
  base_pairs <- c("A", "T", "C", "G")
  all_combinations <- expand.grid(rep(list(base_pairs), number_motif))
  all_combinations$motif <- apply(all_combinations, 1, paste, collapse = "")
  empty.df <- as.data.frame(all_combinations[,-c(1:number_motif)])
  names(empty.df)<-'motif'
  EDM<-merge(empty.df,sedm,by='motif',all.x=TRUE)
  EDM[is.na(EDM)] <- 0

  return(EDM)
}