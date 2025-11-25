#' @title Get OJM
#' @description
#' A function, is used to get OJM.
#'
#' @param bamfile bamfile which you want to read
#' @importFrom Rsamtools ScanBamParam
#' @importFrom Rsamtools scanBamFlag
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicRanges mcols
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges strand
#' @importFrom GenomicAlignments cigar
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom BSgenome getSeq
#' @importFrom tidyr separate
#' @importFrom stringr str_detect
#' @importFrom stringr str_extract_all
#' @importFrom mgsub mgsub
#' @importFrom dplyr "%>%"
#' @return OJM
#' @export
#'
#'
runOJM <- function(bamfile){

  param <- ScanBamParam(flag = scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = NA,
                                           hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
                                           isFirstMateRead = NA, isSecondMateRead = NA,
                                           isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
                                           isDuplicate = NA, isSupplementaryAlignment = NA),
                        what = 'seq',tag = 'SA')#set the parameters for reading bam file 
  differ <- readGAlignments( bamfile, param = param,use.names = T) #read the bam file
  sa <- as.data.frame(differ@elementMetadata$SA)
  
  names(sa) <- 'SA'
  sa <- separate(sa,"SA",into = c('chr','pos','strand','cigar2','mapq','nm'),sep = ",")

  a<-mcols(differ)
  a.seqnames <- seqnames(differ)
  a.start<-start(differ)
  a.end<-end(differ)
  a.strand<-strand(differ)
  a.cigar <- cigar(differ)
  df1<-(a)
  a.id<-rownames(df1)
  df1$p.chr <- a.seqnames
  df1$p.start<-a.start
  df1$p.end<-a.end
  df1$p.strand<-a.strand
  df1$p.cigar <- a.cigar
  df1$p.signal <- NA
  df1$p.signal <- ifelse(grepl("^[0-9]+M.*[0-9]+[SH]$|^[0-9][SH][0-9]+M", df1$p.cigar), "A","B")
  df1$id<-a.id
  df1$cigar2 <- sa$cigar2
  df1$signal2 <- ifelse(grepl("^[0-9]+M.*[0-9]+[SH]$|^[0-9][SH][0-9]+M", df1$cigar2), "A","B")
  df2<-df1[df1[, "p.cigar"] %>% str_detect("^[0-9]+[SHM][0-9]+[MSH]$"), ]
  df2<-df2[df2[, "cigar2"] %>% str_detect("^[0-9]+[SHM][0-9]+[MSH]$"), ]

  hg.38 <- BSgenome.Hsapiens.UCSC.hg38
  # fetch the chromosome names of autosames and chromosome X and Y.  
  non_grch38_chr <- setdiff(seqlevels(hg.38), c(paste0("chr", 1:22), "chrX", "chrY"))
 


  df3 <- df2


  Mlong <- unlist(lapply(str_extract_all(df3$p.cigar, "\\d+(?=M)"),function(x) sum(as.numeric(x))))
  Slong <- unlist(lapply(str_extract_all(df3$cigar2, "\\d+(?=[HS])"),function(x) sum(as.numeric(x))))
  range <- Mlong-Slong
  range <- which(range > 3)#threshold of overlap base 
  df3 <- df3[range,]
  df3$len<- unlist(lapply(str_extract_all(df3$p.cigar, "\\d+(?=M)"),function(x) sum(as.numeric(x))))-unlist(lapply(str_extract_all(df3$cigar2, "\\d+(?=[HS])"),function(x) sum(as.numeric(x))))


  df3$overlap_middle <- NA
  # df3$overlap_middle[df3$p.signal == "A"] <-
  #   with(df3[df3$p.signal == "A", ], floor(p.start + len/2))
  # df3$overlap_middle[df3$p.signal == "B"] <-
  #   with(df3[df3$p.signal == "B", ], floor(p.end - len/2))

  # fetch logical index 
  is_A <- df3$p.signal == "A"
  is_B <- df3$p.signal == "B"

  # process according to the logical index  
  df3$overlap_middle[is_A] <-
    floor(df3$p.start[is_A] + df3$len[is_A]/2)
  df3$overlap_middle[is_B] <-
    floor(df3$p.end[is_B] - df3$len[is_B]/2)

  gr1 <- GRanges(seqnames = df3$p.chr,ranges = IRanges(start = df3$overlap_middle-2,end = df3$overlap_middle+1))
  pa.seq<-getSeq(hg.38,gr1)
  pa.seq.df <- as.character(pa.seq)
  jm <- as.data.frame(prop.table(table(pa.seq.df)))
  names(jm) <- c('motif','frequency')
  number_motif<-"4"
  base_pairs <- c("A", "T", "C", "G")
  all_combinations <- expand.grid(rep(list(base_pairs), number_motif))
  all_combinations$motif <- apply(all_combinations, 1, paste, collapse = "")
  empty.df <- as.data.frame(all_combinations[,-c(1:number_motif)])
  names(empty.df)<-'motif'
  JM<-merge(empty.df,jm,by='motif',all.x=TRUE)
  JM[is.na(JM)] <- 0
  return(JM)
}
