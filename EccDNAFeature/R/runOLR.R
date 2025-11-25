#' @title Get OLR
#' @description
#' A function, is used to get OLR.
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
#' @return OLR
#' @export
#'
#'
runOLR <- function(bamfile){

  param <- ScanBamParam(flag = scanBamFlag(isPaired = NA, isProperPair = NA, isUnmappedQuery = NA,
                                           hasUnmappedMate = NA, isMinusStrand = NA, isMateMinusStrand = NA,
                                           isFirstMateRead = NA, isSecondMateRead = NA,
                                           isSecondaryAlignment = NA, isNotPassingQualityControls = NA,
                                           isDuplicate = NA, isSupplementaryAlignment = NA),
                        what = 'seq',tag = 'SA')
  differ <- readGAlignments( bamfile, param = param,use.names = T)
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

  non_grch38_chr <- setdiff(seqlevels(hg.38), c(paste0("chr", 1:22), "chrX", "chrY"))




  Mlong <- unlist(lapply(str_extract_all(df2$p.cigar, "\\d+(?=M)"),function(x) sum(as.numeric(x))))
  Slong <- unlist(lapply(str_extract_all(df2$cigar2, "\\d+(?=[HS])"),function(x) sum(as.numeric(x))))
  range <- Mlong-Slong
  range <- which(range > -1)
  df2 <- df2[range,]
  df2$len<- unlist(lapply(str_extract_all(df2$p.cigar, "\\d+(?=M)"),function(x) sum(as.numeric(x))))-unlist(lapply(str_extract_all(df2$cigar2, "\\d+(?=[HS])"),function(x) sum(as.numeric(x))))
  frequency <- table(df2$len)
  frequency_df <- data.frame(len = names(frequency), frequency = as.vector(prop.table(frequency)))
  return(frequency_df)
}
