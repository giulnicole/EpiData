
#' @noRd
#'
GRconversion2<- function(cleaned.list){


  nrows <- dim(cleaned.list[["Coverage_matrix"]])[1]
  ncols <- dim(cleaned.list[["Coverage_matrix"]])[2]
  counts <- cleaned.list[["Coverage_matrix"]]
  counts2 <- cleaned.list[["Met_matrix"]]

  names_cpg<-rownames(counts)

  chr<- gsub("-.*", "", rownames(counts))
  chr<- unique(chr)
  freq <- length(chr)
  positions <-  gsub(".*-", "", rownames(counts))

  rowRanges <- GRanges(rep(chr,freq),
                       IRanges(positions),
                       feature_id=names_cpg)



  colData <- DataFrame(Matrix=rep("coverage", ncol(counts)))
  object <- SummarizedExperiment(assays=cleaned.list,
                                 rowRanges=rowRanges, colData=colData)
  return(object)



}


