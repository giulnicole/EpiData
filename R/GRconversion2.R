#' @title GRconversion2
#' @description
#' \code{\link{GRconversion2}} internal function which helps in adjusting methylated matrix and unmethylated matrix when filtering the CpGs according to the already filtered coverage
#'
#'
#' @param input.list with the three experimental matrices from the SummarizedExperiment object: coverage filtered matrix, methylated counts' matrix, unmethylated counts'matrix.
#'
#'
#' @name GRconversion2
#'
#' @return
#' a SummarizedExperiment object with three matrices in the assays:
#'  \item{coverage matrix}{methylated counts matrix}
#'  \item{methylated matrix}{methylated counts matrix}
#'  \item{unmethylated matrix}{unmethylated counts matrix}
#'
#'
#' @examples
#' \dontrun{
#' data("matrices")
#' coverage<- assay(dati,1)
#' met <- assay(dati,2)
#' unmet <- assay(dati,3)
#' input.list= list(Cov=coverage, Met=met, Unmet=unmet)
#' gr.object <- GRconversion2(input.list)
#'  }
#'
#'@export
GRconversion2<- function(input.list){


  nrows <- dim(input.list[[1]])[1]
  ncols <- dim(input.list[[1]])[2]
  counts <- input.list[[1]]
  #counts2 <- input.list[[2]]

  names_cpg<-rownames(counts)

  chr.whole<- gsub("-.*", "", rownames(counts))
  chr<- unique(chr.whole)
  #freq <- length(chr)
  positions <-  gsub(".*-", "", rownames(counts))


  list.new<- lapply(input.list, as.data.frame)


  for (i in 1:length(list.new)){


    rownames(list.new[[i]]) <- names_cpg


  }


  rowRanges <- GRanges(chr.whole,
                       IRanges(positions),
                       feature_id=names_cpg)

  colData <- DataFrame(Matrix=rep("coverage", ncol(counts)))


  names(list.new)<- names(input.list)
  object <- SummarizedExperiment(assays=list.new,
                                 rowRanges=rowRanges, colData=colData)



  return(object)



}

