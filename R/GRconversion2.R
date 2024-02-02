








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

