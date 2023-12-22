


splitting <- function(list.cleaned)   {

  # split the chromosome

  met.cleaned <- list.cleaned[[2]]
  unmet.cleaned <- list.cleaned[[3]]
  cov.cleaned <- list.cleaned[[1]]
  beta.cleaned <- list.cleaned[[4]]
  m.cleaned <- list.cleaned[[5]]

  chr <- rownames(met.cleaned)
  chr <- gsub("-.*", "", chr)
  n <- ncol(met.cleaned)

  met.cleaned <- as.data.frame(cbind(met.cleaned, chr))
  colnames(met.cleaned)[n+1] <- "chr"

  unmet.cleaned <- as.data.frame(cbind(unmet.cleaned, chr))
  colnames(unmet.cleaned)[n+1] <- "chr"

  cov.cleaned <- as.data.frame(cbind(cov.cleaned, chr))
  colnames(cov.cleaned)[n+1] <- "chr"

  beta.cleaned <- as.data.frame(cbind(beta.cleaned, chr))
  colnames(beta.cleaned)[n+1] <- "chr"

  m.cleaned <- as.data.frame(cbind(m.cleaned, chr))
  colnames(m.cleaned)[n+1] <- "chr"

  label.chr <- unique(chr)


  split_met <- list()
  split_unmet <- list()
  split_cov <- list()
  split_beta <- list()
  split_m <- list()

  # Met
  for (i in label.chr) {

    split_data <- met.cleaned[met.cleaned$chr== i, ]
    #split_data <- data.frame(sapply(split_data, as.numeric))
    is_matrix <- is.matrix(split_data)

    if (!is_matrix) {
      split_data <- as.matrix(split_data)

    }

    split_data2 <- apply(split_data, c(1, 2), as.numeric)
    split_met[[i]] <- split_data2
    #split_met2 <- lapply(split_met, as.matrix.data.frame)

    names(split_met)
  }

  # remove chr col
  for (i in label.chr) {

    split_met[[i]] <- split_met[[i]][,-(n+1)]


  }


  # Unmet
  for (i in label.chr) {

    split_data <- unmet.cleaned[unmet.cleaned$chr== i, ]
    #split_data <- data.frame(sapply(split_data, as.numeric))
    is_matrix <- is.matrix(split_data)

    if (!is_matrix) {
      split_data <- as.matrix(split_data)

    }

    split_data2 <- apply(split_data, c(1, 2), as.numeric)
    split_unmet[[i]] <- split_data2
    #split_met2 <- lapply(split_met, as.matrix.data.frame)

    names(split_unmet)
  }

  # remove chr col
  for (i in label.chr) {

    split_unmet[[i]] <- split_unmet[[i]][,-(n+1)]


  }


  # Beta
  for (i in label.chr) {

    split_data <- beta.cleaned[beta.cleaned$chr== i, ]
    #split_data <- data.frame(sapply(split_data, as.numeric))
    is_matrix <- is.matrix(split_data)

    if (!is_matrix) {
      split_data <- as.matrix(split_data)

    }

    split_data2 <- apply(split_data, c(1, 2), as.numeric)
    split_beta[[i]] <- split_data2
    #split_met2 <- lapply(split_met, as.matrix.data.frame)

    names(split_beta)
  }

  # remove chr col
  for (i in label.chr) {

    split_beta[[i]] <- split_beta[[i]][,-(n+1)]


  }



  # M
  for (i in label.chr) {

    split_data <- m.cleaned[m.cleaned$chr== i, ]
    #split_data <- data.frame(sapply(split_data, as.numeric))
    is_matrix <- is.matrix(split_data)

    if (!is_matrix) {
      split_data <- as.matrix(split_data)

    }

    split_data2 <- apply(split_data, c(1, 2), as.numeric)
    split_m[[i]] <- split_data2
    #split_met2 <- lapply(split_met, as.matrix.data.frame)

    names(split_m)
  }

  # remove chr col
  for (i in label.chr) {

    split_m[[i]] <- split_m[[i]][,-(n+1)]


  }



  # Coverage
  for (i in label.chr) {

    split_data <- cov.cleaned[cov.cleaned$chr== i, ]
    #split_data <- data.frame(sapply(split_data, as.numeric))
    is_matrix <- is.matrix(split_data)

    if (!is_matrix) {
      split_data <- as.matrix(split_data)

    }

    split_data2 <- apply(split_data, c(1, 2), as.numeric)
    split_cov[[i]] <- split_data2
    #split_met2 <- lapply(split_met, as.matrix.data.frame)

    names(split_cov)
  }

  # remove chr col
  for (i in label.chr) {

    split_cov[[i]] <- split_cov[[i]][,-(n+1)]


  }


  # From the 5 lists I will unify the 5 lists per each chromosome

  final <- list()

  matrices <- list(Coverage = split_cov, Methylated = split_met, Unmethylated = split_unmet,
                   Beta = split_beta, M = split_m)


  for (i in label.chr){

    a <- as.character(i)
    aa <- list(Coverage_matrix3 = matrices[["Coverage"]][[a]], Meth_matrix3 = matrices[["Methylated"]][[a]],
               Unmeth_matrix3 = matrices[["Unmethylated"]][[a]],
               Beta_matrix3 = matrices[["Beta"]][[a]], M_matrix3 = matrices[["M"]][[a]])

    final[[i]]<- aa

  }


  return(final)

}
