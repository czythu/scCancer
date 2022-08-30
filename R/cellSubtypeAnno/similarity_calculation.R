#' @export
Jaccard <- function(cell.sets, p = 1){
  similarity.mar <- matrix(nrow = length(cell.sets),
                           ncol = length(cell.sets))
  rownames(similarity.mar) <- names(cell.sets)
  colnames(similarity.mar) <- names(cell.sets)
  
  for(i in 1:length(cell.sets)){
    for(j in i:length(cell.sets)){
      m <- intersect(cell.sets[[i]], cell.sets[[j]])
      n <- union(cell.sets[[i]], cell.sets[[j]])
      similarity.mar[i, j] <- length(m)^p / length(n)
      similarity.mar[j, i] <- similarity.mar[i, j]
    }
  }
  return(similarity.mar)
}

#' @export
Spearman <- function(mean.expr){
  similarity.mar <- matrix(nrow = length(mean.expr),
                           ncol = length(mean.expr))
  rownames(similarity.mar) <- names(mean.expr)
  colnames(similarity.mar) <- names(mean.expr)
  for(i in 1:length(mean.expr)){
    for(j in i:length(mean.expr)){
      common.gene <- intersect(names(mean.expr[[i]]), names(mean.expr[[j]]))
      similarity.mar[i, j] <- cor(unname(mean.expr[[i]][common.gene]),
                                  unname(mean.expr[[j]][common.gene]),
                                  method = "spearman")
      similarity.mar[j, i] <- similarity.mar[i, j]
    }
  }
  return(similarity.mar)
}

#' @export
Intergration <- function(all.matrix){
  similarity.mean <- all.matrix[[1]]
  similarity.var <- all.matrix[[1]]
  for(i in 1:nrow(similarity.mean)){
    for(j in 1:ncol(similarity.mean)){
      values <- lapply(all.matrix, function(matrix){
        return(matrix[i, j])
      })
      similarity.mean[i, j] <- mean(unlist(values))
      similarity.var[i, j] <- var(unlist(values))
    }
  }
  return(list(mean = similarity.mean, var = similarity.var))
}
