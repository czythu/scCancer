#' Compute expression entropy.
#' @name Entropy_R
#' @usage Entropy_R(expr, window=120, low = 2000)
#' @param expr The expression dataframe. Rows should be cells and columns should be genes.
#' @param window The window size for expression value discretization.
#' @param low The lower limit for normalizing expression entropy
#' @return A dataframe..
#' @export
Entropy_R <- function(expr, window=120, low = 2000){
  expr <- as.data.frame(expr)

  ent_res <- tibble(
    gene = colnames(expr),
    mean.expr = colMeans(expr)
  ) %>%
    dplyr::filter(mean.expr < 6000)

  expr <- expr[,ent_res$gene]
  #
  #out <- GenEntr(expr, window, n_threads)
  #
  n_cell <- nrow(expr)
  n_gene <- ncol(expr)
  out <- c()
  for(i in 1:n_gene){
    states <- ceiling(1000000.0/window)
    discretize <- vector(mode="numeric", length=states)
    li <- ceiling(expr[,i]/window)+1
    number <- table(li)
    for(j in names(number)){
      a <- number[j]
      discretize[as.numeric(j)] <- a
    }
    sumAll <- sum(discretize)
    discretizeN <- discretize[discretize!=0]
    out <- c(out,-sum(discretizeN/sumAll * log(discretizeN/sumAll)))
  }

  ent_res %>%
    dplyr::mutate(entropy = out) %>%
    dplyr::mutate(fit = 0.18*log(0.03*mean.expr + 1)) -> ent_res

  ent_res %>%
    dplyr::filter(mean.expr > low) %>%
    dplyr::mutate(k = entropy/fit) %>%    #linear normalization of expression entropy
    dplyr::pull(k) %>%
    quantile(0.75) %>%
    as.numeric() -> k

  ent_res <- ent_res %>% dplyr::mutate(norm_ent = entropy/k)

  return(ent_res)
}

#' Calculating the nulltest.
#' @name NullTest_R
#' @usage NullTest_R(ref, query, null_expr, labels,gene_num = 500)
#' @param ref The reference dataset. Rows should be cells, columns should be genes.
#' @param query The query dataset. Rows should be cells and columns should be genes.The genes are intersection of ref , null and query.
#' @param labels The reference dataset's labels of each cell.
#' @param gene_num The number of common markers of reference set used for false positive control.
#' @return A vector of null test.
#' @export
NullTest_R <- function(ref,query,null_expr,labels,gene_num){
   n_ref <- nrow(ref)
   n_query <- nrow(query)
   n_gene <- ncol(ref)
   labell <- levels(labels)
   n_label <- length(labell)
   label_total <- matrix(0,ncol(ref),1)
   for(i in labell){
     label_TPM <- ref[labels==i,]
     label_mean <- colMeans(label_TPM)
     temp <- t(matrix(label_mean,1))
     colnames(temp) <- i
     label_total <- cbind(label_total, temp)
   }
   label_total <- label_total[,-1]
   rownames(label_total)<- colnames(ref)
   e1 <- c()
   ds <- c()
   e1 =  rowSums(label_total)/n_label
   ds = log(e1+1)-log(null_expr+1)
   sum_e1 <- sum(e1)
   sum_null <- sum(null_expr)
   e1 <- log((e1+1)/(sum_e1+n_gene))-log((null_expr+1)/(sum_null+n_gene))
   dss <- sort(ds,decreasing = TRUE)
   prob <- vector(mode = "numeric",length = n_query)
   for(j in 1:n_query){
       prob[j] <- sum(e1[names(dss)[1:gene_num]]*query[j,names(dss)[1:gene_num]])
   }
   return(prob)
}

#' Calculating the confidence score C for false positive control.
#' @name conf_score_R
#' @usage conf_score(ref, query, null_expr, gene_num = 500)
#' @param ref The reference dataset. Rows should be cells, columns should be genes and last column should be "label".
#' @param query The query dataset. Rows should be cells and columns should be genes.
#' @param gene_num The number of common markers of reference set used for false positive control.
#' @return A vector of confidence scores.
#' @export
conf_score_R <- function(ref, query, null_expr, gene_num){
  labels <- factor(ref$label)
  genes <- Reduce(intersect, list(colnames(ref), colnames(query), names(null_expr)))
  ref <- log1p(as.matrix(ref[, genes])) / log(2)
  query <- log1p(as.matrix(query[, genes])) / log(2)
  a <- NullTest_R(ref, query, null_expr[genes], labels, gene_num)
  b <- NullTest_R(ref, ref, null_expr[genes], labels, gene_num)
  prob <- a/max(b)
  prob[prob > 1] <- 1
  return(prob)
}

#' Generate Bet function from a model matrix.
#' @name LoadModel_R
#' @usage LoadModel(x, genes, labels)
#' @param x A SciBet model in the format of a matrix.
#' @param genes (Optional).
#' @param labels (Optional).
#' @return A Bet function.
#' @export
LoadModel_R <- function(x, genes=NULL, labels=NULL){
  prob <- x
  if (is.null(genes))
    genes <- rownames(x)
  if (is.null(labels))
    labels <- colnames(x)
  function(expr, result="list"){
    have_genes <- intersect(genes,colnames(expr))
    expra <- log1p(as.matrix(expr[,have_genes])) / log(2)
    switch(result,
           list = MLEstimate(expra, prob[have_genes,],FALSE),
           table = {
             out <- MLEstimate(expra, prob[have_genes, ], TRUE)
             rownames(out) <- have_genes
             return(out)
           }
    )
  }
}
