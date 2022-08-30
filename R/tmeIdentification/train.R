#' @param expr Rows should be cells and the last column should be "label".
#' @param k Number of genes to select.
#' @return  'out'  selected gene list(class "character)
Entropy_test <- function(expr, k = 1000) {
  labels <- factor(expr$label)
  labels_set <- levels(labels)
  label_total <- matrix(0,ncol(expr)-1,1)
  for(i in labels_set){
    label_TPM <- expr[expr$label==i,][,-ncol(expr)]
    label_mean <- colMeans(label_TPM)
    temp <- t(matrix(label_mean,1))
    colnames(temp) <- i
    label_total <- cbind(label_total, temp)
  }
  label_total <- label_total[,-1]
  #E-test
  log_E <- log2(rowMeans(label_total + 1))
  E_log <- rowMeans(log2(label_total + 1))
  t_scores <- log_E - E_log
  out <- cbind(label_total, t_scores)
  rownames(out) <- colnames(expr)[-ncol(expr)]
  out <- out[order(-out[, "t_scores"]),]
  select <- out[1:k,]
  out <- rownames(select)
  return(out)
}

#' @param expr Rows should be cells and the last column should be "label".
#' @param k Number of genes to select.
#' @return  'out'  selected gene list(class "character)
HRG <- function(counts, k = 1000) {
  object.seurat <- CreateSeuratObject(counts)
  all.genes = rownames(object.seurat)
  object.seurat <- ScaleData(object.seurat, features=all.genes, 
                             verbose=FALSE)
  object.seurat <- RunPCA(object.seurat, features=all.genes, 
                          verbose=FALSE)
  object.seurat <- FindRegionalGenes(object.seurat,dims=1:10,
                                     nfeatures=k, overlap_stop=0.99)
  out <- head(RegionalGenes(object.seurat), k)
  # Replace dashes in Seurat's feature names with underscores
  out <- gsub('[-]', '_', out)
  return(out)
}


# Feature Selection
#' @name SelectGenes
#' @return geneset A user-defined set of genes for prediction.
#' @export
SelectGenes <- function(expr, method="Entropy", k = 2000){
  if(method == "Entropy"){
    geneset <- Entropy_test(expr, k)
  }
  if(method == "HRG"){
    geneset <- HRG(data.frame(t(expr[, 1:dim(expr)[2]-1])), k)
  }
  return(geneset)
}


# Sample Selection
#' @name SelectCells
#' @param object Seurat object of training set.
#' @return representative.index A user-defined set of cells for training.
#' @export
SelectCells <- function(object, marker_file_path, 
                        cutoff = .80, database = org.Hs.eg.db){
  marker_check <- check_markers(object, marker_file_path, 
                                db=database, 
                                cds_gene_id_type = "SYMBOL", 
                                marker_file_gene_id_type = "SYMBOL")
  # print(plot_markers(marker_check))
  
  # select representative samples for training
  set.seed(unclass(Sys.time()))
  training_sample <- select_fine_samples(cds = object,
                                         marker_file = marker_file_path,
                                         db=database,
                                         cds_gene_id_type = "SYMBOL",
                                         cutoff = cutoff,
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL",
                                         return_initial_assign = TRUE)
  representative.index <- which(training_sample$Assignment != "None")
  return(list(index = representative.index,
              markers = marker_check))
}


#' Train SciBet model and return cellType-genes prob-matrix
#' @name Train
#' @usage Train(expr, geneset, k=1000)
#' @param expr The reference expression dataframeThe expression dataframe, with rows being cells, and columns being genes. The last column should be "label".
#' @param geneset A user-defined set of genes for prediction.
#' @return reference model(cellType-genes prob-matrix)
#' @export

Train <- function(expr, geneset=NULL){
  if(is.null(geneset)){
    print("All scRNA-seq features...")
    path <- "D:/FinalProject/scCancer_MicroEnv/single_cell_features.tsv"
    gene.list <- read_tsv(path)
    gene.list <- gene.list$`MIR1302-10`
    index_remain <- which(colnames(expr) %in% intersect(colnames(expr), gene.list))
    geneset <- colnames(expr)[index_remain]
  }
  labels <- factor(expr$label)
  label.level <- levels(labels)
  expr_select <- expr[,geneset]
  label_total <- matrix(0, length(geneset), 1)
  lambda_all <- matrix(0, length(geneset), 1)
  
  # Calculate zero-ratio lambda for every possible cell type.
  dropout.ratio <- c()
  for(j in 1:dim(expr_select)[2]){
    dropout.ratio <- c(dropout.ratio, 
                       length(which(expr_select[,j] == 0)) / dim(expr_select)[1])
  }
  
  # all possible labels
  for(i in label.level){
    label_TPM <- expr_select[labels == i,]
    # Calculate average expression.
    label_mean <- colSums(log2(label_TPM + 1)) # log normalization
    temp <- t(matrix(label_mean, 1))
    colnames(temp) <- i
    label_total <- cbind(label_total, temp)
    temp <- t(matrix(dropout.ratio, 1))
    colnames(temp) <- i
    lambda_all <- cbind(lambda_all, temp)
  }
  label_total <- label_total[,-1]
  rownames(label_total) <- geneset
  lambda_all <- lambda_all[,-1]
  rownames(lambda_all) <- geneset
  # calculate probability and log transform
  label_t <- t(label_total)
  # log2[(X + 1) / sum(X + 1)]
  prob <- log2(label_t + 1) - log2(rowSums(label_t) + length(geneset))
  prob <- t(prob)
  return(cbind(prob, lambda_all))
}
