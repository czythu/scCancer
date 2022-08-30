#' Construct convenient Seurat pipeline
#' @name normalization
#' @usage Seurat_normalization(dataset.scibet)
#' @param dataset The expression dataframe, 
#' with rows being cells, and columns being genes.
#' @return list of Seurat object and Normalized data
normalization <- function(dataset){
  counts <- data.frame(t(dataset))
  object <- CreateSeuratObject(counts = counts)
  object <- NormalizeData(object, verbose = TRUE)
  data.normalized <- t(as.matrix(object@assays[["RNA"]]@data))
  result <- list(object, data.frame(data.normalized))
  return(result)
}

#' Construct convenient Seurat pipeline
#' @name oridata.umap
#' @usage oridata.umap(dataset, label, batch)
#' @param dataset The expression dataframe, 
#' with rows being cells, and columns being genes. 
#' The last column should be "label".
#' @param label groundtruth or output of scibet function Test.
oridata_umap <- function(object, label=NULL, batch=NULL){
  if(is.null(label)){
    label <- rep("unknown cell", time=ncol(object@assays[["RNA"]]@data))
  }
  if(is.null(batch)){
    batch <- rep("unknown batch", time=ncol(object@assays[["RNA"]]@data))
  }
  object$celltype <- label
  object$batchinfo <- batch
  object <- FindVariableFeatures(object, selection.method = "vst",verbose = TRUE)
  object <- ScaleData(object, verbose = FALSE)
  object <- RunPCA(object, npcs = 30, verbose = FALSE)
  object <- RunUMAP(object, reduction = "pca", dims = 1:30)
  p1 <- DimPlot(object, reduction = "umap", 
                group.by = "celltype", repel = TRUE)
  p2 <- DimPlot(object, reduction = "umap", 
                group.by = "batchinfo", repel = TRUE)
  result <- list(object, p1, p2)
  return(result)
}

#' Construct convenient Seurat pipeline
#' @name scibet_visualization
#' @usage Scibet_visualization(dataset, label)
#' @param dataset The expression dataframe, 
#' with rows being cells, and columns being genes. 
#' The last column should be "label".
#' @param label groundtruth or output of scibet function Test.
scibet_visualization <- function(dataset, 
                                 label=NULL,
                                 normalize=TRUE,
                                 reduction="umap",
                                 metacell=FALSE){
  counts <- data.frame(t(dataset[, 1:dim(dataset)[2]-1]))
  object <- CreateSeuratObject(counts = counts)
  if(is.null(label)){
    label <- rep("unknown cell", time=ncol(counts))
  }
  object$celltype <- label
  if(normalize){
    object <- NormalizeData(object, verbose = TRUE)
  }
  object <- FindVariableFeatures(object, selection.method = "vst",verbose = TRUE)
  object <- ScaleData(object, verbose = FALSE)
  # object <- ScaleData(object, do.scale = FALSE, do.center = TRUE, scale.max = 10)
  object <- RunPCA(object, npcs = 30, verbose = FALSE)
  if(reduction == "umap"){
    object <- RunUMAP(object, reduction = "pca", dims = 1:30)
  }
  else{
    object <- RunTSNE(object, reduction = "pca", dims = 1:30)
  }
  p <- DimPlot(object, reduction = reduction, group.by = "celltype", repel = TRUE)
  if(metacell){
    object <- FindNeighbors(object, dims = 1:30)
    object <- FindClusters(object, resolution = 100)
  }
  return(list(object = object, plot = p))
}


#' Expression patterns of informative genes across cell types.
#' @name Marker_heatmap
#' @usage Marker_heatmap(expr, gene)
#' @param expr The expression dataframe. Rows should be cells, columns should be genes and last column should be "label".
#' @param gene A vector of informative genes.
#' @return A figure.
#' @export
Marker_heatmap <- function(expr, gene){
  expr <- expr[,c(gene,'label')]
  type_expr <- expr %>%
    tidyr::nest(-label) %>%
    dplyr::rename(expr = data) %>%
    dplyr::mutate(colmeans = purrr::map(
      .x = expr,
      .f = function(.x){colMeans(.x)}))
  
  type_expr$colmeans %>%
    as.data.frame() %>%
    tibble::remove_rownames() %>%
    t() %>%
    as.data.frame() %>%
    tibble::remove_rownames() -> type_mean_expr
  
  rownames(type_mean_expr) <- type_expr$label
  colnames(type_mean_expr) <- colnames(expr)[-ncol(expr)]
  
  sub_expr <- type_mean_expr
  sub_expr <- sub_expr %>%
    as.tibble() %>%
    dplyr::mutate_all(funs((. - mean(.))/sd(.))) %>%
    t()
  colnames(sub_expr) <- type_expr$label
  get_label <- function(num){
    v <- sub_expr[num,]
    colnames(sub_expr)[which(v == max(v))]
  }
  sub_expr <- sub_expr %>%
    tibble::as.tibble() %>%
    dplyr::mutate(group = purrr::map_chr(1:length(gene), get_label))
  sub_expr <- as.data.frame(sub_expr)
  rownames(sub_expr) <- gene
  sub_expr <- sub_expr %>%
    dplyr::mutate(gene = gene) %>%
    tidyr::gather(key = 'cell_type', value = 'zscore', -group, -gene) %>%
    dplyr::arrange(group, desc(zscore))
  sub_expr %>%
    ggplot(aes(factor(gene, levels = unique(sub_expr$gene)),
               factor(cell_type, levels = sort(unique(sub_expr$cell_type), decreasing = T)))) +
    geom_point(aes(size = zscore, colour = zscore)) +
    theme(
      strip.text.x = element_blank(),
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 13),
      axis.text.y = element_text(color="black"),
      axis.text.x = element_text(color="black", angle = -90, hjust = 0),
      panel.background = element_rect(colour = "black", fill = "white"),
      panel.grid = element_line(colour = "grey", linetype = "dashed"),
      panel.grid.major = element_line(
        colour = "grey",
        linetype = "dashed",
        size = 0.2
      )
    ) +
    facet_grid(. ~ group, scales = "free", space = "free") +
    scale_colour_distiller(palette = "RdYlBu") +
    labs(
      x = '',
      y = ''
    ) -> p
  
  return(p)
}

#' Heatmap of classification result.
#' @name ConfusionMatrix
#' @usage ConfusionMatrix(label.name, label, predict)
#' @param label.name A vector of the sorted unique labels
#' @param label A vector of the original labels for each cell in the test set.
#' @param predict A vector of the predicted labels for each cell in the test set.
#' @return A heatmap for the confusion matrix of the classification result.
#' @export
ConfusionMatrix <- function(name.reference, name.prediction,
                            label, predict,
                            title='Confusion Matrix',
                            xlab='Predicted label',
                            ylab='True label',
                            normalize=F,
                            font.size=20/.pt){
  
  x <- matrix(nrow = length(name.prediction), ncol=length(name.reference))
  x[is.na(x)] <- 0
  colnames(x) <- name.reference
  rownames(x) <- name.prediction
  for (i in 1:length(predict)){
    col <- which(name.reference == label[i])    #truth
    row <- which(name.prediction == predict[i])  #predict
    x[row, col] <- x[row, col] + 1
  }
  
  if(!is.table(x)){
    x = as.table(x)
  }
  if(!is.numeric(x)){
    stop('input should be numeric, not ',mode(x),
         call. = F)
  }
  
  if(normalize){
    # x = round(prop.table(x,1), 2)
    x = round(prop.table(x,2), 2)
    mar = as.data.frame(x)
  }
  else{
    mar = as.data.frame(x)
  }
  mytheme <- theme(plot.title=element_text(
    face="bold.italic", size=12, color="darkblue"), #指定图的标题应该为粗斜体棕色14�?
    axis.title=element_text(face="bold.italic", size=10, color="darkblue"),#轴的标题为粗斜体的棕�?12
    axis.text=element_text(face="bold", size=8, color="darkblue"),#轴标签为粗体的深蓝色8�?
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.background=element_rect(fill="white",color="darkblue"),#图片区域有白色的填充和深蓝色的边框线
    panel.grid.minor.x=element_blank(), #垂直网格不输�?
    legend.position="right") #图例展示在右�?
  ggplot(mar, aes(mar[,2],mar[,1])) +
    geom_tile(aes(fill=Freq),color='black') +
    scale_fill_gradientn(colours = c('gray98','steelblue1','midnightblue'))+
    geom_label(aes(label = Freq), size=font.size) +
    labs(title = title, fill='',x=xlab,y=ylab) +
    ylim(rev(levels(mar[,2]))) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    mytheme
}

#' Heatmap for the confusion matrix of the classification with the false postive control.
#' @name Confusion_heatmap_negctrl
#' @usage Confusion_heatmap_negctrl(res, cutoff = 0.4)
#' @param res Classification result.
#' @param cutoff The cutoff of confifence score C.
#' @return A heatmap for the confusion matrix of the classification result with the false postive control.
#' @export
Confusion_heatmap_negctrl <- function(res, cutoff = 0.4){
  res %>%
    dplyr::mutate(prd = ifelse(c_score < cutoff, 'unassigned', prd)) %>%
    dplyr::count(ori, prd) %>%
    tidyr::spread(key = prd, value = n) -> cla.res
  
  cla.res[is.na(cla.res)] = 0
  cla.res[,-1] <- round(cla.res[,-1]/rowSums(cla.res[,-1]),2)
  cla.res <- cla.res %>% tidyr::gather(key = 'prd', value = 'Prob', -ori)
  label <- cla.res$ori %>% unique()
  cla.res %>%
    ggplot(aes(prd, factor(ori, levels = c(label[-3],'Neg.cell')), fill = Prob)) +
    geom_tile(colour = 'white', lwd = 0.5) +
    theme(axis.title = element_text(size = 12)) +
    theme(axis.text = element_text(size = 12)) +
    theme(legend.title = element_text(size = 0)) +
    theme(legend.text = element_text(size = 12)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank()) +
    theme(axis.text.y = element_text(color="black"),
          axis.text.x = element_text(color="black", angle = 50, hjust = 1)) +
    scale_fill_material('blue')
}

#' Correlation heatmap and hierarchical clustering
#' @export
SimilarityMap <- function(plot.title, reference, similarity.mar, similarity.var = NULL, number.digits = 2, number.cex = 1, tl.cex = 1){
  print(corrplot(similarity.mar, tl.col="black", tl.srt=45, tl.cex=tl.cex, order="hclust",
            number.digits=number.digits, number.cex=number.cex))
  title(main=paste0(plot.title, "-corr heatmap, hclust order"), sub=reference, cex=1)
  print(corrplot(similarity.mar, method="shade", shade.col=NA,
           tl.col="black", tl.srt=45, tl.cex=0.7, addCoef.col="black", cl.pos="n", order="AOE",
           col=colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))(200),
           number.digits=number.digits, number.cex=number.cex))
  title(main=paste0(plot.title, "-corr heatmap, AOE order"), sub=reference, cex=1)
  if(!is.null(similarity.var)){
    print(corrplot(similarity.var, method="shade", shade.col=NA,
             tl.col="black", tl.srt=90, tl.cex=0.7, addCoef.col="black", cl.pos="n", order="AOE",
             col=colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))(200), 
             number.digits=number.digits, number.cex=number.cex))
    title(main=paste0("Heatmap:var(corr)", "-AOE order"), sub=reference, cex=1)
  }
  print(plot(hclust(dist(similarity.mar, method = "euclidean"), method = "ward.D2"),
       main = paste0(plot.title, "-hierarchical clustering\n", reference)))
}
