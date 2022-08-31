#' @export
predSubType <- function(test_set,
                        pretrained.path,
                        savePath,
                        celltype.list,
                        umap.plot = FALSE){
  # Split test dataset with rough labels.
  celltypes <- sort(unique(test_set$rough.labels))
  finelabels.list <- lapply(celltypes, function(celltype){
    if(celltype %in% celltype.list){
      message(celltype)
      testdata <- test_set[which(test_set$rough.labels == celltype),]
      folder.path <- paste0(pretrained.path, "/", celltype, "/")
      file.name <- list.files(folder.path)
      file.path <- paste0(folder.path, file.name)
      # Different classification principles: Several lists of subtype
      pdf(file = paste0(savePath, celltype, ".pdf"), width = 7, height = 7)
      subtypes.predict <- lapply(file.path, function(model.path){
        model.ref <- read.csv(model.path)
        model.ref <- pro.core(model.ref)
        # Return a list of subtype
        label.predict <- CrossTest(model.ref, testdata)
        label.predict <- paste0(label.predict, "[", which(file.path == model.path), "]")
        if(umap.plot){
          print(scibet_visualization(testdata, label.predict)[[2]])
        }
        return(label.predict)
      })
      dev.off()
      message("Classification finished.")
      subtypes.predict <- data.frame(matrix(unlist(subtypes.predict), 
                                            nrow = length(subtypes.predict),
                                            byrow = T))
      names(subtypes.predict) <- rownames(testdata)
      return(subtypes.predict)
    }    
  })
  names(finelabels.list) <- celltypes
  return(finelabels.list)
}

#' @export
Similarity_Calculation <- function(fine.labels, savePath){
  for(i in 1:length(fine.labels)){
    predict <- fine.labels[[i]]
    if(is.null(predict)){
      next
    }
    all.results <- c()
    all.labels <- c()
    for(j in 1:dim(predict)[1]){
      all.results <- c(all.results, predict[j, ])
      all.labels <- c(all.labels, unique(as.list(predict[j, ])))
    }
    celltype <- names(fine.labels)[i]
    cell.sets <- lapply(all.labels, function(label){
      index <- which(all.results == label)
      cells <- names(all.results)[index]
      return(cells)
    })
    names(cell.sets) <- all.labels
    # Jaccard similarity calculation.
    similarity.mar <- Jaccard(cell.sets)
    similarity.mar[is.na(similarity.mar)] <- 0
    # Heatmap and Hierarchical clustering
    plot.title <- paste0("similarity map of ", celltype)
    pdf(file = paste0(savePath, "similarity-", celltype, ".pdf"), width = 12, height = 15)
    SimilarityMap(plot.title, "reference = ...", similarity.mar,
                  number.digits = 1, number.cex = 0.6, tl.cex = 0.7)
    dev.off()
  }
}

#' @export
runCellSubtypeClassify <- function(expr, 
                                   pretrained.path, 
                                   savePath, 
                                   celltype.list,
                                   umap.plot){
  dataset <- data.frame(t(expr@assays$RNA@data))
  dataset$rough.labels <- expr$Cell.Type
  fine.labels <- predSubType(dataset,
                             pretrained.path,
                             savePath,
                             celltype.list,
                             umap.plot)
  Similarity_Calculation(fine.labels, savePath)
  return(fine.labels)
}