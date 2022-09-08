#' @export
trainAnnoModel <- function(expr,
                           input.type = "Seurat",
                           label,
                           markers.Path,
                           model.savePath,
                           gene.method = "Entropy",
                           gene.length = 2000,
                           garnett.cutoff = 0.7,
                           train.ratio = 0.8,
                           repeat.times = 1,
                           dropout.modeling = FALSE,
                           metacell.anno = FALSE){
    if(input.type == "Seurat"){
        data <- data.frame(t(expr@assays$RNA@data))
    }
    data$label <- label
    cnt <- 0
    repeat{
        cnt <- cnt + 1
        if(cnt > repeat.times) {
            break
        }
        set.seed(unclass(Sys.time()))
        ID <- sort(sample(nrow(data), train.ratio * nrow(data)))
        train_set <- data[ID,]      #construct reference set
        test_set <- data[-ID,]      #construct query set
        object <- scibet_visualization(train_set, train_set$label)[["object"]]
        object <- as.CellDataSet(object)
        object <- estimateSizeFactors(object)
        # sample selection
        result <- SelectCells(object,
                              cutoff = garnett.cutoff,
                              marker_file_path = markers.Path)
        representative.index <- result[["index"]]
        markers <- result[["markers"]]
        # feature selection
        geneset <- SelectGenes(train_set[representative.index,],
                               method = gene.method,
                               k = gene.length)
        # training process
        result <- Train(train_set[representative.index,], union(geneset, markers$marker_gene))
        label.name <- sort(unique(data$label))
        prob <- result[,1:length(label.name)]
        lambda <- result[,(length(label.name)+1):dim(result)[2]]
        model.save <- data.frame(t(prob))
        # label prediction
        result <- MarkerScore(test_set,
                              marker_file_path,
                              cutoff=.70,
                              metacell = metacell.anno)
        predict <- Test(prob, lambda, test_set,
                        weighted.markers = result[["markers"]],
                        dropout.modeling = dropout.modeling,
                        average.expr = result[["average"]])
        if(metacell.anno){
            predict <- predict[result[["clustering"]]]
        }
        correct <- 0
        for (i in 1:length(predict)){
            if (test_set$label[i] == predict[i]){
                correct <- correct + 1
            }
        }
        message("Accuracy: ", correct / length(predict))
        saveRDS(model.save, paste0(model.savePath, "model-", cnt ,".rds"))
    }
}

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
similarityCalculation <- function(fine.labels, savePath){
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
    message("[", Sys.time(), "] -----: TME cell subtypes annotation")
    dataset <- data.frame(t(expr@assays$RNA@data))
    dataset$rough.labels <- expr$Cell.Type
    fine.labels <- predSubType(dataset,
                               pretrained.path,
                               savePath,
                               celltype.list,
                               umap.plot)
    similarityCalculation(fine.labels, savePath)
    return(fine.labels)
}
