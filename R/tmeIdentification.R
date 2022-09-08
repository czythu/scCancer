#' @export
trainSubAnnoModel <- function(data,
                              label,
                              markers.Path,
                              model.savePath,
                              gene.method,
                              gene.length,
                              garnett.cutoff,
                              train.ratio,
                              repeat.times,
                              umap.visualization,
                              confusion.matrix,
                              dropout.modeling,
                              metacell.anno){
    data$label <- label
    cnt <- 0
    models <- c()
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
                              cutoff = garnett.cutoff,
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
        predict.unknown <- AssignUnknown(predict, result[["unknown"]])
        # umap visualization
        if(umap.visualization){
            print(scibet_visualization(test_set, predict))
            print(scibet_visualization(test_set, predict.unknown))
        }
        # confusion matrix for analysis
        if(confusion.matrix){
            print(ConfusionMatrix(label.name, label.name,
                                  test_set$label, predict,
                                  xlab='Reference',ylab='Prediction', normalize=F))
            print(ConfusionMatrix(label.name, c(label.name, "unknown"),
                                  test_set$label, predict.unknown,
                                  xlab='Reference',ylab='Prediction', normalize=F))
        }
        models <- c(models, model.save)
        write.csv(model.save, paste0(model.savePath, "model-", cnt ,".csv"))
    }
    return(models)
}

#' @export
trainAnnoModel <- function(expr,
                           rough.label,
                           fine.label,
                           markers.Path,
                           model.savePath,
                           input.type = "Seurat",
                           gene.method = "Entropy",
                           gene.length = 2000,
                           garnett.cutoff = 0.7,
                           train.ratio = 0.8,
                           repeat.times = 1,
                           umap.visualization = TRUE,
                           confusion.matrix = TRUE,
                           dropout.modeling = FALSE,
                           metacell.anno = FALSE){
    if(input.type == "Seurat"){
        data <- data.frame(t(expr@assays$RNA@data))
    }
    celltype <- sort(unique(rough.label))
    allmodels <- list()
    for(i in 1:length(celltype)){
        subtype.label <- fine.label[which(rough.label == celltype[i])]
        save.folder <- paste0(model.savePath, "/", celltype[i], "/")
        models <- trainSubAnnoModel(data = data,
                                    label = subtype.label,
                                    markers.Path = markers.Path,
                                    model.savePath = save.folder,
                                    gene.method = gene.method,
                                    gene.length = gene.length,
                                    garnett.cutoff = garnett.cutoff,
                                    train.ratio = train.ratio,
                                    repeat.times = repeat.times,
                                    umap.visualization = umap.visualization,
                                    confusion.matrix = confusion.matrix,
                                    dropout.modeling = dropout.modeling,
                                    metacell.anno = metacell.anno)
        allmodels[[celltype[i]]] <- models
    }
    return(allmodels)
}

#' @export
predSubType <- function(test_set,
                        pretrained.path,
                        savePath,
                        celltype.list,
                        umap.plot = FALSE){
    # Split test dataset with rough labels.
    finelabels.list <- lapply(celltype.list, function(celltype){
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
