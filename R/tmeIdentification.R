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
                           garnett.cutoff = 0.75,
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
                        submodel.path,
                        markers.path,
                        savePath,
                        celltype.list,
                        unknown.cutoff = 0.3,
                        umap.plot = FALSE){
    # Split test dataset with rough labels.
    finelabels.list <- lapply(celltype.list, function(celltype){
        testdata <- test_set[which(test_set$rough.labels == celltype),]
        folder.path1 <- paste0(submodel.path, "/", celltype, "/")
        folder.path2 <- paste0(markers.path, "/", celltype, "/")
        file.path1 <- paste0(folder.path1, list.files(folder.path1))
        file.path2 <- paste0(folder.path2, list.files(folder.path2))
        # Different classification principles: Several lists of subtype
        pdf(file = paste0(savePath, celltype, ".pdf"), width = 7, height = 7)
        subtypes.predict <- lapply(file.path1, function(model.path){
            index <- which(file.path == model.path)
            model.ref <- read.csv(model.path)
            model.ref <- pro.core(model.ref)
            # Return a list of subtype
            result <- MarkerScore(testdata,
                                  marker_file_path = file.path2[index],
                                  cutoff = unknown.cutoff)
            label.predict <- Test(model.ref, lambda, testdata,
                            weighted.markers = result[["markers"]],
                            dropout.modeling = dropout.modeling,
                            average.expr = result[["average"]])
            label.predict <- paste0(label.predict, "[", index, "]")
            label.predict <- AssignUnknown(label.predict, result[["unknown"]])
            if(umap.plot){
                print(scibet_visualization(testdata, label.predict)[["plot"]])
            }
            return(label.predict)
        })
        dev.off()
        message("[", Sys.time(), "] -----: ", celltype, " subtype annotation finished.")
        subtypes.predict <- data.frame(matrix(unlist(subtypes.predict),
                                              nrow = length(subtypes.predict),
                                              byrow = T))
        names(subtypes.predict) <- rownames(testdata)
        return(subtypes.predict)
    })
    names(finelabels.list) <- celltype.list
    return(finelabels.list)
}

#' similarityCalculation
#' @param fine.labels annotation results (barcode-subtype) of function"predSubType"
#' @inheritParams runScAnnotation
#'
#' @return similarity matrixes of all possible celltypes
#' @export
similarityCalculation <- function(fine.labels,
                                  savePath){
    all.matrix <- list()
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
        # small similarity map
        if(dim(similarity.mar)[1] <= 4^2){
            pdf(file = paste0(savePath, "similarity-", celltype, ".pdf"), width = 6, height = 8)
            SimilarityMap(plot.title, "reference = ...", similarity.mar,
                          number.digits = 2, number.cex = 1, tl.cex = 1)
        }
        # huge similarity map
        else{
            pdf(file = paste0(savePath, "similarity-", celltype, ".pdf"), width = 12, height = 15)
            SimilarityMap(plot.title, "reference = ...", similarity.mar,
                          number.digits = 1, number.cex = 0.6, tl.cex = 0.7)
        }
        dev.off()
        all.matrix[[celltype]] <- similarity.mar
    }
    return(all.matrix)
}

#' runCellSubtypeClassify
#' @param expr A Seurat object.
#' @param cell.annotation A data.frame of cells' annotation.
#' @inheritParams runScAnnotation
#'
#' @return A list of updated Seurat object, cell.annotation.
#' @export
runCellSubtypeClassify <- function(expr,
                                   cell.annotation,
                                   submodel.path,
                                   markers.path,
                                   savePath,
                                   celltype.list,
                                   unknown.cutoff,
                                   umap.plot){
    message("[", Sys.time(), "] -----: TME cell subtypes annotation")
    dataset <- data.frame(t(expr@assays$RNA@data))
    dataset$rough.labels <- expr$Cell.Type
    fine.labels <- predSubType(dataset,
                               submodel.path,
                               markers.path,
                               savePath,
                               celltype.list,
                               unknown.cutoff,
                               umap.plot)
    similarity.matrix <- similarityCalculation(cell.annotation$fine.labels,
                                               savePath)
    expr$cell.subtype <- fine.labels
    cell.annotation$fine.labels <- fine.labels
    cell.annotation$similarity.matrix <- similarity.matrix
    return(list(expr = expr,
                cell.annotation = cell.annotation))
}

#' predMalignantCell
#' @param expr A Seurat object.
#' @param cell.annotation A data.frame of cells' annotation.
#' @param MALIGNANT.THRES A threshold of xgboost score
#' to decide whether a cell is malignant. Default is 0.5.
#' @inheritParams runScAnnotation
#'
#' @return A list of cell.annotation and malignancy plots
#' @export
#' @import xgboost
predMalignantCell <- function(expr,
                              cell.annotation,
                              savePath,
                              coor.names = c("tSNE_1", "tSNE_2"),
                              MALIGNANT.THRES = 0.5,
                              model.path = NULL,
                              genes.path = NULL){
    model.path <- paste0(system.file("txt", package = "scCancer2"), "/sc_xgboost.model")
    genes.path <- paste0(system.file("txt", package = "scCancer2"), "/selectGenesByVar.txt")
    model.ref <- xgb.load(model.path)
    genes.preselected <- read.table(genes.path)$V1
    # expr <- readRDS("D:/scCancer-data/SubtypeAnno/Demo-pipeline/KC-example/result-all/expr.RDS")
    testdata <- t(as.matrix(expr@assays$RNA@data))
    testdata <- testdata[,which(colnames(testdata) %in% genes.preselected)]
    testdata <- xgb.DMatrix(testdata)
    predict.label <- predict(model.ref, testdata)

    # store results
    cell.annotation$Malign.score <- predict.label
    # expr$Malign.score <- predict.label
    predict.label[which(predict.label > MALIGNANT.THRES)] <- "malignant"
    predict.label[which(predict.label <= MALIGNANT.THRES)] <- "nonMalignant"
    cell.annotation$Malign.score <- predict.label
    # expr$Malign.type <- predict.label
    # p1 <- DimPlot(expr, reduction = "tsne", group.by = "Malign.score")
    # p2 <- DimPlot(expr, reduction = "tsne", group.by = "Malign.type")

    # plot
    p.results <- plotMalignancy(cell.annotation = cell.annotation,
                                coor.names = coor.names,
                                savePath = savePath)
    p.results[["p.malignScore"]] <- p.malignScore
    ggsave(filename = file.path(savePath, "figures/malignScore.png"),
           p.malignScore, width = 5, height = 4, dpi = 500)

    return(list(cell.annotation = cell.annotation,
                plot = p.results))
}
