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
        object <- visualization_pipeline(train_set, train_set$label)[["object"]]
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
            print(visualization_pipeline(test_set, predict))
            print(visualization_pipeline(test_set, predict.unknown))
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

#' predSubType_Scoring (Old version, scoring)
#' @param test_set An expression matrix.
#' Rows should be cells and the last column should be "rough label".
#' @param unknown.cutoff A threshold for assignment of unknown label. Default is 0.3.
#' @inheritParams runScAnnotation
#'
#' @return A list of fine.labels containing all possible celltypes

predSubType_Scoring <- function(expr,
                        submodel.path,
                        markers.path,
                        savePath,
                        celltype.list,
                        dropout.modeling,
                        unknown.cutoff,
                        umap.plot){
    test_set <- data.frame(t(expr@assays$RNA@data))
    test_set$rough.labels <- expr$Cell.Type
    finelabels.list <- lapply(celltype.list, function(celltype){
        t.expr <- expr
        message(celltype)
        # Split test dataset with rough labels
        testdata <- test_set[which(test_set$rough.labels == celltype),]
        if(dim(testdata)[1] < 100){
            message(celltype, " not enough for subtype annotation. Skip!")
            return(NULL)
        }
        barcodes <- rownames(testdata)
        folder.path1 <- paste0(submodel.path, "/", celltype, "/")
        folder.path2 <- paste0(markers.path, "/", celltype, "/")
        file.path1 <- paste0(folder.path1, list.files(folder.path1))
        file.path2 <- paste0(folder.path2, list.files(folder.path2))
        # Different classification principles: Several lists of subtype
        pdf(file = file.path(savePath, paste0("umap-", celltype, ".pdf")),
            width = 6, height = 5)
        subtypes.predict <- lapply(file.path1, function(model.path){
            message(model.path)
            index <- which(file.path1 == model.path)
            model.ref <- read.csv(model.path)
            model.ref <- pro.core(model.ref)
            suppressWarnings(result <- MarkerScore(test_set = testdata,
                                  marker_file_path = file.path2[index],
                                  cutoff = unknown.cutoff))
            # Return a list of subtype
            output <- Test(model.ref, lambda, testdata,
                           weighted.markers = result[["markers"]],
                           dropout.modeling = dropout.modeling,
                           average.expr = result[["average"]])
            label.predict <- output[["predict"]]
            likelihoods <- output[["likelihoods"]]
            likelihoods <- t(apply(likelihoods, 1, function(l){return(l - min(l))}))
            probability <- likelihoods / rowSums(likelihoods)
            saveRDS(probability, paste0(savePath, "normalized-likelihood-", index, ".rds"))
            label.predict <- paste0(label.predict, " (", index, ")")
            label.predict <- AssignUnknown(NULL, label.predict, result[["unknown"]])[["predict.unknown"]]
            if(umap.plot){
                names(label.predict) <- barcodes
                tt.expr <- AddMetaData(object = t.expr,
                                      metadata = label.predict,
                                      col.name = "cell.subtype")
                tt.expr$cell.subtype[which(is.na(tt.expr$cell.subtype))] <- "NA"
                colors.assigned <- getDefaultColors(length(unique(label.predict)))
                colors.assigned <- c(colors.assigned, "#c1cdcd")
                names(colors.assigned) <- c(unique(label.predict), "NA")
                if(celltype == "T.cells" | celltype == "Myeloid.cells"){
                    legend.size <- 10
                }
                else{
                    legend.size <- 12
                }
                print(DimPlot(tt.expr, group.by = "cell.subtype",
                              repel = TRUE, label = FALSE,
                              cols = colors.assigned)+
                          theme(legend.position = "right",
                                legend.text = element_text(size = legend.size)))
                rm(tt.expr)
                # print(visualization_pipeline(testdata, label.predict)[["plot"]])
            }
            return(label.predict)
        })
        normalized.likelihood <- lapply(seq_len(length(file.path1)), function(index){
            file.path <- paste0(savePath, "normalized-likelihood-", index, ".rds")
            likelihood <- readRDS(file.path)

            file.remove(file.path)
            return(likelihood)
        })
        dev.off()
        message("[", Sys.time(), "] -----: ", celltype, " subtype annotation finished.")
        subtypes.predict <- data.frame(matrix(unlist(subtypes.predict),
                                              nrow = length(subtypes.predict),
                                              byrow = T))
        colnames(subtypes.predict) <- barcodes
        return(list(label = t(subtypes.predict),
                    normalized.likelihood = normalized.likelihood))
    })
    names(finelabels.list) <- celltype.list
    return(finelabels.list)
}


# adjust testdata for XGBoost(common features)
align_XGBoost <- function(test, barcodes, features){
    temp <- matrix(nrow = nrow(test), ncol = length(features),
                   dimnames = list(barcodes, features))
    current.features <- colnames(test)
    for(j in 1:length(features)){
        if(features[j] %in% current.features){
            temp[,j] <- test[, features[j]]
        }
        else{
            temp[,j] <- rep(0, length(barcodes))
            # temp[,j] <- sample(seq(0,10), length(barcodes), replace=TRUE)
        }

    }
    return(temp)
}

# ensemble learning, vote for cellsubtype labels
ensemble_XGBoost <- function(label.matrix){
    label <- apply(label.matrix, 2, function(labels){
        tab <- table(labels)
        if(max(tab) >= 3){
            return(names(which.max(tab)))
        }
        else{
            return("unknown")
        }
    })
    return(label)
}


#' predSubType_XGBoost (XGBoost version)
#' @param test_set An expression matrix.
#' Rows should be cells and the last column should be "rough label".
#' @inheritParams runScAnnotation
#'
#' @return A list of fine.labels containing all possible celltypes

predSubType_XGBoost <- function(expr,
                        submodel.path,
                        savePath,
                        celltype.list,
                        umap.plot){
    test_set <- data.frame(t(expr@assays$RNA@data))
    models <- readRDS(submodel.path)
    finelabels.list <- lapply(celltype.list, function(celltype){
        t.expr <- expr
        message(celltype)
        # Split test dataset with rough labels
        testdata <- test_set[which(t.expr$Cell.Type == celltype),]
        if(dim(testdata)[1] < 100){
            message(celltype, " not enough for subtype annotation. Skip!")
            return(NULL)
        }
        barcodes <- rownames(testdata)
        testdata <- as.matrix(testdata)
        subtypes.predict <- matrix(nrow = 1, ncol = length(barcodes))
        # Different classification principles: Several lists of subtype
        pdf(file = file.path(savePath, paste0("umap-", celltype, ".pdf")),
            width = 7, height = 7)
        for(index in 1:length(models)){
            model <- models[[index]]
            dataset.name <- names(models)[index]
            model.celltype <- strsplit(dataset.name, split = "_")[[1]][1]
            # "T", "M", "B", "E", "F"
            if(substr(model.celltype, 1, 1) != substr(celltype, 1, 1)){
                next
            }
            message(dataset.name)
            # Boosting(5 models)
            label.predict <- matrix(nrow = length(model[["models"]]), ncol = length(barcodes))
            mapping <- model[["mapping"]]
            for(i in 1:length(model[["models"]])){
                weak.model <- model[["models"]][[i]]
                # Construct testdata
                features <- weak.model[["feature_names"]]
                test <- testdata[,which(colnames(testdata) %in% features)]
                test <- align_XGBoost(test, barcodes, features)
                test <- xgb.DMatrix(test)
                label.predict[i,] <- names(mapping)[1 + predict(weak.model, test)]
                label.predict[i,] <- paste0(label.predict[i,], " (", index, ")")
            }
            label.predict <- ensemble_XGBoost(label.predict)
            if(umap.plot){
                names(label.predict) <- barcodes
                tt.expr <- AddMetaData(object = t.expr,
                                       metadata = label.predict,
                                       col.name = "cell.subtype")
                print(DimPlot(tt.expr, group.by = "cell.subtype",
                              repel = TRUE, label = FALSE, label.size = 3))
                rm(tt.expr)
                # print(visualization_pipeline(testdata, label.predict)[["plot"]])
            }
            subtypes.predict <- rbind(subtypes.predict, label.predict)
        }
        dev.off()
        message("[", Sys.time(), "] -----: ", celltype, " subtype annotation finished.")
        subtypes.predict <- subtypes.predict[-1,]
        colnames(subtypes.predict) <- barcodes
        return(t(subtypes.predict))
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
similarityCalculation <- function(fine.labels, savePath){
    all.matrix <- lapply(names(fine.labels), function(celltype){
        predict <- fine.labels[[celltype]]
        if(!is.null(predict)){
            predict <- t(predict)
            all.results <- c()
            all.labels <- c()
            for(j in 1:dim(predict)[1]){
                all.results <- c(all.results, predict[j, ])
                all.labels <- c(all.labels, unique(as.list(predict[j, ])))
            }
            all.labels <- all.labels[-which(all.labels == "unknown")]
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
                pdf(file = file.path(savePath, paste0("similarity-", celltype, ".pdf")),
                    width = 8, height = 8)
                # SimilarityMap(plot.title, "reference = ...", similarity.mar,
                #               number.digits = 2, number.cex = 1, tl.cex = 1)
                SimilarityHeatmap(similarity.mar, celltype)
            }
            # huge similarity map
            else{
                pdf(file = file.path(savePath, paste0("similarity-", celltype, ".pdf")),
                    width = 15, height = 15)
                # SimilarityMap(plot.title, "reference = ...", similarity.mar,
                #               number.digits = 1, number.cex = 0.6, tl.cex = 0.7)
                SimilarityHeatmap(similarity.mar, celltype)
            }
            dev.off()
            return(similarity.mar)
        }
    })
    names(all.matrix) <- names(fine.labels)
    return(all.matrix)
}

#' runCellSubtypeClassify
#' @param expr A Seurat object.
#' @inheritParams runScAnnotation
#'
#' @return A list of fine.labels and similarity matrix
#' @export
runCellSubtypeClassify <- function(expr,
                                   submodel.path,
                                   markers.path,
                                   savePath,
                                   celltype.list,
                                   dropout.modeling,
                                   unknown.cutoff,
                                   umap.plot){
    message("[", Sys.time(), "] -----: TME cell subtypes annotation")
    fine.labels <- predSubType_Scoring(expr = expr,
                               submodel.path = submodel.path,
                               markers.path = markers.path,
                               savePath = savePath,
                               celltype.list = celltype.list,
                               dropout.modeling = dropout.modeling,
                               unknown.cutoff = unknown.cutoff,
                               umap.plot = umap.plot)

    fine.labels <- lapply(fine.labels, function(label){
        return(label[["label"]])
    })
    # fine.labels <- predSubType_XGBoost(expr,
    #                             submodel.path,
    #                             savePath,
    #                             celltype.list,
    #                             umap.plot)
    similarity.matrix <- similarityCalculation(fine.labels, savePath)

    return(list(fine.labels = fine.labels,
                similarity.matrix = similarity.matrix))
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
                              malignancy.method,
                              savePath,
                              coor.names = c("UMAP_1", "UMAP_2"),
                              MALIGNANT.THRES = 0.5,
                              model.path = NULL,
                              genes.path = NULL){
    # model.path <- paste0(system.file("txt", package = "scCancer2"), "/sc_xgboost.model")
    # genes.path <- paste0(system.file("txt", package = "scCancer2"), "/genes-scRNA-tcga-sorted.txt")
    model.path <- paste0(system.file("txt", package = "scCancer"), "/sc_xgboost.model")
    genes.path <- paste0(system.file("txt", package = "scCancer"), "/genes-scRNA-tcga-sorted.txt")
    model.ref <- xgb.load(model.path)
    # features <- read.table(genes.path)$V1
    features <- as.list(read.table(genes.path))[[1]]
    testdata <- t(as.matrix(expr@assays$RNA@data))
    testdata <- testdata[,which(colnames(testdata) %in% features)]
    testdata <- align_XGBoost(testdata, rownames(testdata), features)
    testdata <- xgb.DMatrix(testdata)
    predict.label <- predict(model.ref, testdata)

    # store results
    cell.annotation$Malign.score <- predict.label
    # expr$Malign.score <- predict.label
    predict.label[which(predict.label > MALIGNANT.THRES)] <- "malignant"
    predict.label[which(predict.label <= MALIGNANT.THRES)] <- "nonMalignant"
    cell.annotation$Malign.type <- predict.label
    # expr$Malign.type <- predict.label
    # p1 <- DimPlot(expr, reduction = "tsne", group.by = "Malign.score")
    # p2 <- DimPlot(expr, reduction = "tsne", group.by = "Malign.type")

    # plot
    p.results <- plotMalignancy(cell.annotation = cell.annotation,
                                malignancy.method = malignancy.method,
                                coor.names = coor.names,
                                savePath = savePath)
    return(list(cell.annotation = cell.annotation,
                plot = p.results))
}
