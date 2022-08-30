# Part1. Methods for gene selection
cds_to_other_id <- function(cds,
                            db,
                            input_file_gene_id_type,
                            new_gene_id_type,
                            verbose = FALSE) {
  matrix <- exprs(cds)
  fdata <- fData(cds)

  new_g <- convert_gene_ids(row.names(fdata),
                            db,
                            input_file_gene_id_type,
                            new_gene_id_type)
  lstart <- length(new_g)
  new_g <- new_g[!is.na(new_g)]
  new_g <- new_g[!duplicated(new_g)]
  lend <- length(new_g)
  if((lstart-lend)/lstart > .7) warning(paste("More than 70% of IDs were lost",
                                              "when converting to",
                                              new_gene_id_type, "IDs. Did you",
                                              "specify the correct gene ID",
                                              "types and the correct db?"))
  if(verbose) message(paste("After converting CDS to", new_gene_id_type,"IDs,",
                            lstart - lend, "IDs were lost"))

  matrix <- matrix[names(new_g),]
  fdata <- fdata[names(new_g),, drop=FALSE]
  row.names(matrix) <- new_g
  row.names(fdata) <- new_g

  pd = new("AnnotatedDataFrame", data = pData(cds))
  fd = new("AnnotatedDataFrame", data = fdata)
  cds = suppressWarnings(newCellDataSet(matrix,
                       phenoData=pd,
                       featureData=fd,
                       expressionFamily=cds@expressionFamily,
                       lowerDetectionLimit=cds@lowerDetectionLimit))

  return(cds)
}

convert_gene_ids <- function(gene_list,
                             db,
                             start_type,
                             end_type) {

  tryCatch({suppressMessages(AnnotationDbi::mapIds(db, keys = gene_list,
                                         column = end_type, start_type))},
           error = function(e) {
             msg <- paste0("Garnett cannot convert the gene IDs using the ",
                         "db and types provided. Please check that your db, ",
                         "cds_gene_id_type and marker_file_gene_id_type ",
                         "parameters are correct. Please note that the ", "
                         cds_gene_id_type refers to the type of the ",
                         "row.names of the feature (gene) table in your cds. ",
                         "Conversion error: ", e)
             stop(msg)
           })
}

#' cell_rules class
#'
#' Representation of cell type derived from marker file.
#'
#' @slot name character. Name of the cell type.
#' @slot gene_names character. A list of all of the genes included in the
#'  definition.
#' @slot expressed character. A list of genes defined as "expressed:".
#' @slot not_expressed character. A list of genes defined as "not expressed:".
#' @slot gene_rules vector of GeneRules-class. A list of genes defined under
#'  specific rules using "expressed below:", "expressed above:", or
#'  "expressed between:".
#' @slot meta data.frame of meta data rules specified in marker file.
#' @slot parenttype character. The name of the parent type - specified by
#'  "subtype of:".
#' @slot references character. A list of references included in the definition.
#'
#' @name cell_rules
#' @rdname cell_rules
#' @aliases cell_rules-class
#' @exportClass cell_rules
setClass("cell_rules", representation(name = "character",
                                      gene_names = "character",
                                      expressed = "character",
                                      not_expressed = "character",
                                      gene_rules = "vector",
                                      meta = "data.frame",
                                      parenttype = "character",
                                      references = "character"))

setGeneric("collect_genes", signature = "x", function(x)
  standardGeneric("collect_genes"))

setMethod("collect_genes", c("x" = "cell_rules"), function(x) {
  n1 <- as.character(c(x@expressed, x@not_expressed))
  n2 <- as.character(unlist(lapply(x@gene_rules, function(y) y@gene_name)))
  x@gene_names <- unique(c(n1, n2))
  x
})

collect_gene_names <- function(cellont_list) {
  genes <- lapply(cellont_list, function(x) {
    if(length(collect_genes(x)@gene_names) != 0) {
      pt <- ifelse(identical(x@parenttype, character(0)), "root", x@parenttype)
      data.frame(genes = collect_genes(x)@gene_names, parent = pt,
                 cell_type = x@name)
    }
  })
  all <- do.call("rbind",genes)
  return(all)
}

# Part2. Methods for garnett classifier
new_garnett_classifier <- function()
{
  garc <- new( "garnett_classifier",
               classification_tree = igraph::graph.empty())
  
  root_node_id <- "root"
  
  garc@classification_tree <- garc@classification_tree +
    igraph::vertex(root_node_id,
                   classify_func=list(function(x) {rep(TRUE, ncol(x))}),
                   model = NULL)
  
  return(garc)
}

add_cell_type <- function(classifier,
                          cell_type_name,
                          classify_func,
                          parent_cell_type_name="root")
{
  if (cell_type_name %in% igraph::V(classifier@classification_tree)$name){
    stop(paste("Error: cell type",cell_type_name, "already exists."))
  }
  
  classifier@classification_tree <- classifier@classification_tree +
    igraph::vertex(cell_type_name, classify_func=list(classify_func),
                   model=NULL)
  
  classifier@classification_tree <- classifier@classification_tree +
    igraph::edge(parent_cell_type_name, cell_type_name)
  return (classifier)
}

make_predictions <- function(cds,
                             classifier,
                             curr_node,
                             rank_prob_ratio,
                             cores = 1,
                             s) {
  cvfit <- igraph::V(classifier@classification_tree)[curr_node]$model[[1]]
  
  predictions <- tryCatch({
    if(is.null(cvfit)) {
      child_cell_types <- igraph::V(classifier@classification_tree)[
        suppressWarnings(outnei(curr_node)) ]$name
      predictions <- matrix(FALSE, nrow=nrow(pData(cds)),
                            ncol=length(child_cell_types),
                            dimnames=list(row.names(pData(cds)),
                                          child_cell_types))
      predictions <- split(predictions, rep(1:ncol(predictions),
                                            each = nrow(predictions)))
      names(predictions) <- child_cell_types
      predictions
    } else {
      candidate_model_genes <- cvfit$glmnet.fit$beta[[1]]@Dimnames[[1]]
      good_genes <- intersect(row.names(exprs(cds)),
                              candidate_model_genes)
      if (length(good_genes) == 0) stop(paste("None of the model genes are in",
                                              "your CDS object. Did you",
                                              "specify the correct",
                                              "cds_gene_id_type and the",
                                              "correct db?"))
      x <- Matrix::t(exprs(cds[intersect(row.names(exprs(cds)),
                                         candidate_model_genes),])) #slow
      
      extra <- as(matrix(0, nrow = nrow(x),
                         ncol = length(setdiff(candidate_model_genes,
                                               colnames(x)))), "sparseMatrix")
      row.names(extra) <- row.names(x)
      colnames(extra) <- setdiff(candidate_model_genes, colnames(x))
      
      x <- cbind(x, extra)
      x <- x[,candidate_model_genes]
      
      # predict probabilities using fitted model
      nonz <- Matrix::rowSums(do.call(cbind,
                                      glmnet:::coef.glmnet(cvfit,
                                                           s="lambda.min")))
      nonz <- nonz[2:length(nonz)]
      nonz <- names(nonz[nonz != 0])
      
      if (sum(!nonz %in% row.names(exprs(cds))) > 0) {
        warning(paste("The following genes used in the classifier are not",
                      "present in the input CDS. Interpret with caution.",
                      nonz[!nonz %in% row.names(exprs(cds))]))
      }
      
      temp <- stats::predict(cvfit, #slow
                             newx = x,
                             s = s,
                             type = "response")
      temp[is.nan(temp)] <- 0
      prediction_probs <- as.matrix(as.data.frame(temp))
      
      # normalize probabilities by dividing by max
      prediction_probs <- prediction_probs/Biobase::rowMax(prediction_probs)
      
      prediction_probs[is.nan(prediction_probs)] <- 0
      
      # find the odds ratio of top prob over second best
      prediction_probs <- apply(prediction_probs, 1, function(x) {
        m <- names(which.max(x))
        s <- sort(x, decreasing = T)
        c(cell_type = m, odds_ratio = s[1]/s[2])
      })
      
      prediction_probs <- as.data.frame(t(prediction_probs))
      prediction_probs$cell_name <- row.names(prediction_probs)
      names(prediction_probs) <- c("cell_type", "odds_ratio", "cell_name")
      prediction_probs$odds_ratio <-
        as.numeric(as.character(prediction_probs$odds_ratio))
      
      # odds ratio has to be larger than rank_prob_ratio
      assignments <- prediction_probs[prediction_probs$odds_ratio >
                                        rank_prob_ratio,]
      
      # odds ratio also must be larger than expected by random guess
      # (1/number of cell types)
      random_guess_thresh <- 1.0 / length(cvfit$glmnet.fit$beta)
      assignments <- assignments[assignments$odds_ratio > random_guess_thresh,]
      
      not_assigned <- row.names(pData(cds))[ !row.names(pData(cds)) %in%
                                               assignments$cell_name]
      if(length(not_assigned) > 0) {
        assignments <- rbind(assignments,
                             data.frame(cell_name = not_assigned,
                                        cell_type = NA, odds_ratio = NA))
      }
      
      assignments$cell_type <- stringr::str_replace_all(assignments$cell_type,
                                                        "\\.1",
                                                        "")
      
      # reformat predictions
      predictions <- reshape2::dcast(assignments, cell_name ~ cell_type,
                                     value.var = "odds_ratio")
      predictions <- predictions[!is.na(predictions$cell_name),]
      row.names(predictions) <- predictions$cell_name
      
      if (ncol(predictions) > 2){
        predictions <- predictions[,setdiff(colnames(predictions), "NA")]
        predictions <- predictions[,-1, drop=FALSE]
        predictions <- predictions[rownames(pData(cds)),,drop=FALSE]
        predictions <- as.matrix(predictions)
        predictions[is.na(predictions)] <- FALSE
        predictions[predictions != 0] <- TRUE
        cell_type_names <- colnames(predictions)
        
        predictions <- split(predictions, rep(1:ncol(predictions),
                                              each = nrow(predictions)))
        names(predictions) <- cell_type_names
        
      } else {
        cell_type_names <- names(cvfit$glmnet.fit$beta)
        one_type <- names(predictions)[2]
        if (one_type == "NA") {
          names(predictions)[2] <- "Unknown"
          one_type <- "Unknown"
        }
        predictions <- matrix(FALSE, nrow=nrow(pData(cds)),
                              ncol=length(cell_type_names),
                              dimnames=list(row.names(pData(cds)),
                                            cell_type_names))
        predictions[,one_type] <- TRUE
        
        predictions <- split(predictions, rep(1:ncol(predictions),
                                              each = nrow(predictions)))
        names(predictions) <- cell_type_names
      }
      predictions
    }
    
  },
  #warning = function(w) print(w),
  error = function(e) {
    if (e$message == paste("None of the model genes are in your CDS object.",
                           "Did you specify the correct cds_gene_id_type and",
                           "the correct db?"))
      stop(e)
    print (e)
    cell_type_names <- names(cvfit$glmnet.fit$beta)
    predictions <- matrix(FALSE, nrow=nrow(pData(cds)),
                          ncol=length(cell_type_names),
                          dimnames=list(row.names(pData(cds)),
                                        cell_type_names))
    predictions <- split(predictions, rep(1:ncol(predictions),
                                          each = nrow(predictions)))
    names(predictions) <- cell_type_names
    predictions
  })
  
  for (i in 1:length(predictions)){
    p <- as(as(predictions[[i]], "sparseVector"), "sparseMatrix")
    row.names(p) <- row.names(pData(cds))
    predictions[[i]] <- p
  }
  
  return(predictions)
}


# Part3. Methods for get training sample
get_training_sample <- function(cds,
                                orig_cds,
                                classifier,
                                tf_idf,
                                gene_table,
                                curr_node,
                                parse_list,
                                name_order,
                                max_training_samples,
                                num_unknown,
                                back_cutoff,
                                training_cutoff,
                                marker_scores,
                                return_initial_assign) {
  
  ##### Find type assignment from expressed/not expressed #####
  
  child_cell_types <- igraph::V(classifier@classification_tree)[
    suppressWarnings(outnei(curr_node)) ]$name
  parent <- igraph::V(classifier@classification_tree)[curr_node]$name
  if (length(child_cell_types) > 0) {
    if (length(intersect(child_cell_types,
                         colnames(marker_scores))) == 0) {
      return(NULL)
    }
    assigns <- assign_type(marker_scores[,intersect(child_cell_types,
                                                    colnames(marker_scores)),
                                         drop=FALSE],
                           training_cutoff, return_initial_assign)
    if(return_initial_assign) {
      return(assigns)
    }
  }
  
  assigns <- as.data.frame(assigns)
  
  names(assigns) <- "assigns"
  pData(orig_cds)$assigns <- assigns[row.names(pData(orig_cds)),"assigns"]
  pData(orig_cds)$assigns <- as.character(pData(orig_cds)$assigns)
  pData(orig_cds)$assigns[is.na(pData(orig_cds)$assigns)] <- "None"
  
  ##### Exclude possibles using other definitions #####
  child_rules <- list()
  for (child in child_cell_types) {
    cell_class_func <-
      igraph::V(classifier@classification_tree)[child]$classify_func[[1]]
    
    parent <- environment(cell_class_func)
    if (is.null(parent))
      parent <- emptyenv()
    e1 <- new.env(parent=parent)
    
    Biobase::multiassign(names(pData(orig_cds)), pData(orig_cds), envir=e1)
    environment(cell_class_func) <- e1
    
    type_res <- cell_class_func(exprs(orig_cds))
    if (length(type_res)!= ncol(orig_cds)){
      message(paste("Error: classification function for",
                    igraph::V(classifier@classification_tree)[child]$name,
                    "returned a malformed result."))
      stop()
    }
    
    type_res <- as(as(type_res,"sparseVector"), "sparseMatrix")
    row.names(type_res) <- row.names(pData(orig_cds))
    colnames(type_res) <- child
    child_rules[[ child ]] <- type_res
  }
  
  
  # Assign cell types by classifier rules and check for conflicts
  ctf_cell_type <- rep("Unknown", length(child_rules[[1]]))
  names(ctf_cell_type) <- row.names(child_rules[[1]])
  for (child in child_cell_types) {
    ctf_cell_type[Matrix::which(as.matrix(child_rules[child][[1]]))] <- child
  }
  
  # Find training sample cells and downsample if necessary
  good_cell_type <- ctf_cell_type[ctf_cell_type %in% child_cell_types]
  obs_counts <- table(good_cell_type)
  training_sample <- c()
  for(i in names(obs_counts)){
    num_obs_for_type_i <- min(max_training_samples, obs_counts[i])
    obs_for_type_i <- sample(which(good_cell_type == i), num_obs_for_type_i)
    training_sample <- append(training_sample, obs_for_type_i)
  }
  training_sample <- good_cell_type[training_sample]
  
  # Find outgroup samples and downsample if necessary
  outgroup_samples <- !ctf_cell_type %in% child_cell_types
  
  if (length(outgroup_samples) > 0){
    outgroup_samples <- ctf_cell_type[outgroup_samples]
    
    out_group_cds <- cds[,names(outgroup_samples)]
    out_group_cds <- out_group_cds[,sample(row.names(pData(out_group_cds)),
                                           min(nrow(pData(out_group_cds)),
                                               num_unknown * 10),
                                           replace = F)]
    
    out_group_cds <- get_communities(out_group_cds)
    
    per_clust <-
      floor(num_unknown/length(unique(pData(out_group_cds)$louv_cluster)))
    
    outg <- lapply(unique(pData(out_group_cds)$louv_cluster), function(x) {
      sub <- pData(out_group_cds)[pData(out_group_cds)$louv_cluster == x,]
      sample(row.names(sub), min(nrow(sub), per_clust))
    })
    
    outg <- unlist(outg)
    if(length(outg) < min(length(outgroup_samples), num_unknown)) {
      to_get <- min(length(outgroup_samples), num_unknown)
      outgroup_samples <- outgroup_samples[!names(outgroup_samples) %in% outg]
      outg <- c(outg, sample(names(outgroup_samples), (to_get - length(outg))))
    }
    outgroup <- rep("Unknown", length(outg))
    names(outgroup) <- outg
    training_sample <- append(training_sample, outgroup)
  }
  
  training_sample <- factor(training_sample)
  training_sample <- droplevels(training_sample)
  
  return(training_sample)
}

assign_type <- function(total_vals,
                        training_cutoff,
                        return_initial_assign) {
  cutoffs <- apply(total_vals, 2, stats::quantile, probs = training_cutoff)
  
  q <- as.data.frame(t(t(total_vals) > cutoffs))
  q$total <- rowSums(q)
  q$assign <- sapply(q$total, function(x) {
    if (x == 0) {
      out <- "None"
    } else if (x > 1) {
      out <- "Ambiguous"
    } else {
      out <- "Assign"
    }
    out
  })
  q$cell <- row.names(total_vals)
  if (return_initial_assign) {
    q$total <- NULL
    types <- colnames(total_vals)
    q$Assignment <- "None"
    for(type in types) {
      q$Assignment[q[,type]] <- type
    }
    q$Assignment[q$assign == "Ambiguous"] <- "Ambiguous"
    q$cell <- NULL
    q$assign <- NULL
    return(q)
  }
  
  out <- q[q$assign != "Assign",]$assign
  names(out) <- q[q$assign != "Assign",]$cell
  sub <- q[q$assign == "Assign",]
  sub$total <- NULL
  sub$assign <- NULL
  sub <- reshape2::melt(sub, id.vars = "cell")
  
  sub <- sub[sub$value,]
  out2 <- as.character(sub$variable)
  names(out2) <- sub$cell
  out <- c(out, out2)
  out <- out[out != "None"]
  out
}

aggregate_positive_markers <- function(cell_type,
                                       tf_idf,
                                       gene_table,
                                       back_cutoff,
                                       agg = TRUE) {
  gene_list <- cell_type@expressed
  bad_genes <- gene_table[!gene_table$in_cds,]$orig_fgenes
  gene_list <- gene_list[!gene_list %in% bad_genes]
  gene_list <- gene_table$fgenes[match(gene_list,gene_table$orig_fgenes)]
  
  rel <- intersect(unlist(gene_list), colnames(tf_idf))
  rel <- unique(rel)
  rel_genes <- tf_idf[,rel, drop=FALSE]
  if(length(rel) == 0) return(NULL)
  
  if(length(unlist(gene_list)) == 1) {
    background_cutoff <- back_cutoff * stats::quantile(rel_genes[,1],
                                                       prob = .95)
    rel_genes[rel_genes < background_cutoff] <- 0
    agg <- rel_genes
  } else {
    background_cutoff <- back_cutoff * apply(rel_genes, 2,
                                             stats::quantile,
                                             prob = .95)
    
    temp <- Matrix::t(rel_genes)
    temp[temp < background_cutoff] <- 0
    rel_genes <- Matrix::t(temp)
    
    if(agg) {
      agg <-  (Matrix::rowSums(rel_genes))
    } else {
      agg <- rel_genes
    }
  }
  agg
}

get_negative_markers <- function(cell_type,
                                 tf_idf,
                                 gene_table,
                                 back_cutoff) {
  not_list <- cell_type@not_expressed
  if (length(not_list) != 0) {
    bad_genes <- gene_table[!gene_table$in_cds,]$orig_fgenes
    not_list <- not_list[!not_list %in% bad_genes]
    not_list <- gene_table$fgenes[match(not_list,gene_table$orig_fgenes)]
    rel <- intersect(unlist(not_list), colnames(tf_idf))
    rel <- unique(rel)
    rel_genes <- tf_idf[,rel, drop=FALSE]
    if(length(rel) == 0) return("")
    
    if(length(unlist(not_list)) == 1) {
      background_cutoff <- back_cutoff * stats::quantile(rel_genes[,1],
                                                         prob = .95)
      bad_cells <- names(rel_genes[rel_genes[,1] > background_cutoff,])
    } else {
      background_cutoff <- back_cutoff * apply(rel_genes, 2,
                                               stats::quantile,
                                               prob = .95)
      
      temp <- Matrix::t(rel_genes)
      temp <- temp > background_cutoff
      temp <- Matrix::t(temp)
      bad_cells <- row.names(temp[Matrix::rowSums(temp) != 0,])
    }
  } else {
    bad_cells <- ""
  }
  return(bad_cells)
}
