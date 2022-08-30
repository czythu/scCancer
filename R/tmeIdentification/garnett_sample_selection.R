#' This function takes single-cell expression data in the form of a CDS object
#' and a cell type definition file (marker file) and trains a multinomial
#' classifier to assign cell types. The resulting \code{garnett_classifier}
#' object can be used to classify the cells in the same dataset, or future
#' datasets from similar tissues/samples.
#'
#' @param cds Input CDS object.
#' @param marker_file A character path to the marker file to define cell types.
#'  See details and documentation for \code{\link{Parser}} by running
#'  \code{?Parser}for more information.
#' @param db Bioconductor AnnotationDb-class package for converting gene IDs.
#'  For example, for humans use org.Hs.eg.db. See available packages at
#'  \href{http://bioconductor.org/packages/3.8/data/annotation/}{Bioconductor}.
#'  If your organism does not have an AnnotationDb-class database available,
#'  you can specify "none", however then Garnett will not check/convert gene
#'  IDs, so your CDS and marker file must have the same gene ID type.
#' @param cds_gene_id_type The type of gene ID used in the CDS. Should be one
#'  of the values in \code{columns(db)}. Default is "ENSEMBL". Ignored if
#'  db = "none".
#' @param marker_file_gene_id_type The type of gene ID used in the marker file.
#'  Should be one of the values in \code{columns(db)}. Default is "SYMBOL".
#'  Ignored if db = "none".
#' @param min_observations An integer. The minimum number of representative
#'  cells per cell type required to include the cell type in the predictive
#'  model. Default is 8.
#' @param max_training_samples An integer. The maximum number of representative
#'  cells per cell type to be included in the model training. Decreasing this
#'  number increases speed, but may hurt performance of the model. Default is
#'  500.
#' @param num_unknown An integer. The number of unknown type cells to use as an
#'  outgroup during classification. Default is 500.
#' @param propogate_markers Logical. Should markers from child nodes of a cell
#'  type be used in finding representatives of the parent type? Should
#'  generally be \code{TRUE}.
#' @param cores An integer. The number of cores to use for computation.
#' @param lambdas \code{NULL} or a numeric vector. Allows the user to pass
#'  their own lambda values to \code{\link[glmnet]{cv.glmnet}}. If \code{NULL},
#'  preset lambda values are used.
#' @param classifier_gene_id_type The type of gene ID that will be used in the
#'  classifier. If possible for your organism, this should be "ENSEMBL", which
#'  is the default. Ignored if db = "none".
#' @param return_initial_assign Logical indicating whether an initial
#'  assignment data frame for the root level should be returned instead of a
#'  classifier. This can be useful while choosing/debugging markers. Please
#'  note that this means that a classifier will not be built, so you will not
#'  be able to move on to the next steps of the workflow until you rerun the
#'  functionwith \code{return_initial_assign = FALSE}. Default is \code{FALSE}.
#'
#' @details This function has three major parts: 1) parsing the marker file 2)
#'  choosing cell representatives and 3) training the classifier. Details on
#'  each of these steps is below:
#'
#'  Parsing the marker file: the first step of this function is to parse the
#'  provided marker file. The marker file is a representation of the cell types
#'  expected in the data and known characteristics about them. Information
#'  about marker file syntax is available in the documentation for the
#'  \code{\link{Parser}} function, and on the
#'  \href{https://cole-trapnell-lab.github.io/garnett}{Garnett website}.
#'
#'  Choosing cell representatives: after parsing the marker file, this function
#'  identifies cells that fit the parameters specified in the file for each cell
#'  type. Depending on how marker genes and other cell type definition
#'  information are specified, expression data is normalized and expression
#'  cutoffs are defined automatically. In addition to the cell types in the
#'  marker file, an outgroup of diverse cells is also chosen.
#'
#'  Training the classifier: lastly, this function trains a multinomial GLMnet
#'  classifier on the chosen representative cells.
#'
#'  Because cell types can be defined hierarchically (i.e. cell types can be
#'  subtypes of other cell types), steps 2 and 3 above are performed iteratively
#'  over all internal nodes in the tree representation of cell types.
#'
#'  See the
#'  \href{https://cole-trapnell-lab.github.io/garnett}{Garnett website} and the
#'  accompanying paper for further details.
#'
#' @export
#'
select_fine_samples <- function(cds,                                  marker_file,
                                db,
                                cds_gene_id_type = "ENSEMBL",
                                marker_file_gene_id_type = "SYMBOL",
                                cutoff = .75,
                                min_observations = 8,
                                max_training_samples = 500,
                                num_unknown = 500,
                                propogate_markers = TRUE,
                                cores=1,
                                lambdas = NULL,
                                classifier_gene_id_type = "ENSEMBL",
                                return_initial_assign = FALSE) {
  
  ##### Check inputs #####
  assertthat::assert_that(is(cds, "CellDataSet"))
  assertthat::assert_that(assertthat::has_name(pData(cds), "Size_Factor"),
                          msg = paste("Must run estimateSizeFactors() on cds",
                                      "before calling train_cell_classifier"))
  assertthat::assert_that(sum(is.na(pData(cds)$Size_Factor)) == 0,
                          msg = paste("Must run estimateSizeFactors() on cds",
                                      "before calling train_cell_classifier"))
  assertthat::assert_that(is.character(marker_file))
  assertthat::is.readable(marker_file)
  if (is(db, "character") && db == "none") {
    cds_gene_id_type <- 'custom'
    classifier_gene_id_type <- 'custom'
    marker_file_gene_id_type <- 'custom'
  } else {
    assertthat::assert_that(is(db, "OrgDb"),
                            msg = paste0("db must be an 'AnnotationDb' object ",
                                         "or 'none' see ",
                                         "http://bioconductor.org/packages/",
                                         "3.8/data/annotation/ for available"))
    assertthat::assert_that(is.character(cds_gene_id_type))
    assertthat::assert_that(is.character(marker_file_gene_id_type))
    assertthat::assert_that(cds_gene_id_type %in% AnnotationDbi::keytypes(db),
                            msg = paste("cds_gene_id_type must be one of",
                                        "keytypes(db)"))
    assertthat::assert_that(classifier_gene_id_type %in% AnnotationDbi::keytypes(db),
                            msg = paste("classifier_gene_id_type must be one of",
                                        "keytypes(db)"))
    assertthat::assert_that(marker_file_gene_id_type %in%
                              AnnotationDbi::keytypes(db),
                            msg = paste("marker_file_gene_id_type must be one of",
                                        "keytypes(db)"))
  }
  assertthat::is.count(num_unknown)
  assertthat::is.count(cores)
  assertthat::assert_that(is.logical(propogate_markers))
  if (!is.null(lambdas)) {
    assertthat::assert_that(is.numeric(lambdas))
  }
  
  ##### Set internal parameters #####
  rel_gene_quantile <- .9 # exclusion criterion for genes expressed at greater
  # than rel_gene_quantile in all training cell subsets
  back_cutoff <- 0.25 # percent of 95th percentile of expression that marks the
  # cutoff between "expressed" and "not expressed"
  perc_cells <- 0.05 # percent of training cells a gene is expressed to be
  # included in glmnet training
  training_cutoff <- cutoff # percentile of marker score required for training
  # assignment
  
  ##### Normalize and rename CDS #####
  if (!is(exprs(cds), "dgCMatrix")) {
    sf <- pData(cds)$Size_Factor
    pd <- new("AnnotatedDataFrame", data = pData(cds))
    fd <- new("AnnotatedDataFrame", data = fData(cds))
    cds <- suppressWarnings(newCellDataSet(as(exprs(cds), "dgCMatrix"),
                                           phenoData = pd,
                                           featureData = fd))
    pData(cds)$Size_Factor <- sf
  }
  
  pData(cds)$num_genes_expressed <- Matrix::colSums(as(exprs(cds),
                                                       "lgCMatrix"))
  cell_totals <-  Matrix::colSums(exprs(cds))
  sf <- pData(cds)$Size_Factor
  
  pd <- new("AnnotatedDataFrame", data = pData(cds))
  fd <- new("AnnotatedDataFrame", data = fData(cds))
  temp <- exprs(cds)
  temp@x <- temp@x / rep.int(pData(cds)$Size_Factor, diff(temp@p))
  norm_cds <- suppressWarnings(newCellDataSet(temp,
                                              phenoData = pd, featureData = fd))
  orig_cds <- cds
  if(cds_gene_id_type != classifier_gene_id_type)  {
    norm_cds <- cds_to_other_id(norm_cds, db=db, cds_gene_id_type,
                                classifier_gene_id_type)
    orig_cds <- cds_to_other_id(cds, db=db, cds_gene_id_type,
                                classifier_gene_id_type)
  }
  pData(norm_cds)$Size_Factor <- sf
  
  
  ##### Parse Marker File #####
  file_str = paste0(readChar(marker_file, file.info(marker_file)$size),"\n")
  
  parse_list <- parse_input(file_str)
  orig_name_order <- unlist(parse_list[["name_order"]])
  rm("name_order", envir=parse_list)
  
  # Check and order subtypes
  ranks <- lapply(orig_name_order, function(i) parse_list[[i]]@parenttype)
  names(ranks) <- orig_name_order
  
  if(length(unlist(unique(ranks[which(!ranks %in% names(ranks) & lengths(ranks) != 0L)])) != 0)) {
    stop(paste("Subtype", unlist(unique(ranks[which(!ranks %in% names(ranks) & lengths(ranks) != 0L)])), "is not defined in marker file."))
  }
  
  if(any(names(ranks) == ranks)) {
    bad <- ranks[names(ranks) == ranks]
    stop(paste0("'", bad,
                "' cannot be a subtype of itself. Please modify marker file."))
  }
  
  name_order <- names(ranks[lengths(ranks) == 0L])
  ranks <- ranks[!names(ranks) %in% name_order]
  while(length(ranks) != 0) {
    name_order <- c(name_order, names(ranks)[ranks %in% name_order])
    ranks <- ranks[!names(ranks) %in% name_order]
  }
  
  
  if(is.null(parse_list)) stop("Parse failed!")
  message(paste("There are", length(parse_list), "cell type definitions"))
  
  # Check gene names and keywords
  gene_table <- make_name_map(parse_list,
                              as.character(row.names(fData(norm_cds))),
                              classifier_gene_id_type,
                              marker_file_gene_id_type,
                              db)
  
  ##### Make garnett_classifier #####
  classifier <- new_garnett_classifier()
  classifier@gene_id_type <- classifier_gene_id_type
  if(is(db, "character") && db == "none") classifier@gene_id_type <- "custom"
  
  for(i in name_order) {
    # check meta data exists
    if (nrow(parse_list[[i]]@meta) != 0) {
      if (!all(parse_list[[i]]@meta$name %in% colnames(pData(norm_cds)))) {
        bad_meta <- parse_list[[i]]@meta$name[!parse_list[[i]]@meta$name %in%
                                                colnames(pData(norm_cds))]
        stop(paste0("Cell type '", parse_list[[i]]@name,
                    "' has a meta data specification '", bad_meta ,
                    "' that's not in the pData table."))
      }
    }
    logic_list <- assemble_logic(parse_list[[i]], gene_table)
    classifier <- add_cell_rule(parse_list[[i]], classifier, logic_list)
  }
  
  classifier@cell_totals <- exp(mean(log(cell_totals)))/
    stats::median(pData(norm_cds)$num_genes_expressed)
  
  ##### Create transformed marker table #####
  if(propogate_markers) {
    root <- propogate_func(curr_node = "root", parse_list, classifier)
  }
  
  tf_idf <- tfidf(norm_cds) #slow
  
  
  ### Aggregate markers ###
  marker_scores <- data.frame(cell = row.names(tf_idf))
  
  for (i in name_order) {
    agg <- aggregate_positive_markers(parse_list[[i]], tf_idf,
                                      gene_table, back_cutoff)
    bad_cells <- get_negative_markers(parse_list[[i]], tf_idf,
                                      gene_table, back_cutoff)
    if(is.null(agg))  {
      warning (paste("Cell type", i, "has no genes that are expressed",
                     "and will be skipped"))
    } else {
      agg[names(agg) %in% bad_cells] <- 0
      marker_scores <- cbind(marker_scores, as.matrix(agg))
      colnames(marker_scores)[ncol(marker_scores)] <- parse_list[[i]]@name
    }
  }
  # return(marker_scores)
  ##### Sample Selection #####
  for (v in igraph::V(classifier@classification_tree)){
    child_cell_types <- igraph::V(classifier@classification_tree)[
      suppressWarnings(outnei(v))]$name
    
    if(length(child_cell_types) > 0) {
      ### Get CDS subset for training ###
      message(igraph::V(classifier@classification_tree) [ v ]$name)
      if(igraph::V(classifier@classification_tree) [ v ]$name == "root") {
        cds_sub <- norm_cds
        orig_sub <- orig_cds
      } else {
        # loosely classify to subset
        new_assign <-
          make_predictions(norm_cds,
                           classifier,
                           igraph::V(classifier@classification_tree)[
                             suppressWarnings(innei(v))]$name,
                           rank_prob_ratio = 1.1,
                           s = "lambda.min")
        if(!igraph::V(classifier@classification_tree)[v]$name %in%
           names(new_assign)) {
          message(paste0("No cells classified as ",
                         igraph::V(classifier@classification_tree) [ v ]$name,
                         ". No subclassification"))
          next
        }
        good_cells <-
          as.matrix(new_assign[
            igraph::V(classifier@classification_tree)[v]$name][[1]])
        good_cells <- names(good_cells[good_cells[,1] != 0,])
        if(length(good_cells) == 0) {
          message(paste0("No cells classified as ",
                         igraph::V(classifier@classification_tree) [ v ]$name,
                         ". No subclassification"))
          next
        }
        cds_sub <- norm_cds[,good_cells]
        orig_sub <- orig_cds[,good_cells]
      }
      
      ### Get training sample ###
      sample <- get_training_sample(cds = cds_sub,
                                    orig_cds = orig_sub,
                                    classifier,
                                    tf_idf,
                                    gene_table,
                                    v,
                                    parse_list,
                                    name_order,
                                    max_training_samples,
                                    num_unknown,
                                    back_cutoff,
                                    training_cutoff,
                                    marker_scores,
                                    return_initial_assign)

      if(return_initial_assign) {
        return(sample)
      }
    }
  }
  return(training_samples)
}

parse_input <- function(file_str,
                        debug = F) {
  # Parse input_file
  lexer  <- rly::lex(Lexer, debug=debug)
  parser <- rly::yacc(Parser, debug=debug)
  parse_list <- parser$parse(file_str, lexer)

  parse_list
}

make_name_map <- function(parse_list,
                          possible_genes,
                          cds_gene_id_type,
                          marker_file_gene_id_type,
                          db) {
  gene_start <- collect_gene_names(parse_list)
  gene_table <- data.frame(fgenes = gene_start[,1], parent = gene_start[,2])
  gene_table$parent <- as.character(gene_table$parent)
  gene_table$fgenes <- as.character(gene_table$fgenes)
  gene_table$orig_fgenes <- gene_table$fgenes
  if(cds_gene_id_type != marker_file_gene_id_type) {
    gene_table$fgenes <- convert_gene_ids(gene_table$orig_fgenes,
                                          db,
                                          marker_file_gene_id_type,
                                          cds_gene_id_type)
    bad_convert <- sum(is.na(gene_table$fgenes))
    if (bad_convert > 0) warning(paste(bad_convert,
                                 "genes could not be converted from",
                                 marker_file_gene_id_type,
                                 "to", cds_gene_id_type, "These genes are",
                                 "listed below:", paste0(gene_table$orig_genes[
                                   is.na(gene_table$fgenes)],
                                   collapse="\n")))
  } else {
    gene_table$cds <- gene_table$fgenes
  }

  if(cds_gene_id_type == "ENSEMBL" | marker_file_gene_id_type == "ENSEMBL") {
    gene_table$cds <- NULL
    possibles <- data.frame(cds = possible_genes,
                            ensembl = as.character(
                              stringr::str_split_fixed(possible_genes,
                                                       "\\.",
                                                       2)[,1]))
    gene_table <- merge(gene_table, possibles, all.x=T,
                        by.x="fgenes", by.y="ensembl")
    gene_table$fgenes <- gene_table$cds
  } else {
    gene_table$cds <- gene_table$fgenes
  }

  gene_table$in_cds <- gene_table$f %in% possible_genes
  gene_table$in_cds[is.na(gene_table$in_cds)] <- FALSE

  bad_genes <- gene_table$orig_fgenes[!gene_table$in_cds]
  if (length(bad_genes) > 0) warning(strwrap("The following genes from
                                             the cell type definition file are
                                             not present in the cell dataset.
                                             Please check these genes for
                                             errors. Cell type determination
                                             will continue, ignoring these
                                             genes."), "\n",
                                     paste0(bad_genes, collapse="\n"))

  gene_table$fgenes <- as.character(gene_table$fgenes)
  gene_table$cds <- as.character(gene_table$cds)

  gene_table
}

add_cell_rule <- function(cell_type,
                          classifier,
                          logic_list) {
  # Set parenttype to root if no parent
  if (length(cell_type@parenttype) == 0) {
    cell_type@parenttype <- "root"
  }
  # subtype of
  if (length(cell_type@parenttype) > 1) stop("only 1 parenttype allowed")
  parent_type <- as.character(cell_type@parenttype)

  # references
  if (length(cell_type@references) > 0) {
    if (length(classifier@references) == 0) {
      classifier@references <- list()
    }

    classifier@references <- c(classifier@references,
                               list(cell_type@references))
    names(classifier@references)[length(classifier@references)] <-
      cell_type@name
  }

  if (length(logic_list) == 0) {
    warning (paste("Cell type", cell_type@name,
                   "has no valid rules and will be skipped"))
    classifier <- add_cell_type(classifier, cell_type@name,
                                classify_func = function(x) {rep(FALSE, ncol(x))},
                                parent_type)
    return(classifier)
  }
  logic <- paste(unlist(logic_list), collapse = ' & ')

  tryCatch(
    if(nchar(logic) == 0) {
      classifier <- add_cell_type(classifier, cell_type@name,
                                  classify_func = function(x) {FALSE},
                                  parent_type)

    } else {
      classifier <- add_cell_type(classifier, cell_type@name,
                                  classify_func = function(x) {
                                    eval(parse(text = logic))
                                  },
                                  parent_type)
    },
    error = function(e) {
      msg <- paste("Cell type rule generation failed on the",
                   "cell definition for ", cell_type@name, ".\nError: ",e)
      stop(msg)
    }
  )
  return(classifier)
}

assemble_logic <- function(cell_type,
                           gene_table) {

  logic = ""
  logic_list = list()
  bad_genes <- gene_table[!gene_table$in_cds,]$orig_fgenes

  # expressed/not expressed
  logic_list <- lapply(cell_type@gene_rules, function(rule) {
    log_piece <- ""
    if (!rule@gene_name %in% bad_genes) {
      paste0("(x['",
             gene_table$fgenes[match(rule@gene_name, gene_table$orig_fgenes)],
             "',] > ", rule@lower,
             ") & (x['",
             gene_table$fgenes[match(rule@gene_name, gene_table$orig_fgenes)],
             "',] < ",
             rule@upper,
             ")")
    }
  })
  if (length(cell_type@expressed) > 0 | length(cell_type@not_expressed) > 0) {
    logic_list <- list(logic_list, paste0("assigns == '", cell_type@name, "'"))
  }

  if(length(logic_list) == 0) warning(paste("Cell type", cell_type@name,
                                            "has no valid expression rules."))

  # meta data
  if (nrow(cell_type@meta) > 0) {
    mlogic <- plyr::dlply(cell_type@meta, plyr::.(name), function(x) {
      if(nrow(x) == 1){
        out <- paste0(x["name"], " %in% c('", x[,"spec"][1],"')")
      } else {
        out <- paste0(x[,"name"][1], " %in% c('", paste(x[,"spec"],
                                                        collapse = "', '"),
                      "')")
      }
      out
    })
    logic_list <- c(logic_list, unname(mlogic))
  }
  logic_list <- logic_list[!is.na(logic_list)]
  logic_list
}

propogate_func <- function(curr_node,
                           parse_list,
                           classifier) {
  children <- igraph::V(classifier@classification_tree)[
    suppressWarnings(outnei(curr_node))]$name

  if(length(children) == 0) {
    return(parse_list[[curr_node]]@expressed)
  } else {
    child_genes <- c()
    if (curr_node != "root") {
      child_genes <- parse_list[[curr_node]]@expressed
    }
    for(child in children) {
      child_genes <- union(child_genes,
                           propogate_func(child, parse_list, classifier))
    }
    if(curr_node != "root") {
      parse_list[[curr_node]]@expressed <- child_genes
    }
    return(child_genes)
  }
}

tfidf <- function(input_cds) {
  ncounts <- exprs(input_cds)
  ncounts <- ncounts[Matrix::rowSums(ncounts) != 0,]
  nfreqs <- ncounts
  nfreqs@x <- ncounts@x / rep.int(Matrix::colSums(ncounts), diff(ncounts@p))
  tf_idf_counts <- nfreqs * log(1 + ncol(ncounts) / Matrix::rowSums(ncounts > 0))
  Matrix::t(tf_idf_counts)
}
