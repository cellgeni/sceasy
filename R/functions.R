# Functions for converting between different objects


#' Regularise dataframe
#'
#' This function checks if certain columns of a dataframe is of a single value
#' and drop them if required
#'
#' @param df Input data frame, usually cell metadata table (data.frame-like
#'   object)
#' @param drop_single_values Drop columns with only a single value (logical)
#'
#' @return Dataframe
.regularise_df <- function(df, drop_single_values = TRUE) {
  if (ncol(df) == 0) df[["name"]] <- rownames(df)
  if (drop_single_values) {
    k_singular <- sapply(df, function(x) length(unique(x)) == 1)
    if (sum(k_singular) > 0) {
      warning(
        paste("Dropping single category variables:"),
        paste(colnames(df)[k_singular], collapse = ", ")
      )
    }
    df <- df[, !k_singular, drop = F]
    if (ncol(df) == 0) df[["name"]] <- rownames(df)
  }
  return(df)
}

#' Convert Seurat object to AnnData
#'
#' This function converts a Seurat object to an Anndata object
#'
#' @param obj Input Seurat object
#' @param outFile Save output AnnData to this file if specified (str or NULL)
#' @param assay Assay to be converted, default "RNA" (str)
#' @param main_layer Slot in `assay` to be converted to AnnData.X, may be
#'   "counts", "data", "scale.data", default "data" (str)
#' @param transfer_layers If specified, convert slots to AnnData.layers[<slot>],
#'   (vector of str)
#' @param drop_single_values Drop single value columns in cell metadata table,
#'   default TRUE (logical)
#'
#' @return AnnData object
#'
#' @import reticulate
#' @import Matrix
seurat2anndata <- function(obj, outFile = NULL, assay = "RNA", main_layer = "data", transfer_layers = NULL, drop_single_values = TRUE) {
  if (!requireNamespace("Seurat")) {
    stop("This function requires the 'Seurat' package.")
  }
  main_layer <- match.arg(main_layer, c("data", "counts", "scale.data"))
  transfer_layers <- transfer_layers[
    transfer_layers %in% c("data", "counts", "scale.data")
  ]
  transfer_layers <- transfer_layers[transfer_layers != main_layer]

  if (compareVersion(as.character(obj@version), "3.0.0") < 0) {
    obj <- Seurat::UpdateSeuratObject(object = obj)
  }

  X <- Seurat::GetAssayData(object = obj, assay = assay, slot = main_layer)

  obs <- .regularise_df(obj@meta.data, drop_single_values = drop_single_values)

  var <- .regularise_df(Seurat::GetAssay(obj, assay = assay)@meta.features, drop_single_values = drop_single_values)

  obsm <- NULL
  reductions <- names(obj@reductions)
  if (length(reductions) > 0) {
    obsm <- sapply(
      reductions,
      function(name) as.matrix(Seurat::Embeddings(obj, reduction = name)),
      simplify = FALSE
    )
    names(obsm) <- paste0("X_", tolower(names(obj@reductions)))
  }

  layers <- list()
  for (layer in transfer_layers) {
    mat <- Seurat::GetAssayData(object = obj, assay = assay, slot = layer)
    if (all(dim(mat) == dim(X))) layers[[layer]] <- Matrix::t(mat)
  }

  anndata <- reticulate::import("anndata", convert = FALSE)

  adata <- anndata$AnnData(
    X = Matrix::t(X),
    obs = obs,
    var = var,
    obsm = obsm,
    layers = layers
  )

  if (!is.null(outFile)) {
    adata$write(outFile, compression = "gzip")
  }

  adata
}

#' Convert SingleCellExperiment object to AnnData
#'
#' This function converts a SingleCellExperiment object to an Anndata object
#'
#' @param obj Input SingleCellExperiment object
#' @param outFile Save output AnnData to this file if specified (str or NULL)
#' @param main_layer Assay to be converted to AnnData.X, default "counts" (str)
#' @param transfer_layers If specified, convert additional assays to
#'   AnnData.layers[<slot>], (vector of str)
#' @param drop_single_values Drop single value columns in cell metadata table,
#'   default TRUE (logical)
#'
#' @return AnnData object
#'
#' @import reticulate
#' @import Matrix
sce2anndata <- function(obj, outFile = NULL, main_layer = "counts", transfer_layers = NULL, drop_single_values = TRUE) {
  if (!requireNamespace("SummarizedExperiment")) {
    stop("This function requires the 'SummarizedExperiment' package.")
  }
  if (!requireNamespace("SingleCellExperiment")) {
    stop("This function requires the 'SingleCellExperiment' package.")
  }
  assay_names <- SummarizedExperiment::assayNames(obj)
  main_layer <- match.arg(main_layer, assay_names)
  transfer_layers <- transfer_layers[transfer_layers %in% assay_names]
  transfer_layers <- transfer_layers[transfer_layers != main_layer]

  X <- SummarizedExperiment::assay(obj, main_layer)

  obs <- .regularise_df(as.data.frame(SummarizedExperiment::colData(obj)), drop_single_values = drop_single_values)

  var <- .regularise_df(as.data.frame(SummarizedExperiment::rowData(obj)), drop_single_values = drop_single_values)

  obsm <- NULL
  reductions <- SingleCellExperiment::reducedDimNames(obj)
  if (length(reductions) > 0) {
    obsm <- sapply(
      reductions,
      function(name) {
        as.matrix(
          SingleCellExperiment::reducedDim(obj, type = name)
        )
      },
      simplify = FALSE
    )
    names(obsm) <- paste0(
      "X_", tolower(SingleCellExperiment::reducedDimNames(obj))
    )
  }

  layers <- list()
  for (layer in transfer_layers) {
    mat <- SummarizedExperiment::assay(obj, layer)
    if (all(dim(mat) == dim(X))) layers[[layer]] <- Matrix::t(mat)
  }

  anndata <- reticulate::import("anndata", convert = FALSE)

  adata <- anndata$AnnData(
    X = Matrix::t(X),
    obs = obs,
    var = var,
    obsm = obsm,
    layers = layers
  )

  if (!is.null(outFile)) {
    adata$write(outFile, compression = "gzip")
  }

  adata
}

#' Convert Loom object to AnnData
#'
#' This function converts a Loom object to an Anndata object
#'
#' @param inFile Path to the input Loom file on disk (str)
#' @param outFile Save output AnnData to this file if specified (str or NULL)
#' @param main_layer Assay to be converted to AnnData.X, can be "spliced" or
#'   "unspliced", default "spliced" (str)
#' @param obs_names Column in cell metadata table that contains cell IDs (str)
#' @param var_names Column in gene metadata table that contains gene names
#'
#' @return AnnData object
#'
#' @import reticulate
loom2anndata <- function(inFile, outFile = NULL, main_layer = c("spliced", "unspliced"),
                         obs_names = "CellID", var_names = "Gene") {
  main_layer <- match.arg(main_layer)

  anndata <- reticulate::import("anndata", convert = FALSE)

  if (compareVersion(as.character(anndata[["__version__"]]), "0.6.20") < 0) {
    message(paste(
      "Warning: anndata <0.6.20 detected.",
      "Upgrade to handle multi-dimensional embeddings."
    ))
  }

  adata <- anndata$read_loom(
    inFile,
    sparse = TRUE, cleanup = TRUE, X_name = main_layer,
    obs_names = obs_names, var_names = var_names
  )

  anndata$AnnData$obs_names_make_unique(adata)
  anndata$AnnData$var_names_make_unique(adata)

  if (!is.null(outFile)) {
    adata$write(outFile, compression = "gzip")
  }

  adata
}

#' Convert Seurat object to SingleCellExperiment
#'
#' This function converts a Seurat object to a SingleCellExperiment object
#'
#' @param obj Input Seurat object
#' @param outFile Save output SingleCellExperiment object to this file if
#'   specified (str or NULL)
#' @param assay Assay to be converted to AnnData.X, default "RNA" (str)
#'
#' @return AnnData object
seurat2sce <- function(obj, outFile = NULL, main_layer = NULL, assay = "RNA", ...) {
  if (!requireNamespace("Seurat")) {
    stop("This function requires the 'Seurat' package.")
  }
  sce <- Seurat::as.SingleCellExperiment(obj, assay = assay, ...)
  if (!is.null(outFile)) {
    saveRDS(sce, outFile)
  }

  sce
}

#' Convert SingleCellExperiment object to Loom file
#'
#' This function converts a SingleCellExperiment object to a LoomExperiment
#' object
#'
#' @param obj Input SingleCellExperiment object
#' @param outFile Save output Loom to this file if specified (str or NULL)
#' @param drop_single_values Drop single value columns in cell metadata table,
#'   default TRUE (logical)
#'
#' @return LoomExperiment object
sce2loom <- function(obj, outFile, main_layer = NULL, drop_single_values = TRUE, ...) {
  if (!requireNamespace("LoomExperiment")) {
    stop("This function requires the 'LoomExperiment' package.")
  }
  scle <- LoomExperiment::SingleCellLoomExperiment(obj)

  if (!is.null(outFile)) {
    saveRDS(scle, outFile)
  }

  scle
}

#' Read a Loom file into a SingleCellExperiment object
#'
#' This function reads a Loom file object into a SingleCellExperiment object
#'
#' @param inFile Path to input Loom file on disk (str)
#' @param outFile Save output SingleCellExperiment object to this file if
#'   specified (str or NULL)
#'
#' @return LoomExperiment object
loom2sce <- function(inFile, outFile = NULL, main_layer = NULL, main_layer_name = NULL, ...) {
  if (!requireNamespace("LoomExperiment")) {
    stop("This function requires the 'LoomExperiment' package.")
  }
  if (!requireNamespace("SingleCellExperiment")) {
    stop("This function requires the 'LoomExperiment' package.")
  }
  scle <- LoomExperiment::import(inFile)
  sce <- as(scle, "SingleCellExperiment")

  if (!is.null(outFile)) {
    saveRDS(sce, outFile)
  }

  sce
}

#' Prepare cell metadata
#'
#' This function prepare cell metadata from AnnData.obs
#'
#' @param obs_pd Input AnnData.obs dataframe
#' @param assay Assay name, default "RNA" (str)
#'
#' @return AnnData object
#'
#' @import reticulate
.obs2metadata <- function(obs_pd, assay = "RNA") {
  obs_df <- .regularise_df(reticulate::py_to_r(obs_pd), drop_single_values = FALSE)
  attr(obs_df, "pandas.index") <- NULL
  colnames(obs_df) <- sub("n_counts", paste0("nCounts_", assay), colnames(obs_df))
  colnames(obs_df) <- sub("n_genes", paste0("nFeaturess_", assay), colnames(obs_df))
  return(obs_df)
}

#' Prepare feature metadata
#'
#' This function prepare feature metadata from AnnData.var
#'
#' @param var_pd Input AnnData.var dataframe
#'
#' @return AnnData object
#'
#' @import reticulate
.var2feature_metadata <- function(var_pd) {
  var_df <- .regularise_df(reticulate::py_to_r(var_pd), drop_single_values = FALSE)
  attr(var_df, "pandas.index") <- NULL
  colnames(var_df) <- sub("dispersions_norm", "mvp.dispersion.scaled", colnames(var_df))
  colnames(var_df) <- sub("dispersions", "mvp.dispersion", colnames(var_df))
  colnames(var_df) <- sub("means", "mvp.mean", colnames(var_df))
  colnames(var_df) <- sub("highly_variable", "highly.variable", colnames(var_df))
  return(var_df)
}

.uns2misc <- function(ad, target_uns_keys = list()) {
  uns_keys <- intersect(target_uns_keys, reticulate::py_to_r(ad$uns_keys()))
  misc <- sapply(uns_keys, function(x) reticulate::py_to_r(ad$uns[x]), simplify = FALSE, USE.NAMES = TRUE)
  return(misc)
}

#' Convert AnnData object to Seurat object
#'
#' This function converts an AnnData object to a Seurat object
#'
#' @param inFile Path to an input AnnData object on disk (str)
#' @param outFile Save output Seurat to this file if specified (str or NULL)
#' @param assay Name of assay in Seurat object to store expression values,
#'   default "RNA" (str)
#' @param main_layer Name of slot in `assay` to store AnnData.X, can be
#'   "counts", "data", "scale.data", default "counts" (str)
#' @param use_seurat Use Seurat::ReadH5AD() to do the conversion, default FALSE (logical)
#' @param lzf Whether AnnData is compressed by `lzf`, default FALSE (logical)
#'
#' @return Seurat object
#'
#' @import reticulate
#' @import Matrix
anndata2seurat <- function(inFile, outFile = NULL, main_layer = "counts", assay = "RNA", use_seurat = FALSE, lzf = FALSE, target_uns_keys = list()) {
  if (!requireNamespace("Seurat")) {
    stop("This function requires the 'Seurat' package.")
  }
  main_layer <- match.arg(main_layer, c("counts", "data", "scale.data"))
  inFile <- path.expand(inFile)

  anndata <- reticulate::import("anndata", convert = FALSE)
  sp <- reticulate::import("scipy.sparse", convert = FALSE)

  if (use_seurat) {
    if (lzf) {
      tmpFile <- paste0(tools::file_path_sans_ext(inFile), ".decompressed.h5ad")
      ad <- anndata$read_h5ad(inFile)
      ad$write(tmpFile)
      tryCatch(
        {
          srt <- Seurat::ReadH5AD(tmpFile)
        },
        finally = {
          file.remove(tmpFile)
        }
      )
    } else {
      srt <- Seurat::ReadH5AD(inFile)
    }
  } else {
    ad <- anndata$read_h5ad(inFile)

    obs_df <- .obs2metadata(ad$obs)
    var_df <- .var2feature_metadata(ad$var)

    if (reticulate::py_to_r(sp$issparse(ad$X))) {
      X <- Matrix::t(reticulate::py_to_r(sp$csc_matrix(ad$X)))
    } else {
      X <- t(reticulate::py_to_r(ad$X))
    }
    colnames(X) <- rownames(obs_df)
    rownames(X) <- rownames(var_df)

    if (!is.null(reticulate::py_to_r(ad$raw))) {
      raw_var_df <- .var2feature_metadata(ad$raw$var)
      raw_X <- Matrix::t(reticulate::py_to_r(sp$csc_matrix(ad$raw$X)))
      colnames(raw_X) <- rownames(obs_df)
      rownames(raw_X) <- rownames(raw_var_df)
    } else {
      raw_var_df <- NULL
      raw_X <- NULL
    }

    if (main_layer == "scale.data" && !is.null(raw_X)) {
      assays <- list(Seurat::CreateAssayObject(data = raw_X))
      assays[[1]] <- Seurat::SetAssayData(assays[[1]], slot = "scale.data", new.data = X)
      message("X -> scale.data; raw.X -> data")
    } else if (main_layer == "data" && !is.null(raw_X)) {
      if (nrow(X) != nrow(raw_X)) {
        message("Raw layer was found with different number of genes than main layer, resizing X and raw.X to match dimensions")
        raw_X <- raw_X[rownames(raw_X) %in% rownames(X), , drop = F]
        X <- X[rownames(raw_X), , drop = F]
      }
      assays <- list(Seurat::CreateAssayObject(counts = raw_X))
      assays[[1]] <- Seurat::SetAssayData(assays[[1]], slot = "data", new.data = X)
      message("X -> data; raw.X -> counts")
    } else if (main_layer == "counts") {
      assays <- list(Seurat::CreateAssayObject(counts = X))
      message("X -> counts")
    } else {
      assays <- list(Seurat::CreateAssayObject(data = X))
      message("X -> data")
    }
    names(assays) <- assay

    if (main_layer == "scale.data" && !is.null(raw_X)) {
      assays[[assay]]@meta.features <- raw_var_df
    } else {
      assays[[assay]]@meta.features <- var_df
    }

    project_name <- sub("\\.h5ad$", "", basename(inFile))
    srt <- new("Seurat", assays = assays, project.name = project_name, version = packageVersion("Seurat"))
    Seurat::DefaultAssay(srt) <- assay
    Seurat::Idents(srt) <- project_name

    srt@meta.data <- obs_df
    embed_names <- unlist(reticulate::py_to_r(ad$obsm_keys()))
    if (length(embed_names) > 0) {
      embeddings <- sapply(embed_names, function(x) reticulate::py_to_r(ad$obsm[x]), simplify = FALSE, USE.NAMES = TRUE)
      names(embeddings) <- embed_names
      for (name in embed_names) {
        rownames(embeddings[[name]]) <- colnames(assays[[assay]])
      }

      dim.reducs <- vector(mode = "list", length = length(embeddings))
      for (i in seq(length(embeddings))) {
        name <- embed_names[i]
        embed <- embeddings[[name]]
        key <- switch(name,
          sub("_(.*)", "\\L\\1", sub("^X_", "", toupper(name)), perl = T),
          "X_pca" = "PC",
          "X_tsne" = "tSNE",
          "X_umap" = "UMAP"
        )
        colnames(embed) <- paste0(key, "_", seq(ncol(embed)))
        dim.reducs[[i]] <- Seurat::CreateDimReducObject(
          embeddings = embed,
          loadings = new("matrix"),
          assay = assay,
          stdev = numeric(0L),
          key = paste0(key, "_")
        )
      }
      names(dim.reducs) <- sub("X_", "", embed_names)

      for (name in names(dim.reducs)) {
        srt[[name]] <- dim.reducs[[name]]
      }
    }
  }

  srt@misc <- .uns2misc(ad, target_uns_keys = target_uns_keys)

  if (!is.null(outFile)) saveRDS(object = srt, file = outFile)

  srt
}

anndata2cds <- function(inFile, outFile = NULL) {
  message("not implemented yet")
}
