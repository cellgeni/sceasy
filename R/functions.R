.regularise_df <- function(df, drop_single_values = TRUE) {
    if (ncol(df) == 0) df[['name']] <- rownames(df)
    if (drop_single_values) {
        k_singular <- sapply(df, function(x) length(unique(x)) == 1)
        if (sum(k_singular) > 0)
            warning(paste('Dropping single category variables:'),
                    paste(colnames(df)[k_singular], collapse=', '))
            df <- df[, !k_singular, drop=F]
        if (ncol(df) == 0) df[['name']] <- rownames(df)
    }
    return(df)
}

seurat2anndata <- function(
    obj, outFile = NULL, assay = 'RNA', main_layer = 'data', transfer_layers = NULL
) {
    main_layer <- match.arg(main_layer, c('data', 'counts', 'scale.data'))
    transfer_layers <- transfer_layers[
        transfer_layers %in% c('data', 'counts', 'scale.data')]
    transfer_layers <- transfer_layers[transfer_layers != main_layer]

    if (compareVersion(as.character(obj@version), '3.0.0') < 0)
        obj <- Seurat::UpdateSeuratObject(object = obj)

    X <- Seurat::GetAssayData(object = obj, assay = assay, slot = main_layer)

    obs <- .regularise_df(obj@meta.data)

    var <- .regularise_df(Seurat::GetAssay(obj, assay = assay)@meta.features)

    obsm <- NULL
    reductions <- names(obj@reductions)
    if (length(reductions) > 0) {
        obsm <- sapply(
            reductions,
            function(name) as.matrix(Seurat::Embeddings(obj, reduction=name)),
            simplify = FALSE
        )
        names(obsm) <- paste0('X_', tolower(names(obj@reductions)))
    }

    layers <- list()
    for (layer in transfer_layers) {
        mat <- Seurat::GetAssayData(object = obj, assay = assay, slot = layer)
        if (all(dim(mat) == dim(X))) layers[[layer]] <- Matrix::t(mat)
    }

    anndata <- reticulate::import('anndata', convert = FALSE)

    adata <- anndata$AnnData(
        X = Matrix::t(X),
        obs = obs,
        var = var,
        obsm = obsm,
        layers = layers
    )

    if (!is.null(outFile))
        adata$write(outFile, compression = 'gzip')

    adata
}

sce2anndata <- function(
    obj, outFile = NULL, main_layer = 'counts', transfer_layers = NULL
) {
    assay_names <- SummarizedExperiment::assayNames(obj)
    main_layer <- match.arg(main_layer, assay_names)
    transfer_layers <- transfer_layers[transfer_layers %in% assay_names]
    transfer_layers <- transfer_layers[transfer_layers != main_layer]

    X <- SummarizedExperiment::assay(obj, main_layer)

    obs <- .regularise_df(as.data.frame(SingleCellExperiment::colData(obj)))

    var <- .regularise_df(as.data.frame(SingleCellExperiment::rowData(obj)))

    obsm <- NULL
    reductions <- SingleCellExperiment::reducedDimNames(obj)
    if (length(reductions) > 0) {
        obsm <- sapply(
            reductions,
            function(name) as.matrix(
                    SingleCellExperiment::reducedDim(obj, type=name)),
            simplify = FALSE
        )
        names(obsm) <- paste0(
            'X_', tolower(SingleCellExperiment::reducedDimNames(obj)))
    }

    layers <- list()
    for (layer in transfer_layers) {
        mat <- SummarizedExperiment::assay(obj, name)
        if (all(dim(mat) == dim(X))) layers[[layer]] <- Matrix::t(mat)
    }

    anndata <- reticulate::import('anndata', convert = FALSE)

    adata <- anndata$AnnData(
        X = Matrix::t(X),
        obs = obs,
        var = var,
        obsm = obsm,
        layers = layers
    )

    if (!is.null(outFile))
        adata$write(outFile, compression = 'gzip')

    adata
}

loom2anndata <- function(
    inFile, outFile = NULL, main_layer = c('spliced', 'unspliced'),
    obs_names = 'CellID', var_names = 'Gene'
) {
    main_layer <- match.arg(main_layer)

    anndata <- reticulate::import('anndata', convert = FALSE)

    if (compareVersion(as.character(anndata[['__version__']]), '0.6.20') < 0)
        message(paste(
            "Warning: anndata <0.6.20 detected.",
            "Upgrade to handle multi-dimensional embeddings."
        ))

    adata <- anndata$read_loom(
        inFile, sparse = TRUE, cleanup = TRUE, X_name = main_layer,
        obs_names = obs_names, var_names = var_names
    )

    anndata$AnnData$obs_names_make_unique(adata)
    anndata$AnnData$var_names_make_unique(adata)

    if (!is.null(outFile))
        adata$write(outFile, compression = 'gzip')

    adata
}

seurat2sce <- function(obj, outFile = NULL, main_layer=NULL, assay='RNA', ...) {
    sce <- Seurat::as.SingleCellExperiment(obj, assay=assay, ...)
    if (!is.null(outFile))
        saveRDS(sce, outFile)

    sce
}

sce2loom <- function(obj, outFile, main_layer = NULL, ...) {
    SingleCellExperiment::colData(obj) <- .regularise_df(SingleCellExperiment::colData(obj))
    SingleCellExperiment::rowData(obj) <- .regularise_df(SingleCellExperiment::rowData(obj))
    writeExchangeableLoom(obj, outFile, main_layer = main_layer, ...)
}

loom2sce <- function(inFile, outFile = NULL, main_layer = NULL, main_layer_name = NULL, ...) {
    sce <- readExchangeableLoom(inFile, backed = FALSE, main_layer_name = main_layer_name, ...)
    if (!is.null(outFile))
        saveRDS(sce, outFile)

    sce
}
