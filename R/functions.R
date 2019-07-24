seurat2anndata <- function(
    obj, outFile = NULL, main_layer = 'data', transfer_layers = NULL
) {
    main_layer <- match.arg(main_layer, c('data', 'counts', 'scale.data'))
    transfer_layers <- transfer_layers[
        transfer_layers %in% c('data', 'counts', 'scale.data')]
    transfer_layers <- transfer_layers[transfer_layers != main_layer]

    if (compareVersion(as.character(obj@version), '3.0.0') < 0)
        obj <- Seurat::UpdateSeuratObject(object = obj)

    X <- Seurat::GetAssayData(object = obj, assay = 'RNA', slot = main_layer)

    var <- Seurat::GetAssay(obj)@meta.features
    if (ncol(var) == 0) var[['feature_name']] <- rownames(var)

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
        mat <- Seurat::GetAssayData(object = obj, assay = 'RNA', slot = layer)
        if (all(dim(mat) == dim(X))) layers[[layer]] <- Matrix::t(mat)
    }

    anndata <- reticulate::import('anndata', convert = FALSE)

    adata <- anndata$AnnData(
        X = Matrix::t(X),
        obs = obj@meta.data,
        var = var,
        obsm = obsm,
        layers = layers
    )

    if (!is.null(outFile))
        anndata$AnnData$write(adata, outFile, compression = 'gzip')

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

    var <- as.data.frame(SingleCellExperiment::rowData(obj))
    if (ncol(var) == 0) var[['feature_name']] <- rownames(var)

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
        obs = as.data.frame(SingleCellExperiment::colData(obj)),
        var = var,
        obsm = obsm,
        layers = layers
    )

    if (!is.null(outFile))
        anndata$AnnData$write(adata, outFile, compression = 'gzip')

    adata
}

loom2anndata <- function(
    inFile, outFile = NULL, main_layer = 'spliced', obs_names = 'CellID',
    var_names = 'Gene'
) {
    anndata <- reticulate::import('anndata', convert = FALSE)

    adata <- anndata$read_loom(
        inFile, sparse = TRUE, cleanup = TRUE, X_name = main_layer,
        obs_names = obs_names, var_names = var_names
    )

    anndata$AnnData$obs_names_make_unique(adata)
    anndata$AnnData$var_names_make_unique(adata)

    if (!is.null(outFile))
        anndata$AnnData$write(adata, outFile, compression = 'gzip')

    adata
}
