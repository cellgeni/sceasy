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
    obj, outFile = NULL, assay = 'RNA', main_layer = 'data', transfer_layers = NULL, drop_single_values = TRUE
) {
    main_layer <- match.arg(main_layer, c('data', 'counts', 'scale.data'))
    transfer_layers <- transfer_layers[
        transfer_layers %in% c('data', 'counts', 'scale.data')]
    transfer_layers <- transfer_layers[transfer_layers != main_layer]

    if (compareVersion(as.character(obj@version), '3.0.0') < 0)
        obj <- Seurat::UpdateSeuratObject(object = obj)

    X <- Seurat::GetAssayData(object = obj, assay = assay, slot = main_layer)

    obs <- .regularise_df(obj@meta.data, drop_single_values = drop_single_values)

    var <- .regularise_df(Seurat::GetAssay(obj, assay = assay)@meta.features, drop_single_values = drop_single_values)

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
    obj, outFile = NULL, main_layer = 'counts', transfer_layers = NULL, drop_single_values = TRUE
) {
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
            function(name) as.matrix(
                    SingleCellExperiment::reducedDim(obj, type=name)),
            simplify = FALSE
        )
        names(obsm) <- paste0(
            'X_', tolower(SingleCellExperiment::reducedDimNames(obj)))
    }

    layers <- list()
    for (layer in transfer_layers) {
        mat <- SummarizedExperiment::assay(obj, layer)
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
    SummarizedExperiment::colData(obj) <- .regularise_df(SummarizedExperiment::colData(obj), drop_single_values = drop_single_values)
    SummarizedExperiment::rowData(obj) <- .regularise_df(SummarizedExperiment::rowData(obj), drop_single_values = drop_single_values)
    writeExchangeableLoom(obj, outFile, main_layer = main_layer, ...)
}

loom2sce <- function(inFile, outFile = NULL, main_layer = NULL, main_layer_name = NULL, ...) {
    sce <- readExchangeableLoom(inFile, backed = FALSE, main_layer_name = main_layer_name, ...)
    if (!is.null(outFile))
        saveRDS(sce, outFile)

    sce
}

.obs2metadata <- function(obs_pd, assay='RNA') {
    obs_df <- .regularise_df(reticulate::py_to_r(obs_pd), drop_single_values = FALSE)
    attr(obs_df, 'pandas.index') <- NULL
    colnames(obs_df) <- sub('n_counts', paste0('nCounts_', assay), colnames(obs_df))
    colnames(obs_df) <- sub('n_genes', paste0('nFeaturess_', assay), colnames(obs_df))
    return(obs_df)
}

.var2feature_metadata <- function(var_pd) {
    var_df <- .regularise_df(reticulate::py_to_r(var_pd), drop_single_values = FALSE)
    attr(var_df, 'pandas.index') <- NULL
    colnames(var_df) <- sub('dispersions_norm', 'mvp.dispersion.scaled', colnames(var_df))
    colnames(var_df) <- sub('dispersions', 'mvp.dispersion', colnames(var_df))
    colnames(var_df) <- sub('means', 'mvp.mean', colnames(var_df))
    colnames(var_df) <- sub('highly_variable', 'highly.variable', colnames(var_df))
    return(var_df)
}

anndata2seurat <- function(inFile, outFile = NULL, main_layer = 'counts', assay = 'RNA', use_seurat = FALSE, lzf = FALSE) {
    main_layer <- match.arg(main_layer, c('counts', 'data', 'scale.data'))
    inFile <- path.expand(inFile)

    anndata <- reticulate::import('anndata', convert = FALSE)
    sp <- reticulate::import('scipy.sparse', convert = FALSE)

    if (use_seurat) {
        if (lzf) {
            tmpFile <- paste0(tools::file_path_sans_ext(inFile), '.decompressed.h5ad')
            ad <- anndata$read(inFile)
            ad$write(tmpFile)
            tryCatch({
                srt <- Seurat::ReadH5AD(tmpFile)
            }, finally = {
                file.remove(tmpFile)
            })
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

        if (main_layer == 'scale.data' && !is.null(raw_X)) {
            assays <- list(Seurat::CreateAssayObject(data = raw_X))
            assays[[1]] <- Seurat::SetAssayData(assays[[1]], slot = 'scale.data', new.data = X)
            message('X -> scale.data; raw.X -> data')
        } else if (main_layer == 'counts') {
            assays <- list(Seurat::CreateAssayObject(counts = X))
            message('X -> counts')
        } else {
            assays <- list(Seurat::CreateAssayObject(data = X))
            message('X -> data')
        }
        names(assays) <- assay

        if (main_layer == 'scale.data' && !is.null(raw_X)) {
            assays[[assay]]@meta.features <- raw_var_df
        } else {
            assays[[assay]]@meta.features <- var_df
        }

        project_name <- sub('\\.h5ad$', '', basename(inFile))
        srt <- new('Seurat', assays = assays, project.name = project_name, version = packageVersion('Seurat'))
        Seurat::DefaultAssay(srt) <- assay
        Seurat::Idents(srt) <- project_name

        srt@meta.data <- obs_df
        embed_names <- unlist(reticulate::py_to_r(ad$obsm$keys()))
        if (length(embed_names) > 0) {
            embeddings <- sapply(embed_names, function(x) reticulate::py_to_r(ad$obsm[x]), simplify = FALSE, USE.NAMES = TRUE)
            names(embeddings) <- embed_names
            for (name in embed_names) {
                rownames(embeddings[[name]]) <- colnames(assays[[assay]])
            }

            dim.reducs <- vector(mode = 'list', length = length(embeddings))
            for (i in seq(length(embeddings))) {
                name <- embed_names[i]
                embed <- embeddings[[name]]
                key <- switch(
                    name,
                    sub('_(.*)', '\\L\\1', sub('^X_', '', toupper(name)), perl=T),
                    'X_pca' = 'PC', 'X_tsne' = 'tSNE', 'X_umap' = 'UMAP'
                )
                colnames(embed) <- paste0(key, '_', seq(ncol(embed)))
                dim.reducs[[i]] <- Seurat::CreateDimReducObject(
                    embeddings = embed,
                    loadings = new('matrix'),
                    assay = assay,
                    stdev = numeric(0L),
                    key = paste0(key, '_')
                )
            }
            names(dim.reducs) <- sub('X_', '', embed_names)

            for (name in names(dim.reducs)) {
                srt[[name]] <- dim.reducs[[name]]
            }
        }
    }

    if (!is.null(outFile)) saveRDS(object = srt, file = outFile)

    srt
}

anndata2cds <- function(inFile, outFile = NULL) {
    message('not implemented yet')
}
