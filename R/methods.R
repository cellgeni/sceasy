convertFormat <- function(
    obj, from = c('anndata', 'seurat', 'sce', 'loom'), to = c('anndata', 'loom', 'sce', 'seurat'), outFile = NULL,
    main_layer = NULL, ...
) {
    from <- match.arg(from)
    to <- match.arg(to)

    tryCatch({
        func <- eval(parse(text = paste(from, to, sep='2')))
    }, error = function(e) {
        stop(paste0('Unsupported conversion from "', from, '" to "', to, '"'), call. = FALSE)
    }, finally = {})

    return(func(obj, outFile = outFile, main_layer = main_layer, ...))
}
