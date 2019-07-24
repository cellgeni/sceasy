convertFormat <- function(
    obj, from = c('seurat', 'sce', 'loom'), to = c('anndata'), outFile = NULL,
    main_layer = NULL, ...
) {
    from <- match.arg(from)
    to <- match.arg(to)

    func <- eval(parse(text=paste(from, to, sep='2')))

    return(func(obj, outFile = outFile, main_layer = main_layer, ...))
}
