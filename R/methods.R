#' Convert between data objects
#'
#' This function converts between data format frequently used to hold single
#' cell data
#'
#' @param obj Input Seurat object
#' @param from Format of input object, e.g. "anndata", "seurat", "sce", "loom",
#'   etc (str)
#' @param to Format of output object (str)
#' @param main_layer Required by some formats, may be "counts", "data",
#'   "scale.data", etc (str)
#'
#' @return Output object
convertFormat <- function(obj, from = c("anndata", "seurat", "sce", "loom"), to = c("anndata", "loom", "sce", "seurat"), outFile = NULL,
                          main_layer = NULL, ...) {
  from <- match.arg(from)
  to <- match.arg(to)

  tryCatch(
    {
      func <- eval(parse(text = paste(from, to, sep = "2")))
    },
    error = function(e) {
      stop(paste0('Unsupported conversion from "', from, '" to "', to, '"'), call. = FALSE)
    },
    finally = {}
  )

  return(func(obj, outFile = outFile, main_layer = main_layer, ...))
}
