#' @title .CalculateMALAT1Seurat
#' Computes an automatic threshold on normalized MALAT1 expression to separate low-quality
#' droplets or empty droplets from real cells, following the MALAT1_threshold approach
#' developed by BaderLab (see references).
#'
#' @description
#' Computes an automatic threshold on normalized MALAT1 expression to separate low-quality
#' droplets or empty droplets from real cells, following the MALAT1_threshold approach
#' developed by BaderLab (see references). This function is designed to work with Seurat
#' objects and assumes that the input Seurat object has been normalized and contains a
#' "MALAT1" feature in its data layers.
#'
#' @details
#' MALAT1 is a nuclear-retained lncRNA whose expression correlates strongly with
#' nuclear (intronic) RNA content and splice ratio, making it a convenient proxy
#' for cell quality in droplet-based scRNA-seq.
#' Low MALAT1 values typically correspond to empty droplets, cytosolic debris,
#' or damaged cells, whereas higher values correspond to intact nuclei.
#' This function assumes you provide a SeuratObject with normalized counts
#' and a "MALAT1" feature in its data layers. It computes a threshold to
#' flag low-quality cells based on the distribution of MALAT1 expression.
#' The computation is performed per provided layer (e.g., "data.pool1", "data.pool2", etc.).
#' The threshold deterined per layer and a boolean value indicating whether each cell passes
#' the threshold are stored in the Seurat object's metadata.
#'
#' @param SeuObj A SeuratObject with normalized counts and a "MALAT1" feature in its data layers.
#' @param assay Character string specifying the assay in the SeuratObject to use for MALAT1 expression. Default is "RNA".
#' @param layers Character vector specifying the layers in the SeuratObject to process.
#'   Default is `NULL`, which processes all layers containing "data" in their names.
#' @param ... Additional arguments passed to the underlying `ComputeMALAT1Threshold` function, such as:
#'  \itemize{
#'    \item \code{bw.bandwidth}: Bandwidth for kernel density estimation (default: 0.01).
#'    \item \code{chosen.min}: Chosen minimum which a peak should be considered the dataset peak (default: 2).
#'    \item \code{smooth.spar}: Smoothing parameter for density estimation (default: 2).
#'    \item \code{abs.min}: Absolute minimum threshold (default: 1).
#'    \item \code{rough.max}: Rough expected position of the MALAT1 expression peak (default: 6).
#' }
#'
#' @import Seurat
#' @import SeuratObject
#' @importFrom data.table rbindlist
#' 
#'
#' @return:
#' \describe{
#'   \item{SeuratObject} {The input SeuratObject with the MALAT1 threshold applied.}
#' }
#'
#'@noRd
#' 
#' @references
#'
#' Clarke, Bader et al. "MALAT1 expression indicates cell quality in
#' single-cell RNA sequencing data." bioRxiv (2024).
#' \url{https://www.biorxiv.org/content/10.1101/2024.07.14.603469v2}
#'
#' @examples
#' \dontrun{
#'   threshold.res <- .CalculateMALAT1Seurat(
#'     SeuratObject = seurat_obj,
#'     Assay = "RNA",
#'     Layers = c("data.pool1", "data.pool2"),
#'     bw.bandwidth = 0.01,
#'     chosen.min = 2,
#'     smooth.spar = 2,
#'     abs.min = 1,
#'     rough.max = 6
#'   )
#' }
#'
#' @seealso
#' \itemize{
#'   \item DropletQC: nuclear fraction-based QC for empty droplets and damaged cells.
#'   \item EmptyDrops: ambient RNA-based empty droplet detection.
#' }

.CalculateMALAT1Seurat <- function(
  SeuratObject,
  assay = "RNA",
  layers = NULL,
  ...
) {
  # Check if the specified assay exists in the Seurat object
  if (!assay %in% names(SeuratObject@assays)) {
    stop(paste("Assay", assay, "not found in the Seurat object."))
  }

  # Fetch all layers present in the Seurat object for the specified assay
  object.layers <- SeuratObject::Layers(
    SeuObj,
    assay = assay,
  )

  # Check if the user provided layers are present in the Seurat object.
  if (!is.null(layers)) {
    if ((!all(layers %in% object.layers))) {
      stop(
        "Some specified layers are not present in the Seurat object. Please check the provided layers argument."
      )
    }
  } else {
    # If no layers are provided, check if any log-normalized data layers are present in the Seurat object.
    if (!any(grepl("data", object.layers))) {
      stop(
        "User did not specified any layers and no log-normalized data layers found in the Seurat object. Please specify layers to process."
      )
    }
    # If no layers are provided, fetch all log-normalized 'data' layers names from the Seurat object.
    layers <- object.layers[grepl("data", object.layers)]
  }

  # The results of the MALAT1 threshold computation will be stored in the Seurat object's metadata.
  MALAT1.metadata <- lapply(
    layers,
    function(layer) {
      # Retrieve the log-normalized MALAT1 expression values for the current layer.
      MALAT1.logcounts <- SeuratObject::FetchData(
        SeuObj,
        assay = assay,
        vars = "MALAT1",
        layer = layer
      )

      # Compute the MALAT1 threshold using the ComputeMALAT1Threshold function, passing any additional arguments.
      threshold.value <- round(
        .ComputeMALAT1Threshold(MALAT1.logcounts$MALAT1, ...),
        digits = 2
      )

      # We will store the results in a data frame with cell IDs, computed threshold,
      # and whether each cell passes the threshold.
      # The cells id are the row names of the MALAT1.logcounts data frame to ensure
      # proper alignment with the Seurat object's metadata.
      cell.id <- row.names(MALAT1.logcounts)
      # Create a vector of the computed threshold value for each cell.
      MALAT1.threshold <- rep(threshold.value, length(cell.id))
      # Create a boolean vector indicating whether each cell's
      # MALAT1 expression exceeds the computed threshold.
      MALAT1.pass <- MALAT1.logcounts$MALAT1 > MALAT1.threshold
      # Combine the cell IDs, threshold values, and pass/fail results into a data frame.
      results.df <- data.frame(
        cell.id = cell.id,
        MALAT1.threshold = MALAT1.threshold,
        MALAT1.pass = MALAT1.pass
      )
      return(results.df)
    }
  ) |> # Combine the results from all layers into a single data frame.
    data.table::rbindlist() |>
    as.data.frame()

  # Set the row names of the MALAT1.metadata data frame to the cell IDs
  # for proper alignment with the Seurat object's metadata.
  row.names(MALAT1.metadata) <- MALAT1.metadata$cell.id
  # Delete the cell.id column from MALAT1.metadata as it is now redundant with the row names.
  MALAT1.metadata$cell.id <- NULL

  # Add the computed MALAT1 threshold and pass/fail results to the Seurat object's metadata.
  SeuObj <- Seurat::AddMetaData(
    object = SeuObj,
    metadata = MALAT1.metadata
  )

  return(SeuObj)
}
