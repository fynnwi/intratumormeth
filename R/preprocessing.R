#' Preprocess Methylation Data
#'
#' Under the hood this function mainly calls methods from the minfi package. It
#' applies different preprocessing methods and saves the resulting methylation
#' data to the disk under the specified directory.
#'
#' @param samplesheet Dataframe containing columns 'sample_id' and 'Basename'.
#' @param outputDir Directory to save results to.
#' @param minfiForce Argument 'force' passed to minfi::read.metharray() call
#'
#' @export
#'
preprocess_methylation_data <- function(samplesheet, outputDir, minfiForce = FALSE) {
  # checks:
  stopifnot(dir.exists(outputDir))
  stopifnot(length(unique(samplesheet[["sample_id"]])) == nrow(samplesheet))

  # rename idat_location column into "Basename"
  samplesheet <- samplesheet %>%
    dplyr::rename("Basename" = "idat_location") %>%
    as.data.frame()

  # create new subdirectory in outputDir
  dir.create(file.path(outputDir, "methylation_preprocessed"))
  outputDir <- file.path(outputDir, "methylation_preprocessed")

  message(paste0("ITM 1/5: Reading idat files for ", nrow(samplesheet), " samples"))
  rgSet <- minfi::read.metharray.exp(targets = samplesheet, force = TRUE) %>%
    `colnames<-`(samplesheet[["sample_id"]])
  rgSetFile <- file.path(outputDir, "rgset.Rdata")
  message(paste0("Saving ", rgSetFile))
  saveRDS(rgSet, file = rgSetFile)

  message("ITM 2/5: Calculating detection p-values")
  detP <- minfi::detectionP(rgSet)
  detPFile <- file.path(outputDir, "detectionpvals.Rdata")
  message(paste0("Saving ", detPFile))
  saveRDS(detP, file = detPFile)

  message("ITM 3/5: Preprocessing raw")
  mSetRaw <- minfi::preprocessRaw(rgSet)
  mSetRawFile <- file.path(outputDir, "mset_raw.Rdata")
  message(paste0("Saving ", mSetRawFile))
  saveRDS(mSetRaw, file = mSetRawFile)

  message("ITM 4/5: Applying normalization (SWAN)")
  mSetSwan <- minfi::preprocessSWAN(rgSet)
  mSetSwanFile <- file.path(outputDir, "mset_swan.Rdata")
  message(paste0("Saving ", mSetSwanFile))
  saveRDS(mSetSwan, file = mSetSwanFile)

  message("ITM 5/5: Applying normalization (Noob)")
  mSetNoob <- minfi::preprocessNoob(rgSet)
  mSetNoobFile <- file.path(outputDir, "mset_noob.Rdata")
  message(paste0("Saving ", mSetNoobFile))
  saveRDS(mSetNoob, file = mSetNoobFile)
}






#' Get Methylation Beta Values
#'
#' Loads the previously generated methylation beta values from disk and retuns
#' them in a matrix. This function allows to specify probe and sample filters.
#' These should be specified as character vectors containing probe/sample ids.
#' Furthermore, the normalization method can be specified. Probes that contain
#' NA values can be removed from all samples by specifying `removeNaProbes =
#' TRUE`.
#'
#' @param outputDir Output directory of the intratumormeth study.
#' @param probeFilter Character vector containing probe ids that should be
#'   filtered out.
#' @param sampleFilter Character vector containing sample ids that should be
#'   filtered out.
#' @param normalization Should be either "raw", "noob", or "swan".
#' @param removeNaProbes If TRUE, beta values of probes that are NA for at least
#'   one sample will be removed.
#'
#' @return Matrix containing methylation beta values with CpG probes as rows and
#'   sample ids as columns.
#' @export
#'
get_betas <- function(outputDir, probeFilter = NULL, sampleFilter = NULL, normalization = "raw", removeNaProbes = FALSE) {
  # find and load minfi object corresponding to specified normalization
  fname <- file.path(outputDir, "methylation_preprocessed", paste0("mset_", normalization, ".Rdata"))
  if (!file.exists(fname)) {
    message(paste0("File ", fname, " not found! Check requested normalization method"))
    return(NULL)
  }
  message(paste0("Loading ", fname, " ..."))
  mSet <- readRDS(fname)

  # extract beta values and filter probes and samples
  betas <- minfi::getBeta(mSet)
  betas <- betas[!rownames(betas) %in% probeFilter, !colnames(betas) %in% sampleFilter]
  if (!is.null(sampleFilter)) {
    message(paste0("Removing ", length(sampleFilter), " samples"))
  }
  if (!is.null(probeFilter)) {
    message(paste0("Removing ", length(probeFilter), " probes"))
  }

  # treat NA values
  if (removeNaProbes) {
    n <- nrow(betas)
    betas <- betas[stats::complete.cases(betas), ]
    message(paste0("Removing ", n-nrow(betas), " probes containing NA for at least one sample"))
  }

  return(betas)
}
