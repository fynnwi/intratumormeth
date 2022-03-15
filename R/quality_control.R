#' Analyse Detection p-Values
#'
#'
#' @param outputDir TODO
#' @param pThreshold TODO
#'
#' @return List
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_hline geom_col theme element_blank geom_histogram labs
#' @importFrom rlang .data
detection_p_analysis <- function(outputDir, pThreshold = 0.01) {
  fname <- file.path(outputDir, "methylation_preprocessed", "detectionpvals.Rdata")
  if (!file.exists(fname)) {
    message(paste0("File ", fname, " not found! Run 'preprocess_methylation_data() to create this file"))
    return(NULL)
  }
  detP <- readRDS(fname)

  # calculate sample-wise mean detection pvals
  meanDetP <- tibble::tibble("sample_id" = colnames(detP),
                             "mean_detp" = colMeans(detP))
  poorSamples <- meanDetP[["sample_id"]][meanDetP[["mean_detp"]] >= pThreshold]

  message(paste(length(poorSamples), "samples have an average detP >", pThreshold))

  exceedsThreshold <- detP > pThreshold
  n <- sum(exceedsThreshold)
  nProbes <- sum(rowSums(exceedsThreshold) > 0)
  message(paste("There are", n , "beta values in the dataset with detP >", pThreshold, "occurring in", nProbes, "different probes"))

  failingProbes <- tibble::tibble("probe_id" = rownames(exceedsThreshold),
                                  "failed_in_n_samples" = rowSums(exceedsThreshold)) %>%
    dplyr::filter(.data$failed_in_n_samples != 0)


  # make a poor sample plot
  poorSamplePlot <- ggplot(meanDetP, aes(x = .data$sample_id, y = .data$mean_detp, fill = .data$mean_detp < pThreshold)) +
    geom_hline(yintercept = pThreshold, linetype = "dashed") +
    geom_col() +
    pdgfra_theme() +
    theme(axis.text.x = element_blank(),
          legend.position = c(0,1),
          legend.justification = c("left", "top")) +
    labs(title = "Sample-wise average detection p-values",
         x = "Samples", y = "Detection p-value, averaged across all probes",
         fill = paste0("Mean det. p-value < ", pThreshold))

  # make failing probes plot
  failingProbesPlot <- ggplot(failingProbes, aes(x = .data$failed_in_n_samples)) +
    geom_histogram(binwidth = 1) +
    pdgfra_theme() +
    labs(title = paste0("Probes exceeding det. p-value (", pThreshold, ") in multiple samples"),
         x = "Number of samples for which a given probe fails",
         y = "Number of probes")

  return(list(meanDetP = meanDetP,
              poorSamplePlot = poorSamplePlot,
              failingProbes = failingProbes,
              failingProbesPlot = failingProbesPlot))
}




#' Predict Sex Based on Methylation Data
#'
#' TODO describe whats happening here (difference between X and Y chr
#' intensities as indicator for male/female samples)
#'
#' @param outputDir TODO
#' @param normalization TODO
#' @param predictionCutoff What should the difference in log2 copynumber be
#'   between males and females. See \code{?minfi::getSet} for details. Default
#'   is -2.
#'
#' @return TODO
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_density labs geom_jitter geom_vline
#' @importFrom rlang .data
predict_sex <- function(outputDir, normalization = "raw", predictionCutoff = -2) {
  # load methylation data
  fname <- file.path(outputDir, "methylation_preprocessed", paste0("mset_", normalization, ".Rdata"))
  if (!file.exists(fname)) {
    message(paste0("File ", fname, " not found! Make sure a correct normalization method was specified or run 'preprocess_methylation_data() to create this file"))
    return(NULL)
  }
  mSet <- readRDS(fname)
  # map to genome to be able to extract X and Y chromosome probes
  gmSet <- minfi::mapToGenome(mSet) %>%
    suppressMessages()
  # free memory
  remove(mSet)
  gc()
  predSex <- minfi::getSex(gmSet, cutoff = predictionCutoff)
  # format
  predSex <- predSex %>%
    tibble::as_tibble(rownames = "sample_id") %>%
    dplyr::rename("predicted_sex" = "predictedSex") %>%
    dplyr::mutate("predicted_sex" = tolower(.data$predicted_sex))

  # extract the beta values of the X and Y chromosomes
  granges <- SummarizedExperiment::rowRanges(gmSet)
  chrXProbes <- names(granges[GenomicRanges::seqnames(granges) == "chrX"])
  chrYProbes <- names(granges[GenomicRanges::seqnames(granges) == "chrY"])

  chrXbetas <- minfi::getBeta(gmSet[chrXProbes]) %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "sample_id", values_to = "beta") %>%
    tidyr::drop_na(beta) %>%
    dplyr::left_join(predSex, by = "sample_id")
  chrYbetas <- minfi::getBeta(gmSet[chrYProbes]) %>%
    tibble::as_tibble() %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "sample_id", values_to = "beta") %>%
    tidyr::drop_na(beta) %>%
    dplyr::left_join(predSex, by = "sample_id")
  # free memory
  remove(gmSet)
  gc()

  # plot
  pX <- ggplot(chrXbetas, aes(x = .data$beta, color = .data$predicted_sex, group = .data$sample_id)) +
    geom_density() +
    labs(title = "Chromosome X beta value densities",
         x = "Beta value",
         y = "Density",
         color = "Predicted Sex") +
    pdgfra_theme()
  pY <- ggplot(chrYbetas, aes(x = .data$beta, color = .data$predicted_sex, group = .data$sample_id)) +
    geom_density() +
    labs(title = "Chromosome Y beta value densities",
         x = "Beta value",
         y = "Density",
         color = "Predicted Sex") +
    pdgfra_theme()

  # cutoff plot
  predSex["xy_gap"] <- predSex[["yMed"]] - predSex[["xMed"]]

  pCutoff <- ggplot(predSex, aes(x = .data$xy_gap, y = "", color = .data$predicted_sex)) +
    geom_jitter(height = 0.05, shape = 18, size = 2) +
    geom_vline(xintercept = predictionCutoff, linetype = "dashed") +
    labs(title = "Cutoff value used for sex prediction",
         subtitle = paste("Samples with signal intensity distance less than",
                          predictionCutoff, "between chrY and chrX are classified as female"),
         x = "chrY median intensity - chrX median intensity",
         y = "", color = "Predicted sex") +
    pdgfra_theme()

  plots <- list()
  plots[["chrX"]] <- pX
  plots[["chrY"]] <- pY
  plots[["cutoff"]] <- pCutoff
  plots[["predicted_sex"]] <- predSex %>%
    dplyr::select(-.data$xMed, -.data$yMed, -.data$xy_gap)

  return(plots)
}
