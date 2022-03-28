



#' Run Conumee Workflow
#'
#' @param outputDir TODO
#' @param samplesheet TODO
#' @param patientsheet TODO
#' @param refF TODO
#' @param refM TODO
#' @param annoF TODO
#' @param annoM TODO
#'
#' @return TODO
#' @export
#'
run_conumee <- function(outputDir, samplesheet, patientsheet, refF, refM, annoF, annoM) {
  # load methylation data
  mSet <- readRDS(file.path(outputDir, "methylation_preprocessed", "mset_raw.Rdata"))
  # create output directory
  outputDir <- file.path(outputDir, "conumee")
  dir.create(outputDir)

  for (sample in samplesheet[["sample_id"]]) {
    message(paste("Processing sample", sample))
    if (!sample %in% colnames(mSet)) {
      message(paste("Sample", sample, "not found in methylation data; skipping"))
      next
    }

    patId <- samplesheet[["patient_id"]][samplesheet[["sample_id"]] == sample]
    sex <- patientsheet[["sex"]][patientsheet[["patient_id"]] == patId]

    # sanity check:
    if (sex == "m") {
      anno <- annoM
      ref <- refM
    } else if (sex == "f") {
      anno <- annoF
      ref <- refF
    } else {
      stop(paste("Patient sex", sex, "not supported; must be 'm' or 'f'"))
    }

    conumeeInput <- mSet[ , sample]

    # load combined intensities values
    cndata <- conumee::CNV.load(conumeeInput)

    # assign sample id as name
    names(cndata) <- sample

    # filter annotation probes (the annotation might contain a few probes that are
    # not present on all EPIC chips)

    # find all probes that are covered in (1) annotation object, (2) query
    # intensities, and (3) reference intensities
    commonProbes <- intersect(names(anno@probes),
                              intersect(rownames(cndata@intensity),
                                        rownames(ref@intensity)))
    # apparently it only matters that anno contains no probe ids that are not
    # present in query or reference intensities...
    keep <- names(anno@probes) %in% commonProbes
    # message(paste("Removing", (length(keep) - sum(keep)), "probes that are not present in both query, and reference data")) # TODO make less verbose?
    anno@probes <- anno@probes[keep]

    # keep <- rownames(cndata@intensity) %in% commonProbes
    # cndata@intensity <- cndata@intensity[keep, ]

    # conumee pipeline
    x <- conumee::CNV.fit(cndata, ref, anno)
    x <- conumee::CNV.bin(x)
    x <- conumee::CNV.detail(x)
    x <- conumee::CNV.segment(x)

    # save as Rdata to hard drive
    fname <- file.path(outputDir, paste0("conumee_", sample, ".Rdata"))
    saveRDS(x, fname)
  }
}









#' Conumee CNV Genomeplots
#'
#' Render a PDF containing conumee CNV genomeplots for every sample organized by
#' patient.
#'
#' @param outputDir Directory that contains all intratumormeth output.
#' @param samplesheet Relates sample id to patient id and more.
#' @param detailRegions Character vector containing names of those detail
#'   regions that should appear as labels on the CNV plot.
#'
#' @export
#'
cnv_genomeplots_per_patient <- function(outputDir, samplesheet, detailRegions = NULL) {

  # look for conume output files
  conumeeFilesList <- list.files(file.path(outputDir, "conumee"), full.names = TRUE)
  availableFiles <- stringr::str_remove(stringr::str_remove(basename(conumeeFilesList), "\\.Rdata$"), "^conumee_")

  # how many cnv plots should be put on one page?
  plotsPerPage <- 5
  dir.create(file.path(outputDir, "copynumber"))
  plotFilepath <- file.path(outputDir, "copynumber", paste0("cnv_genomeplots_", length(availableFiles), "_samples.pdf"))
  message(paste("Writing", length(conumeeFilesList), "genomeplots to", plotFilepath))
  grDevices::pdf(plotFilepath, width = 8.3, height = 11.7)

  for (pat in unique(samplesheet$patient_id)) {

    message(pat)

    graphics::par(mfrow = c(plotsPerPage, 1))

    for (currentSample in samplesheet$sample_id[samplesheet$patient_id == pat]) {
      # get corresponding CNV analysis object

      # check if there exists a corresponding file
      if (sum(availableFiles == currentSample) == 0) {
        message(paste("No conumee file found for sample", currentSample, "- skipping"))
        next
      } else if (sum(availableFiles == currentSample) > 1) {
        message(paste("More than one conumee file found for sample", currentSample, "- skipping"))
        next
      }

      conumeeFile <- conumeeFilesList[which(availableFiles == currentSample)]
      cnv <- readRDS(conumeeFile)

      # The CNV object may contain too many detail regions to be plotted in a meaningful
      # way. Therefore, one can specify the detailRegions argument which filters those

      if (!is.null(detailRegions)) {
        detailProbes <- cnv@detail
        detailAnno <- cnv@anno@detail

        detailProbes$ratio <- detailProbes$ratio[names(detailProbes$ratio) %in% detailRegions]
        detailProbes$probes <- detailProbes$probes[names(detailProbes$probes) %in% detailRegions]

        detailAnno <- detailAnno[detailAnno$name %in% detailRegions, ]

        cnv@detail <- detailProbes
        cnv@anno@detail <- detailAnno
      }

      message(paste0("Creating plots for sample", currentSample, "..."))
      conumee::CNV.genomeplot(cnv, set_par = FALSE, detail = TRUE)
    }
  }
  grDevices::dev.off()
}





#' Postprocess Conumee Files
#'
#' Extract the relevant information from conumee CNV.analysis objects.
#'
#' @param outputDir Directory where output of intratumormeth is stored.
#'
#' @export
postprocess_conumee_files <- function(outputDir) {
  # look for conume output files
  conumeeFilesList <- list.files(file.path(outputDir, "conumee"), full.names = TRUE)
  availableFiles <- stringr::str_remove(stringr::str_remove(basename(conumeeFilesList), "\\.Rdata$"), "^conumee_")
  message(paste0("Postprocessing ", length(availableFiles), " conumee files"))

  segmentData <- NULL
  detailData <- NULL
  # iterate over all conumee files and extract segment and detail information
  pb <- progress::progress_bar$new(total = length(availableFiles))
  pb$tick(0)
  for (i in 1:length(availableFiles)) {
    pb$tick()
    # message(availableFiles[i])

    cnv <- readRDS(conumeeFilesList[i])
    seg <- extract_segments(cnv)
    detail <- extract_details(cnv)

    # store data
    if (is.null(segmentData)) {
      segmentData <- seg
    } else {
      segmentData <- rbind(segmentData, seg)
    }
    if (is.null(detailData)) {
      detailData <- detail
    } else {
      detailData <- rbind(detailData, detail)
    }
  }
  detailData <- detailData %>%
    tibble::as_tibble()
  segmentData <- segmentData %>%
    tibble::as_tibble()

  # write results to disk
  # TODO check if dir exists?
  # dir.create(file.path(outputDir, "copynumber"), showWarnings = TRUE)
  fname <- file.path(outputDir, "copynumber", "cnv_detailregions.Rdata")
  message(paste("Saving file", fname))
  saveRDS(detailData, fname)
  fname <- file.path(outputDir, "copynumber", "cnv_segments.Rdata")
  message(paste("Saving file", fname))
  saveRDS(segmentData, fname)
}





#' Extract Detail Region CNV
#'
#' From conumee CNV.analysis object.
#'
#' @param cnv TODO
#'
#' @return TODO
#' @importFrom rlang .data
extract_details <- function(cnv) {
  currentSample <- names(cnv)
  # extract information about chr, start, and end
  anno <- cnv@anno@detail %>%
    as.data.frame() %>%
    dplyr::select(.data$seqnames, .data$start, .data$end, .data$name) %>%
    `colnames<-`(c("chr", "start", "end", "region")) %>%
    dplyr::mutate("chr" = chr_to_numeric(.data$chr))
  detail <- data.frame("sample_id" = currentSample,
                       "region" = names(cnv@detail$ratio),
                       "log2ratio_median" = cnv@detail$ratio) %>%
    `rownames<-`(c()) %>%
    dplyr::left_join(anno, by = "region")
  return(detail)
}


#' Extract Segment CNV
#'
#' From conumee CNV.analysis object.
#'
#' @param cnv TODO
#'
#' @return TODO
#' @importFrom rlang .data
extract_segments <- function(cnv) {
  currentSample <- names(cnv)

  seg <- cnv@seg$summary
  seg$sample_id <- currentSample
  seg$chr <- chr_to_numeric(seg$chrom)

  seg <- seg %>%
    dplyr::select(.data$sample_id, .data$chr, .data$loc.start, .data$loc.end, -.data$num.mark, .data$seg.mean, .data$seg.sd, .data$seg.median, .data$seg.mad) %>%
    dplyr::arrange(.data$chr, .data$loc.start) %>%
    `colnames<-`(c("sample_id", "chr", "start", "end", "log2ratio_mean", "log2ratio_sd", "log2ratio_median", "log2ratio_mad"))

  return(seg)
}









#' Convert Chromosome Name to Numeric
#'
#' Converts character string of format e.g. chr10 into numeric 10, or chrX into
#' 23
#'
#' @param chrString Character string of convention TODO
#'
#' @return TODO
chr_to_numeric <- function(chrString) {
  # TODO implement style argument
  chr <- substring(chrString, 4)
  chr <- sub("X", "23", chr, fixed = TRUE)
  chr <- sub("Y", "24", chr, fixed = TRUE)
  return(as.numeric(chr))
}


#' Convert numeric to a string
#'
#' E.g. 8 -> chr8
#'
#' @param chrNum TODO
#'
#' @return TODO
#'
chr_to_string <- function(chrNum) {
  # stopifnot(length(chrNum == 1))
  if (chrNum == 24) {
    s <- "Y"
  } else if (chrNum == 23) {
    s <- "X"
  } else if (chrNum >= 1 && chrNum <= 22) {
    s <- chrNum
  } else {
    stop(paste("Unknown chromosome number:", chrNum))
  }

  return(paste0("chr", s))
}















