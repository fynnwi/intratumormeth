



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
























