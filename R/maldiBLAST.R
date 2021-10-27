#'  Search one IDBac database against another IDBac database
#'
#' @param queryPool pool of database to search
#' @param subjectPool pool of database of knowns to search againstr
#' @param searchppm ppm threshold for search
#' @param searchThreshold minimum cosine score of returned matches
#' @importFrom IDBacApp idbac_get_peaks createFuzzyVector
#' @importFrom coop cosine
#' @importFrom data.table as.data.table
#' @inheritParams IDBacApp::idbac_get_peaks
#' @return
#' @export
#'
#' @examples
maldiBLAST <- function(queryPool,
                       subjectPool,
                       minFrequency = 0L,
                       minNumber = 1L,
                       lowerMassCutoff = 4000L,
                       upperMassCutoff = 15000L,
                       minSNR = 6,
                       tolerance = 2,
                       type = "protein",
                       mergeReplicates = TRUE,
                       method = "strict",
                       searchppm = 3000,
                       searchThreshold = 0.7){

  # Retrieve peak data
  subject_peak_data <- IDBacApp::idbac_get_peaks(pool = subjectPool,
                                                 sampleIDs = NULL,
                                                 minFrequency = minFrequency,
                                                 minNumber = minNumber,
                                                 lowerMassCutoff = lowerMassCutoff,
                                                 upperMassCutoff = upperMassCutoff,
                                                 minSNR = minSNR,
                                                 tolerance = tolerance,
                                                 type = type,
                                                 mergeReplicates = mergeReplicates,
                                                 method = method)
  # Retrieve peak data
  query_peak_data <- IDBacApp::idbac_get_peaks(pool = queryPool,
                                               sampleIDs = NULL,
                                               minFrequency = minFrequency,
                                               minNumber = minNumber,
                                               lowerMassCutoff = lowerMassCutoff,
                                               upperMassCutoff = upperMassCutoff,
                                               minSNR = minSNR,
                                               tolerance = tolerance,
                                               type = type,
                                               mergeReplicates = mergeReplicates,
                                               method = method)

  # Represent peaks as probability distributions on higher dimensional vector (alternative to binning)
  subject_peak_matrix <- IDBacApp::createFuzzyVector(massStart = lowerMassCutoff,
                                                     massEnd = upperMassCutoff,
                                                     massList = lapply(subject_peak_data, function(x) x@mass),
                                                     intensityList = lapply(subject_peak_data, function(x) x@intensity),
                                                     ppm = searchppm)

  # Represent peaks as probability distributions on higher dimensional vector (alternative to binning)
  query_peak_matrix <- IDBacApp::createFuzzyVector(massStart = lowerMassCutoff,
                                                    massEnd = upperMassCutoff,
                                                    massList = lapply(query_peak_data, function(x) x@mass),
                                                    intensityList = lapply(query_peak_data, function(x) x@intensity),
                                                    ppm = searchppm)

  # Basically turn normal distributions of peak representation to presence/absence
  # This means we're euqally-weighting variable regions around each peak that varies based on ppm
  subject_peak_matrix[subject_peak_matrix > 0] <- 1
  query_peak_matrix[query_peak_matrix > 0] <- 1

  # combined peak lists from query and subject
  query_length <- ncol(query_peak_matrix)
  subject_length <- ncol(subject_peak_matrix)
  combined <- cbind(query_peak_matrix, subject_peak_matrix)

  # Calculate distances
  combined <- coop::cosine(combined)

  # Remove all values except non-duplicated distances between query and subject
  combined[upper.tri(combined, diag = TRUE)] <- NA
  combined[1:ncol(query_peak_matrix), ] <- NA
  combined[ , seq(1 + query_length, ncol(combined), 1)] <- NA
  combined <- reshape2::melt(combined)
  combined <- as.data.table(combined)
  combined <- combined[!is.na(value), ]
  # Sort, best matches first
  combined <- combined[order(value, decreasing = T), ]
  # Remove matches less than provided threshold
  combined <- combined[value > searchThreshold, ]
  combined$Var1 <- as.character(combined$Var1)
  combined$Var2 <- as.character(combined$Var2)
  colnames(combined) <- c("subject", "query", "cosine")

  return(combined)

}


