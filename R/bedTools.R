
#' Intersect and Extend regions from multiple bed files
#'
#' intersectAndExtendBed is a function which will intersect a given vector
#' of bed files and return a single Granges object of the merged bed file
#'
#' @param bedFiles A character vector of bed files
#' @param minOverlap Min number of files or fraction of file to have a region
#' @param regionFormat bed/raw. If raw the start will be converted to 0 based
#' @param skipLines Number of lines to skip from the bed file
#'
#' @return A granges oject of the merged regions
#' @export
#'
#' @examples
#' \dontrun{
#' intersectAndExtendBed(bedFiles = c("/path/to/A.bed", "/path/to/B.bed", "/path/to/C.bed"),
#' minOverlap = 2, peakformat = "raw", skipLines = 0)
#' }
#'
intersectAndExtendBed <- function(bedFiles, minOverlap, regionFormat = "bed", skipLines = 0){

  # If no minoverlap found assume we need overlap over all files
  if(missing(minOverlap)){
    minOverlap = length(bedFiles)
  }

  # Create a minimum samplesheet needed by diffbind
  samplesheet.df <- data.frame(
    SampleID = paste0("bed_file_", seq(1, length(bedFiles))),
    Peaks = bedFiles,
    PeakCaller = regionFormat,
    ScoreCol = 0
  )

  #
  diff.dba <- DiffBind::dba(minOverlap = minOverlap,
                            sampleSheet = samplesheet.df,
                            skipLines = skipLines)

  # Extract the expanded bed
  returnGrange <- DiffBind::dba.peakset(DBA = diff.dba,
                                        bRetrieve = TRUE)

  GenomicRanges::elementMetadata(returnGrange) <- NULL

  return(returnGrange)
}
