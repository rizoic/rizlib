
#' Intersect and Extend regions from multiple bed files
#'
#' intersectAndExtendBed is a function which will intersect a given vector
#' of bed files and return a single Granges object of the merged bed file
#'
#' @param bedFiles A character vector of bed files
#' @param minOverlap Min number of files or fraction of file to have a region
#'
#' @return A granges oject of the merged regions
#' @export
#'
#' @examples
#' \dontrun{
#' intersectAndExtendBed(bedFiles = c("/path/to/A.bed", "/path/to/B.bed", "/path/to/C.bed"),
#' minOverlap = 2)
#' }
#'
intersectAndExtendBed <- function(bedFiles, minOverlap){

  # If no minoverlap found assume we need overlap over all files
  if(missing(minOverlap)){
    minOverlap = length(bedFiles)
  }

  # Create a minimum samplesheet needed by diffbind
  samplesheet.df <- data.frame(
    SampleID = paste0("bed_file_", seq(1, length(bedFiles))),
    Peaks = bedFiles,
    ScoreCol = 0
  )

  # Create a diffbind obj
  suppressMessages(diff.dba <- DiffBind::dba(minOverlap = as.numeric(minOverlap),
                            sampleSheet = samplesheet.df))

  # Extract the expanded bed
  returnGrange <- DiffBind::dba.peakset(DBA = diff.dba,
                                        bRetrieve = TRUE)

  # Remove metadata that diffbind adds
  GenomicRanges::elementMetadata(returnGrange) <- NULL

  # Return a stripped down granges
  return(returnGrange)
}
