context("Bed Tools")
library(rizlib)
library(GenomicRanges)


bedFiles = list.files(path = "testData/bed/", pattern = "*.bed", full.names = TRUE)

overlapAllOut <- GenomicRanges::GRanges(seqnames = Rle(c("chr1", "chr2", "chr3")),
                                        ranges = IRanges(start = c(60, 260, 460),
                                                         end = c(250, 450, 650)))

overlap2Out <- GenomicRanges::GRanges(seqnames = Rle(c("chr1", "chr2", "chr3", "chr4")),
                                        ranges = IRanges(start = c(60, 260, 460, 700),
                                                         end = c(250, 450, 650, 850)))

testthat::test_that("test merge and extend bed",{
  testthat::expect_equal(sum(rizlib::intersectAndExtendBed(bedFiles = bedFiles) == overlapAllOut), 3)
  testthat::expect_equal(sum(rizlib::intersectAndExtendBed(bedFiles = bedFiles, minOverlap = 2) == overlap2Out), 4)
  })
