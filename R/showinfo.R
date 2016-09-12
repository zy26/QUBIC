#' Show report of biclusters
#'
#' This function can make a report for biclusters.
#' @param matrix microarray matrix
#' @param bic array of biclusters
#'
#' @return Text of report
#' @examples
#' # Load microarray matrix
#' data(BicatYeast)
#' matrix <- BicatYeast[1:50, ];
#' res1 <- biclust::biclust(matrix, method=BCQU(), verbose = FALSE)
#' res2 <- biclust::biclust(matrix, method=BCCC())
#' res3 <- biclust::biclust(matrix, method=BCBimax())
#' # Show the report
#' showinfo(matrix, c(res1, res2, res3))
#' @seealso \code{\link{QUBIC}}
showinfo <- function(matrix, bic) {
  headers <- c("Call and Parameter",
  "number of detected biclusters",
  "nrow of the first bicluster",
  "ncol of the first bicluster",
  "area of the first bicluster",
  "ratio (nrow / ncol) of the first bicluster",
  "ratio (nrow / ncol) of the matrix",
  "max nrow and corresponding bicluster",
  "max ncol and corresponding bicluster",
  "max area and corresponding bicluster",
  "union of rows, (# and %)",
  "union of columns, (# and %)",
  "overlap of first two biclusters (row, col, area)")
  for (i in 1:length(headers)) {
    cat(i)
    cat(": ")
    cat(headers[[i]])
    cat("\n")
  }
  cat("\n")
  cat("\n")
  for (i in 1:length(headers)) {
    cat(i)
    cat("\t")
  }
  cat("\n")
  for (biclust in bic) {
    cat(deparse(biclust@Parameters$Call))
    cat("\t")
    cat(biclust@Number)
    cat("\t")
    bic <- biclust::bicluster(matrix, biclust, 1)[[1]]
    bicall <- biclust::bicluster(matrix, biclust)
    cat(nrow(bic))
    cat("\t")
    cat(ncol(bic))
    cat("\t")
    cat(nrow(bic) * ncol(bic))
    cat("\t")
    cat(nrow(bic) / ncol(bic))
    cat("\t")
    cat(nrow(matrix) / ncol(matrix))
    cat("\t")
    maxnrow <- c(-1, -1)
    maxncol <- c(-1, -1)
    maxarea <- c(-1, -1)
    for (i in 1:biclust@Number) {
      nrow <- nrow(bicall[[i]])
      ncol <- ncol(bicall[[i]])
      area <- nrow * ncol
      if (maxnrow[[1]] < nrow) maxnrow = c(nrow, i)
      if (maxncol[[1]] < ncol) maxncol = c(ncol, i)
      if (maxarea[[1]] < area) maxarea = c(area, i)
    }
    cat(maxnrow)
    cat("\t")
    cat(maxncol)
    cat("\t")
    cat(maxarea)
    cat("\t")
    genes_union <- sum(apply(biclust@RowxNumber, 1, max))
    cat(c(genes_union, genes_union / nrow(biclust@RowxNumber) * 100))
    cat("\t")
    conditions_union <- sum(apply(biclust@NumberxCol, 2, max))
    cat(c(conditions_union, conditions_union / ncol(biclust@NumberxCol) * 100))
    cat("\t")
    if(biclust@Number >= 2) {
      rows_overlap <- sum(apply(biclust@RowxNumber[,1:2], 1, min))
      cols_overlap <- sum(apply(biclust@NumberxCol[1:2,], 2, min))
      cat(c(rows_overlap, cols_overlap, rows_overlap * cols_overlap))
    }
    cat("\n")
  }
}
