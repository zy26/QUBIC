#' Show report of biclusters
#'
#' This function can make a report for biclusters.
#' @param matrix microarray matrix
#' @param bic array of biclusters
#'
#' @return Text of report
#' @examples
#' # Load microarray matrix
#' if (requireNamespace("biclust", quietly = TRUE) && methods::isClass("BCQU")) {
#'   data(BicatYeast)
#'   matrix <- BicatYeast[1:50, ];
#'   res1 <- biclust::biclust(matrix, method=BCQU(), verbose = FALSE)
#'   res2 <- biclust::biclust(matrix, method=BCCC())
#'   res3 <- biclust::biclust(matrix, method=BCBimax())
#'   # Show the report
#'   showinfo(matrix, c(res1, res2, res3))
#' }
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
  for (i in seq_along(headers)) {
    writeLines(paste0(i, ": ", headers[[i]]))
  }
  writeLines("")
  writeLines("")
  writeLines(paste(seq_along(headers), collapse = "\t"))
  matrix_ratio <- nrow(matrix) / ncol(matrix)

  for (result in bic) {
    out_fields <- rep("", length(headers))
    out_fields[[1]] <- paste(deparse(.qubic_result_call(result)), collapse = " ")

    total_biclusters <- .qubic_result_number(result)
    out_fields[[2]] <- as.character(total_biclusters)

    bic_all <- .qubic_bicluster(matrix, result)
    bic_first <- bic_all[[1]]
    out_fields[[3]] <- as.character(nrow(bic_first))
    out_fields[[4]] <- as.character(ncol(bic_first))
    out_fields[[5]] <- as.character(nrow(bic_first) * ncol(bic_first))
    out_fields[[6]] <- as.character(nrow(bic_first) / ncol(bic_first))
    out_fields[[7]] <- as.character(matrix_ratio)

    nrows <- vapply(bic_all, nrow, integer(1L))
    ncols <- vapply(bic_all, ncol, integer(1L))
    areas <- nrows * ncols
    out_fields[[8]] <- paste(c(max(nrows), which.max(nrows)), collapse = " ")
    out_fields[[9]] <- paste(c(max(ncols), which.max(ncols)), collapse = " ")
    out_fields[[10]] <- paste(c(max(areas), which.max(areas)), collapse = " ")

    rowxnumber <- .qubic_result_rowxnumber(result)
    numberxcol <- .qubic_result_numberxcol(result)
    genes_union <- sum(apply(rowxnumber, 1, max))
    out_fields[[11]] <- paste(c(genes_union, genes_union / nrow(rowxnumber) * 100), collapse = " ")

    conditions_union <- sum(apply(numberxcol, 2, max))
    out_fields[[12]] <- paste(c(conditions_union, conditions_union / ncol(numberxcol) * 100), collapse = " ")

    if (total_biclusters >= 2) {
      rows_overlap <- sum(apply(rowxnumber[, seq_len(2), drop = FALSE], 1, min))
      cols_overlap <- sum(apply(numberxcol[seq_len(2), , drop = FALSE], 2, min))
      out_fields[[13]] <- paste(c(rows_overlap, cols_overlap, rows_overlap * cols_overlap), collapse = " ")
    }

    writeLines(paste(out_fields, collapse = "\t"))
  }
}
