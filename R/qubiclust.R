#' @describeIn QUBIC Performs a QUalitative BIClustering for a discret matrix.
#'
#' @usage qubiclust_d(x, c = 0.95, o = 100, f = 1,
#'         k = max(ncol(x) \%/\% 20, 2),
#'         type = 'default', P = FALSE, C = FALSE, verbose = TRUE,
#'         weight = NULL, seedbicluster = NULL)
qubiclust_d <- function(x, c = 0.95, o = 100, f = 1,
                         k = max(ncol(x)%/%20, 2), type = "default",
                         P = FALSE, C = FALSE, verbose = TRUE, weight = NULL, seedbicluster = NULL) {
  MYCALL <- match.call()
  params <- as.list(MYCALL)
  params$Call <- MYCALL
  S <- (type == "area")
  if(!is.null(weight)) {
    w <- Matrix::Matrix(data = 0, nrow = nrow(x), ncol = nrow(x), dimnames = list(rownames(x), rownames(x)))
    weight[] <- rank(weight, ties.method = "average")
    intersect_rows <- intersect(rownames(x), rownames(weight))
    w[intersect_rows, intersect_rows] <- weight[intersect_rows, intersect_rows]
    res <- .qubic_dw(x, c, o, f, k, P, S, C, verbose, w)
  } else if(!is.null(seedbicluster)) {
    if (.qubic_result_number(seedbicluster) >= 1)
      res <- .qubic_de(x, c, verbose, .qubic_result_rowxnumber(seedbicluster), .qubic_result_numberxcol(seedbicluster))
    else
      return(seedbicluster)
  } else res <- .qubic_d(x, c, o, f, k, P, S, C, verbose)
  return(.qubic_new_BiclustResult(params, matrix(unlist(res["RowxNumber"]), ncol = as.numeric(res["Number"]), byrow = FALSE),
                                  matrix(unlist(res["NumberxCol"]), nrow = as.numeric(res["Number"]), byrow = FALSE), as.numeric(res["Number"]),
                                  res["info"]))
}

#' @describeIn QUBIC Performs a QUalitative BIClustering.
#'
#' @usage qubiclust(x, r = 1L, q = 0.06, c = 0.95, o = 100, f = 1,
#'         k = max(ncol(x) \%/\% 20, 2),
#'         type = 'default', P = FALSE, C = FALSE, verbose = TRUE,
#'         weight = NULL, seedbicluster = NULL)
qubiclust <- function(x, r = 1L, q = 0.06, c = 0.95, o = 100, f = 1,
                       k = max(ncol(x)%/%20, 2), type = "default",
                       P = FALSE, C = FALSE, verbose = TRUE, weight = NULL, seedbicluster = NULL) {
  x_d <- qudiscretize(x, r, q)
  res <- qubiclust_d(x_d, c, o, f, k, type, P, C, verbose, weight, seedbicluster)
  .qubic_set_result_call(res, match.call())
}
