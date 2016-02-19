.qubiclust <- function(x, r = 1, q = 0.06, c = 0.95, o = 100, f = 1, k = max(ncol(x)%/%20,
                                                                             2), type = "default", P = FALSE, C = FALSE, verbose = TRUE) {
  MYCALL <- match.call()
  S <- (type == "area")
  res <- .qubic(x, r, q, c, o, f, k, P, S, C, verbose)
  return(BiclustResult(as.list(MYCALL), matrix(unlist(res["RowxNumber"]),
                                               ncol = as.numeric(res["Number"]), byrow = FALSE), matrix(unlist(res["NumberxCol"]),
                                                                                                        nrow = as.numeric(res["Number"]), byrow = FALSE), as.numeric(res["Number"]),
                       res["info"]))
}
.qubiclust_d <- function(x, c = 0.95, o = 100, f = 1, k = max(ncol(x)%/%20,
                                                              2), type = "default", P = FALSE, C = FALSE, verbose = TRUE) {
  MYCALL <- match.call()
  S <- (type == "area")
  res <- .qubic_d(x, c, o, f, k, P, S, C, verbose)
  return(BiclustResult(as.list(MYCALL), matrix(unlist(res["RowxNumber"]),
                                               ncol = as.numeric(res["Number"]), byrow = FALSE), matrix(unlist(res["NumberxCol"]),
                                                                                                        nrow = as.numeric(res["Number"]), byrow = FALSE), as.numeric(res["Number"]),
                       res["info"]))
}
