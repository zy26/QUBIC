#' Construction and visualization of co-expression network
#'
#' This function can automatically create co-expression networks along with their visualization based on identified biclusters in QUBIC.
#' Three correlation methods, Pearson, Kendall and Spearman, are available for a user, facilitating different preferences in practical usage.
#' @aliases Qnetwork network Qunetwork
#' @param x The data matrix
#' @param BicRes Result object returned by \code{\link{qubiclust}} or
#' \code{\link{qubiclust_d}}
#' @param number Which bicluster to be plotted
#' @param groups An object that indicates which nodes belong together.
#' @param method A character string indicating
#' which correlation coefficient (or covariance) is to be computed.
#' One of 'pearson' (default), 'kendall', or 'spearman', can be abbreviated.
#' @return a list contains a weights matrix and groupinfo
#' @examples
#' # Load microarray matrix
#' if (requireNamespace("QUBICdata", quietly = TRUE)) {
#'   data(ecoli, package = "QUBICdata")
#'   res <- qubiclust(ecoli, r = 1, q = 0.06, c = 0.95,
#'                   o = 100, f = 0.25, k = max(ncol(ecoli) %/% 20, 2),
#'                   verbose = FALSE)
#'   # Constructing the network for the 5th identified bicluster.
#'   net <- qunetwork(ecoli, res, number = 5, groups = 5, method = "spearman")
#' }
#' \dontrun{
#' if (requireNamespace('qgraph'))
#'     qgraph::qgraph(net[[1]], groups = net[[2]], layout = 'spring', minimum = 0.6,
#'        color = cbind(rainbow(length(net[[2]]) - 1),'gray'), edge.labels = FALSE)
#'
#' }
#' \dontrun{
#' # Load microarray matrix
#' data(ecoli, package = "QUBICdata")
#' res <- qubiclust(ecoli, r = 1, q = 0.06, c = 0.95,
#'                 o = 100, f = 0.25, k = max(ncol(ecoli) %/% 20, 2),
#'                 verbose = FALSE)
#' # Constructing the networks for the 4th and 8th identified biclusters,
#' #   using the whole network as a background.
#' net <- qunetwork(ecoli, res, number = c(4, 8), groups = c(4, 8), method = "spearman")
#' if (requireNamespace('qgraph'))
#'     qgraph::qgraph(net[[1]], groups = net[[2]], layout = 'spring', minimum = 0.6,
#'        color = c('red', 'blue', 'gold', 'gray'), edge.labels = FALSE)
#' }
#' @seealso \code{\link{qunet2xml}} \code{\link{QUBIC}} \code{\link{cor}}
qunetwork <-  function(x, BicRes, number = NULL,
                       groups = NULL,
                       method = c("pearson", "kendall", "spearman")) {
  if (length(number) < 1)
    stop("at least 1 bicluster needed.")
  if (is.null(rownames(x)))
    stop("can not plot without rownames.")
  if (is.null(colnames(x)))
    stop("can not plot without colnames.")

  total_biclusters <- .qubic_result_number(BicRes)
  if (missing(number)) {
    number <- seq_len(total_biclusters)
  }
  if (is.null(number)) {
    number <- seq_len(total_biclusters)
  }
  if (is.null(groups)) {
    groups <- number[[1]]
  }

  bics <- .qubic_bicluster(x, BicRes, number)
  index <- which(number %in% groups)

  rownamelist <- list()
  colnamelist <- list()
  for (i in seq_along(bics)) {
    rownamelist[[names(bics)[i]]] <- rownames(bics[[i]])
    colnamelist[[names(bics)[i]]] <- colnames(bics[[i]])
  }

  allrownames <- Reduce(union, rownamelist)
  allcolnames <- Reduce(union, colnamelist)

  un <- x[allrownames, allcolnames]
  rowidlist <- list()

  if (length(groups) > 2)
    stop("length(group) > 2")
  if (length(groups) == 1) {
    rowidlist[[names(bics)[index[[1]]]]] <-
      match(rownamelist[[index[[1]]]],
            rownames(un))
  } else if (length(groups) == 2) {
    rowidlist[[paste(names(bics)[index[[1]]], " & ", names(bics)[index[[2]]],
                     sep = "")]] <-
      match(intersect(rownamelist[[index[[1]]]], rownamelist[[index[[2]]]]),
            rownames(un))
    rowidlist[[names(bics)[index[[1]]]]] <-
      match(setdiff(rownamelist[[index[[1]]]],
                    rownamelist[[index[[2]]]]), rownames(un))
    rowidlist[[names(bics)[index[[2]]]]] <-
      match(setdiff(rownamelist[[index[[2]]]],
                    rownamelist[[index[[1]]]]), rownames(un))
    rowidlist[["Others"]] <-
      match(setdiff(allrownames, union(rownamelist[[index[[1]]]],
                                       rownamelist[[index[[2]]]])), rownames(un))
  }

  cort <- stats::cor(t(un), method = method)

  return(list(cort, rowidlist))
}
