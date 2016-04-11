#' Construction and visualization of co-expression network
#'
#' This function can automatically create co-expression networks along with their visualization based on identified biclusters in QUBIC.
#' Three correlation methods, Pearson, Kendall and Spearman, are available for a user, facilitating different preferences in practical usage.
#' @aliases Qnetwork network Qunetwork
#' @param x The data matrix
#' @param BicRes biclust::BiclustResult object
#' @param number Which bicluster to be plotted
#' @param groups An object that indicates which nodes belong together.
#' @param method A character string indicating
#' which correlation coefficient (or covariance) is to be computed.
#' One of 'pearson' (default), 'kendall', or 'spearman', can be abbreviated.
#' @return a list contains a weights matrix and groupinfo
#' @examples
#' # Load microarray matrix
#' data(BicatYeast)
#' res <- biclust::biclust(BicatYeast[1:50, ], method=BCQU(), verbose = FALSE)
#' # Constructing the networks for the 4th and 13th identified biclusters.
#' net <- qunetwork(BicatYeast[1:50, ], res, number = c(4, 13), group = c(4, 13), method = 'spearman')
#' \dontrun{
#' if (requireNamespace('qgraph'))
#'     qgraph::qgraph(net[[1]], groups = net[[2]], layout = 'spring', minimum = 0.6,
#'        color = cbind(rainbow(length(net[[2]]) - 1),'gray'), edge.label = FALSE)
#'
#' }
#' \dontrun{
#' #Load microarray matrix
#' data(BicatYeast)
#' res <- biclust::biclust(BicatYeast[1:50, ], method=BCQU(), verbose = FALSE)
#' # Constructing the networks for the 4th and 13th identified biclusters,
#' #   using the whole network as a background.
#' net <- qunetwork(BicatYeast[1:50, ], res, group = c(4, 13), method = 'spearman')
#' if (requireNamespace('qgraph'))
#'     qgraph::qgraph(net[[1]], groups = net[[2]], layout = 'spring', minimum = 0.6,
#'        color = cbind(rainbow(length(net[[2]]) - 1),'gray'), edge.label = FALSE)
#' }
#' @seealso \code{\link{qunet2xml}} \code{\link{QUBIC}} \code{\link{cor}}
qunetwork <-  function(x, BicRes, number = 1:BicRes@Number,
                       groups = c(number[[1]]),
                       method = c("pearson", "kendall", "spearman")) {
  if (length(number) < 1)
    stop("at least 1 bicluster needed.")
  if (is.null(rownames(x)))
    stop("can not plot without rownames.")
  if (is.null(colnames(x)))
    stop("can not plot without colnames.")
  bics <- biclust::bicluster(x, BicRes, number)
  index <- which(number %in% groups)

  rownamelist <- list()
  colnamelist <- list()
  for (i in 1:length(bics)) {
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
