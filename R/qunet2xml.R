#' Convert newwork to XGMML
#'
#' This function can convert the constructed co-expression networks into XGMML format, which can be used to do further network analysis in Cytoscape, Biomax and JNets.
#' @param net Result of \code{\link{qunetwork}}
#' @param minimum cutoff, default: 0.6
#' @param color default: cbind(grDevices::rainbow(length(net[[2]]) - 1), 'gray')
#' @return Text of XGMML
#' @examples
#' # Load microarray matrix
#' if (requireNamespace("QUBICdata", quietly = TRUE)) {
#'   data(ecoli, package = "QUBICdata")
#'   res <- qubiclust(ecoli, r = 1, q = 0.06, c = 0.95,
#'                   o = 100, f = 0.25, k = max(ncol(ecoli) %/% 20, 2),
#'                   verbose = FALSE)
#'   # Get all biclusters
#'   net <- qunetwork(ecoli, res, number = 5, groups = 5, method = "spearman")
#'   # Save the network to a XGMML file
#'   outfile <- tempfile(fileext = ".gr")
#'   sink(outfile)
#'   qunet2xml(net, minimum = 0.6,
#'            color = cbind(grDevices::rainbow(length(net[[2]]) - 1), "gray"))
#'   sink()
#'   unlink(outfile)
#' }
#' @seealso \code{\link{qunetwork}} \code{\link{QUBIC}}
qunet2xml <- function(net, minimum = 0.6, color = cbind(grDevices::rainbow(length(net[[2]]) -
                                                                             1), "gray")) {
  lines <- c(
    "<?xml version=\"1.0\" encoding=\"utf-8\"?>",
    "<!DOCTYPE graph SYSTEM \"http://www.cs.rpi.edu/~puninj/XGMML/xgmml.dtd\">",
    "<graph directed=\"0\" label=\"QUNET2XML: Hello, I am a graph\" >"
  )

  for (i in seq_along(net[[2]])) {
    for (j in net[[2]][[i]]) {
      lines <- c(lines, paste0(
        "  <node label=\"", rownames(net[[1]])[[j]],
        "\" id=\"", j,
        "\"><graphics type=\"ELLIPSE\" fill=\"", substr(color[[i]], 1, 7),
        "\"/></node>"
      ))
    }
  }

  n_nodes <- nrow(net[[1]])
  if (n_nodes >= 2L) {
    for (i in seq_len(n_nodes - 1L)) {
      for (j in seq.int(i + 1L, n_nodes)) {
      if (net[[1]][i, j] >= minimum) {
        lines <- c(lines, paste0(
          "  <edge source=\"", i,
          "\" target=\"", j,
          "\" ><att name=\"weight\" type=\"real\" value=\"", net[[1]][i, j],
          "\"/></edge>"
        ))
      }
      }
    }
  }

  lines <- c(lines, "</graph>")
  writeLines(lines)
}
