#' Convert newwork to XGMML
#'
#' This function can convert the constructed co-expression networks into XGMML format, which can be used to do further network analysis in Cytoscape, Biomax and JNets.
#' @param net Result of \code{\link{qunetwork}}
#' @param minimum cutoff, default: 0.6
#' @param color default: cbind(grDevices::rainbow(length(net[[2]]) - 1), 'gray')
#' @return Text of XGMML
#' @examples
#' # Load microarray matrix
#' data(BicatYeast)
#' res <- biclust::biclust(BicatYeast[1:50, ], method=BCQU(), verbose = FALSE)
#' # Get all biclusters
#' net <- qunetwork(BicatYeast[1:50, ], res, group = c(4, 13), method = 'spearman')
#' # Save the network to a XGMML file
#' sink('tempnetworkresult.gr')
#' qunet2xml(net, minimum = 0.6, color = cbind(grDevices::rainbow(length(net[[2]]) - 1), 'gray'))
#' sink()
#' # You can use Cytoscape, Biomax or JNets open file named tempnetworkresult.gr
#' @seealso \code{\link{qunetwork}} \code{\link{QUBIC}}
qunet2xml <- function(net, minimum = 0.6, color = cbind(grDevices::rainbow(length(net[[2]]) -
                                                                             1), "gray")) {
  cat("<?xml version=\"1.0\" encoding=\"utf-8\"?>")
  cat("\n")
  cat("<!DOCTYPE graph SYSTEM \"http://www.cs.rpi.edu/~puninj/XGMML/xgmml.dtd\">")
  cat("\n")
  cat("<graph directed=\"0\" label=\"QUNET2XML: Hello, I am a graph\" >")
  cat("\n")
  for (i in 1:length(net[[2]])) {
    for (j in net[[2]][[i]]) {
      cat("  <node label=\"")
      cat(rownames(net[[1]])[[j]])
      cat("\" id=\"")
      cat(j)
      cat("\"><graphics type=\"ELLIPSE\" fill=\"")
      cat(substr(color[[i]], 1, 7))
      cat("\"/></node>")
      cat("\n")
    }
  }
  for (i in 1:(nrow(net[[1]]) - 1)) {
    for (j in (i + 1):nrow(net[[1]])) {
      if (net[[1]][i, j] >= minimum) {
        cat("  <edge source=\"")
        cat(i)
        cat("\" target=\"")
        cat(j)
        cat("\" ><att name=\"weight\" type=\"real\" value=\"")
        cat(net[[1]][i, j])
        cat("\"/></edge>")
        cat("\n")
      }
    }
  }
  cat("</graph>")
  cat("\n")
}
