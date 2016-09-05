######################################################################
#' QUBIC: A Qualitative Biclustering Algorithm for Analyses of Gene Expression Data
#'
#'@description
#' \code{QUBIC} is a biclustering package, with source code upgrading from C code to C++ code.
#' The updated source code can avoid memory allocation error and is much efficient than the original one.
#' Based on our preliminary analysis, it can save 40\% running time on a plant microarray data.
#'
#'@details
#' For a given representing matrix of a microarray data set,
#' we construct a weighted graph G with genes represented as vertices, edges connecting every pair of genes,
#' and the weight of each edge being the similarity level between the two corresponding (entire) rows.
#' Clearly, the higher a weight, the more similar two corresponding rows are.
#' Intuitively, genes in a bicluster should induce a heavier subgraph of G because under a subset of the conditions,
#' these genes have highly similar expression patterns that should make the weight of each involved edge heavier,
#' comparing to the edges in the background.
#' But it should be noted that some heavy subgraph may not necessarily correspond to a bicluster,
#' i.e. genes from a heavy subgraph may not necessarily have similar expression patterns
#' because different edges in a subgraph may have heavier weights under completely different subsets of conditions.
#' It should also be noted that recognizing all heavy subgraphs in a weighted graph itself is
#' computationally intractable because identification of maximum cliques in a graph is a special case of this,
#' and the maximum clique problem is a well known intractable problem (NP-hard).
#' So in our solution, we do not directly solve the problem of finding heavy subgraphs in a graph.
#' Instead, we built our biclustering algorithm based on this graph representation of a microarray gene expression data,
#' and tackle the biclustering problem as follows.
#' We find all feasible biclusters (I,J) in the given data set such that min\{|I|, |J|\} is as large as possible,
#' where I and J are subsets of genes and conditions, respectively.
#'
#' @name QUBIC
#'
#' @aliases QUBIC qubic BCQU bcqu biclust,matrix,BCQU-method
#' @param method \code{BCQU()} or \code{BCQUD()}, to perform QUBIC algorithm
#' @param x the input data matrix, which could be the normalized gene expression matrix or its qualitative representation from Qdiscretization or other discretization ways.
#' (for example: a qualitative representation of gene expression data) \cr
#' For \code{BCQU()}, the data matrix should be real \cr
#' For \code{BCQUD()}, the data matrix should be discretized as integer.
#' Zeros in the matrix will be treated as non-relevant value.
#' @param r Affect the granularity of the biclusters. The range of possible ranks.
#' A user can start with a small value of \code{r}
#' (the default value is \code{1} so the corresponding data matrix consists of values '\code{1}', '\code{-1}' and '\code{0}'),
#' evaluate the results, and then use larger values
#' (should not be larger than half of the number of the columns) to look for fine structures within the identified biclusters.
#' @param q Affect the granularity of the biclusters. The percentage of the regulating conditions for each gene.
#' The choice of \code{q}'s value depends on the specific application goals;
#' that is if the goal is to find genes that are responsive to local regulators,
#' we should use a relatively small \emph{q}-value; otherwise we may want to consider larger \emph{q}-values.
#' The default value of \code{q} is \code{0.06} in QUBIC
#' (this value is selected based on the optimal biclustering results on simulated data).
#' @param c The required consistency level of a bicluster. The default value of \code{c} is \code{0.95}
#' @param o The number of output biclusters. \code{o}'s default value is \code{100}.
#' @param f Control parameter, to control the level of overlaps between to-be-identified biclusters.
#' The filter cut-off for data post-processing. For overlaps among to-be-identified biclusters.
#' Its default value is set to \code{1} to ensure that no two reported biclusters overlap more than \code{f}.
#' @param k The minimum column width of the block, minimum \code{max(ncol(x) \%/\% 20, 2)} columns.
#' @param type The constrain type. \cr
#' If \code{type} is omitted or \code{type='default'}, the original objective function in QUBIC will be used, which is to maximize the minimal value of numbers of rows and columns.
#' If \code{type='area'}, the program tries to identify the bicluster with the maximal value of number of rows multiplied by number of columns.
#' Other types are reserved for future use.
#' @param P The flag to enlarge current bicluster using a \emph{p}-value contrain,
#' which is defined based on its significance of expression consistency  comparing to some simulated submatrix. Default: \code{FALSE}.
#' @param C The flag to set the lower bound of the condition number in a bicluster as 5\% of the total condition number in the input data.
#' Only suggested to use when the input data has a few conditions (e.g. less than \code{20}). Default: \code{FALSE}.
#' @param verbose If '\code{TRUE}', prints extra information on progress.
#' @param weight Alternative weight matrix provided by user, will append to default weight. \code{o}, \code{f}, \code{k}, \code{P}, \code{type}, \code{C} will be ignored if using this parameter.
#' @param seedbicluster Seed provided by user, normally should be a result of function \code{biclust}.
#' @return Returns an Biclust object, which contains bicluster candidates
#'
#' @seealso \code{\link{BCQU-class}} \code{\link{qudiscretize}} \code{\link{qunetwork}} \code{\link{qunet2xml}} \code{\link{biclust}}
#'
#' @references Li G, Ma Q, Tang H, Paterson AH, Xu Y.
#' QUBIC: a qualitative biclustering algorithm for analyses of gene expression data.
#' \emph{Nucleic Acids Research}. 2009;\bold{37(15)}:e101. doi:10.1093/nar/gkp491.
#' @references Zhou F, Ma Q, Li G, Xu Y.
#' QServer: A Biclustering Server for Prediction and Assessment of Co-Expressed Gene Clusters.
#' \emph{PLoS ONE}. 2012;\bold{7(3)}:e32660. doi: 10.1371/journal.pone.0032660
#'
#' @keywords qubic biclust bicluster bi-cluster biclustering bi-clustering
#'
#' @examples
#' # Random matrix with one embedded bicluster
#' test <- matrix(rnorm(5000), 100, 50)
#' test[11:20, 11:20] <- rnorm(100, 3, 0.3)
#' res <- biclust::biclust(test, method = BCQU())
#' summary(res)
#' show(res)
#' names(attributes(res))
#'
#' \dontrun{
#' # Load microarray matrix
#' data(BicatYeast)
#'
#' # Display number of column and row of BicatYeast
#' ncol(BicatYeast)
#' nrow(BicatYeast)
#' # Bicluster on microarray matrix
#' system.time(res <- biclust::biclust(BicatYeast, method = BCQU()))
#'
#' # Show bicluster info
#' res
#' # Show the first bicluster
#' biclust::bicluster(BicatYeast, res, 1)
#' # Get the 4th bicluster
#' bic4 <- biclust::bicluster(BicatYeast, res, 4)[[1]]
#'
#' # or
#' bic4 <- biclust::bicluster(BicatYeast, res)[[4]]
#' # Show rownames of the 4th bicluster
#' rownames(bic4)
#' # Show colnames of the 4th bicluster
#' colnames(bic4)
#'
#' }
#' \dontrun{
#' # Bicluster on selected of genes
#' data(EisenYeast)
#' genes <- c("YHR051W", "YKL181W", "YHR124W", "YHL020C", "YGR072W", "YGR145W",
#'     "YGR218W", "YGL041C", "YOR202W", "YCR005C")
#' # same result as res <- biclust::biclust(EisenYeast[1:10,], method=BCQU())
#' res <- biclust::biclust(EisenYeast[genes, ], method = BCQU())
#' res
#'
#' }
#' \dontrun{
#' # Get bicluster by row name = 249364_at
#' biclust::bicluster(BicatYeast, res, which(res@@RowxNumber[which(rownames(BicatYeast) ==
#'     "249364_at"), ]))
#'
#' }
#' \dontrun{
#' # Get bicluster by col name = cold_roots_6h
#' biclust::bicluster(BicatYeast, res, which(res@@NumberxCol[, which(colnames(BicatYeast) ==
#'     "cold_roots_6h")]))
#'
#' }
#' \dontrun{
#' # Draw a single bicluster using drawHeatmap {bicust}
#' data(BicatYeast)
#' res <- biclust::biclust(BicatYeast, BCQU(), verbose = FALSE)
#' # Draw heatmap of the first cluster
#' biclust::drawHeatmap(BicatYeast, res, 1)
#'
#' }
#' \dontrun{
#' # Draw a single bicluster using heatmap {stats}
#' data(BicatYeast)
#' res <- biclust::biclust(BicatYeast, BCQU(), verbose = FALSE)
#' bic10 <- biclust::bicluster(BicatYeast, res, 10)[[1]]
#'
#' # Draw heatmap of the 10th cluster using heatmap {stats}
#' heatmap(as.matrix(t(bic10)), Rowv = NA, Colv = NA, scale = 'none')
#'
#' # Draw heatmap of the 10th cluster using plot_heatmap {phyloseq}
#' if (requireNamespace('phyloseq'))
#'     phyloseq::plot_heatmap(otu_table(bic10, taxa_are_rows = TRUE))
#'
#' }
#' \dontrun{
#' # Draw a single bicluster with original data background and color options
#' data(BicatYeast)
#' res <- biclust::biclust(BicatYeast, BCQU(), verbose = FALSE)
#' palette <- colorRampPalette(c('red', 'yellow', 'green'))(n = 100)
#' # Draw heatmap of the first cluster with color
#' biclust::drawHeatmap(BicatYeast, res, 1, FALSE, beamercolor = TRUE, paleta = palette)
#'
#' }
#' \dontrun{
#' # Draw some overlapped biclusters
#' data(BicatYeast)
#' res <- biclust::biclust(BicatYeast, BCQU(), verbose = FALSE)
#' biclusternumber(res, 1)
#' biclusternumber(res, 3)
#' # Draw overlapping heatmap
#' biclust::heatmapBC(x = BicatYeast, bicResult = res, number = c(1, 3), local = TRUE)
#'
#' }
#' \dontrun{
#' # Draw all the biclusters
#' data(BicatYeast)
#' res <- biclust::biclust(BicatYeast, BCQU(), verbose = FALSE)
#' # Draw the first bicluster on heatmap
#' biclust::heatmapBC(x = BicatYeast, bicResult = res, number = 1)
#' # Draw all the biclusters, not working well.
#' # Overlap plotting only works for two neighbor bicluster defined by the order in the number slot.
#' biclust::heatmapBC(x = BicatYeast, bicResult = res, number = 0)
#'
#' }
NULL

#' Class BCQU.
#'
#' Class \code{BCQU} define a QUalitative BIClustering calcuator.
#'
#' @name BCQU-class
#' @rdname BCQU-class
#' @seealso \code{\link{BCQU}} \code{\link{qudiscretize}} \code{\link{qunetwork}} \code{\link{qunet2xml}} \code{\link{biclust}}

setClass(Class = "BCQU", contains = "BiclustMethod",
         prototype = prototype(biclustFunction = function(x, ...) {
  qubiclust(x, ...)
}))

#' @describeIn QUBIC Performs a QUalitative BIClustering.
#' @usage \S4method{biclust}{matrix,BCQU}(x, method = BCQU(),
#'         r = 1, q = 0.06,
#'         c = 0.95, o = 100, f = 1,
#'         k = max(ncol(x) \%/\% 20, 2),
#'         type = 'default', P = FALSE, C = FALSE, verbose = TRUE,
#'         weight = NULL, seedbicluster = NULL)
BCQU <- function(x = NULL, r = 1, q = 0.06, c = 0.95, o = 100, f = 1,
                 k = max(ncol(x) %/% 20, 2),
                 type = 'default', P = FALSE, C = FALSE, verbose = TRUE,
                 weight = NULL, seedbicluster = NULL) {
  if (is.null(x)) return(methods::new("BCQU"))
  res <- biclust(x = x, method = BCQU(), r = r, q = q, c = c, o = o, f = f,
                        k = k, type = type, P = P, C = C, verbose = verbose,
                        weight = weight, seedbicluster = seedbicluster)
  res@Parameters$Call = match.call()
  return(res)
}

#' QUBICD
#'
#' \code{BCQUD} performs a QUalitative BIClustering for a discret matrix.
#'
#' @name BCQUD-class
#'
#' @aliases qubic_d QUBICD QUD BCQUD-class biclust,matrix,BCQUD-method
#'
#' @rdname QUBIC
#'
#' @examples
#' # Biclustering of discretized yeast microarray data
#' data(BicatYeast)
#' disc<-qudiscretize(BicatYeast[1:10,1:10])
#' biclust::biclust(disc, method=BCQUD())
setClass("BCQUD", contains = "BiclustMethod",
         prototype = prototype(biclustFunction = function(x, ...) {
  qubiclust_d(x, ...)
}))

#' @describeIn QUBIC Performs a QUalitative BIClustering for a discret matrix.
#'
#' @usage \S4method{biclust}{matrix,BCQUD}(x, method = BCQUD(),
#'         c = 0.95, o = 100, f = 1,
#'         k = max(ncol(x) \%/\% 20, 2),
#'         type = 'default', P = FALSE, C = FALSE, verbose = TRUE,
#'         weight = NULL, seedbicluster = NULL)
BCQUD <- function(x = NULL, c = 0.95, o = 100, f = 1,
                  k = max(ncol(x) %/% 20, 2),
                  type = 'default', P = FALSE, C = FALSE, verbose = TRUE,
                  weight = NULL, seedbicluster = NULL) {
  if (is.null(x)) return(methods::new("BCQUD"))
  res <- biclust(x = x, method = BCQUD(), c = c, o = o, f = f,
                          k = k, type = type, P = P, C = C, verbose = verbose,
                          weight = weight, seedbicluster = seedbicluster)
  res@Parameters$Call = match.call()
  return(res)
}
