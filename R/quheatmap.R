#' Visualization of identified biclusters
#'
#' This function can visualize the identifed biclusters using heatmap in support of overall expression pattern analysis,either for a single bicluster or two biclusters.
#' @aliases quheatmap
#' @param x The data matrix
#' @param bicResult biclust::BiclustResult object
#' @param number which bicluster to be plotted
#' @param col default: c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
#' @param showlabel If TRUE, show the xlabel and ylabel
#' @param ... Additional options in \code{fields::image.plot}
#' @seealso \code{\link{qunet2xml}} \code{\link{QUBIC}} \code{\link{heatmapBC}}
#' @examples
#' # Load microarray matrix
#' data(BicatYeast)
#' res <- biclust::biclust(BicatYeast, method=BCQU(), verbose = FALSE)
#' # Draw heatmap for the 2th identified bicluster
#' par(mar = c(5, 4, 3, 5) + 0.1, mgp = c(0, 1, 0), cex.lab = 1.1, cex.axis = 0.5, cex.main = 1.1)
#' quheatmap(x = BicatYeast, res, number = 2, showlabel = TRUE)
#' # Draw heatmap for the 2th and 3th identified biclusters.
#' par(mar = c(5, 5, 5, 5), cex.lab = 1.1, cex.axis = 0.5, cex.main = 1.1)
#' quheatmap(x = BicatYeast, res, number = c(2, 3), showlabel = TRUE)


quheatmap = function(x, bicResult, number = 1, showlabel = FALSE, col = c("#313695",
    "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61",
    "#F46D43", "#D73027", "#A50026"), ...) {
    if (length(number) > 2)
        stop("Only 1 or 2 numbers are supported.")

    b <- biclust::bicluster(x, bicResult, number)

    if (length(number) == 1) {
        x2 <- b[[1]]
        graphics::image(1:ncol(x2), 1:nrow(x2), t(x2[1:nrow(x2), 1:ncol(x2)])[, nrow(x2):1],
            col = col, axes = FALSE, ylab = "", xlab = "")

        graphics::title(paste("Bicluster ", number, " (size ", nrow(x2), "x", ncol(x2),
            ")"))

        if (showlabel) {
            graphics::axis(side = 1, 1L:ncol(x2), labels = colnames(x2), las = 2,
                line = -0.5, tick = 0)
            graphics::axis(side = 2, nrow(x2):1L, labels = rownames(x2), las = 2,
                line = -0.5, tick = 0)
        }

        graphics::box()

        # add color legend using graphics::image.plot function
        if (requireNamespace("fields")) {
            low <- floor(min(x2[1:nrow(x2), 1:ncol(x2)]))
            up <- ceiling(max(x2[1:nrow(x2), 1:ncol(x2)]))
            fields::image.plot(zlim = c(low, up), legend.only = TRUE, horizontal = FALSE,
                col = col, ...)
        } else {
            warning("We need package fields to display legend")
        }
    } else if (length(number) == 2) {
        biclust::heatmapBC(x = x, bicResult = bicResult, number = number, local = TRUE,
            col = col, ylab = "", xlab = "")

        number1 = number[[1]]
        bicRows1 = which(bicResult@RowxNumber[, number1])
        bicCols1 = which(bicResult@NumberxCol[number1, ])

        number2 = number[[2]]
        bicRows2 = which(bicResult@RowxNumber[, number2])
        bicCols2 = which(bicResult@NumberxCol[number2, ])

        rows = c(setdiff(bicRows1, bicRows2), intersect(bicRows1, bicRows2),
            setdiff(bicRows2, bicRows1))
        cols = c(setdiff(bicCols1, bicCols2), intersect(bicCols1, bicCols2),
            setdiff(bicCols2, bicCols1))

        all = x[rows, cols]

        x1 = b[[1]]
        x2 = b[[2]]

        if (showlabel) {
            columnName <- colnames(all)
            rowName <- rownames(all)
            graphics::axis(1, at = 1:length(columnName), labels = columnName, las = 2,
                tick = 0, line = -0.5)
            graphics::axis(2, at = length(rowName):1, labels = rowName, las = 1,
                tick = 0, line = -0.5)
        }

        graphics::title(paste("Bicluster ", number[[1]], "( size ", nrow(x1), "x",
            ncol(x1), ")", "&", "Bicluster ", number[[2]], "( size ", nrow(x2),
            "x", ncol(x2), ")"), cex.main = 0.8)


        if (requireNamespace("fields")) {
            low <- floor(min(x[1:nrow(x), 1:ncol(x)]))
            up <- ceiling(max(x[1:nrow(x), 1:ncol(x)]))
            fields::image.plot(zlim = c(low, up), legend.only = TRUE, horizontal = FALSE,
                col = col, ...)
        } else {
            warning("We need package fields to display legend")
        }
    }
}
