#' Visualization of identified biclusters
#'
#' This function can visualize the identifed biclusters using heatmap in support of overall expression pattern analysis,either for a single bicluster or two biclusters.
#' @aliases quheatmap
#' @param x The data matrix
#' @param bicResult Result object returned by \code{\link{qubiclust}} or
#' \code{\link{qubiclust_d}}
#' @param number which bicluster to be plotted
#' @param col default: c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
#' @param showlabel If TRUE, show the xlabel and ylabel
#' @param ... Additional options in \code{fields::image.plot}
#' @return Invisibly returns \code{NULL} after drawing the heatmap.
#' @seealso \code{\link{qunet2xml}} \code{\link{QUBIC}} \code{\link{qunetwork}}
#' @examples
#' # Load microarray matrix
#' if (requireNamespace("QUBICdata", quietly = TRUE)) {
#'   data(ecoli, package = "QUBICdata")
#'   res <- qubiclust(ecoli, r = 1, q = 0.06, c = 0.95,
#'                   o = 100, f = 0.25, k = max(ncol(ecoli) %/% 20, 2),
#'                   verbose = FALSE)
#'   # Draw heatmap for the 5th identified bicluster
#'   par(mar = c(5, 4, 3, 5) + 0.1, mgp = c(0, 1, 0), cex.lab = 1.1,
#'       cex.axis = 0.5, cex.main = 1.1)
#'   quheatmap(x = ecoli, res, number = 5, showlabel = TRUE)
#'   # Draw heatmap for the 4th and 8th identified biclusters.
#'   par(mar = c(5, 5, 5, 5), cex.lab = 1.1, cex.axis = 0.5, cex.main = 1.1)
#'   quheatmap(x = ecoli, res, number = c(4, 8), showlabel = TRUE)
#' }


quheatmap <- function(x, bicResult, number = 1, showlabel = FALSE, col = c("#313695",
    "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61",
    "#F46D43", "#D73027", "#A50026"), ...) {
    if (length(number) > 2)
        stop("Only 1 or 2 numbers are supported.")

    b <- .qubic_bicluster(x, bicResult, number)

    if (length(number) == 1) {
        x2 <- b[[1]]
        graphics::image(seq_len(ncol(x2)), seq_len(nrow(x2)),
            t(x2[seq_len(nrow(x2)), seq_len(ncol(x2))])[, rev(seq_len(nrow(x2)))],
            col = col, axes = FALSE, ylab = "", xlab = "")

        graphics::title(paste("Bicluster ", number, " (size ", nrow(x2), "x", ncol(x2),
            ")"))

        if (showlabel) {
            graphics::axis(side = 1, seq_len(ncol(x2)), labels = colnames(x2), las = 2,
                line = -0.5, tick = 0)
            graphics::axis(side = 2, rev(seq_len(nrow(x2))), labels = rownames(x2), las = 2,
                line = -0.5, tick = 0)
        }

        graphics::box()

        # add color legend using graphics::image.plot function
        if (requireNamespace("fields")) {
            low <- floor(min(x2[seq_len(nrow(x2)), seq_len(ncol(x2))]))
            up <- ceiling(max(x2[seq_len(nrow(x2)), seq_len(ncol(x2))]))
            fields::image.plot(zlim = c(low, up), legend.only = TRUE, horizontal = FALSE,
                col = col, ...)
        } else {
            warning("We need package fields to display legend")
        }
    } else if (length(number) == 2) {
        number1 <- number[[1]]
        rowxnumber <- .qubic_result_rowxnumber(bicResult)
        numberxcol <- .qubic_result_numberxcol(bicResult)

        bicRows1 <- which(rowxnumber[, number1])
        bicCols1 <- which(numberxcol[number1, ])

        number2 <- number[[2]]
        bicRows2 <- which(rowxnumber[, number2])
        bicCols2 <- which(numberxcol[number2, ])

        rows <- c(setdiff(bicRows1, bicRows2), intersect(bicRows1, bicRows2),
            setdiff(bicRows2, bicRows1))
        cols <- c(setdiff(bicCols1, bicCols2), intersect(bicCols1, bicCols2),
            setdiff(bicCols2, bicCols1))

        all <- x[rows, cols]

        nr <- nrow(all)
        nc <- ncol(all)
        graphics::image(seq_len(nc), seq_len(nr), t(all)[, rev(seq_len(nr))], col = col, axes = FALSE,
            ylab = "", xlab = "")

        nRow1Only <- length(setdiff(bicRows1, bicRows2))
        nRowBoth <- length(intersect(bicRows1, bicRows2))
        nCol1Only <- length(setdiff(bicCols1, bicCols2))
        nColBoth <- length(intersect(bicCols1, bicCols2))

        # Outline the first bicluster (top-left block in the reordered matrix).
        graphics::rect(xleft = 0.5,
            xright = nCol1Only + nColBoth + 0.5,
            ybottom = nr - (nRow1Only + nRowBoth) + 0.5,
            ytop = nr + 0.5)

        # Outline the second bicluster (bottom-right block in the reordered matrix).
        graphics::rect(xleft = nCol1Only + 0.5,
            xright = nc + 0.5,
            ybottom = 0.5,
            ytop = nr - nRow1Only + 0.5)

        x1 <- b[[1]]
        x2 <- b[[2]]

        if (showlabel) {
            columnName <- colnames(all)
            rowName <- rownames(all)
            graphics::axis(1, at = seq_len(length(columnName)), labels = columnName, las = 2,
                tick = 0, line = -0.5)
            graphics::axis(2, at = rev(seq_len(length(rowName))), labels = rowName, las = 1,
                tick = 0, line = -0.5)
        }

        graphics::box()

        graphics::title(paste("Bicluster ", number[[1]], "( size ", nrow(x1), "x",
            ncol(x1), ")", "&", "Bicluster ", number[[2]], "( size ", nrow(x2),
            "x", ncol(x2), ")"), cex.main = 0.8)


        if (requireNamespace("fields")) {
            low <- floor(min(x[seq_len(nrow(x)), seq_len(ncol(x))]))
            up <- ceiling(max(x[seq_len(nrow(x)), seq_len(ncol(x))]))
            fields::image.plot(zlim = c(low, up), legend.only = TRUE, horizontal = FALSE,
                col = col, ...)
        } else {
            warning("We need package fields to display legend")
        }
    }
}
