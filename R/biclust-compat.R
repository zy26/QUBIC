# Internal compatibility layer for biclust-dependent operations.
# Phase 2 skeleton introduces an internal result object and unified accessors.

if (!methods::isClass("QUBICBiclustResult")) {
  methods::setClass(
    "QUBICBiclustResult",
    slots = c(
      Parameters = "list",
      RowxNumber = "matrix",
      NumberxCol = "matrix",
      Number = "numeric",
      info = "list"
    )
  )
}

.qubic_has_biclust <- function() {
  !isTRUE(getOption("qubic.disable.biclust", FALSE)) &&
    requireNamespace("biclust", quietly = TRUE)
}

.qubic_require_biclust <- function() {
  if (!.qubic_has_biclust()) {
    stop("Package 'biclust' is required for this operation.", call. = FALSE)
  }
}

.qubic_new_internal_result <- function(parameters, rowxnumber, numberxcol, number, info) {
  methods::new(
    "QUBICBiclustResult",
    Parameters = parameters,
    RowxNumber = rowxnumber,
    NumberxCol = numberxcol,
    Number = as.numeric(number),
    info = as.list(info)
  )
}

.qubic_is_internal_result <- function(result) {
  isS4(result) && "QUBICBiclustResult" %in% class(result)
}

.qubic_result_number <- function(result) {
  if (.qubic_is_internal_result(result)) {
    return(as.numeric(result@Number))
  }
  as.numeric(result@Number)
}

.qubic_result_rowxnumber <- function(result) {
  result@RowxNumber
}

.qubic_result_numberxcol <- function(result) {
  result@NumberxCol
}

.qubic_result_call <- function(result) {
  result@Parameters$Call
}

.qubic_set_result_call <- function(result, call) {
  result@Parameters$Call <- call
  result
}

.qubic_extract_biclusters_internal <- function(x, result, number) {
  rowxnumber <- .qubic_result_rowxnumber(result)
  numberxcol <- .qubic_result_numberxcol(result)

  if (missing(number)) {
    number <- seq_len(.qubic_result_number(result))
  }
  number <- as.integer(number)

  out <- lapply(number, function(i) {
    rows <- which(rowxnumber[, i])
    cols <- which(numberxcol[i, ])
    x[rows, cols, drop = FALSE]
  })
  names(out) <- paste0("Bicluster ", number)
  out
}

.qubic_biclust <- function(...) {
  .qubic_require_biclust()
  biclust::biclust(...)
}

.qubic_bicluster <- function(x, result, number = NULL) {
  if (.qubic_has_biclust()) {
    if (is.null(number)) {
      return(biclust::bicluster(x, result))
    }
    return(biclust::bicluster(x, result, number))
  }
  if (is.null(number)) {
    return(.qubic_extract_biclusters_internal(x, result))
  }
  .qubic_extract_biclusters_internal(x, result, number)
}

.qubic_new_BiclustResult <- function(parameters, rowxnumber, numberxcol, number, info) {
  # Keep qubiclust() outputs in QUBIC's own result class regardless of biclust.
  .qubic_new_internal_result(parameters, rowxnumber, numberxcol, number, info)
}

.qubic_print_summary <- function(object) {
  total <- as.integer(.qubic_result_number(object))
  rowxnumber <- .qubic_result_rowxnumber(object)
  numberxcol <- .qubic_result_numberxcol(object)
  class_label <- if (.qubic_is_internal_result(object)) {
    "QUBICBiclustResult (biclust-compatible summary)"
  } else {
    "Biclust"
  }

  call_expr <- .qubic_result_call(object)
  call_text <- if (is.null(call_expr)) "NULL" else paste(deparse(call_expr), collapse = "\n")

  writeLines(c(
    "",
    paste0("An object of class ", class_label),
    "",
    "call:",
    paste0("\t", call_text),
    "",
    paste0("Number of Clusters found: ", total),
    ""
  ))

  if (total <= 0) {
    return(invisible(NULL))
  }

  row_counts <- colSums(rowxnumber != 0)
  col_counts <- rowSums(numberxcol != 0)
  cluster_names <- paste("BC", seq_len(total))

  stats <- rbind(
    "Number of Rows" = row_counts,
    "Number of Columns" = col_counts
  )
  colnames(stats) <- cluster_names

  writeLines("Cluster sizes:")
  writeLines(utils::capture.output(stats))
  writeLines("")
  invisible(NULL)
}

.qubic_print_brief <- function(object) {
  total <- as.integer(.qubic_result_number(object))
  rowxnumber <- .qubic_result_rowxnumber(object)
  numberxcol <- .qubic_result_numberxcol(object)
  class_label <- if (.qubic_is_internal_result(object)) {
    "QUBICBiclustResult"
  } else {
    "Biclust"
  }

  call_expr <- .qubic_result_call(object)
  call_text <- if (is.null(call_expr)) "NULL" else paste(deparse(call_expr), collapse = "\n")

  writeLines(c(
    "",
    paste0("An object of class ", class_label),
    "",
    "call:",
    paste0("\t", call_text),
    "",
    paste0("Number of Clusters found: ", total),
    ""
  ))

  if (total <= 0) {
    return(invisible(NULL))
  }

  n_show <- min(total, 5L)
  idx <- seq_len(n_show)
  row_counts <- colSums(rowxnumber != 0)[idx]
  col_counts <- rowSums(numberxcol != 0)[idx]

  stats <- rbind(
    "Number of Rows:" = row_counts,
    "Number of Columns:" = col_counts
  )
  colnames(stats) <- paste("BC", idx)

  writeLines(paste0("First  ", n_show, "  Cluster sizes:"))
  writeLines(utils::capture.output(stats))
  writeLines("")
  invisible(NULL)
}

methods::setMethod("show", signature(object = "QUBICBiclustResult"), function(object) {
  .qubic_print_brief(object)
})

methods::setMethod("summary", signature(object = "QUBICBiclustResult"), function(object, ...) {
  .qubic_print_summary(object)
})

summary.QUBICBiclustResult <- function(object, ...) {
  .qubic_print_summary(object)
}
