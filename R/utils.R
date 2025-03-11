#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' Check if a directory exists
#'
#' @param dir_path Character string specifying directory path
#' @return Logical indicating if the directory exists
#' @keywords internal
._dirExists <- function(dir_path) {
  if (!is.character(dir_path) || length(dir_path) != 1) {
    stop("dir_path must be a single character string")
  }
  dir.exists(dir_path)
}

#' Create a directory if it doesn't exist
#'
#' @param dir_path Character string specifying directory path
#' @return Logical indicating if the directory exists or was created
#' @keywords internal
._createDirIfNotExists <- function(dir_path) {
  if (!._dirExists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
  }
  return(._dirExists(dir_path))
}

#' Check if a string ends with a pattern
#'
#' @param x Character vector
#' @param pattern Character string with pattern to match at end
#' @return Logical vector indicating which elements end with pattern
#' @keywords internal
#' @importFrom stringr str_detect
._endsWith <- function(x, pattern) {
  str_detect(x, paste0(pattern, "$"))
}

#' Validate that a file exists
#'
#' @param file_path Character string specifying file path
#' @return Logical indicating if the file exists
#' @keywords internal
._fileExists <- function(file_path) {
  if (!is.character(file_path) || length(file_path) != 1) {
    stop("file_path must be a single character string")
  }
  file.exists(file_path)
}