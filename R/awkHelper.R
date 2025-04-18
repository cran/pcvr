#' subset helper function for use reading in large data, called in pcv.sub.read
#'
#'
#'
#' @param inputFile Path to csv file of plantCV output, should be provided internally in read.pcv
#' @param filters filtering conditions, see read.pcv for details. Format as
#' list("trait in area, perimeter", "other contains stringToMatch")
#' @param awk Optional awk command to use instead.
#' @keywords read.csv pcv wide long
#' @details awkHelper attempts to make awk commands from human readable input.
#' Currently when filters are supplied the input file has quotes removed by `sed`
#' then is piped into awk, so an equivalent command line statement may be:
#' \code{sed 's/\"//g' pcvrTest2.csv | awk -F ','  '{ if (NR==1 || $18=="area") { print } }'}
#' @return Returns a character string representing a unix style awk statement
#'   which is typically passed to \code{pipe} or used as a connection in \code{data.table::fread}.
#' @importFrom utils read.csv capture.output
#' @examples
#' tryCatch(
#'   { # in case offline
#'     link1 <- "https://gist.githubusercontent.com/seankross/"
#'     link2 <- "a412dfbd88b3db70b74b/raw/5f23f993cd87c283ce766e7ac6b329ee7cc2e1d1/mtcars.csv"
#'     file <- paste0(link1, link2)
#'     awkHelper(file, list("gear in 4, 3"), awk = NULL)
#'     awkHelper(file, "gear contains 3", awk = NULL)
#'     # note that to be filtered the file has to exist on your local system, this example only shows
#'     # the output of awkHelper, which would then be executed by read.pcv on a unix system
#'     awkHelper(file, list("gear in 4, 3"), awk = "existing_command")
#'   },
#'   error = function(e) {
#'     message(e)
#'   }
#' )
#' @export
awkHelper <- function(inputFile, filters, awk = NULL) {
  if (is.null(awk)) {
    if (!is.list(filters)) {
      filters <- list(filters)
    }
    sed <- paste0("sed 's/\"//g' ", inputFile, " | ")
    awkStart <- "awk -F "
    awkDelim <- "',' "
    awkFiltStart <- "'{ if ("
    awkFiltEnd <- ") { print } }'"
    COLS <- colnames(read.csv(inputFile, nrows = 1))
    awkFilts <- lapply(filters, function(filt) {
      filtCol <- strsplit(filt, " ")[[1]][1]
      affector <- strsplit(filt, " ")[[1]][2]
      values <- trimws(gsub(",", " ", strsplit(filt, " ")[[1]][-c(1:2)]))
      if (affector %in% c("in", "is", "=")) {
        filt_string <- paste(paste0("($", which(COLS == filtCol), '=="', values, '")'),
                             collapse = " || ")
      } else if (affector == "contains") {
        valReg <- paste0(values, collapse = "|")
        filt_string <- paste0("($", which(COLS == filtCol), " ~ /", valReg, "/)")
      }
      return(filt_string)
    })
    awkFilt <- paste(paste("(", awkFilts, ")"), collapse = " && ")
    awkCommand <- capture.output(cat(sed, awkStart, awkDelim, awkFiltStart, awkFilt, awkFiltEnd))
  } else {
    awkCommand <- awk
  }
  return(awkCommand)
}
