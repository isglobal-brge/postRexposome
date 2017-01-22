#' Class \code{EnrichResult}
#'
#' @name EnrichResult-class
#' @aliases EnrichResult-class
#' @rdname EnrichResult-class
#' @exportClass EnrichResult
#' @slot fun_origin \code{list} containing the methods used to create
#' the \code{EnrichResult} object/s.
#' @slot search \code{list} containing the elements used to query the
#' web databses and to create \code{EnrichResult} object/s.
#' @slot database \code{list} containing the database queried to create
#' \code{EnrichResult} object/s.
#' @slot eresult \code{list} containing the results obtained on quering
#' web databases
#' @exportClass EnrichResult
setClass("EnrichResult",
          representation =
              representation( 
                  fun_origin = "list", # paired of methods used
                  search     = "list", # genes used to query
                  database   = "list", # where to search
                  eresult    = "list"  # list of results (one or more)
              ),
          prototype = 
              prototype( 
                  fun_origin = list(),
                  search     = list(),
                  database   = list(),
                  eresult    = list()
              )
)