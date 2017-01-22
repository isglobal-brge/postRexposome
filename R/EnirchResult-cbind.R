#' @param ... Two or more object of class \code{EnrichResult}
#' @describeIn EnrichResult-class 
#' @export join
join <- function(...) {
    n <- sum(sapply(list(...), function(x) class(x) == "EnrichResult"))
    if(n != length(list(...))) {
        stop("At last one of the given arguments is not an 'EnrichedResult'.")
    }
    new("EnrichResult",
        fun_origin=lapply(list(...), function(x) x@fun_origin),
        search=lapply(list(...), function(x) x@search),
        database=lapply(list(...), function(x) x@database),
        eresult=lapply(list(...), function(x) x@eresult)
    )
}