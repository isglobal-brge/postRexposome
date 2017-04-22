#' Function to combine two or more ResultSet objects
#' 
#' This function allows to combine multiple \code{\link{ResultSet}}
#' objects in a single one.
#' 
#' @param ... Two or more object of class \code{EnrichResult}
#' @describeIn EnrichResult-class 
#' @export join
#' @examples 
#' data(air)
#' no2 <- enrichDIS(air, rid="NO2", sel.pval="P.Value", warnings=FALSE,
#'                  sel.feature="gene_assignment", translate.from="REFSEQ")
#' ben <- enrichDIS(air, rid="Ben", sel.pval="P.Value", warnings=FALSE,
#'                  sel.feature="gene_assignment", translate.from="REFSEQ")
#' join(no2, ben)
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