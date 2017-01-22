#' Function to create a \code{ResultSet} from scratch
#' 
#' @param fun_originA Function used to create the \code{ResultSet} (aka. 
#' \code{enirchDIS}, ...).
#' @param fun_originB (optional) Function used to create the base object
#' of the new \code{ResultSet}. Usually \code{rexposome}'s functions (aka.
#' \code{assocGE}, ...).
#' @param search \code{character} used in the database query.
#' @param database \code{character} indicating the databased used.
#' @param eresult \code{list} with a single element, the result obtained from 
#' the query.
#' @describeIn EnrichResult-class 
#' @export EnrichSet
EnrichSet <- function(fun_originA, fun_originB, search, database, eresult) {
    if(missing(fun_originB)) {
        new("EnrichResult",
            fun_origin=list(fun_originA),
            search=list(search),
            database=list(database),
            eresult=eresult
        )
    } else {
        new("EnrichResult",
            fun_origin=list(c(fun_originA, fun_originB)),
            search=list(search),
            database=list(database),
            eresult=eresult
        )
    }
}