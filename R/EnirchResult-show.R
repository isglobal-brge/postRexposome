setMethod(
    f = "show",
    signature = "EnrichResult",
    definition = function(object) {
        cat("Object of class 'EnrichResult'\n")
        cat(" . Number of items: ....", length(object@eresult), "\n")
    }
)