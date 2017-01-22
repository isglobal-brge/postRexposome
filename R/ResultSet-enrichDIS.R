setMethod(
    f = "enrichDIS",
    signature = "ResultSet",
    definition = function(object, rid=1, fData.tag=1,
                          sel.pval="adj.P.Val", th.pval=0.01, 
                          sel.feature="genes", feature.null="", 
                          database="CURATED", verbose=FALSE, 
                          warnings=TRUE) {
        if(length(object) != 1 & missing(rid)) {
            stop("Given 'ResultSet' has more than one result and 'rid' ",
                 "is missing.")
        } else if(length(object) != 1 & !missing(rid) & class(rid) == "character") {
            if(!rid %in% names(object)) {
                stop("Given 'rid' ('", rid, "') not in 'ResultSet'.")
            }
        } else if(length(object) == 1) {
            rid <- 1
        }
        if(class(rid) == "numeric") {
            rid <- names(object@results)[rid]
        }
        
        if(class(fData.tag) == "numeric") {
            fData.tag <- names(Biobase::fData(object))[fData.tag]
        }
        tmet <- grep(fData.tag, names(Biobase::fData(object)))
        if(length(tmet) == 0) {
            stop("No datasets matching '", fData.tag, "' in given ",
                 "'ResultSet'.")
        } else if (length(tmet) > 1) {
            stop("Multiple datasets were used to create this ",
                 "'ResultSet'. There is no option to enrich a 'ResultSet' ",
                 "with multiple annotations for the same type of dataset.")
        }
        tmet <- names(Biobase::fData(object))[tmet]
        pd <- Biobase::fData(object)[[tmet]]
        dta <- object@results[[rid]]$result
        dta$gene <- sapply(pd[rownames(dta), sel.feature], function(x)
            strsplit(x, ";")[[1]][1])
        dta <- dta[!is.na(dta$gene), ]
        dta <- dta[dta$gene != feature.null, ]
        dta <- dta[dta[ , sel.pval] <= th.pval, ]
        
        requireNamespace("disgenet2r", quietly = TRUE)
        
        if(length(unique(dta$gene)) == 0) {
            stop("Filtering result in non probe under requirments. Try againg increasing 'th.pval'.")
        }
        
        
        dis <- disgenet2r::disgenetGene(gene=unique(dta$gene), database=database,
                                 verbose=verbose, warnings=warnings)
        dis <- list(dis)
        names(dis) <- rid
        
        EnrichSet("enrichDIS", object@fun_origin, unique(dta$gene), dis)
    }
)