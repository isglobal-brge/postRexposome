setMethod(
    f = "enrichGWAS",
    signature = "ResultSet",
    definition = function(object, version=c("hg19", "hg38"), 
                          verbose=FALSE, warnings=TRUE) {
        ## CHECK FOR GENOME VERSION
        version <- match.arg(version, choices = c("hg19", "hg38", "current"))
        if(verbose) {
            message("Selected chrmosome version was '", version, "'. ",
                "Loading data from GWAS catalog R package.")
        }
        
        if(version == "hg19") {
            data(ebicat37, package="gwascat")
            dta <- ebicat37
            rm(ebicat37)
        } else if (version == "hg38") {
            data(ebicat38, package="gwascat")
            dta <- ebicat38
            rm(ebicat38)
        } else { # get current versio
            # it raises an error:
            # Error in gwascat::makeCurrentGwascat() : object 'si.hs.38' not found
            dta <- gwascat::makeCurrentGwascat()
        }
        ## /
        
        if(object@fun_origin == "assocSNP") {
            rst <- object@results[[1]]$result
            intr <- intersect(gwascat::getRsids(dta), rst$Name)
            if(verbose) {
                message("Performing intersecion between given 'ResultSet' ",
                    "and GWAS Catalog ('", version, "').")
                message("Obtained ", length(intr), " matches (", 
                        round(length(intr) / nrow(rst) * 100, digits=2), "%).")
            }
            intr <- as.data.frame(dta[intr])
        } else if (object@fun_origin == "crossomics") {
            
        } else {
            stop("This method can only be aplied to 'ResultSet' objects ",
                "generated with 'assocSNP' or 'crossomics'.")
        }
        
        return(intr)
    }
)