setMethod(
    f = "enrichDISf",
    definition = function(object, family, fData.tag=1, fData.tag=1,
                          sel.pval="adj.P.Val", th.pval=0.01, 
                          sel.feature="genes", feature.null="", 
                          database="CURATED", verbose=FALSE, 
                          warnings=TRUE) {
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
        dta <- extract(object)
        ## --TO CHANGE------------------------------------------------------ ##
        dta$gene <- sapply(pd[rownames(dta), sel.feature], function(x)
            strsplit(x, ";")[[1]][1])
        ## ----------------------------------------------------------------- ##
        dta <- dta[!is.na(dta$gene), ]
        dta <- dta[dta$gene != feature.null, ]
        dta <- dta[dta[ , sel.pval] <= th.pval, ]
        
        inst <- requireNamespace("disgenet2r", quietly = TRUE)
        if(!inst) {
            stop("This method requires 'disgenet2r'. Install 'disgenet2r' by running:\n  devtools::
install_bitbucket('ibi_group/disgenet2r')\nThen call again 'enrichDISf'.")
        }
        
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