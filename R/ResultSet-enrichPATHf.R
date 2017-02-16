setMethod(
    f = "enrichPATHf",
    signature = "ResultSet",
    definition = function(object, family, fData.exp=2, fData.omic=1, 
                          version="hg18", database="KEGG", 
                          sel.pval="adj.P.Val", th.pval=0.01, 
                          sel.feature="genes", feature.id="geneSymbol", 
                          feature.null="", plot.fit=FALSE, verbose=FALSE, 
                          warnings=TRUE) {
        
        version <- match.arg(version, choices=c("hg18", "hg19"))
        database <- match.arg(toupper(database), 
                              choices=c("KEGG", "GO:CC", "GO:BP", "GO:MF"))
        
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
        
        ## CREATE VECTOR OF DE GENES
        genes <- sapply(unique(dta$gene), function(gn) {    
            sum(dta[dta$gene == gn, sel.pval] <= th.pval) != 0 
        })
        ## /
        
        ## CREATE VECTOR OF CPGS PER GENE
        countbias=as.data.frame(table(dta$gene))
        n <- countbias[ , 1]
        countbias <- countbias[ , 2]
        names(countbias) <- n; rm(n)
        ## /
        
        ## ORDER BOTH VECTORS
        or <- names(genes)[order(names(genes))]
        countbias <- countbias[or]
        genes <- genes[or]; rm(or)
        ## /
        
        if(verbose | warnings) {
            if(sum(dta=="") != 0 ) {
                warning("Detected ", sum(dta==""), " probes (CpG) with no gene ",
                        "identifier assigned.")
            }
        }
        
        
        if(verbose) {
            message("Probability Weighting Function was called with biased data.")
        }
        pwd <- goseq::nullp(genes, genome="hg19", id=feature.id, 
                            bias.data=NULL, plot.fit=plot.fit)
        
        if(verbose) {
            message("Performing query to database ('", database, "').")
        }
        wall <- goseq::goseq(pwd, genome=version, id=feature.id, 
                             test.cats=database)
        wall <- list(wall)
        names(wall) <- rid
        
        EnrichSet("enrichPATH", object@fun_origin, unique(dta$gene), database, wall)
    }
)