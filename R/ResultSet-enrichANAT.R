setMethod(
    f = "enrichANAT",
    signature = "ResultSet",
    definition = function(object, rid=1, fData.tag=1,
                          sel.pval="adj.P.Val", th.pval=0.01, 
                          sel.feature="genes", feature.null="", 
                          bgee.data="affymetrix", statistic="fisher",
                          verbose=FALSE, warnings=TRUE) {
        stop("- TODO -")
        if(length(object) != 1 & missing(rid)) {
            stop("Given 'ResultSet' has more than one result and 'rid' is missing.")
        } else if(length(object) != 1 & !missing(rid) & class(rid) == "character") {
            if(!rid %in% names(object)) {
                stop("Given 'rid' ('", rid, "') not in 'ResultSet'.")
            }
        } else if(length(object) == 1) {
            rid <- 1
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
        
        if(length(unique(dta$gene)) == 0) {
            stop("Filtering result in non probe under requirments. Try ",
                "againg increasing 'th.pval'.")
        }
        
        
        bgee <- BgeeDB::Bgee$new(species="Homo_sapiens", dataType=bgee.data)
        myTopAnatData <- BgeeDB::loadTopAnatData(bgee)
        
        
        conv <- read.delim(system.file(paste0("extdata", .Platform$file.sep,
                                              "ensemble2gene_symbol.tsv"), 
                                       package = "postRexposome"), 
                           stringsAsFactors = FALSE)
        
        myGenes <- conv[conv$hgnc_symbol %in% dta$gene, 1, drop = FALSE]
        colnames(myGenes) <- "ensembl_gene_id"
        
        ## The universe is composed by all the genes with ensemble id
        universe <- conv[ , 1, drop = FALSE]
        ## /
        
        ## Prepares the gene list vector 
        geneList <- factor(as.integer(universe[,1] %in% myGenes[,1]))
        names(geneList) <- universe[,1]
        ## /
        
        ## Perform enrichment in Bgee (TopAnat result)
        myTopAnatObject <-  BgeeDB::topAnat(myTopAnatData, geneList)
        ## /
        
        ## Prepares the topGO object
        results <- topGO::runTest(myTopAnatObject, 
                                  algorithm = 'classic', 
                                  statistic = statistic)
        ## 7
        
        ## Create a table with the full result without filtering for FDR 
        ## (setting the cutOff to 1 we include all results)
        tableOver <- BgeeDB::makeTable(myTopAnatData, 
                                       myTopAnatObject, 
                                       results, 
                                       cutOff=1)
        ## /
        
        EnrichSet("enrichANAT", object@fun_origin, unique(dta$gene), tableOver)
        
    }
)