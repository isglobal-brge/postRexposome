setMethod(
    f = "enrichANAT",
    signature = "ResultSet",
    definition = function(object, rid=1, coef=2, contrast=1,
        fData.tag=2, sel.pval="adj.P.Val", th.pval=0.01, sel.feature="genes", 
        feature.null=c("", "---"), translate.from=NA, bgee.data="affymetrix", 
        statistic="fisher", verbose=FALSE, warnings=TRUE) {
        
        if(length(object) == 1) {
            rid <- 1
        }
        
        if(class(rid) == "character") {
            if(!rid %in% omicRexposome::rid(object)) {
                stop("Given 'rid' ('", rid, "') not in 'ResultSet'.")
            }
        } else { # if(class(rid) == "numeric") {
            rid <- omicRexposome::rid(object)[rid]
        }
        
        if(class(fData.tag) != "numeric") {
            stop("Content of 'fData.tag' must be numeric.")
        }
        
        pd <- Biobase::fData(object)[[fData.tag]]
        dta <- omicRexposome::topTable(object, rid=rid, coef=coef, contrast=contrast)
        dta <- dta[!is.na(dta$gene), ]
        dta <- dta[!dta$gene %in% feature.null, ]
        
        ## ----------------------------------------------------------------- ##
        dta2 <- data.frame(transcrit="", gene="", p.value=0, stringsAsFactors=FALSE)
        for(ii in rownames(dta)) {
            for(it in strsplit(pd[ii, sel.feature], ";")[[1]]) {
                dta2 <- rbind(dta2, c(ii, it, dta[ii, sel.pval]))
            }
        }
        dta <- dta2[-1, ];
        dta$p.value <- as.numeric(dta$p.value)
        rm(dta2, pd)
        ## ----------------------------------------------------------------- ##
        
        
        if(length(unique(dta$gene)) == 0) {
            stop("Filtering result in non probe under requirments. Try ",
                "againg increasing 'th.pval'.")
        }
        
        
        bgee <- BgeeDB::Bgee$new(species="Homo_sapiens", dataType=bgee.data)
        myTopAnatData <- BgeeDB::loadTopAnatData(bgee)
        
        
        if(!is.na(translate.from)) {
            xx <- clusterProfiler::bitr(unique(dta$gene), 
                                           fromType = translate.from, toType = "ENSEMBLTRANS", 
                                           OrgDb="org.Hs.eg.db")
            xx$p.value <- sapply(xx[ , 1], function(it) {
                min(dta$p.value[dta$gene == it])
            })
            dta <- xx; rm(xx)
            colnames(dta) <- c("gene", "_trs", "p.value")
        } else {
            dta$`_trs` <- dta[ , sel.feature]
        }
        
        filtered <- dta$`_trs`[dta[ , "p.value"] <= th.pval]
        if(length(filtered) == 0) {
            stop("Any features is under given significant threashold (p.value <= ", th.pval, ")")
        }
        geneList <- factor(as.integer(dta$`_trs` %in% filtered))
        names(geneList) <- dta$`_trs`
        
        ## Perform enrichment in Bgee (TopAnat result)
        myTopAnatObject <- tryCatch({
            BgeeDB::topAnat(myTopAnatData, geneList, nodeSize=1)
        }, error = function(e) {
            stop("Any of selected features are in BgeeDB.\n\n", e)
        })
        ## /
        
        ## Prepares the topGO object
        results <- topGO::runTest(myTopAnatObject, 
                                  algorithm = 'classic', 
                                  statistic = statistic)
        ## /
        
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