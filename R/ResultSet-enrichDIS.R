setMethod(
    f = "enrichDIS",
    signature = "ResultSet",
    definition = function(object, rid=1, coef=2, contrast=1,
        fData.tag=2, sel.pval="adj.P.Val", th.pval=0.01, sel.feature="gene", 
        feature.null=c("", "---"), translate.from=NA, database="CURATED", 
        verbose=FALSE, warnings=TRUE) {
        
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
        dta <- dta[dta[ , sel.pval] <= th.pval, ]
        
        ## ----------------------------------------------------------------- ##
        dta2 <- data.frame(transcrit="", gene="", p.value=0, stringsAsFactors=FALSE)
        for(ii in rownames(dta)) {
            for(it in strsplit(pd[ii, sel.feature], ";")[[1]]) {
                dta2 <- rbind(dta2, c(ii, it, dta[ii, sel.pval]))
            }
        }
        dta <- dta2[-1, ];
        rm(dta2, pd)
        ## ----------------------------------------------------------------- ##
        
        
        inst <- requireNamespace("disgenet2r", quietly = TRUE)
        if(!inst) {
            stop("This method requires 'disgenet2r'. Install 'disgenet2r' by running:\n  devtools::
install_bitbucket('ibi_group/disgenet2r')\nThen call again 'enrichDIS'.")
        }
        
        if(length(unique(dta$gene)) == 0) {
            stop("Filtering result in non probe under requirments. Try againg increasing 'th.pval'.")
        }
        
        
        if(!is.na(translate.from)) {
            genes <- clusterProfiler::bitr(unique(dta$gene), 
                fromType = translate.from, toType = "SYMBOL", 
                OrgDb="org.Hs.eg.db")
            genes <- unique(genes$SYMBOL)
        } else {
            genes <- unique(dta$gene)
        }
        
        dis <- disgenet2r::disgenetGene(gene=genes, database=database, 
                                        verbose=verbose, warnings=warnings)
        dis <- list(dis)
        names(dis) <- rid
        
        EnrichSet("enrichDIS", object@fun_origin, unique(dta$gene), database, dis)
    }
)