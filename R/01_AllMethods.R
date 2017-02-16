#' Function to enrich results from differential studies on DisGeNET database
#' 
#' This method takes the results from a \code{ResultSet} and uses the genes
#' associated to the probes to perform a query in DisGeNET and obtain the
#' diseases associated to the given genes.
#' 
#' @param object Object of class \code{ResultSet}.
#' @param rid (default \code{1}) If given \code{ResultSet} hat has more than 
#' one result, id to select the result where the enrichment will be performed.
#' @param fData.tag (default \code{1}) Identifier of the annotation to be
#' used to mach the probes of \code{ResultSet} with gene identifier. In other
#' words, index of the element in \code{fData} that matches the omic dataset 
#' to be used in the enrichment.
#' @param sel.pval (default \code{"adj.P.Val"}) Name of the P.Value columns to 
#' filter datset's probes given \code{th.pval}.
#' @param th.pval (default \code{0.01}) Threshold used to include a probe (CpG)
#' into the enrichment analysis.
#' @param sel.feature (default \code{"genes"}) Name containing gene identifier
#' in \code{ResultSet}'s \code{fData}.
#' @param feature.null (default \code{""}) String identifier for non 
#' gene-assigned probes.
#' @param database (default \code{"CURATED"}) Name of the version of the
#' DisGeNET's databse where the query will be computed. Check 
#' \code{?disgenetGene} for a description of each name.
#' @param verbose (default \code{FALSE}) If set to \code{TRUE} messages
#' indicating the steps done by the method are raised.
#' @param warnings (default \code{TRUE}) If set to \code{FALSE} important
#' notes raised by the method are hide.
#' @export enrichDIS
#' @importClassesFrom rexposome ResultSet
#' @note 
#' Piñero J, Queralt-Rosinach N, Bravo A, Deu-Pons J, Bauer-Mehren A, Baron M, 
#' Sanz F, Furlong LI;
#' DisGeNET: a discovery platform for the dynamical exploration of human 
#' diseases and their genes.; Database (2015).
#' 
#' Queralt-Rosinach N, Piñero J, Bravo A, Sanz F, Furlong LI;
#' DisGeNET-RDF: harnessing the innovative power of the Semantic Web to explore
#' the genetic basis of diseases Bioinformatics; Bioinformatics (2016).

setGeneric("enrichDIS", function(object, rid=1, fData.tag=1,
      sel.pval=c("P.Value", "adj.P.Val"), th.pval=0.01, sel.feature="genes",
      feature.null="", database=c("CTD_human", "UNIPROT", "CLINVAR", "GWASCAT",
      "ORPHANET", "CURATED", "RGD", "MGD", "CTD_rat", "CTD_mouse",
      "PREDICTED", "ALL"), verbose=FALSE, warnings=TRUE)
    standardGeneric("enrichDIS")
)

#' Function to enrich results from differential studies on DisGeNET database
#' 
#' This method takes the results from a \code{ResultSet} and uses the genes
#' associated to the probes to perform a query in DisGeNET and obtain the
#' diseases associated to the given genes.
#' 
#' @param object Object of class \code{ResultSet}.
#' @param family Family of exposures used to perform the enrichment analysis.
#' @param fData.exp (default \code{2}) Identifier of the exposures decrition 
#' to be used to mach the exposures of \code{ResultSet} with their families.
#' @param fData.omic (default \code{1}) Identifier of the annotation to be
#' used to mach the probes of \code{ResultSet} with gene identifier. In other
#' words, index of the element in \code{fData} that matches the omic dataset 
#' to be used in the enrichment.
#' @param sel.pval (default \code{"adj.P.Val"}) Name of the P.Value columns to 
#' filter datset's probes given \code{th.pval}.
#' @param th.pval (default \code{0.01}) Threshold used to include a probe (CpG)
#' into the enrichment analysis.
#' @param sel.feature (default \code{"genes"}) Name containing gene identifier
#' in \code{ResultSet}'s \code{fData}.
#' @param feature.null (default \code{""}) String identifier for non 
#' gene-assigned probes.
#' @param database (default \code{"CURATED"}) Name of the version of the
#' DisGeNET's databse where the query will be computed. Check 
#' \code{?disgenetGene} for a description of each name.
#' @param verbose (default \code{FALSE}) If set to \code{TRUE} messages
#' indicating the steps done by the method are raised.
#' @param warnings (default \code{TRUE}) If set to \code{FALSE} important
#' notes raised by the method are hide.
#' @export enrichDISf
#' @importClassesFrom rexposome ResultSet
#' @note 
#' Piñero J, Queralt-Rosinach N, Bravo A, Deu-Pons J, Bauer-Mehren A, Baron M, 
#' Sanz F, Furlong LI;
#' DisGeNET: a discovery platform for the dynamical exploration of human 
#' diseases and their genes.; Database (2015).
#' 
#' Queralt-Rosinach N, Piñero J, Bravo A, Sanz F, Furlong LI;
#' DisGeNET-RDF: harnessing the innovative power of the Semantic Web to explore
#' the genetic basis of diseases Bioinformatics; Bioinformatics (2016).
setGeneric("enrichDISf", function(object, family, fData.exp=2,  fData.omic=1,
        sel.pval=c("P.Value", "adj.P.Val"), th.pval=0.01, sel.feature="genes",
        feature.null="", database=c("CTD_human", "UNIPROT", "CLINVAR", "GWASCAT",
        "ORPHANET", "CURATED", "RGD", "MGD", "CTD_rat", "CTD_mouse",
        "PREDICTED", "ALL"), verbose=FALSE, warnings=TRUE)
    standardGeneric("enrichDISf")
)

#' Function to enrich in anatomical terms provided by Bgee Database
#' 
#' TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO TODO
#' 
#' @param object Object of class \code{ResultSet}.
#' @param rid (default \code{1}) If given \code{ResultSet} hat has more than 
#' one result, id to select the result where the enrichment will be performed.
#' @param fData.tag (default \code{1}) Identifier of the annotation to be
#' used to mach the probes of \code{ResultSet} with gene identifier. In other
#' words, index of the element in \code{fData} that matches the omic dataset 
#' to be used in the enrichment.
#' @param sel.pval (default \code{"adj.P.Val"}) Name of the P.Value columns to 
#' filter datset's probes given \code{th.pval}.
#' @param th.pval (default \code{0.01}) Threshold used to include a probe (CpG)
#' into the enrichment analysis.
#' @param sel.feature (default \code{"genes"}) Name containing gene identifier
#' in \code{ResultSet}'s \code{fData}.
#' @param feature.null (default \code{""}) String identifier for non 
#' gene-assigned probes.
#' @param bgee.data (default \code{"affymetrix"}) Type of data requested to 
#' Bgee to perform the functinal enrichment. Can take any of them or more than
#' one of the possibles options.
#' @param statistic (default \code{"fisher"}) Statistic used to test the
#' enrichment in Bgee anatomical categories throught \code{topGO} R package.
#' Check \code{whichTests()} for available tests.
#' @param verbose (default \code{FALSE}) If set to \code{TRUE} messages
#' indicating the steps done by the method are raised.
#' @param warnings (default \code{TRUE}) If set to \code{FALSE} important
#' notes raised by the method are hide.
#' @export enrichANAT
#' @importClassesFrom rexposome ResultSet
#' @note 
#' Komljenovic A, Roux J, Robinson-Rechavi M, Bastian F;
#' BgeeDB, an R package for retrieval of curated expression datasets and for 
#' gene list enrichment tests. 
#' F1000Research (2016).
#' 
#' Bastian F, Parmentier G, Roux J, Moretti S, Laudet V, Robinson-Rechavi M;
#' Bgee: Integrating and Comparing Heterogeneous Transcriptome Data Among 
#' Species; Data Integration Life Sci. Lecture Notes in Computer Science (2008)
setGeneric("enrichANAT", function(object, rid=1, fData.tag=1,
      sel.pval=c("P.Value", "adj.P.Val"), th.pval=0.01, sel.feature="genes",
      feature.null="", bgee.data=c("rna_seq","affymetrix","est","in_situ"), 
      statistic="fisher", verbose=FALSE, warnings=TRUE)
    standardGeneric("enrichANAT")
)

#' Function to enrich results from differential studies on KEGG and GO 
#' databases
#' 
#' This method takes a \code{ResultSet} and performs an enrichment on KEGG
#' pathways or GO terms in base to a given set of genes. Argument 
#' \code{database} allows to select the database where the query will be 
#' performed. \code{sel.pval} allows to select the column of the resulting
#' table in \code{ResultSet} with the p-value to be used in the process.
#' Argument \code{th.pval} is used to exclude features using \code{sel.pval}
#' column. Argument \code{sel.feature} indicates the column that has the
#' gene associated to the feature (row in the table of results).
#' \code{feature.id} indicates the coding system for the \code{sel.feature}
#' column.
#' 
#' @param object Object of class \code{ResultSet}.
#' @param rid (default \code{1}) If given \code{ResultSet} that has more than 
#' one result, id to select the result where the enrichment will be performed.
#' @param fData.tag (default \code{1}) Identifier of the annotation to be
#' used to mach the probes of \code{ResultSet} with gene identifier. In other
#' words, index of the element in \code{fData} that matches the omic dataset 
#' to be used in the enrichment.
#' @param version (default \code{"hg18"}) Can takes \code{"hg18"} and 
#' \code{"hg19"} and corresponds to the genome version to use for enrichment 
#' in KEGG/GO.
#' @param database (default \code{"KEGG"}) Can takes values \code{"KEGG"}, 
#' \code{"GO:CC"}, \code{"GO:BP"} and \code{"GO:MF"}; corresponding to 
#' the available enrichments.
#' @param sel.pval (default \code{"adj.P.Val"}) Name of the P.Value columns to 
#' filter datset's probes given \code{th.pval}.
#' @param th.pval (default \code{0.01}) Threshold used to include a probe (CpG)
#' into the enrichment analysis.
#' @param sel.feature (default \code{"genes"}) Name containing gene identifier
#'  in \code{ResultSet}'s \code{fData}.
#' @param feature.id (default \code{"geneSymbol"}) Name of the format of the 
#' gene identifier. Run \link{supportedGeneIDs} from \code{geneLenDataBase} for 
#' more information. This must match with the identifier in \code{ResultSet}'s 
#' \code{fData}.
#' @param feature.null (default \code{""}) String identifier for non 
#' gene-assigned probes.
#' @param plot.fit (default \code{FALSE}) If set to \code{TRUE} calls plot 
#' \link{PWF} with default values.
#' @param verbose (default \code{FALSE}) If set to \code{TRUE} messages
#' indicating the steps done by the method are raised.
#' @param warnings (default \code{TRUE}) If set to \code{FALSE} important
#' notes raised by the method are hide.
#' @export enrichPATH
#' @importClassesFrom rexposome ResultSet
#' @note 
#' The Gene Ontology Consortium; Gene Ontology Consortium: going forward;
#' Nucl Acids Res (2015).
#' 
#' Ashburner et al.; Gene ontology: tool for the unification of biology;
#' Nat Genet (2000).
#' 
#' Kanehisa M, Furumichi M, Tanabe M, Sato Y, Morishima K; 
#' KEGG: new perspectives on genomes, pathways, diseases and drugs. 
#' Nucleic Acids Res (2017).
#' 
#' Kanehisa M, Sato Y, Kawashima M, Furumichi M, Tanabe M; 
#' KEGG as a reference resource for gene and protein annotation. 
#' Nucleic Acids Res (2016).
#' 
#' Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. 
#' Nucleic Acids Res. 28, 27-30 (2000).
setGeneric("enrichPATH", function(object, rid=1, fData.tag=1,
                                  version=c("hg18", "hg19"),
                                  database=c("KEGG", "GO:CC", "GO:BP", "GO:MF"), 
                                  sel.pval=c("P.Value", "adj.P.Val"), 
                                  th.pval=0.01, sel.feature="genes",
                                  feature.id="geneSymbol", feature.null="",
                                  plot.fit=FALSE, verbose=FALSE, 
                                  warnings=TRUE)
    standardGeneric("enrichPATH")
)

#' Function to enrich results from differential studies on KEGG and GO 
#' databases
#' 
#' This method takes a \code{ResultSet} and performs an enrichment on KEGG
#' pathways or GO terms in base to a given set of genes. Argument 
#' \code{database} allows to select the database where the query will be 
#' performed. \code{sel.pval} allows to select the column of the resulting
#' table in \code{ResultSet} with the p-value to be used in the process.
#' Argument \code{th.pval} is used to exclude features using \code{sel.pval}
#' column. Argument \code{sel.feature} indicates the column that has the
#' gene associated to the feature (row in the table of results).
#' \code{feature.id} indicates the coding system for the \code{sel.feature}
#' column.
#' 
#' @param object Object of class \code{ResultSet}.
#' @param family Family of exposures used to perform the enrichment analysis.
#' @param fData.exp (default \code{2}) Identifier of the exposures decrition 
#' to be used to mach the exposures of \code{ResultSet} with their families.
#' @param fData.omic (default \code{1}) Identifier of the annotation to be
#' used to mach the probes of \code{ResultSet} with gene identifier. In other
#' words, index of the element in \code{fData} that matches the omic dataset 
#' to be used in the enrichment.
#' @param version (default \code{"hg18"}) Can takes \code{"hg18"} and 
#' \code{"hg19"} and corresponds to the genome version to use for enrichment 
#' in KEGG/GO.
#' @param database (default \code{"KEGG"}) Can takes values \code{"KEGG"}, 
#' \code{"GO:CC"}, \code{"GO:BP"} and \code{"GO:MF"}; corresponding to 
#' the available enrichments.
#' @param sel.pval (default \code{"adj.P.Val"}) Name of the P.Value columns to 
#' filter datset's probes given \code{th.pval}.
#' @param th.pval (default \code{0.01}) Threshold used to include a probe (CpG)
#' into the enrichment analysis.
#' @param sel.feature (default \code{"genes"}) Name containing gene identifier
#'  in \code{ResultSet}'s \code{fData}.
#' @param feature.id (default \code{"geneSymbol"}) Name of the format of the 
#' gene identifier. Run \link{supportedGeneIDs} from \code{geneLenDataBase} for 
#' more information. This must match with the identifier in \code{ResultSet}'s 
#' \code{fData}.
#' @param feature.null (default \code{""}) String identifier for non 
#' gene-assigned probes.
#' @param plot.fit (default \code{FALSE}) If set to \code{TRUE} calls plot 
#' \link{PWF} with default values.
#' @param verbose (default \code{FALSE}) If set to \code{TRUE} messages
#' indicating the steps done by the method are raised.
#' @param warnings (default \code{TRUE}) If set to \code{FALSE} important
#' notes raised by the method are hide.
#' @export enrichPATHf
#' @importClassesFrom rexposome ResultSet
#' @note 
#' The Gene Ontology Consortium; Gene Ontology Consortium: going forward;
#' Nucl Acids Res (2015).
#' 
#' Ashburner et al.; Gene ontology: tool for the unification of biology;
#' Nat Genet (2000).
#' 
#' Kanehisa M, Furumichi M, Tanabe M, Sato Y, Morishima K; 
#' KEGG: new perspectives on genomes, pathways, diseases and drugs. 
#' Nucleic Acids Res (2017).
#' 
#' Kanehisa M, Sato Y, Kawashima M, Furumichi M, Tanabe M; 
#' KEGG as a reference resource for gene and protein annotation. 
#' Nucleic Acids Res (2016).
#' 
#' Kanehisa, M. and Goto, S.; KEGG: Kyoto Encyclopedia of Genes and Genomes. 
#' Nucleic Acids Res. 28, 27-30 (2000).
setGeneric("enrichPATHf", function(object, family, fData.exp=2, fData.omic=1,
                                  version=c("hg18", "hg19"),
                                  database=c("KEGG", "GO:CC", "GO:BP", "GO:MF"), 
                                  sel.pval=c("P.Value", "adj.P.Val"), 
                                  th.pval=0.01, sel.feature="genes",
                                  feature.id="geneSymbol", feature.null="",
                                  plot.fit=FALSE, verbose=FALSE, 
                                  warnings=TRUE)
    standardGeneric("enrichPATHf")
)

#' Function to enrich GWAS results on GWAS Catalog
#' 
#' This method takes a \code{ResultSet} and performs a query to GWAS Catalog
#' data giving a set of SNPs and retrieving the associated results from the
#' database.
#' 
#' @param object Object of class \code{ResultSet}.
#' @param version Version of the GWAS Catalog to be used. This argument can
#' take \code{"hg19"}, \code{"hg38"} and \code{"current"}. If it is set to 
#' \code{"current"} an update snapshot from GWAS Catalog will be downloaded,
#' otherwise stored version in \code{gwascat} will be used.
#' @param verbose (default \code{FALSE}) If set to \code{TRUE} messages
#' indicating the steps done by the method are raised.
#' @param warnings (default \code{TRUE}) If set to \code{FALSE} important
#' notes raised by the method are hide.
#' @export enrichGWAS
#' @importClassesFrom rexposome ResultSet
#' @seealso See \link{makeCurrentGwascat} for more information on updated
#' version of GWAS Catalog.
#' @note 
#' Welter D, MacArthur J, Morales J, Burdett T, Hall P, Junkins H, Klemm A, 
#' Flicek P, Manolio T, Hindorff L, Parkinson H;
#' The NHGRI GWAS Catalog, a curated resource of SNP-trait associations.
#' Nucleic Acids Research (2014)
setGeneric("enrichGWAS", function(object, version=c("hg19", "hg38", "current"), 
                                  verbose=FALSE, warnings=TRUE)
    standardGeneric("enrichGWAS")
)