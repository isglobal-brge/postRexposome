% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_AllMethods.R
\name{enrichDIS}
\alias{enrichDIS}
\title{Method to enrich results from differential studies on DisGeNET database}
\usage{
enrichDIS(object, rid = 1, coef = 2, contrast = 1, fData.tag = 2,
  sel.pval = "adj.P.Val", th.pval = 0.01, sel.feature = "gene",
  feature.null = c("", "---"), translate.from = NA, database = "CURATED",
  verbose = FALSE, warnings = TRUE)
}
\arguments{
\item{object}{Object of class \code{ResultSet}.}

\item{rid}{(default \code{1}) If given \code{ResultSet} hat has more than 
one result, id to select the result where the enrichment will be performed.}

\item{fData.tag}{(default \code{1}) Identifier of the annotation to be
used to mach the probes of \code{ResultSet} with gene identifier. In other
words, index of the element in \code{fData} that matches the omic dataset 
to be used in the enrichment.}

\item{sel.pval}{(default \code{"adj.P.Val"}) Name of the P.Value columns to 
filter datset's probes given \code{th.pval}.}

\item{th.pval}{(default \code{0.01}) Threshold used to include a probe (CpG)
into the enrichment analysis.}

\item{sel.feature}{(default \code{"genes"}) Name containing gene identifier
in \code{ResultSet}'s \code{fData}.}

\item{feature.null}{(default \code{""}) String identifier for non 
gene-assigned probes.}

\item{database}{(default \code{"CURATED"}) Name of the version of the
DisGeNET's databse where the query will be computed. Check 
\code{?disgenetGene} for a description of each name.}

\item{verbose}{(default \code{FALSE}) If set to \code{TRUE} messages
indicating the steps done by the method are raised.}

\item{warnings}{(default \code{TRUE}) If set to \code{FALSE} important
notes raised by the method are hide.}
}
\description{
This method takes the results from a \code{ResultSet} and uses the genes
associated to the probes to perform a query in DisGeNET and obtain the
diseases associated to the given genes.
}
\note{
Pinero J, Queralt-Rosinach N, Bravo A, Deu-Pons J, Bauer-Mehren A, Baron M, 
Sanz F, Furlong LI;
DisGeNET: a discovery platform for the dynamical exploration of human 
diseases and their genes.; Database (2015).

Queralt-Rosinach N, Pinero J, Bravo A, Sanz F, Furlong LI;
DisGeNET-RDF: harnessing the innovative power of the Semantic Web to explore
the genetic basis of diseases Bioinformatics; Bioinformatics (2016).
}
\examples{
data(air)
ans <- enrichDIS(air, rid="NO2", sel.pval="P.Value", warnings=FALSE,
    sel.feature="gene_assignment", translate.from="REFSEQ")
ans
}
