% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/Nickase-class.R
\docType{class}
\name{nickaseName}
\alias{nickaseName}
\alias{nickaseName<-}
\alias{nickingStrand}
\alias{nickingStrand<-}
\alias{Nickase-class}
\alias{Nickase}
\alias{show,Nickase-method}
\alias{nickaseName,Nickase-method}
\alias{nickaseName<-,Nickase-method}
\alias{nickingStrand,Nickase-method}
\alias{nickingStrand<-,Nickase-method}
\alias{weights,Nickase-method}
\alias{weights<-,Nickase-method}
\alias{isCutting,Nickase-method}
\alias{motifs,Nickase-method}
\alias{motifLength,Nickase-method}
\alias{cutSites,Nickase-method}
\title{An S4 class to represent a nickase}
\usage{
nickaseName(object)

nickaseName(object) <- value

nickingStrand(object)

nickingStrand(object) <- value

Nickase(
  nickaseName,
  nickingStrand = c("original", "opposite"),
  motifs = NULL,
  cutSites = NULL,
  weights = rep(1, length(motifs)),
  metadata = list()
)

\S4method{show}{Nickase}(object)

\S4method{nickaseName}{Nickase}(object)

\S4method{nickaseName}{Nickase}(object) <- value

\S4method{nickingStrand}{Nickase}(object)

\S4method{nickingStrand}{Nickase}(object) <- value

\S4method{weights}{Nickase}(object, expand = FALSE)

\S4method{weights}{Nickase}(object) <- value

\S4method{isCutting}{Nickase}(object)

\S4method{motifs}{Nickase}(
  object,
  primary = FALSE,
  strand = c("+", "-"),
  expand = FALSE,
  as.character = FALSE
)

\S4method{motifLength}{Nickase}(object)

\S4method{cutSites}{Nickase}(object, combine = TRUE)
}
\arguments{
\item{object}{\linkS4class{Nickase} object.}

\item{value}{New value to pass to the setter functions.}

\item{nickaseName}{Name of the nickase.}

\item{nickingStrand}{String specifying with strand with respect
to the motif sequence (5' to 3') is nicked.
Must be either "original" (default) or "opposite".}

\item{motifs}{Character vector of recognition sequence motifs
written from 5' to 3' written in Rebase convention.
If the point of cleavage has been determined, the
precise site is marked with ^. Only letters in the
IUPAC code are accepted. For nickases that cleave
away from their recognition sequence, the cleavage
sites are indicated in parentheses. See details for
more information.}

\item{cutSites}{Vector specifying the cleavage coordinates relative
to the first nucleotide of the motif sequence. Each column corresponds
to a motif specified in the \code{motifs} slot.}

\item{weights}{Optional numeric vector specifying relative weights
for the recognition motifs to specify cleavage probabilities.}

\item{metadata}{Optional list providing global metadata information.}

\item{expand}{Should sequences be expanded to only contain ATCG nucleotides?
FALSE by default.}

\item{primary}{Should only the motif with the highest weight be returned?
FALSE by default. Only relevant if weights are stored in the 
\linkS4class{Nickase} object.}

\item{strand}{Strand to allow reverse complementation of the motif.
"+" by default.}

\item{as.character}{Should the motif sequences be returned as a 
character vector? FALSE by default.}

\item{combine}{Should only unique values be considered?
TRUE by default.}
}
\value{
A Nickase object
}
\description{
An S4 class to represent a nickase
}
\section{Functions}{
\itemize{
\item \code{Nickase()}: Create a \linkS4class{Nickase} object

}}
\section{Slots}{

\describe{
\item{\code{nickaseName}}{Name of the nickase}

\item{\code{motifs}}{DNAStringSet of recognition sequence motifs
written from 5' to 3'.}

\item{\code{nickingStrand}}{String specifying with strand with respect
to the motif sequence (5' to 3') is nicked.
Must be either "original" (default) or "opposite".}

\item{\code{cutSites}}{Vector specifying the cleavage coordinates relative
to the first nucleotide of the motif sequence.
Each column corresponds to a motif specified in the \code{motifs} slot.}

\item{\code{weights}}{Optional numeric vector specifying relative weights
for the motifs corresponding to cleavage probabilities.}

\item{\code{metadata}}{Optional string providing a description of the nickase.}
}}

\section{Constructors}{

    Use the constructor \code{link{Nickase}} to create a Nickase object.
}

\section{Accessors}{

\describe{
    \item{\code{nickaseName}:}{To get the name of the nickase.} 
    \item{\code{nickingStrand}:}{To get the nicking strand.} 
    \item{\code{metadata}:}{To get the metadata list of the nickase} 
    \item{\code{motifs}:}{To get the recognition mofif
         nucleotide sequences.}
    \item{\code{weights}:}{To get nickase weights.} 
    \item{\code{cutSites}:}{To get nickase cut sites.} 
}
}

\examples{
Nb.BsmI <- Nickase("Nb.BsmI",
                   motifs=c("GAATG^C"),
                   nickingStrand="opposite",
                   metadata=list(description="Nb.BsmI nicking enzyme."))

}
\seealso{
See the \linkS4class{CrisprNickase} for CRISPR-specific nickases.
}
