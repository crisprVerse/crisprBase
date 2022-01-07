#' An S4 class to represent a nickase
#' 
#' @slot nickaseName Name of the nickase
#' @slot motifs DNAStringSet of recognition sequence motifs
#'           written from 5' to 3'.
#' @slot nickingStrand String specifying with strand with respect
#'     to the motif sequence (5' to 3') is nicked.
#'     Must be either "original" (default) or "opposite".
#' @slot cutSites Vector specifying the cleavage coordinates relative
#'     to the first nucleotide of the motif sequence.
#'     Each column corresponds to a motif specified in the \code{motifs} slot. 
#' @slot weights Optional numeric vector specifying relative weights
#'           for the motifs corresponding to cleavage probabilities. 
#' @slot metadata Optional string providing a description of the nickase.
#' 
#' @seealso See the \linkS4class{CrisprNickase} for CRISPR-specific nickases.
#' 
#' @section Constructors:
#'     Use the constructor \code{link{Nickase}} to create a Nickase object.
#' 
#' @section Accessors:
#' \describe{
#'     \item{\code{nickaseName}:}{To get the name of the nickase.} 
#'     \item{\code{nickingStrand}:}{To get the nicking strand.} 
#'     \item{\code{metadata}:}{To get the metadata list of the nickase} 
#'     \item{\code{motifs}:}{To get the recognition mofif
#'          nucleotide sequences.}
#'     \item{\code{weights}:}{To get nickase weights.} 
#'     \item{\code{cutSites}:}{To get nickase cut sites.} 
#' }
#' 
#' @examples
#' Nb.BsmI <- Nickase("Nb.BsmI",
#'                    motifs=c("GAATG^C"),
#'                    nickingStrand="original",
#'                    metadata=list(description="Nb.BsmI nicking enzyme."))
#' 
#' @return A Nickase object
#' @export
#' @importFrom Biostrings DNAStringSet
#' @importClassesFrom S4Vectors Annotated
setClass("Nickase",
    contains = "Annotated",
    slots = c(
        nickaseName = "character", 
        nickingStrand = "character", 
        motifs = "DNAStringSet",
        cutSites = "vector",
        weights = "numeric"),
    prototype = list(
        nickaseName = NA_character_,
        nickingStrand = NA_character_,
        motifs = NULL,
        cutSites = NULL,
        weights = NA_real_)
)




#' @describeIn Nickase Create a \linkS4class{Nickase} object
#' @param nickaseName Name of the nickase.
#' @param nickingStrand String specifying with strand with respect
#'     to the motif sequence (5' to 3') is nicked.
#'     Must be either "original" (default) or "opposite".
#' @param motifs Character vector of recognition sequence motifs
#'           written from 5' to 3' written in Rebase convention.
#'           If the point of cleavage has been determined, the
#'           precise site is marked with ^. Only letters in the
#'           IUPAC code are accepted. For nickases that cleave
#'           away from their recognition sequence, the cleavage
#'           sites are indicated in parentheses. See details for
#'           more information. 
#' @param cutSites Vector specifying the cleavage coordinates relative
#'     to the first nucleotide of the motif sequence. Each column corresponds
#'     to a motif specified in the \code{motifs} slot. 
#' @param weights Optional numeric vector specifying relative weights
#'           for the recognition motifs to specify cleavage probabilities. 
#' @param metadata Optional list providing global metadata information.
#' @export
#' @importFrom S4Vectors metadata metadata<-
Nickase <- function(nickaseName,
                    nickingStrand=c("original", "opposite"),
                    motifs = NULL,
                    cutSites = NULL,
                    weights = rep(1, length(motifs)),
                    metadata = list()
){
    nickingStrand <- match.arg(nickingStrand)
    nickaseName <- as.character(nickaseName)
    weights <- as.numeric(weights)

    if (is(motifs, "DNAStringSet")){    
        if (is.null(cutSites)){
            stop("When 'motifs' is provided as a DNAStringset, ",
                 "'cutSites' must be provided. ")
        }
        sequences <- motifs
    } else if (is.character(motifs) & !is.null(cutSites)){
        .checkDNAAlphabet(motifs)
        sequences <- DNAStringSet(motifs)
    } else if (is.character(motifs) & is.null(cutSites)){
        .checkRebaseMotifs(motifs, type="nickase")
        sequences <- .extractSequencesFromRebaseMotifs(motifs)
        rebaseStrand <- ifelse(nickingStrand=="original", "fwd", "rev")
        cutSites  <- .extractCutSitesFromRebaseMotifs(motifs,
                                                      type="nickase",
                                                      strand=rebaseStrand)
    } else {
        stop("'motifs' must be either a character vector or a ",
             "DNAStringSet.")
    }
    
    if (!is.list(metadata)){
        stop("metadata must be a list.")
    }
    out <- new("Nickase",
               nickaseName=nickaseName,
               nickingStrand=nickingStrand,
               motifs=sequences,
               cutSites=cutSites,
               weights=weights)
    metadata(out) <- metadata
    return(out)
}





setValidity("Nickase", function(object){
    out <- TRUE
    if (length(nickaseName(object))!=1){
       out <- "Slot nickaseName must be a character vector of length 1."
    } 
  
    if (length(weights(object))>1){
        if (length(weights(object)) != length(motifs(object))){
            out <- "Slots weights and motifs must be of the same length."
        } 
    } 

    sites <- cutSites(object, combine=FALSE)
    if (length(sites) != length(motifs(object))){
        out <- "Length in cutSites must be equal to the length of motifs."
    }

    return(out)
})




#' @rdname Nickase-class
#' @param object \linkS4class{Nickase} object.
setMethod("show", "Nickase", function(object){
    len <- length(metadata(object))
    cat(paste0("Class: ", is(object)[[1]]), "\n",
      "  Name: ", nickaseName(object), "\n",
      "  Nicking strand: ", nickingStrand(object), "\n",
      "  Metadata: list of length ", len, "\n",
      "  Motifs: ", .printVectorNicely(motifs(object)), "\n",
      "  Weights: ", .printVectorNicely(weights(object)), "\n",
      sep = "")
})



#' @rdname Nickase-class
#' @export
setMethod("nickaseName", "Nickase", 
    function(object){
    return(object@nickaseName)
})

#' @rdname Nickase-class
#' @param value New value to pass to the setter functions.
#' @export
setMethod("nickaseName<-", "Nickase", 
    function(object, value){
    object@nickaseName <- as.character(value)
    return(object)
})


#' @rdname Nickase-class
#' @export
setMethod("nickingStrand", "Nickase", 
    function(object){
    return(object@nickingStrand)
})

#' @rdname Nickase-class
#' @param value New value to pass to the setter functions.
#' @export
setMethod("nickingStrand<-", "Nickase", 
    function(object, value){
    if (!value %in% c("original", "opposite")){
        stop("value must be either 'original' or 'opposite'.")
    }
    object@nickingStrand <- as.character(value)
    return(object)
})



#' @rdname Nickase-class
#' @export
setMethod("weights", "Nickase", 
    function(object,
             expand=FALSE
){
    ws <- object@weights
    names(ws) <- motifs(object)
    if (expand){
        motifs <- .expandMotifs(names(ws))
        ws <- ws[names(motifs)]
        names(ws) <- as.character(motifs)
    }
    return(ws)
})


#' @rdname Nickase-class
#' @export
setMethod("weights<-", "Nickase", 
    function(object, value){
    seqs <- object@motifs
    if (is.null(value)){
        object@motifs <- NA
    } else if (length(value)==1){
        object@motifs <- value
    } else if (length(value)==length(seqs)){
        object@motifs <- value
    } 
    object@weights <- value
    return(object)
})



#' @rdname Nickase-class
#' @export
setMethod("isCutting", "Nickase", 
    function(object){
    cuts <- cutSites(object)
    !all(is.na(c(cuts)))
})



#' @rdname Nickase-class
#' @param primary Should only the motif with the highest weight be returned?
#'     FALSE by default. Only relevant if weights are stored in the 
#'     \linkS4class{Nickase} object.
#' @param expand Should sequences be expanded to only contain ATCG nucleotides?
#'     FALSE by default.
#' @param strand Strand to allow reverse complementation of the motif.
#'     "+" by default. 
#' @param as.character Should the motif sequences be returned as a 
#'     character vector? FALSE by default.
#' @importFrom Biostrings reverseComplement
#' @export
setMethod("motifs", "Nickase", 
    function(object,
             primary=FALSE,
             strand=c("+", "-"),
             expand=FALSE,
             as.character=FALSE
){
    strand <- match.arg(strand)
    motifs  <- object@motifs
    weights <- object@weights

    if (length(weights)>1 & primary){
        w <- max(weights, na.rm=TRUE)
        motifs <- motifs[which(weights==w)]
    }

    if (strand=="-"){
        motifs <- reverseComplement(motifs)
    }
    if (expand){
        motifs <- .expandMotifs(motifs)
    }
    if (as.character){
        motifs <- as.character(motifs)
    }
    return(motifs)
})





#' @rdname Nickase-class
#' @param object \linkS4class{Nickase} object.
#' @importFrom BiocGenerics width
#' @export
setMethod("motifLength",
          "Nickase",function(object) {
    seqs <- motifs(object, primary=TRUE)[1]
    names(seqs) <- NULL
    return(BiocGenerics::width(seqs))
})



#' @rdname Nickase-class
#' @param combine Should only unique values be considered?
#'     TRUE by default. 
#' @export
setMethod("cutSites",
          "Nickase", function(object, combine=TRUE)
{
    sites  <- object@cutSites
    motifs <- object@motifs
    names(sites) <- motifs
    if (combine){
        choices <- unique(sites)
        if (length(choices)==1){
            sites <- choices
            names(sites) <- NULL
        }
    }
    return(sites)
})








