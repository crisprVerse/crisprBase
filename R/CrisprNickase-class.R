#' An S4 class to represent a CRISPR nickase.
#' 
#' @slot pam_side String specifying the side of the PAM sequence
#'     with respect to the protospacer sequence. Must be either 
#'     '3prime' (e.g. SpCas9) or '5prime' (e.g. AsCas12a)
#' @slot spacer_length Integer specifying the
#'     length of the spacer sequence
#' @slot spacer_gap Integer specifying the length (in nucleotides)
#'     between the spacer sequence and the PAM sequence
#'     (e.g. 0 for SpCas9 and AsCas12a).
#' 
#' @section Constructors:
#'     Use the constructor \code{link{CrisprNickase}} to create
#'         a CrisprNickase object.
#' 
#' @section Accessors:
#' \describe{
#'     \item{\code{nickaseName}:}{To get the name of the CRISPR nickase.} 
#'     \item{\code{spacerLength}:}{To return the length of the
#'         spacer sequence.}
#'     \item{\code{targetLength}:}{To return the length of the
#'         target sequence (protospacer + pam).}
#'     \item{\code{pamLength}:}{To return the length of the PAM sequence.}
#'     \item{\code{pamSide}:}{To return the side of the PAM sequence with
#'         respect to the spacer sequence.}
#'     \item{\code{spacerGap}:}{To return the length of the gap between the
#'         PAM and spacer sequences.}
#'     \item{\code{pams}:}{To return the list of PAM sequences.}
#' }
#' 
#' @section Setters:
#' \describe{
#'     \item{\code{spacerGap<-}:}{To change the length of the gap between the
#'         PAM and spacer sequences.}
#'     \item{\code{pamSide<-}:}{To change the side of the PAM sequence
#'         with respect to the protospacer sequence.}
#'     \item{\code{spacerLength<-}:}{To change the length of the spacer
#'         sequence.}
#' }
#' 
#' @section Utility functions for genomic arithmetics:
#' \describe{
#'     \item{\code{pamIndices}:}{To return the relative coordinates of the
#'         PAM sequence within the protospacer sequence.}
#'     \item{\code{spacerIndices}:}{To return the relatiive coordinates of
#'         the spacer sequence within the protospacer sequence.}
#' }
#' @examples
#' Cas9D10A <- CrisprNickase("Cas9D10A",
#'                           nickingStrand="opposite",
#'                           pams=c("(3)NGG", "(3)NAG", "(3)NGA"),
#'                           weights=c(1, 0.2593, 0.0694),
#'                           metadata=list(description="D10A-mutated  Streptococcus
#'                                         pyogenes Cas9 (SpCas9) nickase"),
#'                           pam_side="3prime",
#'                           spacer_length=20)
#' 
#' Cas9H840A <- CrisprNickase("Cas9H840A",
#'                            nickingStrand="original",
#'                            pams=c("(3)NGG", "(3)NAG", "(3)NGA"),
#'                            weights=c(1, 0.2593, 0.0694),
#'                            metadata=list(description="H840A-mutated  Streptococcus
#'                                          pyogenes Cas9 (SpCas9) nickase"),
#'                            pam_side="3prime",
#'                            spacer_length=20)
#' 
#' @return A CrisprNickase object
#' 
#' @export
setClass("CrisprNickase", 
    contains = "Nickase", 
    slots = c(
        pam_side = "character",
        spacer_gap = "integer",
        spacer_length = "integer"),
    prototype = list(
        pam_side = NA_character_,
        spacer_gap = 0L,
        spacer_length =NA_integer_
    )
)



#' @describeIn CrisprNickase Create a \linkS4class{CrisprNickase} object
#' @param nickaseName Name of the CRISPR nickase.
#' @param nickingStrand String specifying with strand with respect
#'     to the motif sequence (5' to 3') is nicked.
#'     Must be either "original" (default) or "opposite".
#' @param pams Character vector of PAM sequence motifs
#'           written from 5' to 3. If the point of cleavage has
#'           been determined, the precise site is marked with ^.
#'           Only letters in the IUPAC code are accepted.
#'           For nickases that cleave away from their
#'           recognition sequence, the cleavage sites are indicated
#'           in parentheses. See details for more information.
#' @param weights Optional numeric vector specifying relative weights
#'           of the PAM sequences to specify cleavage probabilities. 
#' @param metadata Optional list providing global metadata information.
#' @param pam_side String specifying the side of the PAM sequence
#'     sequence with respect to the protospacer sequence. Must be either 
#'     '3prime' (e.g. Cas9) or '5prime' (e.g. Cas12a)
#' @param spacer_length Integer specifying the length of the spacer sequence
#' @param spacer_gap Integer specifying the length (in nucleotides) between
#'     the spacer sequence and the PAM sequence (e.g. 0 for Cas9 and Cas12a).
#' @export
CrisprNickase <- function(nickaseName,
                          nickingStrand = c("original", "opposite"),
                          pams = NA_character_,
                          weights = rep(1, length(pams)),
                          metadata = list(),
                          pam_side = NA_character_,
                          spacer_gap = 0L,
                          spacer_length = NA_integer_
){
    nickingStrand <- match.arg(nickingStrand)
    nick <- Nickase(nickaseName=nickaseName,
                    nickingStrand=nickingStrand,
                    motifs=pams,
                    weights=weights,
                    metadata=metadata)
    new("CrisprNickase",
        nick,
        pam_side = as.character(pam_side),
        spacer_gap = as.integer(spacer_gap),
        spacer_length = as.integer(spacer_length))
}



#' @rdname CrisprNickase-class
#' @export
setMethod("show", "CrisprNickase", function(object) {
    len <- length(metadata(object))
    pams.line <- ("  PAMs: ")
    pams.line.distance <- "    Distance from PAM: "
    pam.side.line <- "  PAM side: "
    
    cat(paste0("Class: ", is(object)[[1]]), "\n",
        "  Name: ", nickaseName(object), "\n",
        "  Nicking strand: ", nickingStrand(object), "\n",
        "  Metadata: list of length ", len, "\n",
        pams.line, .printVectorNicely(motifs(object)), "\n",
        "  Weights: ", .printVectorNicely(weights(object)), "\n",
        "  Spacer length: ",  spacerLength(object), "\n",
        pam.side.line,  pamSide(object), "\n",
        pams.line.distance, spacerGap(object), "\n",
        "  Prototype protospacers: ",
        .printVectorNicely(prototypeSequence(object, primary=FALSE)),
        "\n",
        sep = "")
})






setValidity("CrisprNickase", function(object) {
    out <- TRUE
    if (length(pamSide(object))!=1){
        out <- "Slot pam_side must be a character vector of length 1."
    } 
    if (!pamSide(object) %in% c("5prime", "3prime")){
        out <- "Slot pam_side must be either '5prime' or '3prime'."
    } 
    if (length(spacerLength(object))!=1){
        out <- "Slot spacer_length must be an integer vector of length 1."
    } 
    if (spacerLength(object)<0){
        out <- "Slot spacer_length must be a positive integer."
    } 
    if (length(spacerGap(object))!=1){
        out <- "Slot spacer_gap must be an integer vector of length 1."
    } 
    if (spacerGap(object)<0){
        out <- "Slot spacer_gap must be either 0 or a positive integer."
    } 
    return(out)
})











############################################################
##############       Getter and setters   ##################
#' @rdname CrisprNickase-class
#' @param object \linkS4class{CrisprNickase} object.
#' @export
setMethod("pamLength",
          "CrisprNickase",function(object){
    len <- motifLength(object)[1]
    names(len) <- NULL
    return(len)
})

#' @rdname CrisprNickase-class
#' @export
setMethod("spacerLength", "CrisprNickase", function(object){
    return(object@spacer_length)
})


#' @rdname CrisprNickase-class
#' @param value For \code{spacerLength<-} and \code{gapLength<-}, must be 
#'     a non-negative integer. For \code{pamSide}, must be either
#'     '5prime' or '3prime'.
#' 
#' @export
setMethod("spacerLength<-", "CrisprNickase", function(object, value){
    value <- .validateNonNegativeInteger(value)
    object@spacer_length <- value
    return(object)
})


#' @rdname CrisprNickase-class
#' @export
setMethod("pamSide",
          "CrisprNickase",function(object){
    return(object@pam_side)
})

#' @rdname CrisprNickase-class
#' @export
setMethod("pamSide<-", "CrisprNickase", function(object, value){
    if (!value %in% c("5prime", "3prime")){
        stop("value must be either '5prime' or '3prime'")
    }
    object@pam_side <- value
    return(object)
})





#' @rdname CrisprNickase-class
#' @export
setMethod("spacerGap",
          "CrisprNickase",function(object){
    return(object@spacer_gap)
})


#' @rdname CrisprNickase-class
#' @export
setMethod("spacerGap<-", "CrisprNickase", function(object, value){
    value <- .validateNonNegativeInteger(value)
    object@spacer_gap <- value
    return(object)
})

#' @rdname CrisprNickase-class
#' @export
setMethod("hasSpacerGap",
          "CrisprNickase",function(object){
    object@spacer_gap>0
})




#' @rdname CrisprNickase-class
#' @export
setMethod("targetLength",
          "CrisprNickase",function(object){
    pamLength(object) + spacerLength(object) + spacerGap(object)
})






#' @rdname CrisprNickase-class
#' @param primary Should only the PAM sequence with the heighest weight
#'     be returned? If no cleavage weights are stored in the
#'     \linkS4class{CrisprNickase} object, all sequences are returned.
#'     TRUE by default.
#' @param ignore_pam Should all possible k-mer sequences for a given PAM length
#'     be returned, irrespetively of the PAM sequence motifs stored in the
#'     \linkS4class{CrisprNickase} object? FALSE by default.
#' @param as.character Should the PAM sequences be returned as a 
#'     character vector? FALSE by default.
#' @export
setMethod("pams", "CrisprNickase",
    function(object,
             primary=TRUE,
             ignore_pam=FALSE,
             as.character=FALSE
){
    pam_length <- pamLength(object)
    if (ignore_pam){
        motif <- paste0(rep("N", pam_length), collapse="")
        pams <- .expandMotifs(motif)
    } else {
        pams <- motifs(object,
                       expand=TRUE,
                       primary=primary)
    } 
    if (as.character){
        pams <- as.character(pams)
    }
    return(pams)
})









############################################################
##############        Coordinates         ##################
#' @rdname CrisprNickase-class
#' @export
setMethod("pamIndices", "CrisprNickase",
    function(object){
    indices  <- seq_len(targetLength(object))
    pam_len  <- pamLength(object)
    pam_side <- pamSide(object)
    if (pam_side=="5prime"){
        indices <- indices[seq_len(pam_len)]
    } else {
        indices <- rev(rev(indices)[seq_len(pam_len)])
    }
    return(indices)
})


#' @rdname CrisprNickase-class
#' @export
setMethod("spacerIndices", "CrisprNickase", 
    function(object){
    indices <- seq_len(targetLength(object))
    spacer_len <- spacerLength(object)
    pam_side <- pamSide(object)
    if (pam_side=="3prime"){
        indices <- indices[seq_len(spacer_len)]
    } else {
        indices <- rev(rev(indices)[seq_len(spacer_len)])
    }
    return(indices)
})

#' @rdname CrisprNickase-class
#' @export
setMethod("prototypeSequence",
          "CrisprNickase",
          function(object, primary=TRUE){
    pam_len    <- pamLength(object)
    spacer_len <- spacerLength(object)
    gap_len    <- spacerGap(object)
    pam_side   <- pamSide(object)
    # Building sequence:
    spacer <- paste0(rep("S", spacer_len), collapse="")
    gap <- paste0(rep("-", gap_len), collapse="")
    pam <- motifs(object, primary=primary)
    if (pam_side=="3prime"){
        seq <- paste0(spacer, gap, "[" , pam, "]")
    } else {
        seq <- paste0("[", pam, "]", gap, spacer)
    }
    seq <- paste0("5'--", seq, "--3'")
    return(seq)
})  
