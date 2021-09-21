#' An S4 class to represent a CRISPR nuclease.
#' 
#' @slot spacer_side String specifying the side of the gRNA spacer
#'     sequence with respect to the PAM motif. Must be either 
#'     '5prime' (e.g. Cas9) or '3prime' (e.g. Cas12a)
#' @slot spacer_length Integer specifying the
#'     length of the spacer sequence
#' @slot spacer_gap Integer specifying the length (in nucleotides)
#'     between the spacer sequence and the PAM sequence
#'     (e.g. 0 for Cas9 and Cas12a).
#' 
#' @section Constructors:
#'     Use the constructor \code{link{CrisprNuclease}} to create
#'         a CrisprNuclease object.
#' 
#' @section Accessors:
#' \describe{
#'     \item{\code{name}:}{To get the name of the CRISPR nuclease.} 
#'     \item{\code{spacerLength}:}{To return the length of the
#'         spacer sequence.}
#'     \item{\code{protospacerLength}:}{To return the length of the
#'         protospacer sequence.}
#'     \item{\code{pamLength}:}{To return the length of the PAM sequence.}
#'     \item{\code{spacerSide}:}{To return the side of the spacer sequence
#'         with respect to the PAM sequence.}
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
#'     \item{\code{spacerSide<-}:}{To change the side of the spacer sequence
#'         with respect to the PAM sequence.}
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
#' SpCas9 <- CrisprNuclease("SpCas9",
#'                          pams=c("(3/3)NGG", "(3/3)NAG", "(3/3)NGA"),
#'                          weights=c(1, 0.2593, 0.0694),
#'                          info="Wildtype Streptococcus pyogenes Cas9 
#'                               (SpCas9) nuclease",
#'                          spacer_side="5prime",
#'                          spacer_length=20)
#' 
#' @export
setClass("CrisprNuclease", 
    contains = "Nuclease", 
    slots = c(
        spacer_side = "character",
        spacer_gap = "integer",
        spacer_length = "integer"),
    prototype = list(
        spacer_side = NA_character_,
        spacer_gap = 0L,
        spacer_length =NA_integer_
    )
)



#' @describeIn CrisprNuclease Create a \linkS4class{CrisprNuclease} object
#' @param name Name of the CRISPR nuclease.
#' @param pams Character vector of PAM sequence motifs
#'           written from 5' to 3. If the point of cleavage has
#'           been determined, the precise site is marked with ^.
#'           Only letters in the IUPAC code are accepted.
#'           For nucleases that cleave away from their
#'           recognition sequence, the cleavage sites are indicated
#'           in parentheses. See details for more information.
#' @param weights Optional numeric vector specifying relative weights
#'           of the PAM sequences to specify cleavage probabilities. 
#' @param info Optional string providing a description of the CRISPR nuclease.
#' @param spacer_side String specifying the side of the gRNA spacer
#'     sequence with respect to the PAM motif. Must be either 
#'     '5prime' (e.g. Cas9) or '3prime' (e.g. Cas12a)
#' @param spacer_length Integer specifying the length of the spacer sequence
#' @param spacer_gap Integer specifying the length (in nucleotides) between
#'     the spacer sequence and the PAM sequence (e.g. 0 for Cas9 and Cas12a).
#' @export
#' @export
CrisprNuclease <- function(name,
                           pams = NA_character_,
                           weights = rep(1, length(pams)),
                           info = NA_character_,
                           spacer_side = NA_character_,
                           spacer_gap = 0L,
                           spacer_length = NA_integer_
){

    nuc <- Nuclease(name=name,
                    motifs=pams,
                    weights=weights,
                    info=info)
    new("CrisprNuclease",
        nuc,
        spacer_side = as.character(spacer_side),
        spacer_gap = as.integer(spacer_gap),
        spacer_length = as.integer(spacer_length))
}



#' @rdname CrisprNuclease-class
#' @export
setMethod("show", "CrisprNuclease", function(object) {
    cat(paste0("Class: ", is(object)[[1]]), "\n",
        "  Name: ", object@name, "\n",
        "  Info: ", object@info, "\n",
        "  Motifs: ", .printVectorNicely(object@motifs), "\n",
        "  Weights: ", .printVectorNicely(object@weights), "\n",
        "  Spacer: \n",
        "    Side: ", object@spacer_side, "\n",
        "    Length: ", object@spacer_length, "\n",
        "    Distance from PAM: ", object@spacer_gap, "\n",
        sep = "")
})




setValidity("CrisprNuclease", function(object) {
    out <- TRUE
    if (length(object@spacer_side)!=1){
        out <- "@spacer_side must be of length 1"
    } 
    if (!object@spacer_side %in% c("5prime", "3prime")){
        out <- "@spacer_side must be either 5prime or 3 prime"
    } 
    if (length(object@spacer_length)!=1){
        out <- "@spacer_length must be of length 1"
    } 
    if (object@spacer_length<0){
        out <- "@spacer_length must be a positive integer"
    } 
    if (length(object@spacer_gap)!=1){
        out <- "@spacer_gap must be of length 1"
    } 
    if (object@spacer_gap<0){
        out <- "@spacer_gap must be either 0 or a positive integer"
    } 
    return(out)
})







############################################################
##############       Getter and setters   ##################
#' @rdname CrisprNuclease-class
#' @param object \linkS4class{CrisprNuclease} object.
#' @export
setMethod("pamLength",
          "CrisprNuclease",function(object){
    len <- motifLength(object)[1]
    names(len) <- NULL
    return(len)
})

#' @rdname CrisprNuclease-class
#' @export
setMethod("spacerLength", "CrisprNuclease", function(object){
    return(object@spacer_length)
})


#' @rdname CrisprNuclease-class
#' @param value For \code{spacerLength<-} and \code{gapLength<-}, must be 
#'     a non-negative integer. For \code{spacerSide}, must be either
#'     '5prime' or '3prime'.
#' 
#' @export
setMethod("spacerLength<-", "CrisprNuclease", function(object, value){
    value <- .validateNonNegativeInteger(value)
    object@spacer_length <- value
    return(object)
})


#' @rdname CrisprNuclease-class
#' @export
setMethod("spacerSide",
          "CrisprNuclease",function(object){
    return(object@spacer_side)
})

#' @rdname CrisprNuclease-class
#' @export
setMethod("spacerSide<-", "CrisprNuclease", function(object, value){
    if (!value %in% c("5prime", "3prime")){
        stop("value must be either '5prime' or '3prime'")
    }
    object@spacer_side <- value
    return(object)
})


#' @rdname CrisprNuclease-class
#' @export
setMethod("spacerGap",
          "CrisprNuclease",function(object){
    return(object@spacer_gap)
})


#' @rdname CrisprNuclease-class
#' @export
setMethod("pamSide",
          "CrisprNuclease",function(object){
    side <- object@spacer_side
    ifelse(side=="3prime", "5prime", "3prime")
})

#' @rdname CrisprNuclease-class
#' @export
setMethod("spacerGap<-", "CrisprNuclease", function(object, value){
    value <- .validateNonNegativeInteger(value)
    object@spacer_gap <- value
    return(object)
})

#' @rdname CrisprNuclease-class
#' @export
setMethod("protospacerLength",
          "CrisprNuclease",function(object){
    pamLength(object) + spacerLength(object) + spacerGap(object)
})



#' @rdname CrisprNuclease-class
#' @param primary Should only the PAM sequence with the heighest weight
#'     be returned? If no cleavage weights are stored in the
#'     \linkS4class{CrisprNuclease} object, all sequences are returned.
#'     TRUE by default.
#' @param ignore_pam Should all possible k-mer sequences for a given PAM length
#'     be returned, irrespetively of the PAM sequence motifs stored in the
#'     \linkS4class{CrisprNuclease} object? FALSE by default.
#' @param as.character Should the PAM sequences be returned as a 
#'     character vector? FALSE by default.
#' @export
setMethod("pams", "CrisprNuclease",
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
#' @rdname CrisprNuclease-class
#' @export
setMethod("pamIndices", "CrisprNuclease",
    function(object){
    indices  <- seq_len(protospacerLength(object))
    pam_len  <- pamLength(object)
    pam_side <- pamSide(object)
    if (pam_side=="5prime"){
        indices <- indices[seq_len(pam_len)]
    } else {
        indices <- rev(rev(indices)[seq_len(pam_len)])
    }
    return(indices)
})


#' @rdname CrisprNuclease-class
#' @export
setMethod("spacerIndices", "CrisprNuclease", 
    function(object){
    indices <- seq_len(protospacerLength(object))
    spacer_len <- spacerLength(object)
    spacer_side <- spacerSide(object)
    if (spacer_side=="5prime"){
        indices <- indices[seq_len(spacer_len)]
    } else {
        indices <- rev(rev(indices)[seq_len(spacer_len)])
    }
    return(indices)
})

#' @rdname CrisprNuclease-class
#' @export
setMethod("prototypeSequence",
          "CrisprNuclease",
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
        seq <- paste0(spacer, gap, pam)
    } else {
        seq <- paste0(pam, gap, spacer)
    }
    return(seq)
})  
