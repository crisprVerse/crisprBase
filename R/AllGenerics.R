#### Generics for Nuclease object

#' @description Return motif string representations of recognition sites.
#' @rdname Nuclease-class
#' @param object Object of class \linkS4class{Nuclease} or
#'     \linkS4class{CrisprNuclease}
#' @param x Object of class \linkS4class{Nuclease} or
#'     \linkS4class{CrisprNuclease}
#' @param ... Additional arguments for class-specific methods
#' 
#' @export
setGeneric("motifs", function(object, ...) standardGeneric("motifs"))


#' @description Return length of the recognition sites sequences.
#' @rdname Nuclease-class
#' @export
setGeneric("motifLength",
           function(object, ...) standardGeneric("motifLength"))


#' @rdname Nuclease-class
#' @export
setGeneric("nucleaseName", function(object) standardGeneric("nucleaseName"))



#' @rdname Nuclease-class
#' @export
setGeneric("weights", function(object) standardGeneric("weights"))



#' @rdname Nuclease-class
#' @export
setGeneric("nucleaseName<-", function(object, value) standardGeneric("nucleaseName<-"))


#' @rdname Nuclease-class
#' @export
setGeneric("weights<-", function(object, value) standardGeneric("weights<-"))


#' @rdname Nuclease-class
#' @export
setGeneric("cutSites", function(object, ...) standardGeneric("cutSites"))


#' @rdname Nuclease-class
#' @export
setGeneric("isCutting", function(object) standardGeneric("isCutting"))






#### Generics for CrisprNuclease
#' @rdname CrisprNuclease-class
#' @export
setGeneric("spacerLength", function(object) standardGeneric("spacerLength"))

#' @rdname CrisprNuclease-class
#' @export
setGeneric("protospacerLength",
           function(object) standardGeneric("protospacerLength"))

#' @rdname CrisprNuclease-class
#' @export
setGeneric("pamLength", function(object) standardGeneric("pamLength"))

#' @rdname CrisprNuclease-class
#' @export
setGeneric("spacerSide", function(object) standardGeneric("spacerSide"))


#' @rdname CrisprNuclease-class
#' @export
setGeneric("spacerSide<-",
           function(object, value) standardGeneric("spacerSide<-"))


#' @rdname CrisprNuclease-class
#' @export
setGeneric("spacerGap", function(object) standardGeneric("spacerGap"))


#' @rdname CrisprNuclease-class
#' @export
setGeneric("hasSpacerGap", function(object) standardGeneric("hasSpacerGap"))



#' @rdname CrisprNuclease-class
#' @export
setGeneric("spacerGap<-",
           function(object, value) standardGeneric("spacerGap<-"))


#' @rdname CrisprNuclease-class
#' @export
setGeneric("spacerLength<-",
           function(object, value) standardGeneric("spacerLength<-"))


#' @rdname CrisprNuclease-class
#' @export
setGeneric("pamSide", function(object) standardGeneric("pamSide"))

#' @rdname CrisprNuclease-class
#' @param ... Additional arguments for class-specific methods
#' @export
setGeneric("pams", function(object, ...) standardGeneric("pams"))

#' @rdname CrisprNuclease-class
#' @export
setGeneric("pamIndices", function(object, ...) standardGeneric("pamIndices"))

#' @rdname CrisprNuclease-class
#' @export
setGeneric("spacerIndices",
           function(object, ...) standardGeneric("spacerIndices"))

#' @rdname CrisprNuclease-class
#' @export
setGeneric("prototypeSequence",
           function(object, ...) standardGeneric("prototypeSequence"))






