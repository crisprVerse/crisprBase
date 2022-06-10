#### Generics for Nuclease object

#' @description Return motif string representations of recognition sites.
#' @rdname Nuclease-class
#' @param object Object of class \linkS4class{Nuclease} or
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
setGeneric("targetType", function(object) standardGeneric("targetType"))





#' @rdname Nuclease-class
#' @export
setGeneric("weights", function(object) standardGeneric("weights"))



#' @rdname Nuclease-class
#' @export
setGeneric("nucleaseName<-", function(object, value) standardGeneric("nucleaseName<-"))


#' @rdname Nuclease-class
#' @export
setGeneric("targetType<-", function(object, value) standardGeneric("targetType<-"))



#' @rdname Nuclease-class
#' @export
setGeneric("weights<-", function(object, value) standardGeneric("weights<-"))


#' @rdname Nuclease-class
#' @export
setGeneric("cutSites", function(object, ...) standardGeneric("cutSites"))


#' @rdname Nuclease-class
#' @export
setGeneric("isCutting", function(object) standardGeneric("isCutting"))


#' @rdname Nuclease-class
#' @export
setGeneric("isRnase", function(object) standardGeneric("isRnase"))


#' @rdname Nuclease-class
#' @export
setGeneric("isDnase", function(object) standardGeneric("isDnase"))







#### Generics for CrisprNuclease
#' @rdname CrisprNuclease-class
#' @export
setGeneric("spacerLength", function(object, ...) standardGeneric("spacerLength"))

#' @rdname CrisprNuclease-class
#' @export
setGeneric("targetLength",
           function(object, ...) standardGeneric("targetLength"))

#' @rdname CrisprNuclease-class
#' @export
setGeneric("pamLength", function(object, ...) standardGeneric("pamLength"))




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
setGeneric("pamSide", function(object, ...) standardGeneric("pamSide"))


#' @rdname CrisprNuclease-class
#' @export
setGeneric("pamSide<-", function(object, value) standardGeneric("pamSide<-"))



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



#' @rdname BaseEditor-class
#' @export
setGeneric("baseEditorName", function(object) standardGeneric("baseEditorName"))


#' @rdname BaseEditor-class
#' @export
setGeneric("baseEditorName<-", function(object, value) standardGeneric("baseEditorName<-"))


#' @rdname BaseEditor-class
#' @param ... Additional arguments for class-specific methods
#' @export
setGeneric("editingWeights", function(object, ...) standardGeneric("editingWeights"))


#' @rdname BaseEditor-class
#' @export
setGeneric("editingWeights<-", function(object, value) standardGeneric("editingWeights<-"))



#' @rdname BaseEditor-class
#' @export
setGeneric("editingStrand", function(object, ...) standardGeneric("editingStrand"))


#' @rdname BaseEditor-class
#' @export
setGeneric("editingStrand<-", function(object, value) standardGeneric("editingStrand<-"))







