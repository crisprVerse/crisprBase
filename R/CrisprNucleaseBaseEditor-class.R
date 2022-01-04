#' An S4 class to represent a CRISPR nuclease with base editor.
#' 
#' @slot baseEditorName Name of the base editor enzyme.
#' @slot editingWeights Matrix of editing weights.
#' @slot editingStrand String indicating which strand with
#'     respect to the target protospacer sequence will be 
#'     edited. Must be either "original" or "opposite".
#'     "original" by default.
#' 
#' @section Constructors:
#'     Use the constructor \code{link{CrisprNucleaseBaseEditor}} to create
#'     a CrisprNucleaseBaseEditor object.
#' 
#' @return A CrisprNucleaseWithBaseEditor object
#' @export
setClass("CrisprNucleaseBaseEditor",
    contains = "CrisprNuclease",
    slots = c(
        baseEditorName = "character", 
        editingWeights = "matrix",
        editingStrand = "character"),
    prototype = list(
        baseEditorName = NA_character_,
        editingWeights = NULL,
        editingStrand = "original"
    )
)



#' @describeIn CrisprNucleaseBaseEditor Create a \linkS4class{CrisprNucleaseBaseEditor} object
#' @param CrisprNuclease A \linkS4class{CrisprNuclease} object.
#' @param baseEditorName String specifying base editor name.
#' @param editingStrand String indicating which strand with
#'     respect to the target protospacer sequence will be 
#'     edited. Must be either "original" or "opposite".
#'     "original" by default.
#' @param editingWeights Numeric matrix of editing weights.
#'     Column names must be indicating relative position to the PAM site.
#'     Row names must be of the form "X2Y" where "X" represents the origin
#'     base, and "Y" represents the subtituted base. For instance, "C2T"
#'     indicates the row corresponding to C to T editing. 
#' @export
CrisprNucleaseBaseEditor <- function(CrisprNuclease,
                                     baseEditorName = NA_character_,
                                     editingStrand = c("original", "opposite"),
                                     editingWeights = NULL
){
    editingStrand <- match.arg(editingStrand)
    new("CrisprNucleaseBaseEditor",
        CrisprNuclease,
        baseEditorName = as.character(baseEditorName),
        editingStrand = editingStrand,
        editingWeights = .buildEditingWeightsMatrix(editingWeights)
    )
}



#' @rdname CrisprNucleaseBaseEditor-class
#' @param object \linkS4class{CrisprNucleaseBaseEditor} object.
#' @export
setMethod("baseEditorName",
          "CrisprNucleaseBaseEditor",function(object){
    object@baseEditorName
})


#' @rdname CrisprNucleaseBaseEditor-class
#' @param value Value to replaced with.
#' @export
setMethod("baseEditorName<-",
          "CrisprNucleaseBaseEditor",function(object, value){
    value <- as.character(value)
    object@baseEditorName <- value
    return(object)
})


#' @rdname CrisprNucleaseBaseEditor-class
#' @export
setMethod("editingWeights",
          "CrisprNucleaseBaseEditor",function(object){
    object@editingWeights
})


#' @rdname CrisprNucleaseBaseEditor-class
#' @export
setMethod("editingWeights<-",
          "CrisprNucleaseBaseEditor",function(object, value){
    value <- .buildEditingWeightsMatrix(value)
    object@editingWeights <- value
    return(object)
})


#' @rdname CrisprNucleaseBaseEditor-class
#' @export
setMethod("editingStrand",
          "CrisprNucleaseBaseEditor",function(object){
    object@editingStrand
})


#' @rdname CrisprNucleaseBaseEditor-class
#' @export
setMethod("editingStrand<-",
          "CrisprNucleaseBaseEditor",function(object, value){
    value <- as.character(value)
    if (value %in% c("original", "oppposite")){
        stop("value must be either 'original' or 'opposite'.")
    }
    object@editingStrand <- value
    return(object)
})









#x <- matrix(1, nrow=1, ncol=5)
#rownames(x) <- "C2T"
#colnames(x) <- -18:-14
.buildEditingWeightsMatrix <- function(x=NULL){
    editChoices <- .getComboNames()
    posChoices <- -1000:1000
    if (is.null(x)){
        out <- matrix(0,
                      nrow=length(editChoices),
                      ncol=1)
        rownames(out) <- editChoices
        colnames(out) <- "0"
    } else if (is.matrix(x)){
        if (!all(rownames(x) %in% editChoices)){
            diff <- setdiff(rownames(x), editChoices)
            diff <- paste0(diff, collapse=";")
            stop("The following rownames are not valid",
                 " names: ", diff)
        }
        if (!all(colnames(x) %in% posChoices)){
            diff <- setdiff(colnames(x), posChoices)
            diff <- paste0(diff, collapse=";")
            stop("The following colnames are not valid",
                 " names: ", diff)
        }
        seq_start <- min(as.integer(colnames(x)))
        seq_end <- max(as.integer(colnames(x)))
        seq <- seq_start:seq_end
        out <- matrix(0,
                      nrow=length(editChoices),
                      ncol=length(seq))
        rownames(out) <- editChoices
        colnames(out) <- seq
        out[rownames(x), colnames(x)] <- x

    } else {
        stop("x must be NULL or a matrix.")
    }
    return(out)
} 



.getComboNames <- function(){
    dnaLetters <- c("A", "C", "G", "T")
    combs <- expand.grid(dnaLetters, dnaLetters)
    combs <- combs[combs[,1]!=combs[,2],]
    combs <- paste0(combs[,1], "2", combs[,2])
    return(combs)
}


.getOriginBaseFromRownames <- function(x){
    substr(x,1,1)
}

.getTargetBaseFromRownames <- function(x){
    substr(x,3,3)
}





#combs <- .getComboNames()
#.getOriginBase(combs)
#.getTargetBase(combs)


#data(SpCas9, package="crisprBase")
#nuc <- CrisprNucleaseBaseEditor(SpCas9)








