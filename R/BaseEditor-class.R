#' An S4 class to represent a base editor
#' 
#' @slot baseEditorName Name of the base editor. 
#' @slot editingWeights Matrix of editing weights.
#' @slot editingStrand String indicating which strand with
#'     respect to the target protospacer sequence will be 
#'     edited. Must be either "original" or "opposite".
#'     "original" by default.
#' 
#' @section Constructors:
#'     Use the constructor \code{link{BaseEditor}} to create
#'     a BaseEditor object.
#' 
#' @section Accessors:
#' \describe{
#'     \item{\code{baseEditorName}:}{To get the name of the base editor.}
#'     \item{\code{editingWeights}:}{To return the matrix of editing weights.}
#'     \item{\code{editingStrand}:}{To return the editing strand.}
#' }
#' 
#' @section Setters:
#' \describe{
#'     \item{\code{baseEditorName<-}:}{To change the name of the base editor.}
#'     \item{\code{editingWeights<-}:}{To change the matrix of editing weights.}
#'     \item{\code{editingStrand<-}:}{To change the editing strand.}
#' }
#' 
#' @examples
#' # Creating an object for BE4max (C to T editor)
#' # based on experimental weights
#' 
#' ws <- c(0.7, 0.7, 0.8, 1.8, 1, 2, 1.4, 1.2, 2.3, 1.3, 2.4, 2.2, 3.4, 
#'       2.2, 2.1, 3.5, 5.8, 16.2, 31.8, 63.2, 90.3, 100, 87, 62, 31.4, 
#'       16.3, 10, 5.6, 3.3, 1.9, 1.8, 2.4, 1.7, 0.5, 0.2, 0.1)
#' ws <- matrix(ws, nrow=1, ncol=length(ws))
#' rownames(ws) <- "C2T"
#' colnames(ws) <- -36:-1
#' data(SpCas9, package="crisprBase")
#' BE4max <- BaseEditor(SpCas9,
#'                      baseEditorName="BE4max",
#'                      editingStrand="original",
#'                      editingWeights=ws)
#' metadata(BE4max)$description_base_editor <- "BE4max cytosine base editor."
#' 
#' @return A BaseEditor object
#' 
#' @export
setClass("BaseEditor",
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



#' @describeIn BaseEditor Create a \linkS4class{BaseEditor} object
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
BaseEditor <- function(CrisprNuclease,
                       baseEditorName = NA_character_,
                       editingStrand = c("original", "opposite"),
                       editingWeights = NULL
){
    editingStrand <- match.arg(editingStrand)
    new("BaseEditor",
        CrisprNuclease,
        baseEditorName = as.character(baseEditorName),
        editingStrand = editingStrand,
        editingWeights = .buildEditingWeightsMatrix(editingWeights)
    )
}





#' @rdname BaseEditor-class
#' @export
setMethod("show", "BaseEditor", function(object) {
    len <- length(metadata(object))
    dnase <- isDnase(object)
    if (dnase){
        pams.line <- ("      PAMs: ")
        pams.line.distance <- "        Distance from PAM: "
        pam.side.line <- "      PAM side: "
    } else {
        pams.line <- ("      PFS: ")
        pams.line.distance <- "        Distance from PFS: "
        pam.side.line <- "      PFS side: "
    }
    cat(paste0("Class: ", is(object)[[1]]), "\n",
        "  CRISPR Nuclease name: ", nucleaseName(object), "\n",
        "      Target type: ", targetType(object), "\n",
        "      Metadata: list of length ", len, "\n",
        pams.line, .printVectorNicely(motifs(object)), "\n",
        "      Weights: ", .printVectorNicely(weights(object)), "\n",
        "      Spacer length: ",  spacerLength(object), "\n",
        pam.side.line,  pamSide(object), "\n",
        pams.line.distance, spacerGap(object), "\n",
        "      Prototype protospacers: ",
        .printVectorNicely(prototypeSequence(object, primary=FALSE)),
        "\n",
        "  Base editor name: ", baseEditorName(object), "\n",
        "      Editing strand: ", editingStrand(object), "\n",
        "      Maximum editing weight: ", .getMaxEditingWeight(object), "\n",
        sep = "")
})



.getMaxEditingWeight <- function(object){
    x <- editingWeights(object)
    ind <- which.max(x)
    ind <- arrayInd(ind, .dim=c(nrow(x), ncol(x)))
    sub <- rownames(x)[ind[1]]
    pos <- colnames(x)[ind[2]]
    paste0(sub, " at position ", pos, "")
}



#' @rdname BaseEditor-class
#' @param object \linkS4class{BaseEditor} object.
#' @export
setMethod("baseEditorName",
          "BaseEditor",function(object){
    object@baseEditorName
})


#' @rdname BaseEditor-class
#' @param value Value to replaced with.
#' @export
setMethod("baseEditorName<-",
          "BaseEditor",function(object, value){
    value <- as.character(value)
    object@baseEditorName <- value
    return(object)
})


#' @rdname BaseEditor-class
#' @param substitutions Character vector indicating which substitutions
#'     should be returned. 
#' @export
setMethod("editingWeights",
          "BaseEditor",function(object,
                                substitutions=NULL
){
    ws <- object@editingWeights
    if (is.null(substitutions)){
        substitutions <- rownames(ws)
    }
    diff <- setdiff(substitutions, rownames(ws))
    if (length(diff)!=0){
        stop("Some of the substitutions are not found.")
    }
    out <- ws[substitutions,,drop=FALSE]
    return(out)
})


#' @rdname BaseEditor-class
#' @export
setMethod("editingWeights<-",
          "BaseEditor",function(object, value){
    value <- .buildEditingWeightsMatrix(value)
    object@editingWeights <- value
    return(object)
})


#' @rdname BaseEditor-class
#' @export
setMethod("editingStrand",
          "BaseEditor",function(object){
    object@editingStrand
})


#' @rdname BaseEditor-class
#' @export
setMethod("editingStrand<-",
          "BaseEditor",function(object, value){
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
.buildEditingWeightsMatrix <- function(x=NULL,
                                       scale=TRUE
){
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
        if (scale){
            out <- out/max(out, na.rm=TRUE)
        }
    } else {
        stop("x must be NULL or a matrix.")
    }
    return(out)
} 


.getReducedEditingMatrix <- function(ws){
    good <- rowSums(ws==0, na.rm=TRUE)!=ncol(ws)
    ws[good,,drop=FALSE]
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




#' Quick plot to visualize editing weights 
#' 
#' Quick plot to visualize editing weights from a 
#' BaseEditor object. 
#' 
#' @param baseEditor A \code{\linkS4class{BaseEditor}} object. 
#' @param discardEmptyRows Should rows that have all weight equal
#'     to 0 be discarded? TRUE by default. 
#' @param substitutions Character vector specifying substitutions
#'     to be plotted. If NULL (default), all substitutions are shown.
#' @param ... Additional arguments to be passed to \code{plot}
#' 
#' @return Nothing. A plot is generated as a side effect. 
#' 
#' @examples
#' if (interactive()){
#'     data(BE4max, package="crisprBase")
#'     plotEditingWeights(BE4max)
#' }
#' @importFrom graphics legend lines
#' @export
plotEditingWeights <- function(baseEditor,
                               discardEmptyRows=TRUE,
                               substitutions=NULL,
                               ...
){
    .isBaseEditorOrStop(baseEditor)
    choices <- .getComboNames()
    if (is.null(substitutions)){
        substitutions <- choices
    } else {
        diff <- setdiff(substitutions, choices)
        if (length(diff)!=0){
            diff <- paste0(diff, collapse=",")
            stop("The following substitutions are not valid: ",
                 diff, ".")
        }
    }
    ws <- editingWeights(baseEditor,
                         substitutions=substitutions)
    if (discardEmptyRows){
        ws <- .getReducedEditingMatrix(ws)
    }
    x <- as.numeric(colnames(ws))
    top <- max(ws, na.rm=TRUE)
    ylim <- c(0,top)
    plot(x, ws[1,], col="white",
         xlab="Position relative to PAM site",
         ylab="Relative weight",
         ylim=ylim,
         ...)
    ns <- nrow(ws)
    col <- seq_len(ns)
    for (k in seq_len(ns)){
        lines(x, ws[k,], col=col[k])
    }
    legend("topleft",
           legend=rownames(ws),
           col=col,
           lty=1, cex=0.75)
}











