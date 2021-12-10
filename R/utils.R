#' @import methods
NULL




#' @export
#' @importFrom S4Vectors metadata
S4Vectors::metadata


#' @export
#' @importFrom S4Vectors metadata<-
S4Vectors::`metadata<-`




.checkDNAAlphabet <- function(x){
    chars <- unlist(strsplit(x, split=""))
    choices <- names(Biostrings::IUPAC_CODE_MAP)
    if (!all(chars %in% choices)){
        stop("Some characters are not part of the extended ",
             "nucleotide alphabet.")
    }
}



.validateStrand <- function(strand,
                            values=c("+","-")
){
    if (!all(strand %in% values) | is.null(strand)){
        choices <- paste0(values, collapse=", ")
        stop("strand must have one of the following values: ",
             choices)
    }
    return(strand)
}


.isCrisprNuclease <- function(object){
    is(object, "CrisprNuclease")
}

.isNuclease <- function(object){
    is(object, "Nuclease")
}

.isNucleaseOrStop <- function(object) {
    if (!is(object, "Nuclease")) {
        stop("object is of class '", class(object), "', but needs to be of ",
             "class 'Nuclease'")
    }
}

.isCrisprNucleaseOrStop <- function(object) {
    if (!is(object, "CrisprNuclease")) {
        stop("object is of class '", class(object), "', but needs to be of ",
             "class 'CrisprNuclease'")
    } 
}



# convert values in scientific notation to integers
.makeLongInteger <- function(n){
    as.integer(format(n, scientific=FALSE))
}



.validateNonNegativeInteger <- function(x){
    x  <- .validateInteger(x)
    if (x < 0){
        stop("x must be 0 or positive")
    }
    return(as.integer(x))
}
.validateNonNegativeInteger <- Vectorize(.validateNonNegativeInteger)

.validateInteger <- function(x){
    tol = .Machine$double.eps^0.5
    check <- abs(x - as.integer(x)) < tol
    if (!check){
        stop("x must be an integer")
    } 
    return(as.integer(x))
}



.printVectorNicely <- function(x){
    if (length(x)<=3){
        x <- paste0(x, collapse=", ")
    } else {
        x.last <- x[length(x)]
        x <- paste0(x[seq_len(3)], collapse=", ")
        x <- paste0(x, ",..., ", x.last)
    }
    return(x)
}




.validateGROrNull <- function(gr){
    if (is.null(gr)){
        return(gr)
    }
    if (!is(gr, "GRanges")){
        stop("gr must be a GenomicRanges object")
    }
    return(gr)
}

#' @importFrom BiocGenerics width
.validatePosGrOrNull <- function(gr){
    .validateGROrNull(gr)
    if (is.null(gr)){
        return(gr)
    }
    ws <- width(gr)
    if (!all(ws==1)){
        stop("Supplied GRanges must have ranges of width 1.")
    } 
    return(gr)
}




#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
.buildGRFromPamSite <- function(seqnames,
                                pam_site,
                                strand
){
    len1  <- length(seqnames)
    len2  <- length(pam_site)
    len3  <- length(strand)
    lens  <- unique(c(len1, len2, len3))
    cond1 <- length(lens)==3
    cond2 <- length(lens)==2 & !(1 %in% lens)

    if (cond1 | cond2){
        stop("seqnames, pam_site, strand must all have",
             " the same length if not of length 1")
    } 
    pam_site <- .validateNonNegativeInteger(pam_site)
    strand   <- .validateStrand(strand)
    seqnames <- as.character(seqnames)

    gr <- GRanges(seqnames,
                  IRanges(start=pam_site, width=1),
                  strand=strand)
    return(gr)
}



.resetGRCoordinates <- function(gr){
    start(gr) <- rep(1, length(gr))
    end(gr)   <- rep(1, length(gr))
    return(gr)
}


#' @title Return list of available CrisprNuclease objects in crisprBase
#' @description Return list of available CrisprNuclease objects in crisprBase. 
#'     
#' @return Character vector of available \code{CrisprNuclease} objects found
#'     in \code{crisprBase}.
#' 
#' @author Jean-Philippe Fortin
#'
#' @examples
#' getAvailableCrisprNucleases()
#' 
#' @importFrom utils data
#' @export
getAvailableCrisprNucleases <- function(){
    results <- data(package="crisprBase")[["results"]]
    results <- as.data.frame(results)
    results <- results[, c("Item", "Title")]
    results <- results[grepl("CrisprNuclease", results[["Title"]]),]
    as.character(results[["Item"]])
}


