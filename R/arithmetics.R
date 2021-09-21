#' @title Extract PAM sequences from protospacer sequences
#' @description Extract PAM sequences from protospacer sequences using
#'     information stored in a CrisprNuclease object. 
#'     
#' 
#' @param protospacers Character vector of protospacer sequences. 
#' @param object \code{CrisprNuclease} corresponding to the protospacer 
#'     sequences.
#' 
#' @return Character vector of PAM sequences of length equal to that of the 
#'     \code{protospacers} character vector. 
#' 
#' @author Jean-Philippe Fortin
#'
#' @examples
#' # Extracting PAM sequences from Cas9 protospacers:
#' protospacers <- c("AGGTGCTGATTGTAGTGCTGCGG",
#'                   "AGGTGCTGATTGTAGTGCTGAGG")
#' extractPamFromProtospacer(protospacers, SpCas9)
#' # Extracting PAM sequences from Cas12a protospacers:
#' protospacers <- c("TTTAAGGTGCTGATTGTAGTGCTGTGT",
#'                   "TTTCAGGTGCTGATTGTAGTGCTGAAA")
#' extractPamFromProtospacer(protospacers, AsCas12a)
#' 
#' @export
extractPamFromProtospacer <- function(protospacers,
                                      object
){
    .isCrisprNucleaseOrStop(object)
    .validateProtospacers(protospacers, object)
    wh <- pamIndices(object)
    pams <- substr(protospacers,
                   wh[1],
                   wh[length(wh)])
    return(pams)
}



#' @title Extract spacer sequences from protospacer sequences
#' @description Extract spacer sequences from protospacer sequences using
#'     information stored in a CrisprNuclease object. 
#'     
#' 
#' @param protospacers Character vector of protospacer sequences. 
#' @param object \code{CrisprNuclease} corresponding to the protospacer 
#'     sequences.
#' 
#' @return Character vector of spacer sequences of length equal to that of the 
#'     \code{protospacers} character vector. 
#' 
#' @author Jean-Philippe Fortin
#'
#' @examples
#' # Extracting spacer sequences from Cas9 protospacers:
#' protospacers <- c("AGGTGCTGATTGTAGTGCTGCGG",
#'                   "AGGTGCTGATTGTAGTGCTGAGG")
#' extractSpacerFromProtospacer(protospacers, SpCas9)
#' # Extracting spacer sequences from Cas12a protospacers:
#' protospacers <- c("TTTAAGGTGCTGATTGTAGTGCTGTGT",
#'                   "TTTCAGGTGCTGATTGTAGTGCTGAAA")
#' extractSpacerFromProtospacer(protospacers, AsCas12a)
#' 
#' @export
extractSpacerFromProtospacer <- function(protospacers,
                                         object
){
    .isCrisprNucleaseOrStop(object)
    .validateProtospacers(protospacers, object)
    wh <- spacerIndices(object)
    spacers <- substr(protospacers,
                      wh[1],
                      wh[length(wh)])
    return(spacers)
}





#' @title Construct a protospacer GRanges from a list of PAM sites 
#' 
#' @description Construct a protospacer GRanges from a list of PAM sites
#'     using information stored in a CrisprNuclease object. 
#'     
#' @param gr GRanges object of width 1 specifying the coordinates
#'     of the first nucleotide of the PAM sequences. 
#' @param seqnames Character vector of genomic sequence names.
#'     Ignored if \code{gr} is not NULL.
#' @param pam_site Numeric vector specifying the coordinates of the
#'     first nucleotide of the PAM sequences corresponding to the
#'     protospacers. Ignored if \code{gr} is not NULL.
#' @param strand Character vector specifying the strand of the protospacer. 
#'     Ignored if \code{gr} is not NULL.
#' @param nuclease CrisprNuclease object.
#' @param spacer_len Non-negative integer to overwrite the default spacer
#'     length stored in the CrisprNuclease object.
#' 
#' @return GRanges object representing genomic coordinates of spacer sequences.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples 
#' if (require(GenomicRanges)){
#' gr <- GRanges("chr10",
#'               IRanges(start=c(100,120), width=1),
#'               strand=c("+","-"))
#' getProtospacerRanges(gr, nuclease=SpCas9)
#' getProtospacerRanges(gr, nuclease=AsCas12a)
#' }
#' 
#' @export
#' @importFrom BiocGenerics start end start<- end<- strand
getProtospacerRanges <- function(gr=NULL,
                                 seqnames=NULL,
                                 pam_site=NULL,
                                 strand=NULL,
                                 nuclease=NULL,
                                 spacer_len=NULL
){
    .isCrisprNucleaseOrStop(nuclease)
    if (!is.null(spacer_len)){
        if (length(spacer_len)>1){
            stop("spacer_len must be either NULL or of length 1.")
        }
        message("Overwriting spacerLength(nuclease) with spacer_len")
        spacerLength(nuclease) <- spacer_len
    }
    gr <- .validatePosGrOrNull(gr)
    if (is.null(gr) & 
        (is.null(seqnames) | is.null(pam_site) | is.null(strand))){
        stop("seqnames, pam_site, and strand must be provided if gr=NULL")
    }
    if (is.null(gr)){
        gr <- .buildGRFromPamSite(seqnames=seqnames,
                                  pam_site=pam_site,
                                  strand=strand)
    }
    pam_len    <- pamLength(nuclease)
    spacer_len <- spacerLength(nuclease)
    gap_len    <- spacerGap(nuclease)
    pam_side   <- pamSide(nuclease)

    r <-  as.character(BiocGenerics::strand(gr))=='-'
    if (pam_side=="3prime"){
        start    <- start(gr) - gap_len - spacer_len
        end      <- start(gr) + pam_len - 1
        start[r] <- start(gr)[r] - pam_len +1
        end[r]   <- start(gr)[r] + spacer_len + gap_len
    } else {
        start    <- start(gr)
        end      <- start(gr) + pam_len + spacer_len + gap_len - 1
        start[r] <- start(gr)[r] - spacer_len - pam_len - gap_len + 1
        end[r]   <- start(gr)[r]
    }
    gr.new <- .resetGRCoordinates(gr)
    end(gr.new)   <- end
    start(gr.new) <- start
    return(gr.new)
}






#' @title Construct a spacer GRanges from a list of PAM sites 
#' 
#' @description Construct a spacer GRanges from a list of PAM sites
#'     using information stored in a CrisprNuclease object. 
#'     
#' @param gr GRanges object of width 1 specifying the coordinates
#'     of the first nucleotide of the PAM sequences. 
#' @param seqnames Character vector of genomic sequence names.
#'     Ignored if \code{gr} is not NULL.
#' @param pam_site Numeric vector specifying the coordinates of the
#'     first nucleotide of the PAM sequences corresponding to the
#'     protospacers. Ignored if \code{gr} is not NULL.
#' @param strand Character vector specifying the strand of the protospacer. 
#'     Ignored if \code{gr} is not NULL.
#' @param nuclease CrisprNuclease object.
#' @param spacer_len Non-negative integer to overwrite the default spacer
#'     length stored in the CrisprNuclease object.
#' s
#' @return GRanges object representing genomic coordinates of spacer sequences.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples 
#' if (require(GenomicRanges)){
#' gr <- GRanges("chr10",
#'               IRanges(start=c(100,120), width=1),
#'               strand=c("+","-"))
#' getSpacerRanges(gr, nuclease=SpCas9)
#' getSpacerRanges(gr, nuclease=AsCas12a)
#' }
#' 
#' @export
getSpacerRanges <- function(gr=NULL,
                            seqnames=NULL,
                            pam_site=NULL,
                            strand=NULL,
                            nuclease=NULL,
                            spacer_len=NULL
){
    .isCrisprNucleaseOrStop(nuclease)
    if (!is.null(spacer_len)){
        if (length(spacer_len)>1){
            stop("spacer_len must be either NULL or of length 1.")
        }
        message("Overwriting spacerLength(nuclease) with spacer_len")
        spacerLength(nuclease) <- spacer_len
    }
    gr <- .validatePosGrOrNull(gr)
    if (is.null(gr) & 
        (is.null(seqnames) | is.null(pam_site) | is.null(strand))){
        stop("seqnames, pam_site, and strand must be provided if gr=NULL")
    }
    if (is.null(gr)){
        gr <- .buildGRFromPamSite(seqnames=seqnames,
                                  pam_site=pam_site,
                                  strand=strand)
    }
    pam_len    <- pamLength(nuclease)
    spacer_len <- spacerLength(nuclease)
    gap_len    <- spacerGap(nuclease)
    pam_side   <- pamSide(nuclease)

    r <-  as.character(BiocGenerics::strand(gr))=='-'
    if (pam_side=="3prime"){
        start    <- start(gr) - gap_len - spacer_len
        end      <- start(gr) - gap_len - 1
        start[r] <- start(gr)[r] + gap_len +1
        end[r]   <- start(gr)[r] + spacer_len + gap_len
    } else {
        start    <- start(gr) + pam_len + gap_len
        end      <- start(gr) + pam_len + spacer_len + gap_len - 1
        start[r] <- start(gr)[r] - spacer_len - pam_len - gap_len + 1
        end[r]   <- start(gr)[r] - pam_len - gap_len
    }
    gr.new <- .resetGRCoordinates(gr)
    end(gr.new)   <- end
    start(gr.new) <- start
    return(gr.new)
}


#' @title Construct a PAM GRanges from a list of PAM sites 
#' 
#' @description Construct a PAM GRanges from a list of PAM sites
#'     using information stored in a CrisprNuclease object. 
#'     
#' @param gr GRanges object of width 1 specifying the coordinates
#'     of the first nucleotide of the PAM sequences. 
#' @param seqnames Character vector of genomic sequence names.
#'     Ignored if \code{gr} is not NULL.
#' @param pam_site Numeric vector specifying the coordinates of the
#'     first nucleotide of the PAM sequences corresponding to the
#'     PAM sequences. Ignored if \code{gr} is not NULL.
#' @param strand Character vector specifying the strand of the PAM. 
#'     Ignored if \code{gr} is not NULL.
#' @param nuclease CrisprNuclease object.
#' 
#' @return GRanges object representing genomic coordinates of PAM sequences.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' if (require(GenomicRanges)){
#' gr <- GRanges("chr10",
#'               IRanges(start=c(100,120), width=1),
#'               strand=c("+","-"))
#' getPamRanges(gr, nuclease=SpCas9)
#' getPamRanges(gr, nuclease=AsCas12a)
#' }
#'
#' @export
getPamRanges <- function(gr=NULL,
                         seqnames=NULL,
                         pam_site=NULL,
                         strand=NULL,
                         nuclease=NULL
){
    .isCrisprNucleaseOrStop(nuclease)
    gr <- .validatePosGrOrNull(gr)
    if (is.null(gr) & 
        (is.null(seqnames) | is.null(pam_site) | is.null(strand))){
        stop("seqnames, pam_site, and strand must be provided if gr=NULL")
    }
    if (is.null(gr)){
        gr <- .buildGRFromPamSite(seqnames=seqnames,
                                  pam_site=pam_site,
                                  strand=strand)
    }
    r <-  as.character(BiocGenerics::strand(gr))=='-'
    pam_len  <- pamLength(nuclease)
    start    <- start(gr)     
    end      <- start(gr) + pam_len - 1             
    start[r] <- start(gr)[r] - pam_len + 1           
    end[r]   <- start(gr)[r]
  
    gr.new <- .resetGRCoordinates(gr)
    end(gr.new)   <- end
    start(gr.new) <- start
    return(gr.new)
}



.validateProtospacers <- function(protospacers,
                                  object
){
    .isCrisprNucleaseOrStop(object)
    protospacer.len <- protospacerLength(object)
    n <- unique(nchar(protospacers))
    if (n!=protospacer.len){
        stop("provided protospacers are of length ",n,
             ", but it should be ", protospacer.len, 
             " for the provided ", name(object))  
    } 
    return(protospacers)
}





