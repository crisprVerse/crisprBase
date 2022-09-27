#' @title Extract PAM sequences from target sequences
#' @description Extract PAM sequences from target sequences
#'     (protospacer + PAM) using information stored in a
#'     CrisprNuclease object.
#'  
#'     
#' 
#' @param targets Character vector of target sequences. 
#' @param object \code{CrisprNuclease} corresponding to the target 
#'     sequences.
#' 
#' @return Character vector of PAM sequences of length equal to that of the 
#'     \code{targets} character vector. 
#' 
#' @author Jean-Philippe Fortin
#'
#' @examples
#' data(SpCas9, AsCas12a, package="crisprBase")
#' # Extracting PAM sequences from Cas9 protospacers:
#' targets <- c("AGGTGCTGATTGTAGTGCTGCGG",
#'               "AGGTGCTGATTGTAGTGCTGAGG")
#' extractPamFromTarget(targets, SpCas9)
#' # Extracting PAM sequences from Cas12a targets:
#' targets <- c("TTTAAGGTGCTGATTGTAGTGCTGTGT",
#'              "TTTCAGGTGCTGATTGTAGTGCTGAAA")
#' extractPamFromTarget(targets, AsCas12a)
#' 
#' @export
extractPamFromTarget <- function(targets,
                                 object
){
    .isCrisprNucleaseOrStop(object)
    .validateTargets(targets, object)
    wh <- pamIndices(object)
    pams <- substr(targets,
                   wh[1],
                   wh[length(wh)])
    return(pams)
}



#' @title Extract protospacer sequences from target sequences
#' @description Extract protospacer sequences from target sequences
#'     (protospacer + PAM) using information stored in a
#'     CrisprNuclease object. 
#'     
#' 
#' @param targets Character vector of targets sequences. 
#' @param object \code{CrisprNuclease} corresponding to the targets 
#'     sequences.
#' 
#' @return Character vector of protospacer sequences of 
#'     length equal to that of the 
#'     \code{targets} character vector. 
#' 
#' @author Jean-Philippe Fortin
#'
#' @examples
#' data(SpCas9, AsCas12a, package="crisprBase")
#' # Extracting protospacer sequences from Cas9 targets:
#' targets <- c("AGGTGCTGATTGTAGTGCTGCGG",
#'              "AGGTGCTGATTGTAGTGCTGAGG")
#' extractProtospacerFromTarget(targets, SpCas9)
#' # Extracting protospacer sequences from Cas12a targets:
#' targets <- c("TTTAAGGTGCTGATTGTAGTGCTGTGT",
#'              "TTTCAGGTGCTGATTGTAGTGCTGAAA")
#' extractProtospacerFromTarget(targets, AsCas12a)
#' 
#' @export
extractProtospacerFromTarget <- function(targets,
                                         object
){
    .isCrisprNucleaseOrStop(object)
    .validateTargets(targets, object)
    wh <- spacerIndices(object)
    out <- substr(targets,
                  wh[1],
                  wh[length(wh)])
    return(out)
}





#' @title Construct a target GRanges from a list of PAM sites 
#' 
#' @description Construct a target (protospacer + PAM) GRanges
#'     from a list of PAM sites using information stored in a
#'     CrisprNuclease object. 
#'     
#' @param gr GRanges object of width 1 specifying the coordinates
#'     of the first nucleotide of the PAM sequences. 
#' @param seqnames Character vector of genomic sequence names.
#'     Ignored if \code{gr} is not NULL.
#' @param pam_site Numeric vector specifying the coordinates of the
#'     first nucleotide of the PAM sequences corresponding to the
#'     targets. Ignored if \code{gr} is not NULL.
#' @param strand Character vector specifying the strand of the target. 
#'     Ignored if \code{gr} is not NULL.
#' @param nuclease CrisprNuclease object.
#' @param spacer_len Non-negative integer to overwrite the default spacer
#'     length stored in the CrisprNuclease object.
#' 
#' @return GRanges object representing genomic coordinates of
#'     the target sequences.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples 
#' data(SpCas9, AsCas12a, package="crisprBase")
#' library(GenomicRanges)
#' gr <- GRanges("chr10",
#'               IRanges(start=c(100,120), width=1),
#'               strand=c("+","-"))
#' getTargetRanges(gr, nuclease=SpCas9)
#' getTargetRanges(gr, nuclease=AsCas12a)
#' 
#' 
#' @export
#' @importFrom BiocGenerics start end start<- end<- strand
getTargetRanges <- function(gr=NULL,
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
#' s
#' @return GRanges object representing genomic coordinates of
#'     protospacer sequences.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples 
#' data(SpCas9, AsCas12a, package="crisprBase")
#' library(GenomicRanges)
#' gr <- GRanges("chr10",
#'               IRanges(start=c(100,120), width=1),
#'               strand=c("+","-"))
#' getProtospacerRanges(gr, nuclease=SpCas9)
#' getProtospacerRanges(gr, nuclease=AsCas12a)
#' 
#' @export
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
#' data(SpCas9, AsCas12a, package="crisprBase")
#' library(GenomicRanges)
#' gr <- GRanges("chr10",
#'               IRanges(start=c(100,120), width=1),
#'               strand=c("+","-"))
#' getPamRanges(gr, nuclease=SpCas9)
#' getPamRanges(gr, nuclease=AsCas12a)
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




#' @title Construct a cut site GRanges from a list of PAM sites 
#' 
#' @description Construct a cut site GRanges from a list of PAM sites 
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
#' @return GRanges object representing genomic coordinates of the cut sites.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' data(SpCas9, AsCas12a, package="crisprBase")
#' library(GenomicRanges)
#' gr <- GRanges("chr10",
#'               IRanges(start=c(100,120), width=1),
#'               strand=c("+","-"))
#' getCutSiteRanges(gr, nuclease=SpCas9)
#' getCutSiteRanges(gr, nuclease=AsCas12a)
#'
#' @export
getCutSiteRanges <- function(gr=NULL,
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
    pamSites <- BiocGenerics::start(gr)
    strand <- as.character(strand(gr))
    cutSites <- getCutSiteFromPamSite(pam_site=pamSites,
                                      strand=strand,
                                      nuclease=nuclease)

    gr.new <- .resetGRCoordinates(gr)
    BiocGenerics::start(gr.new) <- BiocGenerics::end(gr.new) <- cutSites
    return(gr.new)
}








#' @title Return cut site coordinates from PAM site coordinates
#' 
#' @description Return cut site coordinates from PAM site coordinates.
#'     
#' @param pam_site Coordinate of the first nucleotide of the PAM sequence.
#' @param strand Either "+" or "-".
#' @param nuclease A \linkS4class{CrisprNuclease} object.
#' 
#' @return numeric vector of cut sites
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' data(SpCas9, package="crisprBase")
#' getCutSiteFromPamSite(pam_site=100, strand="+", nuclease=SpCas9)
#' getCutSiteFromPamSite(pam_site=100, strand="-", nuclease=SpCas9)
#' 
#' @export
getCutSiteFromPamSite <- function(pam_site,
                                  strand,
                                  nuclease=NULL
){
    .isCrisprNucleaseOrStop(nuclease)
    strand <- .validateStrand(strand)
    stopifnot(length(pam_site)==length(strand))
    cut_offset <- cutSites(nuclease, middle=TRUE)
    cut_site <- pam_site
    cut_site[strand=='+'] <- pam_site[strand=='+'] + cut_offset
    cut_site[strand=='-'] <- pam_site[strand=='-'] - cut_offset
    return(cut_site)
}




#' @title Return optimal editing site coordinates from PAM site coordinates
#' 
#' @description Return optimal editing site coordinates from PAM 
#'     site coordinates. 
#'     
#' @param pam_site Coordinate of the first nucleotide of the PAM sequence.
#' @param strand Either "+" or "-".
#' @param baseEditor A \linkS4class{BaseEditor} object.
#' @param substitution String indicating which substitution
#'     should be used to estimate the optimal editing position. 
#'     E.g. "C2T" will return the optimal editing position for C to T
#'     editing. 
#' 
#' @return numeric vector of editing sites.
#' 
#' @author Jean-Philippe Fortin
#' 
#' @examples
#' data(BE4max, package="crisprBase")
#' getEditingSiteFromPamSite(pam_site=100, strand="+", baseEditor=BE4max, "C2T")
#' 
#' @export
getEditingSiteFromPamSite <- function(pam_site,
                                      strand,
                                      baseEditor=NULL,
                                      substitution=NULL
){
    .isBaseEditorOrStop(baseEditor)
    strand <- .validateStrand(strand)
    stopifnot(length(pam_site)==length(strand))
    if (is.null(substitution)){
        stop("[getEditingSiteFromPamSite] 'substitution' must be specified. ")
    }
    if (length(substitution)!=1){
        stop("[getEditingSiteFromPamSite] 'substitution' must have length 1.")
    }
    ws <- editingWeights(baseEditor)
    if (!substitution %in% rownames(ws)){
        stop("[getEditingSiteFromPamSite] 'substitution' not found in the",
             " matrix of editing weights.")
    }
    ws <- ws[substitution,,drop=FALSE]
    maxEditing <- colnames(ws)[which.max(ws)]
    editing_offset <- as.integer(maxEditing)
    editing_site <- pam_site
    editing_site[strand=='+'] <- editing_site[strand=='+'] + editing_offset
    editing_site[strand=='-'] <- editing_site[strand=='-'] - editing_offset
    return(editing_site)
}







.validateTargets <- function(targets,
                             object
){
    .isCrisprNucleaseOrStop(object)
    target.len <- targetLength(object)
    n <- unique(nchar(targets))
    if (n!=target.len){
        stop("provided targets are of length ",n,
             ", but it should be ", target.len, 
             " for the provided ", nucleaseName(object))  
    } 
    return(targets)
}





