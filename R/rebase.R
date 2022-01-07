#' @importFrom stringr str_match
.checkRebaseMotif <- function(rebaseMotif,
                              type=c("nuclease", "nickase")
){
    type <- match.arg(type)
    choices_nickase <- c(names(Biostrings::IUPAC_CODE_MAP),
                         "^", "(", ")",
                         0:9,
                         "-")
    choices_nuclease <- c(choices_nickase, "/")
    if (type=="nuclease"){
        choices <- choices_nuclease
    } else {
        choices <- choices_nickase
    }

    symbols <- strsplit(rebaseMotif, split="")[[1]]
    if(!all(symbols %in% choices)){
        offenders <- unique(setdiff(symbols, choices))
        offenders <- paste0(offenders, collapse=",")
        stop("The following characters are not allowed: ",
             offenders)
    }

    #Checking if Rebase motif follows expected pattern:
    if (type=="nuclease"){
        pattern="(\\([0-9]+/[0-9]+\\))?([A-Z]*\\^?[A-Z]*)(\\([0-9]+/[0-9]+\\))?"
    } else {
        pattern="(\\([0-9]+\\))?([A-Z]*\\^?[A-Z]*)(\\([0-9]+\\))?"
    }
    
    if (str_match(rebaseMotif, pattern)[1]!=rebaseMotif){
        stop("Motif does not follow the notation needed for ",
             "a recognition site. See details.")
    }

    if (.isCutUpstream(rebaseMotif, type=type) & 
        .isCutDownstream(rebaseMotif, type=type)){
        stop("For nucleases cleaving outside of the recognition site, 
             only one set of parentheses should be provided. 
             Example: (9/10)ACCTG or ACCTG(9/10) are valid, 
             but not (9/10)ACCTG(9/10).")
    }
    if ((.isCutUpstream(rebaseMotif , type=type)|
         .isCutDownstream(rebaseMotif,  type=type)) &
         .isCutWithin(rebaseMotif)){
        stop("Concurrent cleavage outside of the recognition site and
             within the recognition sequence are not supported at the moment.
             Example: (9/10)ACCTG is valid, but (9/10)AC^CTG is invalid.")
    }
}

.checkRebaseMotifs <- function(rebaseMotifs,
                               type=c("nuclease", "nickase")
){
    type <- match.arg(type)
    x=lapply(rebaseMotifs, .checkRebaseMotif, type=type)
}


.isCutUpstream <- function(rebaseMotif,
                           type=c("nuclease", "nickase")
){  
    type <- match.arg(type)
    if (type=="nuclease"){
        out <- grepl("^\\([0-9]+/[0-9]+\\)", rebaseMotif)
    } else {
        out <- grepl("^\\([0-9]+\\)", rebaseMotif)
    }
    return(out)
}

.isCutDownstream <- function(rebaseMotif,
                             type=c("nuclease", "nickase")
){
    type <- match.arg(type)
    if (type=="nuclease"){
        out <- grepl("\\([0-9]+/[0-9]+\\)$", rebaseMotif)
    } else {
        out <- grepl("\\([0-9]+\\)$", rebaseMotif)
    }
    return(out)
}

.isCutWithin <- function(rebaseMotif){
    out <- grepl("\\^", rebaseMotif)
    return(out)
}

.isCutNotSpecified <- function(rebaseMotif,
                               type=c("nuclease", "nickase")
){
    type <- match.arg(type)
    !.isCutUpstream(rebaseMotif, type=type) & 
    !.isCutDownstream(rebaseMotif, type=type) &
    !.isCutWithin(rebaseMotif) 
}



#' @importFrom stringr str_match
.extractSequencesFromRebaseMotifs <- function(rebaseMotifs){
    #pattern <- "([A-Z]+\\^?[A-Z]+)"
    pattern <- "([A-Z]+)"
    seqs <- vapply(rebaseMotifs, function(xx){
        xx <- gsub("\\^", "", xx)
        return(str_match(xx, pattern)[1])
    }, FUN.VALUE="a")
    names(seqs) <- NULL
    seqs <- DNAStringSet(seqs)
    return(seqs)
}


.extractCutSitesFromRebaseMotifs <- function(rebaseMotifs,
                                             type=c("nuclease", "nickase"),
                                             strand=c("both", "fwd", "rev")
){
    strand <- match.arg(strand)
    cutSites <- lapply(rebaseMotifs, .getCutSites, type=type)
    cutSites <- do.call("cbind", cutSites)
    cutSites <- as.matrix(cutSites)
    rownames(cutSites) <- c("fwd", "rev")
    motifs <- .extractSequencesFromRebaseMotifs(rebaseMotifs)
    colnames(cutSites) <- as.character(motifs)
    if (strand=="fwd"){
        cutSites <- cutSites["fwd",,drop=FALSE]
    } else if (strand=="rev"){
        cutSites <- cutSites["rev",,drop=FALSE]
    }
    return(cutSites)
}




.getCutSites <- function(rebaseMotif,
                         type=c("nuclease", "nickase")
){
    type <- match.arg(type)
    if (.isCutWithin(rebaseMotif)){
        cutSites <- .getCutSites_within(rebaseMotif,
                                        type=type)
    } else if (.isCutUpstream(rebaseMotif,type=type)){
        cutSites <- .getCutSites_upstream(rebaseMotif,
                                          type=type)
    } else if (.isCutDownstream(rebaseMotif, type=type)){
        cutSites <- .getCutSites_downstream(rebaseMotif,
                                            type=type)
    } else if (.isCutNotSpecified(rebaseMotif, type=type)){
        cutSites <- .getCutSites_unspecified(rebaseMotif)
    }
    return(cutSites)
}


.getCutSites_unspecified <- function(rebaseMotif){
    out <- c(fwd=NA, rev=NA)
    return(out)
}

.getCutSites_upstream <- function(rebaseMotif,
                                  type=c("nuclease", "nickase")
){
    type <- match.arg(type)
    stopifnot(.isCutUpstream(rebaseMotif, type=type))
    if (type=="nuclease"){
        pattern <- "^\\([0-9]+/[0-9]+\\)"
    } else {
        pattern <- "^\\([0-9]\\)"
    }
    pos <- regexpr(pattern, rebaseMotif)
    cut_info <- substr(rebaseMotif, pos, pos+attr(pos, "match.length")-1)
    cut_info <- gsub("\\(|\\)", "", cut_info)
    if (type=="nuclease"){
        cut_pos <- as.integer(strsplit(cut_info, split="/")[[1]][[1]])
        cut_neg <- as.integer(strsplit(cut_info, split="/")[[1]][[2]])
    } else {
        cut_pos <- cut_neg <- as.integer(cut_info)
    }
    out <- c(fwd=-cut_pos, rev=-cut_neg)
    return(out)
}

.getCutSites_downstream <- function(rebaseMotif,
                                    type=c("nuclease", "nickase")
){
    type <- match.arg(type)
    stopifnot(.isCutDownstream(rebaseMotif, type=type))
    seq <- .extractSequencesFromRebaseMotifs(rebaseMotif)
    motifLength <- BiocGenerics::width(seq)
    if (type=="nuclease"){
        pattern <- "\\([0-9]+/[0-9]+\\)$"
    } else {
        pattern <- "\\([0-9]+\\)$"
    }
    pos <- regexpr(pattern, rebaseMotif)
    cut_info <- substr(rebaseMotif,
                       pos,
                       pos+attr(pos, "match.length")-1)
    cut_info <- gsub("\\(|\\)", "", cut_info)
    if (type=="nuclease"){
        cut_pos <- as.integer(strsplit(cut_info, split="/")[[1]][[1]])
        cut_neg <- as.integer(strsplit(cut_info, split="/")[[1]][[2]])
    } else {
        cut_pos <- cut_neg <- as.integer(cut_info)
    }
    out <- c(fwd=cut_pos + motifLength,
             rev=cut_neg + motifLength)
    return(out)
}


.getCutSites_within <- function(rebaseMotif,
                                type=c("nuclease", "nickase")
){
    type <- match.arg(type)
    stopifnot(.isCutWithin(rebaseMotif))
    seq <- .extractSequencesFromRebaseMotifs(rebaseMotif)
    motifLength <- BiocGenerics::width(seq)
    chars <- strsplit(rebaseMotif, split="")[[1]]
    pos <- which(chars=="^")
    if (type=="nuclease"){
        cut_pos <- pos-1
        cut_neg <- motifLength - cut_pos
    } else {
        cut_pos <- cut_neg <- pos-1
    }
    out <- c(fwd=cut_pos, rev=cut_neg)
    return(out)
}