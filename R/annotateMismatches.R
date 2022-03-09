#' @title Annotate mismatches between spacer and protospacer sequences
#' @description Annotate mismatches between spacer and protospacer sequences.
#' 
#' @param spacers A character vector specifying spacer sequences (gRNA).
#' @param protospacers A character vector specifying protospacer sequences
#'     (target DNA).
#' @param rnase Is it for an RNAse? FALSE by default.
#'     If TRUE, spacers and protospacers are expected to be
#'     the reverse complement of each other. 
#' 
#' @return A data.frame storing spacer and protospacer columns,
#'     as well as number of mismatches, and positions for 
#'     the different mismatches, if any.
#'     Positions are relative to the 5' end of the 
#'     spacer sequences. For RNAses (e.g. CasRx), this means
#'     that a mismatch at position 1 corresponds to the last
#'     nucleotide of the protospacer sequence. 
#' 
#' 
#' @examples
#' spacers <- c("CCGGAGCGAGTTGCAGTAAGCAG",
#'     "GCCGGAGCGAGTTGCAGTAAGCA", 
#'     "GGCCGGAGCGAGTTGCAGTAAGC")
#' 
#' protospacers=c("CTGCTTACTGCAACTCGCTCTGG",
#'     "TGCTTAATGCAACCCGCTCCGGC", 
#'     "GCTTACTGCAACTCGCTCCGGCC")
#' 
#' ann <- annotateMismatches(spacers,
#'     protospacers,
#'     rnase=TRUE)
#' 
#' @author Jean-Philippe Fortin
#' 
#' @export
#' @importFrom Biostrings DNAStringSet reverseComplement
annotateMismatches <- function(spacers,
                               protospacers,
                               rnase=FALSE
){

    if (is(spacers, "DNAStringSet")){
        spacers <- as.character(spacers)
    }
    if (is(protospacers, "DNAStringSet")){
        protospacers <- as.character(protospacers)
    }
    if (length(spacers)!=length(protospacers)){
        stop("spacers and protospacers must have the same length.")
    }
    if (!all.equal(nchar(spacers),nchar(protospacers))){
        stop("spacers and protospacers must have words of the same length.")
    }

    max_mismatches <- 3
    results <- data.frame(spacer=spacers,
                          protospacer=protospacers)
    words1 <- spacers
    words2 <- protospacers
    if (rnase){
        words2 <- reverseComplement(DNAStringSet(words2))
        words2 <- as.character(words2)
    }
    results$n_mismatches <- vapply(seq_along(words1), function(i){
        cost <- c(sub=1, del=1000, ins=1000)
        utils::adist(words1[i],
                     words2[i],
                     cost=cost)[[1]]
    }, FUN.VALUE=1)


    max_mismatches <- max(max(results$n_mismatches),
                          max_mismatches)
    words1 <- Biostrings::DNAStringSet(words1)
    words2 <- Biostrings::DNAStringSet(words2)
    x1 <- as.matrix(words1)
    x2 <- as.matrix(words2)
    whs <- apply(x1!=x2,1,which, simplify=FALSE)
    mm <- lapply(whs, function(wh){
        wh <- c(wh, rep(NA, max_mismatches-length(wh)))
        return(wh)
    }) 
    mm <- do.call(rbind, mm)
    cols <- paste0("mm", seq_len(max_mismatches))
    colnames(mm) <- cols
    results <- cbind(results, mm)
    return(results)
}
