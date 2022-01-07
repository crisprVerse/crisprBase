#' @importFrom methods validObject
.nickaseToNuclease <- function(x){
    out <- new("Nuclease")
    targetType(out) <- "DNA"
    nucleaseName(out) <- paste0(nickaseName(x),
                                " converted to nuclease.")
    metadata(out) <- metadata(x)
    out@motifs <- motifs(x)
    weights(out) <- weights(x)
    cuts <- rbind(cutSites(x),cutSites(x))
    rownames(cuts) <- c("fwd", 'rev')
    colnames(cuts) <- motifs(x, as.character=TRUE)
    out@cutSites <- cuts
    stopifnot(validObject(out))
    return(out)
}


#' @importFrom methods validObject
.crisprNickaseToCrisprNuclease <- function(x){
    out <- new("CrisprNuclease")
    targetType(out) <- "DNA"
    nucleaseName(out) <- paste0(nickaseName(x),
                                " converted to a CrisprNuclease.")
    metadata(out) <- metadata(x)
    out@motifs <- motifs(x)
    weights(out) <- weights(x)
    cuts <- rbind(cutSites(x),cutSites(x))
    rownames(cuts) <- c("fwd", 'rev')
    colnames(cuts) <- motifs(x, as.character=TRUE)
    out@cutSites <- cuts

    pamSide(out) <- pamSide(x)
    spacerGap(out) <- spacerGap(x)
    spacerLength(out) <- spacerLength(x)
    return(out)
}



