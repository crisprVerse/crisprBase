library(stringr)
library(crisprBase)
data <- readLines("allenz.txt")
wh   <- which(grepl("References", data))
data <- data[1:wh]
wh1 <- which(grepl("<1>", data))
wh5 <- which(grepl("<5>", data))
enzymes <- data.frame(name=data[wh1])
enzymes$motif <- data[wh5] 
enzymes$name  <- gsub("<1>","", enzymes$name)
enzymes$motif <- gsub("<5>","", enzymes$motif)
enzymes <- enzymes[!grepl("\\?",enzymes$motif),]
enzymes <- enzymes[!grepl("\\-",enzymes$motif),]
good <- str_count(pattern="\\(",enzymes$motif)<=1
enzymes <- enzymes[good,]
enzymes <- enzymes[!grepl(",", enzymes$motif),]


checks <- lapply(enzymes$motif, function(x){
    crisprBase:::.checkRebaseMotif(x)
})
#valid <- vapply(checks, function(x) x$valid, FUN.VALUE=TRUE)
#enzymes <- enzymes[valid,]
enzymes$motifSequence <- gsub("\\^", "", enzymes$motif)
enzymes$motifSequence <- gsub("\\([0-9]+/[0-9]+\\)", "",
                              enzymes$motifSequence)
resEnzymes <- enzymes


nucs <- lapply(1:nrow(resEnzymes), function(i){
    Nuclease(name=resEnzymes$name[i],
             info=paste0(resEnzymes$name[i], " from REBASE."),
             motifs=resEnzymes$motif[i])
})
names(nucs) <- resEnzymes$name
restrictionEnzymes <- nucs
save(restrictionEnzymes,
     file="../../data/restrictionEnzymes.rda")

