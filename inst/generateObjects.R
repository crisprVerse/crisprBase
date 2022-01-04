library(crisprBase)
library(devtools)


SpCas9 <- CrisprNuclease("SpCas9",
                         pams=c("(3/3)NGG", "(3/3)NAG", "(3/3)NGA"),
                         weights=c(1, 0.2593, 0.0694),
                         metadata=list(description="Wildtype Streptococcus pyogenes Cas9 (SpCas9) nuclease"),
                         pam_side="3prime",
                         spacer_length=20)

SaCas9 <- CrisprNuclease("SaCas9",
                         pams=c("(3/3)NNGRRT"),
                         metadata=list(description="Wildtype Staphylococcus aureus Cas9 (SaCas9) nuclease"),
                         pam_side="3prime",
                         spacer_length=21)

SpGCas9 <- CrisprNuclease("SpGCas9",
                          pams=c("(3/3)NGN"),
                          metadata=list(description="Engineered Streptococcus pyogenes Cas9 (SpCas9) nuclease named SpG"),
                          pam_side="3prime",
                          spacer_length=20)

AsCas12a <- CrisprNuclease("AsCas12a",
                           pams="TTTV(18/23)",
                           metadata=list(description="Wildtype Acidaminococcus Cas12a (AsCas12a) nuclease."),
                           pam_side="5prime",
                           spacer_length=23)

CasRx <- CrisprNuclease("CasRx",
                        targetType="RNA",
                        pams="N",
                        metadata=list(description="Cas13d-NLS from Ruminococcus flavefaciens strain XPD3002.",
                                      doi="10.1016/j.cell.2018.02.033"),
                        pam_side="3prime",
                        spacer_length=22)



load("pams/cas12a/cas12a.pams.rda")
pams <- cas12a.pams[, c("PAM", "Score_Doench")]
pams <- pams[!is.na(pams$Score_Doench),]
pams$Score_Doench[pams$PAM %in% c("TTTC","TTTA", "TTTG")] <- 1
pams <- pams[order(-pams$Score_Doench),,]
motifs <- paste0(pams$PAM, "(18/23)")
enAsCas12a <- CrisprNuclease("enAsCas12a",
                             pams=motifs,
                             weights=pams$Score_Doench,
                             metadata=list(description="Enhanced Acidaminococcus Cas12a (AsCas12a) nuclease."),
                             pam_side="5prime",
                             spacer_length=23)


#Generate base editor
ws <- t(read.csv("../inst/be/b4max.csv"))
pos <- ws["Position",]-21
colnames(ws) <- pos
ws <- ws[-c(match("Position", rownames(ws))),]
load("../data/SpCas9.rda")
BE4max <- CrisprNucleaseBaseEditor(SpCas9,
                                   baseEditorName="BE4max",
                                   editingStrand="original",
                                   editingWeights=ws)
metadata(BE4max)$description_base_editor <- "BE4max cytosine base editor."


#dir.create("../data")
use_data(SpCas9,
         SaCas9,
         SpGCas9,
         AsCas12a,
         enAsCas12a,
         CasRx,
         BE4max,
         compress="xz", internal=FALSE, overwrite=TRUE)





# Enzymes:
# EcoRI <- Nuclease("EcoRI",
#                   motifs=c("G^AATTC"),
#                   metadata="EcoRI restriction enzyme")

# SmaI <- Nuclease("SmaI",
#                   motifs=c("CCC^GGG"),
#                   metadata="SmaI restriction enzyme")

# HgaI <- Nuclease("HgaI",
#                  motifs=c("GACGC(5/10)"),
#                  metadata="HgaI restriction enzyme")

# PfaAI <- Nuclease("PfaAI",
#                   motifs=c("G^GYRCC"),
#                   metadata="PfaAI restriction enzyme")








