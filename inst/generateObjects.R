library(crisprBase)
library(devtools)


SpCas9 <- CrisprNuclease("SpCas9",
                         pams=c("(3/3)NGG", "(3/3)NAG", "(3/3)NGA"),
                         weights=c(1, 0.2593, 0.0694),
                         metadata="Wildtype Streptococcus pyogenes Cas9 (SpCas9) nuclease",
                         spacer_side="5prime",
                         spacer_length=20)

SaCas9 <- CrisprNuclease("SaCas9",
                         pams=c("(3/3)NNGRRT"),
                         metadata="Wildtype Staphylococcus aureus Cas9 (SaCas9) nuclease",
                         spacer_side="5prime",
                         spacer_length=21)

SpGCas9 <- CrisprNuclease("SpGCas9",
                          pams=c("(3/3)NGN"),
                          metadata="Engineered Streptococcus pyogenes Cas9 (SpCas9) nuclease named SpG",
                          spacer_side="5prime",
                          spacer_length=20)

AsCas12a <- CrisprNuclease("AsCas12a",
                           pams="TTTV(18/23)",
                           metadata="Wildtype Acidaminococcus Cas12a (AsCas12a) nuclease.",
                           spacer_side="3prime",
                           spacer_length=23)


load("pams/cas12a/cas12a.pams.rda")
pams <- cas12a.pams[, c("PAM", "Score_Doench")]
pams <- pams[!is.na(pams$Score_Doench),]
pams$Score_Doench[pams$PAM %in% c("TTTC","TTTA", "TTTG")] <- 1
pams <- pams[order(-pams$Score_Doench),,]
motifs <- paste0(pams$PAM, "(18/23)")
enAsCas12a <- CrisprNuclease("enAsCas12a",
                             pams=motifs,
                             weights=pams$Score_Doench,
                             metadata="Enhanced Acidaminococcus Cas12a (AsCas12a) nuclease.",
                             spacer_side="3prime",
                             spacer_length=23)


#dir.create("../data")
use_data(SpCas9,
         SaCas9,
         SpGCas9,
         AsCas12a,
         enAsCas12a,
         compress="xz", internal=FALSE, overwrite=TRUE)


# Enzymes:
EcoRI <- Nuclease("EcoRI",
                  motifs=c("G^AATTC"),
                  metadata="EcoRI restriction enzyme")

SmaI <- Nuclease("SmaI",
                  motifs=c("CCC^GGG"),
                  metadata="SmaI restriction enzyme")

HgaI <- Nuclease("HgaI",
                 motifs=c("GACGC(5/10)"),
                 metadata="HgaI restriction enzyme")

PfaAI <- Nuclease("PfaAI",
                  motifs=c("G^GYRCC"),
                  metadata="PfaAI restriction enzyme")








