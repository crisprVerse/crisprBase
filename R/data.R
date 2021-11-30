#' List of Nuclease objects representing common restriction enzymes
#'
#' List of Nuclease objects representing common restriction enzymes
#'     from REBASE database. 
#' 
#' @format List of Nuclease objects. 
#' 
#' @details List of Nuclease objects representing common restriction enzymes
#'     from REBASE database. 
#' @usage data(restrictionEnzymes, package="crisprBase")
"restrictionEnzymes"


#' SpCas9 CrisprNuclease object
#'
#' CrisprNuclease object for the wildtype Streptococcus pyogenes Cas9 (SpCas9)
#'     nuclease.
#' 
#' @format CrisprNuclease object.
#' 
#' @details The SpCas9 nuclease recognizes NGG PAM sequences. Spacer
#'     sequences must be located upstream of PAM sequences.
#' @usage data(SpCas9, package="crisprBase")
"SpCas9"




#' SaCas9 CrisprNuclease object
#'
#' CrisprNuclease object for the wildtype Staphylococcus aureus Cas9 (SaCas9)
#'     nuclease.
#' 
#' @format CrisprNuclease object.
#' 
#' @details The AsCas9 nuclease recognizes NNGRRT PAM sequences. Spacer
#'     sequences must be located upstream of PAM sequences.
#' @usage data(SaCas9, package="crisprBase")
"SaCas9"


#' SpGCas9 CrisprNuclease object
#'
#' CrisprNuclease object for the engineered Streptococcus pyogenes Cas9
#'     SpG nuclease.
#' 
#' @format CrisprNuclease object.
#' 
#' @details The SpGCas9 nuclease recognizes NGN PAM sequences. Spacer
#'     sequences must be located upstream of PAM sequences.
#' @usage data(SpGCas9, package="crisprBase")
"SpGCas9"


#' AsCas12a CrisprNuclease object
#'
#' CrisprNuclease object for the Wildtype Acidaminococcus Cas12a (AsCas12a)
#'     nuclease.
#' 
#' @format CrisprNuclease object.
#' @details The AsCas12a nuclease recognizes TTTV PAM sequences. Spacer
#'     sequences must be located downstream of PAM sequences.
#' @usage data(AsCas12a, package="crisprBase")
"AsCas12a"


 

#' enAsCas12a CrisprNuclease object
#'
#' CrisprNuclease object for the Enhanced Acidaminococcus Cas12a (AsCas12a)
#'     nuclease.
#' 
#' @format CrisprNuclease object.
#' @details The enAsCas12a nuclease recognizes an extended set of PAM sequences
#'     beyong the canonical TTTV sequence for AsCas12a. Spacer sequences must 
#'     be located downstream of PAM sequences.
#' @usage data(enAsCas12a, package="crisprBase")
"enAsCas12a"



