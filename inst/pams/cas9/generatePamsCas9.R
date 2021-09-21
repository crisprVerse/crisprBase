#Script to create Cas9 PAMs:
pams <- c("AGG", "AAG", "AGA", "CGG",
          "CAG", "CGA", "GGG", "GAG", 
          "GGA", "TGG", "TAG", "TGA")
cas9.pams <- data.frame(PAM=pams,
                        Tier="Tier1",
                        stringsAsFactors=FALSE)
# Now let's generate the second tier:
choices <- expand.grid(rep(list(c('A', 'G', 'T', 'C')), 3))
choices <- do.call(paste0, choices)
choices <- setdiff(choices, pams)
temp <- data.frame(PAM=choices,
                   Tier="Tier2",
                   stringsAsFactors=FALSE)
cas9.pams <- rbind(cas9.pams, temp)
save(cas9.pams,
     file="cas9.pams.rda")
