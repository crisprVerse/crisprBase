load("cas12a.pams.rda")
load("doench/encas12a.scores.doench.rda")
wh <- match(cas12a.pams$PAM, encas12a.scores.doench$PAM)
cas12a.pams$Score_Doench <- encas12a.scores.doench$Score[wh]

# Rescaling Doench score:
#common <- which(!is.na(cas12a.pams$Score) & !is.na(cas12a.pams$Score_Doench))
#top    <- max(1-cas12a.pams$Score_Doench[common], na.rm=TRUE)
#bottom <- max(1-cas12a.pams$Score[common], na.rm=TRUE)
#factor <- top/bottom
#temp=1-cas12a.pams$Score_Doench
#temp <- temp/factor
#temp <- 1-temp
#cas12a.pams$Score_Doench <- temp


#Filling the gaps:
#wh <- which(is.na(cas12a.pams$Score_Doench))
#cas12a.pams$Score_Doench[wh] <- cas12a.pams$Score[wh]
save(cas12a.pams, file="cas12a.pams.rda")

#col <- as.numeric(as.factor(cas12a.pams$Tier))
#plot(cas12a.pams$Score,cas12a.pams$Score_Doench,col=col, pch=20)
#abline(a=0,b=1)
