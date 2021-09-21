# Script to generate Cas12a PAM sequences from PMID:30742127
library(readxl)
library(tidyverse)
data <- read_excel("41587_2018_11_MOESM4_ESM.xlsx") %>% as.data.frame
data <- data[,c(1,2,8)]
colnames(data) <- c("PAM","Score", "Tier")
data <- data[-1,]
data <- data[1:256,]
data$Score <- as.numeric(data$Score)
data$Score <- data$Score/max(data$Score)
data$Tier[data$Tier=="-"] <- "Tier3"
data$Tier <- gsub("tier ", "Tier",data$Tier)
data <- data[order(data$Tier,-data$Score, data$PAM),]
rownames(data) <- NULL
cas12a.pams <- data
#save(cas12a.pams, file="../../../data/cas12a.pams.rda")
save(cas12a.pams,
     file="cas12a.pams.rda")
