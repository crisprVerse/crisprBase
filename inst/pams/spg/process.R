library(readxl)
data <- read_excel("SpG activity scores by PAM.xlsx")
data <- as.data.frame(data)
colnames(data) <- c("PAM", "Score")
spg.pams <- data
save(spg.pams, file="spg.pams.rda")
