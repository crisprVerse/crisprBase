###Create all possibilities for 2 genes:
prefix1 <- "ACCACC"
prefix2 <- "ATGACC"
suffix1 <- "GCTGATTGCTAGGACAGCTAGAGGTCG"
suffix2 <- "GCTGATTGCTAGGAGGGCTAAGAGCCG"
pams <- expand.grid(rep(list(c('A', 'G', 'T', 'C')), 4))
pams <- do.call(paste0, pams)
seqs1 <- paste0(prefix1, pams,suffix1)
seqs2 <- paste0(prefix2, pams,suffix2)
names1 <- paste0(">MYL1", 1:length(seqs1))
names2 <- paste0(">MYC1_", 1:length(seqs2))
df1 <- cbind(names1,seqs1)
df2 <- cbind(names2,seqs2)
df1 <- as.vector(t(df1))
df2 <- as.vector(t(df2))
writeLines(df1[1:256], "files/myl1_myc1_11.fasta")
writeLines(df1[1:256+256], "files/myl1_myc1_12.fasta")
writeLines(df2[1:256], "files/myl1_myc1_21.fasta")
writeLines(df2[1:256+256], "files/myl1_myc1_22.fasta")

###Reading results in:
cols <- c("spacer", "score", "pam")
results11 <- read.csv("files/results11.txt", sep="\t")
results12 <- read.csv("files/results12.txt", sep="\t")
results1 <- rbind(results11, results12)
results1 <- results1[,c("sgRNA.Sequence","On.Target.Efficacy.Score","PAM.Sequence")]
colnames(results1) <- cols
results1 <- results1[results1$spacer=="GCTGATTGCTAGGACAGCTAGAG",]
results1 <- results1[order(results1$pam),]

results21 <- read.csv("files/results21.txt", sep="\t")
results22 <- read.csv("files/results22.txt", sep="\t")
results2 <- rbind(results21, results22)
results2 <- results2[,c("sgRNA.Sequence","On.Target.Efficacy.Score","PAM.Sequence")]
colnames(results2) <- cols
results2 <- results2[results2$spacer=="GCTGATTGCTAGGAGGGCTAAGA",]
results2 <- results2[order(results2$pam),]

results1$score <- results1$score/max(results1$score)
results2$score <- results2$score/max(results2$score)
#plot(results1$score, results2$score, ylim=c(0,1), xlim=c(0,1))
#abline(a=0,b=1, lty=3)
results1 <- results1[,c("score", "pam")]
results2 <- results2[,c("score", "pam")]
results <- results1
results$score <- (results1$score+results2$score)/2 
colnames(results) <- c("Score", "PAM")
results <- results[,2:1]
results <- results[order(-results$Score),]
rownames(results) <- NULL
encas12a.scores.doench <- results
save(encas12a.scores.doench, file="encas12a.scores.doench.rda")
