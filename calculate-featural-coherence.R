dist.align <- read.table("distinctive-alignment.txt", sep="\t", header = TRUE)

# The variable names used here correspond to terminology used in our paper as follows:
# dist.align = "featural coherence"
# Atr = "A (agent-like)" 
# P = "P (patient-like)"
# prime.subj = "Persistence"
# crmv.crct = "Craner's V (with correction)"
# Stock = "language family"

dist.align$Atr <- as.character(dist.align$Atr)
dist.align$P <- as.character(dist.align$P)
dist.align <- subset(dist.align, !(Language=='Chintang' | Language=='Bantawa'))

library(rcompanion)
dist.align$crmv <- rep(NA, nrow(dist.align))
dist.align$crmv.crct <- rep(NA, nrow(dist.align))
for (i in 1:nrow(dist.align)) {
  AtrVec <- as.integer(unlist(strsplit(dist.align$Atr[i], ",")))
  PVec <- as.integer(unlist(strsplit(dist.align$P[i], ",")))
  # if there is only one slot, then there is no independence whatsoever
  if ( length(AtrVec)==1 ) {
    dist.align$cramerv[i] <- NA
  } else {
    mat <- rbind(AtrVec, PVec)
    # For chisq, and therefore for CramerV, we need to remove all columns (slots) in the table that have zero for both A and P
    matPos <- mat[,apply(mat,2,function(mat) !all(mat==0))]
    dist.align$crmv[i] <- round(cramerV(matPos, bias.correct=FALSE),digits=3)
    dist.align$crmv.crct[i] <- round(cramerV(matPos, bias.correct=TRUE),digits=3)
  }
}

write.table(dist.align, file = "category-clustering/distinctive-alignment-crmv.txt", sep="\t", row.names = FALSE, quote = FALSE)

library(ggplot2)
ggplot(dist.align, aes(x=crmv.crct)) + geom_histogram(binwidth = 0.1) + xlab("Featural Coherence: Cramér's V") + ylab("Count")+ theme_bw(base_size = 24)

# Investigate by language family
dist.align$Stock.grouped <- as.character(dist.align$Stock)
dist.align$Stock.grouped[dist.align$Stock %in% names(table(dist.align$Stock))[table(dist.align$Stock) < 10]] <- "Other"
dist.align$Stock.grouped[dist.align$Stock.grouped=="Sino-Tibetan"] <- "Kiranti"# all the STs in this data are Kiranti
dist.align$Stock.grouped[dist.align$Stock.grouped=="Algic"] <- "Algonquian"
dist.align$Stock.grouped <- factor(dist.align$Stock.grouped)
dist.align$Stock.grouped <- factor(dist.align$Stock.grouped, c("Algonquian", "Kiranti", "Other" ))
ggplot(dist.align, aes(x=crmv.crct)) + facet_wrap(~Stock.grouped) + geom_histogram(binwidth = 0.2) + xlab("Featural Coherence: Cramér's V") + ylab("Count" ) + theme_bw(base_size = 24)

# Statistical test: beta regression on Cramer's V corrected
# first need this correction to transform 0,1 vals
dist.align$corr.crmv <- (dist.align$crmv.crct * (length(dist.align$crmv.crct) - 1) + .5) / length(dist.align$crmv.crct)
library(glmmTMB)
glm.model <- glmmTMB(corr.crmv ~ (1|Stock), data=dist.align, family=list(family="beta",link="logit"))
summary(glm.model)
confints <- confint(glm.model)
