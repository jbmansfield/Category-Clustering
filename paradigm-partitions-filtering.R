# Some filtering and plotting of paradigm partitions
paradigm.partitions <- read.table("paradigm-partition-probs.txt", header = TRUE, sep = "\t")

# The variable names used here correspond to terminology used in our paper as follows:
# Atr, A = "A (agent-like)" 
# P = "P (patient-like)"
# SLOTS = "Positions available"
# SLOTS.USED = "Positions.used"
# Category = "grammatical role (A or P)"
# CLUST.INDEX = "paradigmatic alignment index"

ppp <- paradigm.partitions
ppp <- subset(ppp, !(Language=='Chintang' | Language=='Bantawa'))
ppp <- subset(ppp, (Category=='Atr' | Category=='P'))
levels(ppp$Category)[levels(ppp$Category)=="Atr"] <- "A"
ppp <- droplevels(ppp)

library(dplyr)
langs <- ppp %>%
  select(Language, SLOTS, Category) %>%
  group_by(Language) %>% summarise(SLOTS = max(SLOTS), CATS=paste(Category, collapse=" "))

library(gridExtra)
library(ggplot2)
p <- ggplot(data=langs, aes(x=SLOTS)) + geom_histogram(binwidth=1) + xlab("(a) Positions available") + ylab("Count") + scale_x_continuous(breaks=seq(1, 13, by=2), limits=c(0,14)) + theme_bw(base_size = 24)
q <- ggplot(data=ppp, aes(x=SLOTS.USED)) + geom_histogram(binwidth=1) + xlab("(b) Positions used") + ylab("Count") + scale_x_continuous(breaks=seq(1, 13, by=2), limits=c(0,14)) + theme_bw(base_size = 24)
grid.arrange(p, q, ncol=1)

# Now for alignment probabilities, we only look at languages with multiple known slots
ppp <- subset(ppp, SLOTS > 1)
ppp <- droplevels(ppp)

# And drop a couple of extraneous columns
ppp$HOST.CAT <- NULL
ppp$Position.binned5 <- NULL

write.table(ppp, file="category-clustering/paradigm-partitions-probs_filt.txt", sep="\t", row.names = FALSE, quote = FALSE)

library(ggplot2)
ggplot(ppp, aes(x=CLUST.INDEX)) + geom_histogram(binwidth = 0.1) + xlab("Alignment Index") + ylab("Count")
ggplot(data=subset(ppp, !is.na(Category)), aes(x=CLUST.INDEX)) + facet_wrap(~Category) + geom_histogram(binwidth = 0.2) + xlab("Paradigmatic alignment index") + ylab("Count") + theme_bw(base_size = 24)
