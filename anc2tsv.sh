#!/usr/bin/env Rscript
arguments <- commandArgs(trailingOnly=TRUE)
fname <- arguments[1]
mrca <- read.csv(fname)
root2tip <- round(mrca$value[1])
beast <- round(mrca$value[2])
beastl <- round(mrca$lower.CI[2])
beasth <- round(mrca$upper.CI[2])
baseDate <- as.Date('2000-01-01')
r <- as.character(baseDate+root2tip)
b <- as.character(baseDate + beast)
bl <- as.character(baseDate + beastl)
bh <- as.character(baseDate + beasth)
mtx <- matrix(c(r,b,bl,bh),ncol=4,byrow=TRUE)
colnames(mtx) <- c("root2tip","beast","beast_low","beast_high")
ancTab <- as.table(mtx)
rownames(ancTab) <- NULL
write.table(ancTab,file="",sep="\t",row.names=FALSE,quote=FALSE)
