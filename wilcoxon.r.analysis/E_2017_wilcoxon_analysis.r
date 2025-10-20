#######################################
## Wilcoxon results with GLM results ##
#######################################

setwd("~/Documents/josh_r_analysis/fst/Github_E.2017.analysis/wilcoxon.r.analysis/")
library(ggplot2)
library(dplyr)
library(ggpubr)

#load SNP glm results with metadata
snp.gff.fulloverlap.glm <- read.delim("snp.gff.fulloverlap.slim.glm.txt", sep = "\t", header=TRUE)
#load candidate gene and SNP filtering tables
cand.filt.table <- read.table("cand.filt.table.txt", header=TRUE)
candSNP.filt.table <- read.table("candSNP.filt.table.txt", header=TRUE)

## For SNPs within candidate genes, as determined by locations of candidate SNPs
# Run the Wilcox Rank Sum test between candidates and matched sets for n iterations
# Initialize results tables
wilcox.glm.t1_2.p <- data.frame(matrix(NA, nrow = 1000, ncol = 4))
wilcox.glm.t2_3.p <- data.frame(matrix(NA, nrow = 1000, ncol = 4))
wilcox.glm.t3_4.p <- data.frame(matrix(NA, nrow = 1000, ncol = 4))
wilcox.glm.t1_4.p <- data.frame(matrix(NA, nrow = 1000, ncol = 4))
wilcox.glm.t2_4.p <- data.frame(matrix(NA, nrow = 1000, ncol = 4))

## Run nested loops
for(l in c(12,15,17,14,16)) { #set number of iterations based on pairwise comparisons
	glm.temp <- data.frame(matrix(NA, nrow = 1000, ncol = 4))
	for(j in 1:4) { #set number of iterations based on candidate lists
  		#Establish Rank Sum test sets for CMH p values
  		#candidate set
  		candidates <- snp.gff.fulloverlap.glm[cand.filt.table[,j]==1 & snp.gff.fulloverlap.glm$genelength>0,l]
  		#background non-candidate set
  		bg <- snp.gff.fulloverlap.glm[cand.filt.table[,j]!=1,l]  
  		for(k in 1:1000) { #set number of iterations for test
    		bg.samp.list <- c()
    		#This loop finds a random set of background genes of roughly the same genomic length
    		#as the list of candidate genes
    		for(i in unique(snp.gff.fulloverlap.glm[cand.filt.table[,j]==1 & +
                                            snp.gff.fulloverlap.glm$genelength>0,]$genelength)) {
      			bg.samp.temp <- sample(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$V1 == head(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$genelength == i,],1)$V1 & +
                                                              snp.gff.fulloverlap.glm$genelength < i*1.25 & +
                                                              snp.gff.fulloverlap.glm$genelength > i*0.75 & +
                                                              cand.filt.table[,j]!=1,]$V9),1)
      			bg.samp.list <- append(bg.samp.list,bg.samp.temp)
    			}
    			bg.samp <- subset(snp.gff.fulloverlap.glm, V9 %in% bg.samp.list)[,l]
    			candidates.samp <- candidates
    			#this will subsample the SNPs contained in one list to the number of SNPs in the shorter list
    			if(length(candidates.samp) > length(bg.samp)){
      			candidates.samp <- sample(candidates.samp, length(bg.samp))
    			} else {
      			bg.samp <- sample(bg.samp, length(candidates.samp))
    		}
    		#finally, run the tests and append results
    		glm.temp[k,j] <- wilcox.test(candidates.samp, bg.samp, alternative = "less")$p.value
  		}
	}
  assign(paste("wilcox.glm", names(snp.gff.fulloverlap.glm)[l], sep="."), glm.temp)
}

## For single candidate SNPs only (not SNPs within candidate genes)
# Run the Wilcox Rank Sum test between candidates and matched sets for n iterations
# Initialize results tables
wilcox.glm.t1_2.p.snp <- data.frame(matrix(NA, nrow = 1000, ncol = 3))
wilcox.glm.t2_3.p.snp <- data.frame(matrix(NA, nrow = 1000, ncol = 3))
wilcox.glm.t3_4.p.snp <- data.frame(matrix(NA, nrow = 1000, ncol = 3))
wilcox.glm.t1_4.p.snp <- data.frame(matrix(NA, nrow = 1000, ncol = 3))
wilcox.glm.t2_4.p.snp <- data.frame(matrix(NA, nrow = 1000, ncol = 3))

## Run nested loops
for(l in c(12,15,17,14,16)) { #set number of iterations based on pairwise comparisons
  glm.temp <- data.frame(matrix(NA, nrow = 1000, ncol = 3))
	for(j in 1:3) { #set number of iterations
  		#Establish Rank Sum test sets for CMH p values
  		#candidate set
  		candidates <- snp.gff.fulloverlap.glm[candSNP.filt.table[,j]==1 & snp.gff.fulloverlap.glm$genelength>0,l]
  		#background non-candidate set
  		bg <- snp.gff.fulloverlap.glm[candSNP.filt.table[,j]!=1,l]
  
  		for(k in 1:1000) { #set number of iterations
    		bg.samp.list <- c()
    		#This loop finds a random set of background genes of roughly the same genomic length
    		#as the list of candidate genes
    		for(i in unique(snp.gff.fulloverlap.glm[candSNP.filt.table[,j]==1 & +
                                            snp.gff.fulloverlap.glm$genelength>0,]$genelength)) {
      			bg.samp.temp <- sample(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$V1 == head(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$genelength == i,],1)$V1 & +
                                                              snp.gff.fulloverlap.glm$genelength < i*1.25 & +
                                                              snp.gff.fulloverlap.glm$genelength > i*0.75 & +
                                                              candSNP.filt.table[,j]!=1,]$V9),1)
      			bg.samp.list <- append(bg.samp.list,bg.samp.temp)
    			}
    			bg.samp <- subset(snp.gff.fulloverlap.glm, V9 %in% bg.samp.list)[,l]
    			candidates.samp <- candidates
    			#this will subsample the SNPs contained in one list to the number of SNPs in the shorter list
    			if(length(candidates.samp) > length(bg.samp)){
      			candidates.samp <- sample(candidates.samp, length(bg.samp))
    			} else {
      			bg.samp <- sample(bg.samp, length(candidates.samp))
    		}
    		#finally, run the tests and append results
    		glm.temp[k,j] <- wilcox.test(candidates.samp, bg.samp, alternative = "less")$p.value
  		}
	}
  	assign(paste("wilcox.glm", names(snp.gff.fulloverlap.glm)[l], "snp", sep="."), glm.temp)
}


## save results  
write.table(wilcox.glm.t1_2.p, file="./results/wilcox.glm.t1_2.p.txt", sep="\t", quote = FALSE, row.names = F)
write.table(wilcox.glm.t2_3.p, file="./results/wilcox.glm.t2_3.p.txt", sep="\t", quote = FALSE, row.names = F)
write.table(wilcox.glm.t3_4.p, file="./results/wilcox.glm.t3_4.p.txt", sep="\t", quote = FALSE, row.names = F)
write.table(wilcox.glm.t1_4.p, file="./results/wilcox.glm.t1_4.p.txt", sep="\t", quote = FALSE, row.names = F)
write.table(wilcox.glm.t2_4.p, file="./results/wilcox.glm.t2_4.p.txt", sep="\t", quote = FALSE, row.names = F)

write.table(wilcox.glm.t1_2.p.snp, file="./results/wilcox.glm.t1_2.p.snp.txt", sep="\t", quote = FALSE, row.names = F)
write.table(wilcox.glm.t2_3.p.snp, file="./results/wilcox.glm.t2_3.p.snp.txt", sep="\t", quote = FALSE, row.names = F)
write.table(wilcox.glm.t3_4.p.snp, file="./results/wilcox.glm.t3_4.p.snp.txt", sep="\t", quote = FALSE, row.names = F)
write.table(wilcox.glm.t1_4.p.snp, file="./results/wilcox.glm.t1_4.p.snp.txt", sep="\t", quote = FALSE, row.names = F)
write.table(wilcox.glm.t2_4.p.snp, file="./results/wilcox.glm.t2_4.p.snp.txt", sep="\t", quote = FALSE, row.names = F)
