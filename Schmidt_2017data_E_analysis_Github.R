########################################################################
###  Selection analysis for Rudman 2017 Drosphila data: E treatment  ###
########################################################################

# Script developed starting 12/12/2023 by Josh Faber-Hammond

#set up working environment
setwd("~/Documents/josh_r_analysis/fst/Github_E.2017.analysis")
library(ggfortify)
library(ggpubr)
library(matrixStats)
library(stats)
library(dplyr)
library(plyr)

#Load allele frequency datasets and metadata
load("af.orch17.Rdata")
afmat -> afmat.17 #allele frequency table
rm(afmat)
rd -> rd.17 #read depth table
rm(rd)
samps -> samps.17 #sample metadata
rm(samps)
sites -> sites.17 #loci
rm(sites)
#Turn sites object into a bed file that can be used to look for candidate SNP/gene overlap
sites.17_bed <- data.frame(cbind('chrom'=sites.17$chrom, 'start'=sites.17$pos,
                                 'stop'=sites.17$pos + 1,'ref'=sites.17$ref, 'alt'=sites.17$alt))
sites.17_bed <- transform(sites.17_bed, start = as.integer(start))
sites.17_bed <- transform(sites.17_bed, stop = as.integer(stop))
sapply(sites.17_bed, class)
#write.table(sites.17_bed, file="sites.17.bed", sep="\t", quote = FALSE, row.names = F)
#write.table(sites.17, file="sites.17.txt", sep="\t", quote = FALSE, row.names = F)

## Load files associated with functional candidate analysis

# First lets load SNPs in ion-transport candidate genes. They were found using 
# the bash script "ion_candidate_retrieval.sh".
# These gene intervals include 500bp both upstream and downstream to capture 
# SNPs in potentially regulatory regions. 
ioncandidate_gene_overlap <- read.table("./functional.analysis/ion_candidate_dm5_fixed_snp_overlap_nr.txt", header=TRUE)
nrow(sites.17_bed) #expected number of rows
nrow(ioncandidate_gene_overlap)
all.equal(ioncandidate_gene_overlap$start,sites.17$pos) #TRUE
#if FALSE look into this and specifically what SNPs, if any, are duplicated and why
nrow(ioncandidate_gene_overlap)
nrow(ioncandidate_gene_overlap[duplicated(ioncandidate_gene_overlap[,1:2])==FALSE,])
# originally answer was false but I fixed by identifying duplicated candidate gene entries
# that were listed under different names, plus some TEH gene intervals overlapped after 
# adding the 500bp margins to capture upstream and downstream elements. Check column
# "duplicated", which indicates if a SNP had multiple candidate gene overlaps. 

## Seasonal candidates from Rudman 2022 Science paper
## Pull in seasonal-related genic SNPs from bedtools results ##
rudcandidateSNPs_bed <- read.table("./functional.analysis/rudman2022_candidateSNPs_snp_overlap.txt", header=FALSE)
rudcandidategenes500_bed <- read.table("./functional.analysis/rudman2022_candidategenes+-500_snp_overlap.txt", header=FALSE)

#Check the rows match
nrow(rudcandidateSNPs_bed) == nrow(sites.17_bed) #TRUE
all.equal(rudcandidateSNPs_bed[,2],sites.17$pos) #TRUE
nrow(rudcandidategenes500_bed) == nrow(sites.17_bed) #FALSE
all.equal(rudcandidategenes500_bed[,2],sites.17$pos) #FALSE

#if FALSE look into this and specifically what SNPs, if any, are duplicated and why
nrow(sites.17_bed) #expected number of rows
nrow(rudcandidategenes500_bed)
nrow(rudcandidategenes500_bed[duplicated(rudcandidategenes500_bed[,1:2])==FALSE,])

## Removing duplicates and convert to data frame
rudcandidateSNPs_overlap <- data.frame(rudcandidateSNPs_bed)
rudcandidategenes_overlap <- data.frame(rudcandidategenes500_bed[duplicated(rudcandidategenes500_bed[,1:2])==FALSE,])

## Fecundity candidates from https://pubmed.ncbi.nlm.nih.gov/25000897/
## Pull in fecundity-related genic SNPs from bedtools results ##
fecundcandidateSNPs_bed <- read.table("./functional.analysis/DGRP_fecundity_candidateSNPs_snp_overlap.txt", header=FALSE)
fecundcandidategenes500_bed <- read.table("./functional.analysis/DGRP_fecundity_candidategenes+-500_snp_overlap.txt", header=FALSE)
#Check the rows match
nrow(fecundcandidateSNPs_bed) == nrow(sites.17_bed) #TRUE
all.equal(fecundcandidateSNPs_bed[,2],sites.17$pos) #TRUE
nrow(fecundcandidategenes500_bed) == nrow(sites.17_bed) #FALSE
all.equal(fecundcandidategenes500_bed[,2],sites.17$pos) #FALSE

#if FALSE look into this and specifically what SNPs, if any, are duplicated and why
nrow(sites.17_bed) #expected number of rows
nrow(fecundcandidategenes500_bed)
nrow(fecundcandidategenes500_bed[duplicated(fecundcandidategenes500_bed[,1:2])==FALSE,])

## Removing duplicates and convert to data frame
fecundcandidateSNPs_overlap <- data.frame(fecundcandidateSNPs_bed)
fecundcandidategenes_overlap <- data.frame(fecundcandidategenes500_bed[duplicated(fecundcandidategenes500_bed[,1:2])==FALSE,])

## Starvation candidates from https://pmc.ncbi.nlm.nih.gov/articles/PMC3683990/#SD4
## Pull in starvation-related genic SNPs from bedtools results ##
starvcandidateSNPs_bed <- read.table("./functional.analysis/starvation_candidateSNPs_snp_overlap.txt", header=FALSE)
starvcandidategenes500_bed <- read.table("./functional.analysis/starvation_candidategenes+-500_snp_overlap.txt", header=FALSE)
#Check the rows match
nrow(starvcandidateSNPs_bed) == nrow(sites.17) #FALSE
all.equal(starvcandidateSNPs_bed[,2],sites.17$pos) #FALSE
nrow(starvcandidategenes500_bed) == nrow(sites.17) #TRUE
all.equal(starvcandidategenes500_bed[,2],sites.17$pos) #TRUE

#if FALSE look into this and specifically what SNPs, if any, are duplicated and why
nrow(sites.17_bed) #expected number of rows
nrow(starvcandidateSNPs_bed)
nrow(starvcandidateSNPs_bed[duplicated(starvcandidateSNPs_bed[,1:2])==FALSE,])

## Removing duplicates and convert to data frame
starvcandidateSNPs_overlap <- data.frame(starvcandidateSNPs_bed[duplicated(starvcandidateSNPs_bed[,1:2])==FALSE,])
starvcandidategenes_overlap <- data.frame(starvcandidategenes500_bed[duplicated(starvcandidategenes500_bed[,1:2])==FALSE,])

## Seasonal candidate SNPs from Bitter et al: https://pubmed.ncbi.nlm.nih.gov/39143223/
bitter.start2end <- read.table("./functional.analysis/bitter.start2end.candidateSNPs.bed", header=FALSE)
bitter.bestmatch <- read.table("./functional.analysis/bitter.bestmatch.candidateSNPs.bed", header=FALSE)

#Check the rows match
nrow(bitter.start2end) == nrow(sites.17) #TRUE
all.equal(bitter.start2end[,2],sites.17$pos) #TRUE
nrow(bitter.bestmatch) == nrow(sites.17) #FALSE
all.equal(bitter.bestmatch[,2],sites.17$pos) #"Numeric: lengths (1893190, 1885069) differ"

#if FALSE look into this and specifically what SNPs, if any, are duplicated and why
nrow(sites.17_bed) #expected number of rows
nrow(bitter.bestmatch)
nrow(bitter.bestmatch[duplicated(bitter.bestmatch[,1:2])==FALSE,])

## Removing duplicates and convert to data frame
bitter.start2end_overlap <- data.frame(bitter.start2end)
bitter.bestmatch_overlap <- data.frame(bitter.bestmatch[duplicated(bitter.bestmatch[,1:2])==FALSE,])

## Add column that indicates significant SNPs in 2017 data
bitter.start2end_overlap$V12 <- "1" 
bitter.start2end_overlap[bitter.start2end_overlap$V8=="-1",]$V12 <- "0"
bitter.bestmatch_overlap$V12 <- "1" 
bitter.bestmatch_overlap[bitter.bestmatch_overlap$V8=="-1",]$V12 <- "0"

## Now load all SNP-gene overlaps for overall functional analysis
## Gene coordinates are based on genome build dm5.39, with +/-500bp added to
## account for potential upstream or downstream cis-regulatory elements.
snp.gff.fulloverlap <- read.delim("./functional.analysis/dm5.39_snp_overlap.500bpbuffer.txt", sep="\t", comment.char = "", header=FALSE)
#check for row duplication
nrow(snp.gff.fulloverlap) == nrow(sites.17_bed)
all.equal(snp.gff.fulloverlap$V2,sites.17$pos)
#if FALSE look into this and specifically what SNPs, if any, are duplicated and why
dim(sites.17)
dim(snp.gff.fulloverlap)
dim(snp.gff.fulloverlap[duplicated(snp.gff.fulloverlap[,1:2])==FALSE,])
#hopefully this will make it TRUE
snp.gff.fulloverlap <- snp.gff.fulloverlap[duplicated(snp.gff.fulloverlap[,1:2])==FALSE,]
nrow(snp.gff.fulloverlap) == nrow(sites.17_bed) #success

#### Filter data for E samples only ####
afmat.17.E <- afmat.17[,samps.17$treatment == "E"]
rd.17.E <- rd.17[,samps.17$treatment == "E"]
dim(afmat.17.E)
dim(rd.17.E)
dim(samps.17[samps.17$treatment == "E",])
samps.17[samps.17$treatment == "E",]$sampID == cbind.data.frame(colnames(afmat.17.E))
samps.17.E <- transform(samps.17[samps.17$treatment == "E",],rep=as.numeric(factor(cage)))

### PCA to see differences in round before and after filtering out T1 ###
#Format objects so they properly work with prcomp
samps.17.E$cage = as.character(samps.17.E$cage)
samps.17.E$tpt = as.character(samps.17.E$tpt)
#Include mean AF of replicate founder samples for reference
# plot candidate by subbing in line number from above results
afmat.17.E.F <- afmat.17[,samps.17$treatment == "E" | (samps.17$treatment == "Founder")]
samps.17.E.F <- rbind(samps.17[1:2,],samps.17.E[,-7])
samps.17.E.F[1,4] <- "founder"
samps.17.E.F[2,4] <- "founder"
samps.17.E.F[1,5] <- "founder"
samps.17.E.F[2,5] <- "founder"
samps.17.E.F$tpt = as.character(samps.17.E.F$tpt)
samps.17.E.F.filt <- samps.17.E.F[samps.17.E.F$tpt!=1,]

# Make plots
pca_res.founder <- prcomp(t(cbind(afmat.17[,1:2],afmat.17.E)), scale. = TRUE)
#Remove tpt1 due to noise from combined sequencing runs
pca_res.founder.filt <- prcomp(t(cbind(afmat.17[,1:2],afmat.17.E)[,which(samps.17.E.F$tpt != 1)]), scale. = TRUE)

autoplot(pca_res.founder, data = samps.17.E.F, colour = 'sequencingRd')
autoplot(pca_res.founder, data = samps.17.E.F, colour = 'cage', shape = 'sequencingRd')
autoplot(pca_res.founder, data = samps.17.E.F, colour = 'cage', shape = 'tpt', size = 3) + theme_bw()
#no apparent separation by time-point or sequencing round, strong clustering by cage
#with some more divergent from found pop than others

#Are different patterns seen when dropping tpt 1 (for variation in sequencing round)?
autoplot(pca_res.founder.filt, data = samps.17.E.F.filt, colour = 'sequencingRd')
autoplot(pca_res.founder.filt, data = samps.17.E.F.filt, colour = 'cage', shape = 'sequencingRd')
autoplot(pca_res.founder.filt, data = samps.17.E.F.filt, colour = 'cage', shape = 'tpt', size = 3) + theme_bw()
#Generally the same patterns with or without tpt 1
autoplot(pca_res.founder.filt, data = samps.17.E.F.filt,
         colour = 'tpt', frame = TRUE, frame.type = 'norm', size = 3) + theme_bw()
#When looking at distribution by time-point, centroids are all similar suggesting
#they are not explained well by PC1 or PC2, however the overall distributions of
#samples grows wider over time, possibly due to the compounding impact of drift/selection

### Examine GLM results run on kamiak using script "glm2017.E.CAS.r" ###
#This is a 6-column table containing raw p-values for each time-point contrast
# 1 vs. 2, 1 vs. 3, 1 vs. 4, 2 vs. 3, 2 vs. 4, 3 vs. 4
#Each row is a separate SNP corresponding to rows in afmat.17.E
E.glm.contrastp <- read.table("./glm.r.analysis/E.2017.contrastTPT.txt", header=TRUE)
#rename column headers
names(E.glm.contrastp)[1] <- "t1_2.p"
names(E.glm.contrastp)[2] <- "t1_3.p"
names(E.glm.contrastp)[3] <- "t1_4.p"
names(E.glm.contrastp)[4] <- "t2_3.p"
names(E.glm.contrastp)[5] <- "t2_4.p"
names(E.glm.contrastp)[6] <- "t3_4.p"
#p-value correction
E.glm.contrastp$t1_2.fdr <- as.numeric(p.adjust(E.glm.contrastp[,1], method = "fdr"))
E.glm.contrastp$t1_3.fdr <- as.numeric(p.adjust(E.glm.contrastp[,2], method = "fdr"))
E.glm.contrastp$t1_4.fdr <- as.numeric(p.adjust(E.glm.contrastp[,3], method = "fdr"))
E.glm.contrastp$t2_3.fdr <- as.numeric(p.adjust(E.glm.contrastp[,4], method = "fdr"))
E.glm.contrastp$t2_4.fdr <- as.numeric(p.adjust(E.glm.contrastp[,5], method = "fdr"))
E.glm.contrastp$t3_4.fdr <- as.numeric(p.adjust(E.glm.contrastp[,6], method = "fdr"))
#how many pairwise comparisons significant per row at several thresholds
E.glm.contrastp$SelCount <- rowSums(E.glm.contrastp[,1:6] < 0.01)
E.glm.contrastp$SelCount.noT1 <- rowSums(E.glm.contrastp[,4:6] < 0.01) #skip t1 contrasts
E.glm.contrastp$fdrCount <- rowSums(E.glm.contrastp[,10:12] < 0.05)
E.glm.contrastp$fdrCount.noT1 <- rowSums(E.glm.contrastp[,10:12] < 0.05) #skip t1 contrasts

## Let's look at stats for each GLM pairwise correlations
#FDR<0.05
nrow(E.glm.contrastp[E.glm.contrastp$t1_2.fdr<0.05,])
nrow(E.glm.contrastp[E.glm.contrastp$t2_3.fdr<0.05,])
nrow(E.glm.contrastp[E.glm.contrastp$t3_4.fdr<0.05,])
nrow(E.glm.contrastp[E.glm.contrastp$t1_4.fdr<0.05,])
nrow(E.glm.contrastp[E.glm.contrastp$t2_4.fdr<0.05,])
#less strict FDR
nrow(E.glm.contrastp[E.glm.contrastp$t1_2.fdr<0.25,])
nrow(E.glm.contrastp[E.glm.contrastp$t2_3.fdr<0.25,])
nrow(E.glm.contrastp[E.glm.contrastp$t3_4.fdr<0.25,])
nrow(E.glm.contrastp[E.glm.contrastp$t1_4.fdr<0.25,])
nrow(E.glm.contrastp[E.glm.contrastp$t2_4.fdr<0.25,])
#unadjusted pvals
#p<0.01
nrow(E.glm.contrastp[E.glm.contrastp$t1_2.p<0.01,])
nrow(E.glm.contrastp[E.glm.contrastp$t2_3.p<0.01,])
nrow(E.glm.contrastp[E.glm.contrastp$t3_4.p<0.01,])
nrow(E.glm.contrastp[E.glm.contrastp$t1_4.p<0.01,])
nrow(E.glm.contrastp[E.glm.contrastp$t2_4.p<0.01,])
#p<0.05
nrow(E.glm.contrastp[E.glm.contrastp$t1_2.p<0.05,])
nrow(E.glm.contrastp[E.glm.contrastp$t2_3.p<0.05,])
nrow(E.glm.contrastp[E.glm.contrastp$t3_4.p<0.05,])
nrow(E.glm.contrastp[E.glm.contrastp$t1_4.p<0.05,])
nrow(E.glm.contrastp[E.glm.contrastp$t2_4.p<0.05,])

### Make Manhattan Plot while highlighting significant genes
#Calculate LOD scores
E.glm.contrastp$logp1 <- -log10(E.glm.contrastp[,1])
E.glm.contrastp$logp2 <- -log10(E.glm.contrastp[,2])
E.glm.contrastp$logp3 <- -log10(E.glm.contrastp[,3])
E.glm.contrastp$logp4 <- -log10(E.glm.contrastp[,4])
E.glm.contrastp$logp5 <- -log10(E.glm.contrastp[,5])
E.glm.contrastp$logp6 <- -log10(E.glm.contrastp[,6])

#Create new table from all SNPs as input for 'manhattan' function
manh <-data.frame(cbind(chrom=ioncandidate_gene_overlap$chrom,
                        bp=ioncandidate_gene_overlap$start,
                        snp=paste("chr",ioncandidate_gene_overlap$chrom,
                        ioncandidate_gene_overlap$start, sep=":"),
                        E.glm.contrastp))

#Format columns so they plot correctly
manh$chr <- manh$chrom
manh$chr[manh$chr == '2L'] <- '1'
manh$chr[manh$chr == '2R'] <- '2'
manh$chr[manh$chr == '3L'] <- '3'
manh$chr[manh$chr == '3R'] <- '4'
manh$chr[manh$chr == 'X'] <- '5'
manh$chr <- as.numeric(manh$chr)
manh$bp <- as.numeric(manh$bp)

#Joining candidate lists together into single column to plot in different colors
manh.gene.cand <- cbind(manh,
                        ion=ioncandidate_gene_overlap$overlap,
                        seasonal=rudcandidategenes_overlap$V11,
                        fecundity=fecundcandidategenes_overlap$V12,
                        starvation=starvcandidategenes_overlap$V12)
templabel <- ifelse(manh.gene.cand$seasonal==1,"seasonal","noncand");
templabel2 <- ifelse(manh.gene.cand$fecundity==1,"fecundity",templabel);
templabel <- ifelse(manh.gene.cand$starvation==1,"starvation",templabel2);
templabel2 <- ifelse(manh.gene.cand$ion==1,"ion-transport",templabel);
manh.gene.cand$label <- templabel2
unique(manh.gene.cand$label) #make sure all functional categories listed

## Making plots ##
# Note: Thanks to some analysis not shown here, we identified that large portion 
# of outliers or significant SNPs in time-point 1 contrasts were actually driven
# by differences between sequencing round within time-point 1. This was done by
# running simple t-tests on allele frequencies between sequencing rounds 
# within time-point 1, then looking at overlap of significant SNPs with sets of 
# significant SNPs from time-point contrasts including tpt 1, and found that
# this overlap is higher than expected by chance. From here, we'll drop 
# contrasts that include tpt 1.

#T2 vs T3
#this block of code will be use to highlight any highly significant SNPs
manh.gene.cand.sig <- manh.gene.cand[manh.gene.cand$t2_3.fdr<0.05,]
manh.gene.cand.sig$label <- "FDR<0.05"
manh.gene.cand2 <- rbind(manh.gene.cand.sig,manh.gene.cand) #prioritize highlighting sig. candidate genes
#Plot
glmmanh1 <- ggplot(manh.gene.cand2, aes(bp, logp4)) + geom_point(colour = "grey") +
  geom_hline(yintercept=min(manh.gene.cand[manh.gene.cand$t2_3.fdr<0.25,]$logp4),color="black",linetype="dashed",linewidth=.25) +
  geom_hline(yintercept=min(manh.gene.cand[manh.gene.cand$t2_3.fdr<0.1,]$logp4),color="black",linetype="dashed",linewidth=.25) +
  geom_hline(yintercept=min(manh.gene.cand[manh.gene.cand$t2_3.fdr<0.05,]$logp4),color="black",linetype="dashed",linewidth=.25) +
  geom_jitter(data=manh.gene.cand2[manh.gene.cand2$label!="noncand",],aes(color = label),position = "jitter", size = 2.5) +
  scale_color_manual(values=c("black","#56B4E9","#009E73","#D81B60","#FFC107")) +
  facet_grid(~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(n.breaks = 4,guide = guide_axis(angle = 45)) + 
  labs(col="candidate\ngene\ncategory") +
  scale_y_continuous(limits=c(2, NA)) +
  xlab("chromosome position") +
  ylab("-log10(p)") +
  ggtitle("Sept 22nd vs Oct 18th") +
  theme_classic()

#T3 vs T4
#this block of code will be use to highlight any highly significant SNPs
manh.gene.cand.sig <- manh.gene.cand[manh.gene.cand$t3_4.fdr<0.05,]
manh.gene.cand.sig$label <- "FDR<0.05"
manh.gene.cand2 <- rbind(manh.gene.cand.sig,manh.gene.cand) #prioritize highlighting sig. candidate genes
#Plot
glmmanh2 <- ggplot(manh.gene.cand2, aes(bp, logp6)) + geom_point(colour = "grey") +
  geom_hline(yintercept=min(manh.gene.cand[manh.gene.cand$t3_4.fdr<0.25,]$logp6),color="black",linetype="dashed",linewidth=.25) +
  geom_hline(yintercept=min(manh.gene.cand[manh.gene.cand$t3_4.fdr<0.1,]$logp6),color="black",linetype="dashed",linewidth=.25) +
  geom_hline(yintercept=min(manh.gene.cand[manh.gene.cand$t3_4.fdr<0.05,]$logp6),color="black",linetype="dashed",linewidth=.25) +
  geom_jitter(data=manh.gene.cand2[manh.gene.cand2$label!="noncand",],aes(color = label),position = "jitter", size = 2.5) +
  scale_color_manual(values=c("black","#56B4E9","#009E73","#D81B60","#FFC107")) +
  facet_grid(~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(n.breaks = 4,guide = guide_axis(angle = 45)) + 
  labs(col="candidate\ngene\ncategory") +
  scale_y_continuous(limits=c(2, NA)) +
  xlab("chromosome position") +
  ylab("-log10(p)") +
  ggtitle("Oct 18th vs Nov 12th") +
  theme_classic()


#T2 vs T4
#this block of code will be use to highlight any highly significant SNPs
manh.gene.cand.sig <- manh.gene.cand[manh.gene.cand$t2_4.fdr<0.05,]
manh.gene.cand.sig$label <- "FDR<0.05"
manh.gene.cand2 <- rbind(manh.gene.cand.sig,manh.gene.cand) #prioritize highlighting sig. candidate genes
#Plot
glmmanh3 <- ggplot(manh.gene.cand2, aes(bp, logp5)) + geom_point(colour = "grey") +
  geom_hline(yintercept=min(manh.gene.cand[manh.gene.cand$t2_4.fdr<0.25,]$logp5),color="black",linetype="dashed",linewidth=.25) +
  geom_hline(yintercept=min(manh.gene.cand[manh.gene.cand$t2_4.fdr<0.1,]$logp5),color="black",linetype="dashed",linewidth=.25) +
  geom_hline(yintercept=min(manh.gene.cand[manh.gene.cand$t2_4.fdr<0.05,]$logp5),color="black",linetype="dashed",linewidth=.25) +
  geom_jitter(data=manh.gene.cand2[manh.gene.cand2$label!="noncand",],aes(color = label),position = "jitter", size = 2.5) +
  scale_color_manual(values=c("black","#56B4E9","#009E73","#D81B60","#FFC107")) +
  facet_grid(~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(n.breaks = 4,guide = guide_axis(angle = 45)) + 
  labs(col="candidate\ngene\ncategory") +
  scale_y_continuous(limits=c(2, NA)) +
  xlab("chromosome position") +
  ylab("-log10(p)") +
  ggtitle("Sept 22nd vs Nov 12th") +
  theme_classic()

#Put them together
png("/Users/rudmanlab/Documents/josh_r_analysis/fst/E.2017.glm.cand.genes.not1.png",
    width = 1000, height = 750)
glmmanhAll.gene <- ggarrange(glmmanh1, glmmanh2, glmmanh3,
                        ncol = 1, nrow = 3)
glmmanhAll.gene
dev.off()


## This time, let's try SNPs only ##
#Joining candidate lists together into single column to plot in different colors
manh.snp.cand <- cbind(manh,
                       seasonal=rudcandidateSNPs_overlap$V10,
                       fecundity=fecundcandidateSNPs_overlap$V11,
                       starvation=starvcandidateSNPs_overlap$V10,
                       start2end=bitter.start2end_overlap$V12,
                       bestmatch=bitter.bestmatch_overlap$V12)
templabel <- ifelse(manh.snp.cand$seasonal==1,"seasonal","noncand");
templabel2 <- ifelse(manh.snp.cand$fecundity==1,"fecundity",templabel);
templabel <- ifelse(manh.snp.cand$starvation==1,"starvation",templabel2);
templabel2 <- ifelse(manh.snp.cand$start2end==1,"start2end",templabel);
templabel <- ifelse(manh.snp.cand$bestmatch==1,"bestmatch",templabel2);
manh.snp.cand$label <- templabel
unique(manh.snp.cand$label) #make sure all functional categories listed

#T2 vs T3
#this block of code will be use to highlight any highly significant SNPs
manh.snp.cand.sig <- manh.snp.cand[manh.snp.cand$t2_3.fdr<0.05,]
manh.snp.cand.sig$label <- "FDR<0.05"
manh.snp.cand2 <- rbind(manh.snp.cand.sig,manh.snp.cand) #prioritize highlighting sig. candidate genes
#Plot
glmmanh4 <- ggplot(manh.snp.cand2, aes(bp, logp4)) + geom_point(colour = "grey") +
  geom_hline(yintercept=min(manh.snp.cand[manh.snp.cand$t2_3.fdr<0.25,]$logp4),color="black",linetype="dashed",linewidth=.25) +
  geom_hline(yintercept=min(manh.snp.cand[manh.snp.cand$t2_3.fdr<0.1,]$logp4),color="black",linetype="dashed",linewidth=.25) +
  geom_hline(yintercept=min(manh.snp.cand[manh.snp.cand$t2_3.fdr<0.05,]$logp4),color="black",linetype="dashed",linewidth=.25) +
  geom_jitter(data=manh.snp.cand2[manh.snp.cand2$label!="noncand",],aes(color = label),position = "jitter", size = 2.5) +
  scale_color_manual(values=c("black","#56B4E9","#D81B60","#FFC107","green2",'purple')) +
  facet_grid(~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(n.breaks = 4,guide = guide_axis(angle = 45)) + 
  labs(col="candidate\nSNP\ncategory") +
  xlab("chromosome position") +
  ylab("-log10(p)") +
  ggtitle("Sept 22nd vs Oct 18th") +
  theme_classic()

#T3 vs T4
#this block of code will be use to highlight any highly significant SNPs
manh.snp.cand.sig <- manh.snp.cand[manh.snp.cand$t3_4.fdr<0.05,]
manh.snp.cand.sig$label <- "FDR<0.05"
manh.snp.cand2 <- rbind(manh.snp.cand.sig,manh.snp.cand) #prioritize highlighting sig. candidate genes
#Plot
glmmanh5 <- ggplot(manh.snp.cand2, aes(bp, logp6)) + geom_point(colour = "grey") +
  geom_hline(yintercept=min(manh.snp.cand[manh.snp.cand$t3_4.fdr<0.25,]$logp6),color="black",linetype="dashed",linewidth=.25) +
  geom_hline(yintercept=min(manh.snp.cand[manh.snp.cand$t3_4.fdr<0.1,]$logp6),color="black",linetype="dashed",linewidth=.25) +
  geom_hline(yintercept=min(manh.snp.cand[manh.snp.cand$t3_4.fdr<0.05,]$logp6),color="black",linetype="dashed",linewidth=.25) +
  geom_jitter(data=manh.snp.cand2[manh.snp.cand2$label!="noncand",],aes(color = label),position = "jitter", size = 2.5) +
  scale_color_manual(values=c("black","#56B4E9","#D81B60","#FFC107","green2",'purple')) +
  facet_grid(~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(n.breaks = 4,guide = guide_axis(angle = 45)) + 
  labs(col="candidate\nSNP\ncategory") +
  xlab("chromosome position") +
  ylab("-log10(p)") +
  ggtitle("Oct 18th vs Nov 12th") +
  theme_classic()

#T2 vs T4
#this block of code will be use to highlight any highly significant SNPs
manh.snp.cand.sig <- manh.snp.cand[manh.snp.cand$t2_4.fdr<0.05,]
manh.snp.cand.sig$label <- "FDR<0.05"
manh.snp.cand2 <- rbind(manh.snp.cand.sig,manh.snp.cand) #prioritize highlighting sig. candidate genes
#Plot
glmmanh6 <- ggplot(manh.snp.cand2, aes(bp, logp5)) + geom_point(colour = "grey") +
  geom_hline(yintercept=min(manh.snp.cand[manh.snp.cand$t2_4.fdr<0.25,]$logp5),color="black",linetype="dashed",linewidth=.25) +
  geom_hline(yintercept=min(manh.snp.cand[manh.snp.cand$t2_4.fdr<0.1,]$logp5),color="black",linetype="dashed",linewidth=.25) +
  geom_hline(yintercept=min(manh.snp.cand[manh.snp.cand$t2_4.fdr<0.05,]$logp5),color="black",linetype="dashed",linewidth=.25) +
  geom_jitter(data=manh.snp.cand2[manh.snp.cand2$label!="noncand",],aes(color = label),position = "jitter", size = 2.5) +
  scale_color_manual(values=c("black","#56B4E9","#D81B60","#FFC107","green2",'purple')) +
  facet_grid(~ chrom, scales = "free_x", space = "free_x") +
  scale_x_continuous(n.breaks = 4,guide = guide_axis(angle = 45)) + 
  labs(col="candidate\nSNP\ncategory") +
  xlab("chromosome position") +
  ylab("-log10(p)") +
  ggtitle("Sept 22nd vs Nov 12th") +
  theme_classic()

#Put them together
png("/Users/rudmanlab/Documents/josh_r_analysis/fst/E.2017.glm.cand.SNPs.not1.wBitter.png",
    width = 1000, height = 750)
glmmanhAll.snp <- ggarrange(glmmanh4, glmmanh5, glmmanh6,
                        ncol = 1, nrow = 3)
glmmanhAll.snp
dev.off()


######################################
## GLM Hypergeometric Overlap Tests ##
######################################
## How many significant SNPs in functional candidate genes?
# Make a candidate filtering list
gene.cand <- cbind(ion=ioncandidate_gene_overlap$overlap,
                   seasonal=rudcandidategenes_overlap$V11,
                   fecundity=fecundcandidategenes_overlap$V12,
                   starvation=starvcandidategenes_overlap$V12)
#write.table(gene.cand, file="cand.filt.table.txt", sep="\t", quote = FALSE, row.names = F)
# Nested loops for making separate count lists per functional category
for(h in c(1:4)) { #loop through each functional candidate list
  #initialize count tables
  fdr05 <- c()
  fdr25 <- c()
  p01 <- c()
  p05 <- c()
  #loop through columns representing fdr for tpt 2 vs. 3, 3 vs. 4, 2 vs. 4
  for(i in c(10,12,11)) { 
    fdr05 <- append(fdr05,nrow(E.glm.contrastp[E.glm.contrastp[,i]<0.05 &
                           gene.cand[,h]==1,]))
    fdr25 <- append(fdr25,nrow(E.glm.contrastp[E.glm.contrastp[,i]<0.25 &
                           gene.cand[,h]==1,]))
  }
  #loop through columns representing p-values for tpt 2 vs. 3, 3 vs. 4, 2 vs. 4
  for(j in c(4,6,5)) { 
    p01 <- append(p01,nrow(E.glm.contrastp[E.glm.contrastp[,j]<0.01 &
                           gene.cand[,h]==1,]))
    p05 <- append(p05,nrow(E.glm.contrastp[E.glm.contrastp[,j]<0.05 &
                           gene.cand[,h]==1,]))    
  }
  cand.sig.table <- rbind(fdr05,fdr25,p01,p05)
  colnames(cand.sig.table) <- c("t2.vs.t3","t3.vs.t4","t2.vs.t4")
  assign(paste(colnames(gene.cand)[h],"cand.gene.sig.table", sep="."), cand.sig.table)
}
#Check resulting tables
ion.cand.gene.sig.table
seasonal.cand.gene.sig.table
fecundity.cand.gene.sig.table
starvation.cand.gene.sig.table

## How many significant SNPs are specifically functional candidate SNPs?
# Make a candidate filtering list
snp.cand <- cbind(seasonal=rudcandidateSNPs_overlap$V10,
                  fecundity=fecundcandidateSNPs_overlap$V11,
                  starvation=starvcandidateSNPs_overlap$V10,
                  start2end=bitter.start2end_overlap$V12,
                  bestmatch=bitter.bestmatch_overlap$V12)
#write.table(snp.cand, file="candSNP.filt.table.txt", sep="\t", quote = FALSE, row.names = F)

# Nested loops for making separate count lists per functional category
for(h in c(1:5)) { #loop through each functional candidate list
  #initialize count tables
  fdr05 <- c()
  fdr25 <- c()
  p01 <- c()
  p05 <- c()
  #loop through columns representing fdr for tpt 2 vs. 3, 3 vs. 4, 2 vs. 4
  for(i in c(10,12,11)) { 
    fdr05 <- append(fdr05,nrow(E.glm.contrastp[E.glm.contrastp[,i]<0.05 &
                                                 snp.cand[,h]==1,]))
    fdr25 <- append(fdr25,nrow(E.glm.contrastp[E.glm.contrastp[,i]<0.25 &
                                                 snp.cand[,h]==1,]))
  }
  #loop through columns representing p-values for tpt 2 vs. 3, 3 vs. 4, 2 vs. 4
  for(j in c(4,6,5)) { 
    p01 <- append(p01,nrow(E.glm.contrastp[E.glm.contrastp[,j]<0.01 &
                                             snp.cand[,h]==1,]))
    p05 <- append(p05,nrow(E.glm.contrastp[E.glm.contrastp[,j]<0.05 &
                                             snp.cand[,h]==1,]))    
  }
  cand.sig.table <- rbind(fdr05,fdr25,p01,p05)
  colnames(cand.sig.table) <- c("t2.vs.t3","t3.vs.t4","t2.vs.t4")
  assign(paste(colnames(snp.cand)[h],"cand.snp.sig.table", sep="."), cand.sig.table)
}
#Check resulting tables
seasonal.cand.snp.sig.table
fecundity.cand.snp.sig.table
starvation.cand.snp.sig.table
start2end.cand.snp.sig.table
bestmatch.cand.snp.sig.table


## Hypergeometric overlap time ##
# Run hypergeometric overlap test per functional category, GLM contrast, and 
# significance threshold via nested loops.
## Start with SNPs within candidate genes (+/-500bp) 
for(h in c(1:4)) { #loop through each functional candidate list
  #initialize count tables
  phyper.fdr05 <- c()
  phyper.fdr25 <- c()
  phyper.p01 <- c()
  phyper.p05 <- c()
  #loop through columns representing fdr for tpt 2 vs. 3, 3 vs. 4, 2 vs. 4
  for(i in c(10,12,11)) { 
    #fdr<0.05
    q <- nrow(E.glm.contrastp[gene.cand[,h]==1 & E.glm.contrastp[,i]<0.05,])
    m <- nrow(E.glm.contrastp[gene.cand[,h]==1,])
    n <- nrow(E.glm.contrastp[gene.cand[,h]==0,])
    k <- nrow(E.glm.contrastp[E.glm.contrastp[,i]<0.05,])
    ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
    phyper.fdr05 <- append(phyper.fdr05,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
    #fdr<0.25
    q <- nrow(E.glm.contrastp[gene.cand[,h]==1 & E.glm.contrastp[,i]<0.25,])
    m <- nrow(E.glm.contrastp[gene.cand[,h]==1,])
    n <- nrow(E.glm.contrastp[gene.cand[,h]==0,])
    k <- nrow(E.glm.contrastp[E.glm.contrastp[,i]<0.25,])
    ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
    phyper.fdr25 <- append(phyper.fdr25,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  }
  #loop through columns representing p-values for tpt 2 vs. 3, 3 vs. 4, 2 vs. 4
  for(j in c(4,6,5)) { 
    #fdr<0.05
    q <- nrow(E.glm.contrastp[gene.cand[,h]==1 & E.glm.contrastp[,j]<0.01,])
    m <- nrow(E.glm.contrastp[gene.cand[,h]==1,])
    n <- nrow(E.glm.contrastp[gene.cand[,h]==0,])
    k <- nrow(E.glm.contrastp[E.glm.contrastp[,j]<0.01,])
    ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
    phyper.p01 <- append(phyper.p01,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
    #fdr<0.25
    q <- nrow(E.glm.contrastp[gene.cand[,h]==1 & E.glm.contrastp[,j]<0.05,])
    m <- nrow(E.glm.contrastp[gene.cand[,h]==1,])
    n <- nrow(E.glm.contrastp[gene.cand[,h]==0,])
    k <- nrow(E.glm.contrastp[E.glm.contrastp[,j]<0.05,])
    ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
    phyper.p05 <- append(phyper.p05,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  }
  cand.hyper.table <- rbind(phyper.fdr05,phyper.fdr25,phyper.p01,phyper.p05)
  colnames(cand.hyper.table) <- c("t2.vs.t3","t3.vs.t4","t2.vs.t4")
  assign(paste(colnames(gene.cand)[h],"cand.gene.hyper.table", sep="."), cand.hyper.table)
}
#Check resulting tables
ion.cand.gene.hyper.table
seasonal.cand.gene.hyper.table
fecundity.cand.gene.hyper.table
starvation.cand.gene.hyper.table

## Next, look whether specific functional candidate SNPs are significant 
for(h in c(1:5)) { #loop through each functional candidate list
  #initialize count tables
  phyper.fdr05 <- c()
  phyper.fdr25 <- c()
  phyper.p01 <- c()
  phyper.p05 <- c()
  #loop through columns representing fdr for tpt 2 vs. 3, 3 vs. 4, 2 vs. 4
  for(i in c(10,12,11)) { 
    #fdr<0.05
    q <- nrow(E.glm.contrastp[snp.cand[,h]==1 & E.glm.contrastp[,i]<0.05,])
    m <- nrow(E.glm.contrastp[snp.cand[,h]==1,])
    n <- nrow(E.glm.contrastp[snp.cand[,h]==0,])
    k <- nrow(E.glm.contrastp[E.glm.contrastp[,i]<0.05,])
    ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
    phyper.fdr05 <- append(phyper.fdr05,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
    #fdr<0.25
    q <- nrow(E.glm.contrastp[snp.cand[,h]==1 & E.glm.contrastp[,i]<0.25,])
    m <- nrow(E.glm.contrastp[snp.cand[,h]==1,])
    n <- nrow(E.glm.contrastp[snp.cand[,h]==0,])
    k <- nrow(E.glm.contrastp[E.glm.contrastp[,i]<0.25,])
    ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
    phyper.fdr25 <- append(phyper.fdr25,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  }
  #loop through columns representing p-values for tpt 2 vs. 3, 3 vs. 4, 2 vs. 4
  for(j in c(4,6,5)) { 
    #fdr<0.05
    q <- nrow(E.glm.contrastp[snp.cand[,h]==1 & E.glm.contrastp[,j]<0.01,])
    m <- nrow(E.glm.contrastp[snp.cand[,h]==1,])
    n <- nrow(E.glm.contrastp[snp.cand[,h]==0,])
    k <- nrow(E.glm.contrastp[E.glm.contrastp[,j]<0.01,])
    ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
    phyper.p01 <- append(phyper.p01,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
    #fdr<0.25
    q <- nrow(E.glm.contrastp[snp.cand[,h]==1 & E.glm.contrastp[,j]<0.05,])
    m <- nrow(E.glm.contrastp[snp.cand[,h]==1,])
    n <- nrow(E.glm.contrastp[snp.cand[,h]==0,])
    k <- nrow(E.glm.contrastp[E.glm.contrastp[,j]<0.05,])
    ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
    phyper.p05 <- append(phyper.p05,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
  }
  cand.hyper.table <- rbind(phyper.fdr05,phyper.fdr25,phyper.p01,phyper.p05)
  colnames(cand.hyper.table) <- c("t2.vs.t3","t3.vs.t4","t2.vs.t4")
  assign(paste(colnames(snp.cand)[h],"cand.snp.hyper.table", sep="."), cand.hyper.table)
}
#Check resulting tables
seasonal.cand.snp.hyper.table
fecundity.cand.snp.hyper.table
starvation.cand.snp.hyper.table
start2end.cand.snp.hyper.table
bestmatch.cand.snp.hyper.table

###################################################################################
# Run Wilcoxon Rank Sum test between candidates and matched sets for n iterations #
###################################################################################
# The following input files were generated by running the R-script:
# "E_2017_wilcoxon_analysis.r" on our HPCC for speed. See that script for more
# details.

### Stats from ouput files from kamiak
# Columns in each of these files refer to SNPs associated with 3 functional categories:
# Seasonal, Fecundity, Starvation (plus the two Bitter et al seasonal datasets)
wilcox.glm.t2_3.p.snp <- read.delim("./wilcoxon.r.analysis/results/wilcox.glm.t2_3.p.snp.txt", 
                                    sep = "\t", header=TRUE)
wilcox.glm.t3_4.p.snp <- read.delim("./wilcoxon.r.analysis/results/wilcox.glm.t3_4.p.snp.txt", 
                                    sep = "\t", header=TRUE)
wilcox.glm.t2_4.p.snp <- read.delim("./wilcoxon.r.analysis/results/wilcox.glm.t2_4.p.snp.txt", 
                                    sep = "\t", header=TRUE)
# Columns in each of these files refer to SNPs within 4 functional categories of genes:
# Ion-transport, Seasonal, Fecundity, Starvation
wilcox.glm.t2_3.p.gene <- read.delim("./wilcoxon.r.analysis/results/wilcox.glm.t2_3.p.txt", 
                                     sep = "\t", header=TRUE)
wilcox.glm.t3_4.p.gene <- read.delim("./wilcoxon.r.analysis/results/wilcox.glm.t3_4.p.txt", 
                                     sep = "\t", header=TRUE)
wilcox.glm.t2_4.p.gene <- read.delim("./wilcoxon.r.analysis/results/wilcox.glm.t2_4.p.txt", 
                                     sep = "\t", header=TRUE)

# Median values
#SNPs-only
snp.wilcox.table <- rbind(t2_3=colMedians(as.matrix(wilcox.glm.t2_3.p.snp)),
                          t3_4=colMedians(as.matrix(wilcox.glm.t3_4.p.snp)),
                          t2_4=colMedians(as.matrix(wilcox.glm.t2_4.p.snp)))
colnames(snp.wilcox.table) <- c("Seasonal", "Fecundity", "Starvation")
snp.wilcox.table

#SNPs within genes
gene.wilcox.table <- rbind(t2_3=colMedians(as.matrix(wilcox.glm.t2_3.p.gene)),
                           t3_4=colMedians(as.matrix(wilcox.glm.t3_4.p.gene)),
                           t2_4=colMedians(as.matrix(wilcox.glm.t2_4.p.gene)))
colnames(gene.wilcox.table) <- c("Ion-transport", "Seasonal", "Fecundity", "Starvation")
gene.wilcox.table

# Percent significant
#SNPs-only
snp.wilcox.table.perc <- rbind(t2_3=c(nrow(wilcox.glm.t2_3.p.snp[wilcox.glm.t2_3.p.snp$X1<0.05,]),
                                      nrow(wilcox.glm.t2_3.p.snp[wilcox.glm.t2_3.p.snp$X2<0.05,]),
                                      nrow(wilcox.glm.t2_3.p.snp[wilcox.glm.t2_3.p.snp$X3<0.05,])),
                               t3_4=c(nrow(wilcox.glm.t3_4.p.snp[wilcox.glm.t3_4.p.snp$X1<0.05,]),
                                      nrow(wilcox.glm.t3_4.p.snp[wilcox.glm.t3_4.p.snp$X2<0.05,]),
                                      nrow(wilcox.glm.t3_4.p.snp[wilcox.glm.t3_4.p.snp$X3<0.05,])),
                               t2_4=c(nrow(wilcox.glm.t2_4.p.snp[wilcox.glm.t2_4.p.snp$X1<0.05,]),
                                      nrow(wilcox.glm.t2_4.p.snp[wilcox.glm.t2_4.p.snp$X2<0.05,]),
                                      nrow(wilcox.glm.t2_4.p.snp[wilcox.glm.t2_4.p.snp$X3<0.05,])))
colnames(snp.wilcox.table.perc) <- c("Seasonal", "Fecundity", "Starvation")
snp.wilcox.table.perc

#SNPs within genes
gene.wilcox.table.perc <- rbind(t2_3=c(nrow(wilcox.glm.t2_3.p.gene[wilcox.glm.t2_3.p.gene$X1<0.05,]),
                                       nrow(wilcox.glm.t2_3.p.gene[wilcox.glm.t2_3.p.gene$X2<0.05,]),
                                       nrow(wilcox.glm.t2_3.p.gene[wilcox.glm.t2_3.p.gene$X3<0.05,]),
                                       nrow(wilcox.glm.t2_3.p.gene[wilcox.glm.t2_3.p.gene$X4<0.05,])),
                                t3_4=c(nrow(wilcox.glm.t3_4.p.gene[wilcox.glm.t3_4.p.gene$X1<0.05,]),
                                       nrow(wilcox.glm.t3_4.p.gene[wilcox.glm.t3_4.p.gene$X2<0.05,]),
                                       nrow(wilcox.glm.t3_4.p.gene[wilcox.glm.t3_4.p.gene$X3<0.05,]),
                                       nrow(wilcox.glm.t3_4.p.gene[wilcox.glm.t3_4.p.gene$X4<0.05,])),
                                t2_3=c(nrow(wilcox.glm.t2_4.p.gene[wilcox.glm.t2_4.p.gene$X1<0.05,]),
                                       nrow(wilcox.glm.t2_4.p.gene[wilcox.glm.t2_4.p.gene$X2<0.05,]),
                                       nrow(wilcox.glm.t2_4.p.gene[wilcox.glm.t2_4.p.gene$X3<0.05,]),
                                       nrow(wilcox.glm.t2_4.p.gene[wilcox.glm.t2_4.p.gene$X4<0.05,])))
colnames(gene.wilcox.table.perc) <- c("Ion-transport", "Seasonal", "Fecundity", "Starvation")
gene.wilcox.table.perc



#################################
## Additional GLM stats for MS ##
#################################
glm.rowsums <- cbind(p01=rowSums(E.glm.contrastp[,c(4:6)] < 0.01),
                     fdr05=rowSums(E.glm.contrastp[,c(10:12)] < 0.05))

# total number of SNPs significant at p<0.01 in any contrast (excluding t1)
nrow(glm.rowsums[glm.rowsums[,1]>0,]) 
#[1] 13684

# total number of SNPs significant at fdr<0.05 in any contrast (excluding t1)
nrow(glm.rowsums[glm.rowsums[,2]>0,]) 
#[1] 23

#next run hypergeometric test across each chromosome, possible evidence of sweep
#P<0.01
glm_phyper_chrom.p01 <- c() #store phyper test p-values
glm_perc_chrom.p01 <- c() #percent of significant SNPs per chrom at p<0.01
glm_perc_chrom <- c() #percent of SNPs per chrom
for(chr in unique(sites.17$chrom)) { #loop through each chrom
  q <- nrow(glm.rowsums[sites.17$chrom==chr & glm.rowsums[,1]>0,])
  m <- nrow(glm.rowsums[sites.17$chrom==chr,])
  n <- nrow(glm.rowsums[sites.17$chrom!=chr,])
  k <- nrow(glm.rowsums[glm.rowsums[,1]>0,])
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  glm_phyper_chrom.p01 <- append(glm_phyper_chrom.p01,phyper(q-1, m, n, k, 
                                                             lower.tail = FALSE, 
                                                             log.p = FALSE))
  ## Calculate significant SNPs per chrom and SNPs per chrom for comparison
  glm_perc_chrom.p01 <- append(glm_perc_chrom.p01,(q/k)*100)
  glm_perc_chrom <- append(glm_perc_chrom,m/(m+n)*100)
}
unique(sites.17$chrom)
#[1] "2L" "2R" "3L" "3R" "X" 
glm_phyper_chrom.p01 < 0.05
#[1] FALSE FALSE FALSE  TRUE FALSE
glm_phyper_chrom.p01
#[1] 1.000000e+00 1.000000e+00 6.890061e-02 1.968629e-85 1.000000e+00
glm_perc_chrom.p01
#[1] 20.60801 16.52295 22.36188 30.49547 10.01169
glm_perc_chrom
#[1] 23.40615 18.74244 21.83522 23.23337 12.78282
glm_perc_chrom.p01/glm_perc_chrom #proportional enrichment
#[1] 0.8804529 0.8815791 1.0241198 1.3125721 0.7832146

# Looks like significant SNPs are enriched on 3R, as expected from manhattan 
# plots. 3R contains 30.5% of all significant SNPs but 23.2% of SNPs overall

##Repeat for separate time-point contrasts
glm_phyper_chrom.p01.all.contrast <- c() #store phyper test p-values
glm_perc_chrom.p01.all.contrast <- c() #percent of significant SNPs per chrom at p<0.01
glm_perc_chrom.all.contrast <- c() #percent of SNPs per chrom

for(x in c(4:6)) { #loop through each time-point contrast
  #make some objects to store temporary stats before adding them to tables
  temp <- c()
  temp2 <- c()
  temp3 <- c()
  for(chr in unique(sites.17$chrom)) { #loop through each chrom
    q <- nrow(E.glm.contrastp[sites.17$chrom==chr & E.glm.contrastp[,x]<0.01,])
    m <- nrow(E.glm.contrastp[sites.17$chrom==chr,])
    n <- nrow(E.glm.contrastp[sites.17$chrom!=chr,])
    k <- nrow(E.glm.contrastp[E.glm.contrastp[,x]<0.01,])
    ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
    temp <- append(temp,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
    ## Calculate significant SNPs per chrom and SNPs per chrom for comparison
    temp2 <- append(temp2,(q/k)*100)
    temp3 <- append(temp3,m/(m+n)*100)
  }
  # Store all stats in tables
  glm_phyper_chrom.p01.all.contrast <- rbind(glm_phyper_chrom.p01.all.contrast,temp)
  glm_perc_chrom.p01.all.contrast <- rbind(glm_perc_chrom.p01.all.contrast,temp2)
  glm_perc_chrom.all.contrast <- rbind(glm_perc_chrom.all.contrast,temp3)
  # Rename rows based on GLM contrast
  rownames(glm_phyper_chrom.p01.all.contrast)[x-3] <- colnames(E.glm.contrastp)[x]
  rownames(glm_perc_chrom.p01.all.contrast)[x-3] <- colnames(E.glm.contrastp)[x]
  rownames(glm_perc_chrom.all.contrast)[x-3] <- colnames(E.glm.contrastp)[x]
}
# Fix table headers
colnames(glm_phyper_chrom.p01.all.contrast) <- unique(sites.17$chrom)
colnames(glm_perc_chrom.p01.all.contrast) <- unique(sites.17$chrom)
colnames(glm_perc_chrom.all.contrast) <- unique(sites.17$chrom)

# Which chroms are significantly enriched in each contrast?
glm_phyper_chrom.p01.all.contrast < 0.05
#          2L    2R    3L   3R     X
#t2_3.p FALSE FALSE FALSE TRUE FALSE
#t2_4.p FALSE FALSE FALSE TRUE FALSE
#t3_4.p FALSE  TRUE FALSE TRUE FALSE

# What are the p-values
glm_phyper_chrom.p01.all.contrast
#              2L        2R       3L       3R         X
#t2_3.p 1.0277854 0.8486865 1.016899 1.168968 0.8350087
#t2_4.p 0.8052396 0.7549845 1.028190 1.544598 0.6778791
#t3_4.p 0.8418064 1.1277059 1.047393 1.056838 0.9181555

# Calculate percent proportional enrichment
glm_perc_chrom.p01.all.contrast/glm_perc_chrom.all.contrast
#              2L        2R       3L       3R         X
#t2_3.p 1.0277854 0.8486865 1.016899 1.168968 0.8350087
#t2_4.p 0.8052396 0.7549845 1.028190 1.544598 0.6778791
#t3_4.p 0.8418064 1.1277059 1.047393 1.056838 0.9181555

### Change significance threshold to FDR<0.05 and repeat analysis ###
glm_phyper_chrom.fdr05 <- c() #store phyper test p-values
glm_perc_chrom.fdr05 <- c() #percent of significant SNPs per chrom at p<0.01
glm_perc_chrom <- c() #percent of SNPs per chrom
for(chr in unique(sites.17$chrom)) { #loop through each chrom
  q <- nrow(glm.rowsums[sites.17$chrom==chr & glm.rowsums[,2]>0,])
  m <- nrow(glm.rowsums[sites.17$chrom==chr,])
  n <- nrow(glm.rowsums[sites.17$chrom!=chr,])
  k <- nrow(glm.rowsums[glm.rowsums[,2]>0,])
  ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
  glm_phyper_chrom.fdr05 <- append(glm_phyper_chrom.fdr05,phyper(q-1, m, n, k, 
                                                                 lower.tail = FALSE, 
                                                                 log.p = FALSE))
  ## Calculate significant SNPs per chrom and SNPs per chrom for comparison
  glm_perc_chrom.fdr05 <- append(glm_perc_chrom.fdr05,(q/k)*100)
  glm_perc_chrom <- append(glm_perc_chrom,m/(m+n)*100)
}
#only 2L, 3L, and 3R analyses run due only single SNPs on 2R and X being significant
glm_phyper_chrom.fdr05 < 0.05
#[1] FALSE FALSE FALSE
glm_phyper_chrom.fdr05
#[1] 0.4594848 0.2214947 0.1439677
glm_perc_chrom.fdr05
#[1] 26.08696 30.43478 34.78261
glm_perc_chrom[c(1,3,4)] #only 2L, 3L, 3R
#[1] 23.40615 21.83522 23.23337
glm_perc_chrom.fdr05/glm_perc_chrom[c(1,3,4)] #proportional enrichment
#[1] 1.114534 1.393839 1.497097

# Looks like significant SNPs are elevated on 2L, 3L, and 3R, but none reach 
# level of significant enrichment from phyper tests. This is likely related to
# the small sample size of significant SNPs at this threshold.

##Repeat for separate time-point contrasts
glm_phyper_chrom.fdr05.all.contrast <- c() #store phyper test p-values
glm_perc_chrom.fdr05.all.contrast <- c() #percent of significant SNPs per chrom at p<0.01
glm_perc_chrom.all.contrast <- c() #percent of SNPs per chrom

for(x in c(10:12)) { #loop through each time-point contrast
  #make some objects to store temporary stats before adding them to tables
  temp <- c()
  temp2 <- c()
  temp3 <- c()
  for(chr in unique(sites.17$chrom)) { #loop through each chrom
    q <- nrow(E.glm.contrastp[sites.17$chrom==chr & E.glm.contrastp[,x]<0.01,])
    m <- nrow(E.glm.contrastp[sites.17$chrom==chr,])
    n <- nrow(E.glm.contrastp[sites.17$chrom!=chr,])
    k <- nrow(E.glm.contrastp[E.glm.contrastp[,x]<0.01,])
    ## Fill in Values below, and make sure to subtract 1 from overlap for phyper
    temp <- append(temp,phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE))
    ## Calculate significant SNPs per chrom and SNPs per chrom for comparison
    temp2 <- append(temp2,(q/k)*100)
    temp3 <- append(temp3,m/(m+n)*100)
  }
  # Store all stats in tables
  glm_phyper_chrom.fdr05.all.contrast <- rbind(glm_phyper_chrom.fdr05.all.contrast,temp)
  glm_perc_chrom.fdr05.all.contrast <- rbind(glm_perc_chrom.fdr05.all.contrast,temp2)
  glm_perc_chrom.all.contrast <- rbind(glm_perc_chrom.all.contrast,temp3)
  # Rename rows based on GLM contrast
  rownames(glm_phyper_chrom.fdr05.all.contrast)[x-9] <- colnames(E.glm.contrastp)[x]
  rownames(glm_perc_chrom.fdr05.all.contrast)[x-9] <- colnames(E.glm.contrastp)[x]
  rownames(glm_perc_chrom.all.contrast)[x-9] <- colnames(E.glm.contrastp)[x]
}
# Fix table headers
colnames(glm_phyper_chrom.fdr05.all.contrast) <- unique(sites.17$chrom)
colnames(glm_perc_chrom.fdr05.all.contrast) <- unique(sites.17$chrom)
colnames(glm_perc_chrom.all.contrast) <- unique(sites.17$chrom)

# Which chroms are significantly enriched in each contrast?
glm_phyper_chrom.fdr05.all.contrast < 0.05
#            2L    2R    3L    3R     X
#t2_3.fdr FALSE FALSE FALSE FALSE FALSE
#t2_4.fdr FALSE FALSE FALSE FALSE FALSE
#t3_4.fdr FALSE FALSE FALSE FALSE FALSE

# What are the p-values
glm_phyper_chrom.fdr05.all.contrast
#                2L 2R        3L        3R         X
#t2_3.fdr 1.0000000  1 1.0000000 1.0000000 1.0000000
#t2_4.fdr 0.5506533  1 1.0000000 0.1368544 1.0000000
#t3_4.fdr 0.5506533  1 0.5224343 1.0000000 0.3365534

# Calculate percent proportional enrichment
glm_perc_chrom.fdr05.all.contrast/glm_perc_chrom.all.contrast
#               2L  2R       3L       3R        X
#t2_3.fdr      NaN NaN      NaN      NaN      NaN
#t2_4.fdr 1.424127   0 0.000000 2.869436 0.000000
#t3_4.fdr 1.424127   0 1.526586 0.000000 2.607666

## Nothing significant, since very few, if any, SNPs show significant AF changes
## in an individual time-point contrast

################################################################################
### Look into the possibility that well known inversions on 3R could explain ###
### any enrichment patterns on 3R, as opposed to a putative selective sweep. ###
################################################################################

## Use GLM p<0.01 as a low stringency filter
inv.markers <- read.table("./functional.analysis/Inversion_markers.bed", header=TRUE) #inversion loci
colnames(inv.markers)[2] <- "pos"
#double checked and confirmed that these loci correspond to dm5 genome build
#join GLM results (significant in any timepoint contrast) with inversion bed
inv.markers.pvals <- inv.markers %>% inner_join(cbind(sites.17,glm.rowsums), 
                                                by=c('chrom','pos'))
## Check biggest inversion "In(3R)Mo"
nrow(inv.markers[inv.markers$Inversion == "In(3R)Mo",]) #how many inversion markers
nrow(inv.markers.pvals[inv.markers.pvals$Inversion == "In(3R)Mo",]) #how many markers in 2017 dataset
#how many inversion markers significant at unadjusted p
nrow(inv.markers.pvals[inv.markers.pvals$Inversion == "In(3R)Mo" & inv.markers.pvals$p01 > 0,])

##Only inversion with any significant results, so might as well check overlap significance
q <- nrow(inv.markers.pvals[inv.markers.pvals$Inversion == "In(3R)Mo" & inv.markers.pvals$p01 > 0,])
m <- nrow(inv.markers.pvals[inv.markers.pvals$Inversion == "In(3R)Mo",])
n <- nrow(glm.rowsums) - m
k <- nrow(glm.rowsums[glm.rowsums[,1] > 0,])
## Fill in Values below, and make sure to subtract 1 from overlap for phyper
phyper(q-1, m, n, k, lower.tail = FALSE, log.p = FALSE)
#[1] 0.002909615

## Technically yes, In(3R)Mo and GLM significant SNPs, but only 5/13684 significant
## SNPs are found within inversion. Inversion SNPs don't appear to contribute much 
## to other enrichment on 3R.

## Check contribution of SNPs within other 3R inversions
#Check another 3R inversion "In(3R)Payne"
nrow(inv.markers[inv.markers$Inversion == "In(3R)Payne",]) #how many inversion markers
nrow(inv.markers.pvals[inv.markers.pvals$Inversion == "In(3R)Payne",]) #how many markers in 2017 dataset
#how many inversion markers significant at unadjusted p
nrow(inv.markers.pvals[inv.markers.pvals$Inversion == "In(3R)Payne" & inv.markers.pvals$p01 > 0,])

#Check another 3R inversion "In(3R)C"
nrow(inv.markers[inv.markers$Inversion == "In(3R)C",]) #how many inversion markers
nrow(inv.markers.pvals[inv.markers.pvals$Inversion == "In(3R)C",]) #how many markers in 2017 dataset
#how many inversion markers significant at unadjusted p
nrow(inv.markers.pvals[inv.markers.pvals$Inversion == "In(3R)C" & inv.markers.pvals$p01 > 0,])

#Check another 3R inversion "In(3R)K"
nrow(inv.markers[inv.markers$Inversion == "In(3R)K",]) #how many inversion markers
nrow(inv.markers.pvals[inv.markers.pvals$Inversion == "In(3R)K",]) #how many markers in 2017 dataset
#how many inversion markers significant at unadjusted p
nrow(inv.markers.pvals[inv.markers.pvals$Inversion == "In(3R)K" & inv.markers.pvals$p01 > 0,])

## No other significant SNPs within inversions

## How many significant GLM SNPs from any contrast overlap with functional candidate lists?
#GLM p<0.01 within functional gene categories (+/-500bp to account for regulatory elements)
nrow(manh.gene.cand[glm.rowsums[,1]>0 & manh.gene.cand$ion>0,])
nrow(manh.gene.cand[glm.rowsums[,1]>0 & manh.gene.cand$seasonal>0,])
nrow(manh.gene.cand[glm.rowsums[,1]>0 & manh.gene.cand$fecundity>0,])
nrow(manh.gene.cand[glm.rowsums[,1]>0 & manh.gene.cand$starvation>0,])
#GLM p<0.01 within functional SNPs
nrow(manh.gene.cand[glm.rowsums[,1]>0 & manh.snp.cand$seasonal>0,])
nrow(manh.gene.cand[glm.rowsums[,1]>0 & manh.snp.cand$fecundity>0,])
nrow(manh.gene.cand[glm.rowsums[,1]>0 & manh.snp.cand$starvation>0,])
#GLM fdr<0.05 within functional gene categories (+/-500bp to account for regulatory elements)
nrow(manh.gene.cand[glm.rowsums[,2]>0 & manh.gene.cand$ion>0,])
nrow(manh.gene.cand[glm.rowsums[,2]>0 & manh.gene.cand$seasonal>0,])
nrow(manh.gene.cand[glm.rowsums[,2]>0 & manh.gene.cand$fecundity>0,])
nrow(manh.gene.cand[glm.rowsums[,2]>0 & manh.gene.cand$starvation>0,])
#GLM fdr<0.05 within functional SNPs
nrow(manh.gene.cand[glm.rowsums[,2]>0 & manh.snp.cand$seasonal>0,])
nrow(manh.gene.cand[glm.rowsums[,2]>0 & manh.snp.cand$fecundity>0,])
nrow(manh.gene.cand[glm.rowsums[,2]>0 & manh.snp.cand$starvation>0,])

#Only one SNP significant within any functional category at fdr<0.05
manh.gene.cand[glm.rowsums[,2]>0 & manh.gene.cand$ion>0,]
ioncandidate_gene_overlap[581837,] #row 581837 is a SNP within "Dh44-R1"

### Make some plot for SNP in candidate gene Dh44-R1, as seen in t3 vs t4 contrast ###
# Modify sample metadata table for plotting
samps.17.E.F.mod <- transform(samps.17[samps.17$treatment == "E" | samps.17$treatment == "Founder",],
                              rep=as.numeric(factor(cage)))
samps.17.E.F.mod[1,3] <- "E"
samps.17.E.F.mod[2,3] <- "E"
#Combine allele frequency stats with metadata 
af.plot <- cbind(data.frame(af=afmat.17.E.F[581837,]),samps.17.E.F.mod,alpha="0.35")[,c(5,4,6,1,9)]
#summarize allele frequencies by timepoint
af.plot.mean <- af.plot %>%
  group_by(treatment, tpt) %>% 
  summarise_at(vars("af"), mean)
af.plot.mean <- as.data.frame(cbind(cage=c("founder","emean","emean","emean","emean"),af.plot.mean,alpha="1"))
ggplot(af.plot.mean,aes(x=tpt,y=af)) + geom_line(alpha = as.numeric(af.plot.mean$alpha)) + theme_classic2()
#add tpt summary stats to plotting table
af.plot <- rbind(af.plot,af.plot.mean)
af.plot2 <- af.plot #store temp version
#set starting point for all samples as the founder mean
af.plot2[af.plot2$tpt=="1" & af.plot2$alpha==0.35,]$af <- af.plot.mean[1,4] 
af.plot2[af.plot$tpt=="1" & af.plot$alpha==0.35,]$tpt <- 0 
af.plot <- af.plot[-c(1,2),] #remove founder replicates
af.plot2 <- af.plot2[-c(1,2),]
#combine starting points with all other time-point data
af.plot <- rbind(af.plot2[af.plot2$tpt=="0" & af.plot2$alpha==0.35,],af.plot)
af.plot <- af.plot[af.plot$tpt!="1",] #remove tpt1
af.plot$tpt <- as.factor(af.plot$tpt)
#remove tpt1 from tpt mean table
af.plot.mean2 <- af.plot.mean[c(1,3:5),]
af.plot.mean2$tpt <- as.factor(af.plot.mean2$tpt)

# Make primary af plot, which includes tpt means
dh44 <- ggplot(af.plot,aes(x=tpt,y=af, group=cage)) + 
  geom_line(alpha = as.numeric(af.plot$alpha), size=as.numeric(af.plot$alpha)) +
  scale_x_discrete(labels= c("June 15th (founder)", "Sep 23rd", "Oct 18th", "Nov 12th")) +
  #geom_point() +
  xlab("Sample Point") + ylab("Allele Frequency") +
  labs(title = "2R:10273261 (Dh44-R1)") + 
  theme_classic2() + theme(legend.position="none")

#calculate means and se for adding points and error bars per time-point
af.plot.stats <- ddply(af.plot, "tpt", summarise,
                       N    = sum(!is.na(af)),
                       mean = mean(af, na.rm=TRUE),
                       sd   = sd(af, na.rm=TRUE),
                       se   = sd / sqrt(N))
#make line darker
af.plot.stats$cage <- af.plot[af.plot$alpha==1,]$cage

#combine primary plot with more pronounces tpt means and SE bars
dh44.points <- dh44+geom_point(data=af.plot.stats, aes(x=tpt, y=mean),size=5)+
  geom_line(data=af.plot.stats, aes(x=tpt, y=mean), size=2.5)+
  geom_linerange(data=af.plot.stats, aes(x=tpt, y=mean,ymin=mean-se,ymax=mean+se), size=1.25)

#final plot
dh44.points

## Which genes contain significant GLM SNPs ##
# First join SNP metadata with GLM results
snp.gff.fulloverlap.glm <- cbind(snp.gff.fulloverlap,E.glm.contrastp)
# Lets check how many genes contain significant SNPs within each contrast at p<0.01
length(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t2_3.p<0.01,]$V9))
length(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t3_4.p<0.01,]$V9))
length(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t2_4.p<0.01,]$V9))
# Lets check how many genes contain significant SNPs within any contrast at p<0.01
length(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$SelCount.noT1>0,]$V9))

# Lets check how many genes contain significant SNPs within each contrast at fdr<0.05
length(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t2_3.fdr<0.05,]$V9))
length(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t3_4.fdr<0.05,]$V9))
length(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t2_4.fdr<0.05,]$V9))
# Lets check how many genes contain significant SNPs within any contrast at fdr<0.05
length(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$fdrCount.noT1>0,]$V9))

## Export gene lists so that we can run Gene Ontology enrichment in Cytoscape's BiNGO app
write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t2_3.p<0.01,]$V9), 
            file="./functional.analysis/glm.gene.lists/glm.t2_3.p01.genelist.txt", 
            sep="\t", quote = FALSE, row.names = F)
write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t3_4.p<0.01,]$V9), 
            file="./functional.analysis/glm.gene.lists/glm.t3_4.p01.genelist.txt", 
            sep="\t", quote = FALSE, row.names = F)
write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t2_4.p<0.01,]$V9), 
            file="./functional.analysis/glm.gene.lists/glm.t2_4.p01.genelist.txt", 
            sep="\t", quote = FALSE, row.names = F)
write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$SelCount.noT1>0,]$V9), 
            file="./functional.analysis/glm.gene.lists/glm.noT1.p01.fullgenelist.txt", 
            sep="\t", quote = FALSE, row.names = F)

write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t2_3.fdr<0.05,]$V9),
            file="./functional.analysis/glm.gene.lists/glm.t2_3.fdr05.genelist.txt", 
            sep="\t", quote = FALSE, row.names = F)
write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t3_4.fdr<0.05,]$V9),
            file="./functional.analysis/glm.gene.lists/glm.t3_4.fdr05.genelist.txt", 
            sep="\t", quote = FALSE, row.names = F)
write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t2_4.fdr<0.05,]$V9), 
            file="./functional.analysis/glm.gene.lists/glm.t2_4.fdr05.genelist.txt", 
            sep="\t", quote = FALSE, row.names = F)
write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$fdrCount.noT1>0,]$V9), 
            file="./functional.analysis/glm.gene.lists/glm.noT1.fdr05.fullgenelist.txt", 
            sep="\t", quote = FALSE, row.names = F)


## Cross reference SNPeff results with candidate sites
# This SNPeff table contains predicted function impact of SNPs
snpeff.results <- read.delim("./functional.analysis/sites.17.SNPeff_results.mod.txt", header=FALSE, sep = "\t")
dim(snpeff.results)

# A little messy, but we want to re-export gene tables, but remove SNPs that are likely silent
# and have no impact on resulting protein product.
#P<0.01
write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t2_3.p<0.01 & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "intergenic_variant" & 
                                             rowMeans(afmat.17.E) > 0.2 & 
                                             rowMeans(afmat.17.E) < 0.8,]$V9),
            file="./functional.analysis/glm.gene.lists/glm.t2_3.p01.filtgenelist.txt", 
            sep="\t", quote = FALSE, row.names = F)

write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t3_4.p<0.01 & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "intergenic_variant" & 
                                             rowMeans(afmat.17.E) > 0.2 & 
                                             rowMeans(afmat.17.E) < 0.8,]$V9), 
            file="./functional.analysis/glm.gene.lists/glm.t3_4.p01.filtgenelist.txt", 
            sep="\t", quote = FALSE, row.names = F)

write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t2_4.p<0.01 & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "intergenic_variant" & 
                                             rowMeans(afmat.17.E) > 0.2 & 
                                             rowMeans(afmat.17.E) < 0.8,]$V9), 
            file="./functional.analysis/glm.gene.lists/glm.t2_4.p01.filtgenelist.txt", 
            sep="\t", quote = FALSE, row.names = F)

write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$SelCount.noT1>0 & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "intergenic_variant" & 
                                             rowMeans(afmat.17.E) > 0.2 & 
                                             rowMeans(afmat.17.E) < 0.8,]$V9), 
            file="./functional.analysis/glm.gene.lists/glm.noT1.p01.filtgenelist.txt", 
            sep="\t", quote = FALSE, row.names = F)

#FDR<0.05
write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t2_3.fdr<0.05 & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "intergenic_variant" & 
                                             rowMeans(afmat.17.E) > 0.2 & 
                                             rowMeans(afmat.17.E) < 0.8,]$V9), 
            file="./functional.analysis/glm.gene.lists/glm.t2_3.fdr05.filtgenelist.txt", 
            sep="\t", quote = FALSE, row.names = F)

write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t3_4.fdr<0.05 & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "intergenic_variant" & 
                                             rowMeans(afmat.17.E) > 0.2 & 
                                             rowMeans(afmat.17.E) < 0.8,]$V9), 
            file="./functional.analysis/glm.gene.lists/glm.t3_4.fdr05.filtgenelist.txt", 
            sep="\t", quote = FALSE, row.names = F)

write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$t2_4.fdr<0.05 & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "intergenic_variant" & 
                                             rowMeans(afmat.17.E) > 0.2 & 
                                             rowMeans(afmat.17.E) < 0.8,]$V9), 
            file="./functional.analysis/glm.gene.lists/glm.t2_4.fdr05.filtgenelist.txt", 
            sep="\t", quote = FALSE, row.names = F)

write.table(unique(snp.gff.fulloverlap.glm[snp.gff.fulloverlap.glm$fdrCount.noT1>0 & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "synonymous_variant" & 
                                             snpeff.results$V5 != "intergenic_variant" & 
                                             rowMeans(afmat.17.E) > 0.2 & 
                                             rowMeans(afmat.17.E) < 0.8,]$V9), 
            file="./functional.analysis/glm.gene.lists/glm.noT1.fdr05.filtgenelist.txt", 
            sep="\t", quote = FALSE, row.names = F)

