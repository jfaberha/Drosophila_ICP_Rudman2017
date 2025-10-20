########################################################################
## This R script will run GLMs for all ~1.9m SNPs in Dmel E analysis  ##
## from 2017 data. Specifically, the model examines all pairwise      ##
## timepoint contrasts and extracts p-values with emmeans.            ##
########################################################################

#configure r environment
library(emmeans)
setwd("~/Documents/josh_r_analysis/fst/Github_E.2017.analysis/glm.r.analysis")

#load required input files
samps.17.E <- read.table("samps.17.E.txt", header=TRUE)
afmat.17.E <- read.table("afmat.17.E.txt", header=TRUE)
rd.17.E <- read.table("rd.17.E.txt", header=TRUE)

#tables need a bit of reformatting for glm to run
samps.17.E$sampID == cbind.data.frame(colnames(afmat.17.E))
samps.17.E$treatment <- as.factor(samps.17.E$treatment)
samps.17.E$tpt <- as.factor(samps.17.E$tpt)


glmres <- c() #initialize glm object
contrastout <- c() #initialize p-val summary table

#loop through all rows of allele frequency table
for(g in 1:nrow(afmat.17.E)) {
  tryCatch({
    glmdf <- cbind(samps.17.E,t(round(afmat.17.E[g,]*100)))
    colnames(glmdf)[8] <- "count"
    #sapply(glmdf, class)
    glmres[[g]] <- glm(count~tpt, data = glmdf, family = quasipoisson, na.action = na.exclude)
    emtemp <- emmeans(glmres[[g]], pairwise ~ tpt)$contrasts
    #summary(emtemp)$p.value
    contrastout <- rbind(contrastout,summary(emtemp)$p.value)
  }, error=function(e){})
}

#save results so we don't have to do this again
write.table(contrastout, file="E.2017.contrastTPT.txt", sep="\t", quote = FALSE, row.names = F)