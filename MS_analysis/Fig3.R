require(RColorBrewer)
require(gplots)
#install.packages("Rtsne")
require(ggplot2)
require(plyr)
require(rstudioapi)
require(ggrepel)
require(dplyr)
#install.packages("multcomp")
require(multcomp)
#might be OS dependent
#set working dir
rstudioapi::getActiveDocumentContext
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#/settings

#install.packages("tidyr")
require(tidyr)
require(stringr)
require(data.table)
require(plotly)
require(gridExtra)
require(glmnet)
require(selectiveInference)
require(Rmpfr)
require(igraph)
require(biomaRt)
require(ggrepel)

#importing the data
df.raw <- as.data.frame(fread(".\\proteinGroups.txt", header=TRUE, sep="\t"))
colnames(df.raw) <- gsub(" ",".",colnames(df.raw))

#selecting main columns
df <- df.raw %>% dplyr::select(id, Majority.protein.IDs, Gene.names, `Razor.+.unique.peptides`, Unique.peptides,
                               LFQ.intensity.216_MadForAnja_IP_OnBeadDigest_45S_Rep1:LFQ.intensity.802_MadForAnja_OnBeadDigest_ITS_5_3_Sample_18,
                               Only.identified.by.site, Reverse, Potential.contaminant, Fasta.headers)

#Annotating nucleolar proteins
require(biomaRt)
mart <- useEnsembl("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
listAttributes(mart)
bmIDs = getBM(attributes=c('uniprotswissprot','go_id', "namespace_1003"),mart = mart)
bmIDs <- bmIDs %>% dplyr::filter(!uniprotswissprot == "" & !go_id == "" & !namespace_1003 == "" &
                                   !is.na(uniprotswissprot) & !is.na(go_id) & !is.na(namespace_1003))
head(bmIDs)

bmIDs.sel <- bmIDs %>% dplyr::filter(go_id == "GO:0005730" & namespace_1003 == "cellular_component")

df <- df %>% rowwise() %>% do({
  temp <- .
  temp$is.nol = any(strsplit(temp$Majority.protein.IDs, ";", T)[[1]] %in% bmIDs.sel$uniprotswissprot)
  as_tibble(temp)
}) %>% ungroup() %>% as.data.frame()

table(df$is.nol)

#Renaming columns
names(df)[names(df) == "LFQ.intensity.216_MadForAnja_IP_OnBeadDigest_45S_Rep1"] <- "S45_3min_1"
names(df)[names(df) == "LFQ.intensity.216_MadForAnja_IP_OnBeadDigest_45S_Rep2"] <- "S45_3min_2"
names(df)[names(df) == "LFQ.intensity.216_MadForAnja_IP_OnBeadDigest_45S_Rep3"] <- "S45_3min_3"
names(df)[names(df) == "LFQ.intensity.216_MadForAnja_IP_OnBeadDigest_45S_Rep4"] <- "S45_3min_4"

names(df)[names(df) == "LFQ.intensity.609_MadforAnja_RNAproxLabel_45S_Rep1_1"] <- "S45_5min_1"
names(df)[names(df) == "LFQ.intensity.609_MadforAnja_RNAproxLabel_45S_Rep2_2"] <- "S45_5min_2"
names(df)[names(df) == "LFQ.intensity.609_MadforAnja_RNAproxLabel_45S_Rep3_3"] <- "S45_5min_3"
names(df)[names(df) == "LFQ.intensity.609_MadforAnja_RNAproxLabel_45S_Rep4_4"] <- "S45_5min_4"

names(df)[names(df) == "LFQ.intensity.609_MadforAnja_RNAproxLabel_47S_Rep1_5"] <- "S47_5min_1"
names(df)[names(df) == "LFQ.intensity.609_MadforAnja_RNAproxLabel_47S_Rep2_6"] <- "S47_5min_2"
names(df)[names(df) == "LFQ.intensity.609_MadforAnja_RNAproxLabel_47S_Rep3_7"] <- "S47_5min_3"
names(df)[names(df) == "LFQ.intensity.609_MadforAnja_RNAproxLabel_47S_Rep4_8"] <- "S47_5min_4"

names(df)[names(df) == "LFQ.intensity.609_MadforAnja_RNAproxLabel_ITS2_Rep1_13"] <- "ITS2_5min_1"
names(df)[names(df) == "LFQ.intensity.609_MadforAnja_RNAproxLabel_ITS2_Rep2_14"] <- "ITS2_5min_2"
names(df)[names(df) == "LFQ.intensity.609_MadforAnja_RNAproxLabel_ITS2_Rep3_15"] <- "ITS2_5min_3"
names(df)[names(df) == "LFQ.intensity.609_MadforAnja_RNAproxLabel_ITS2_Rep4_16"] <- "ITS2_5min_4"

names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_45_1_1_Sample_7"] <- "S45_1min_1"
names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_45_1_2_Sample_8"] <- "S45_1min_2"
names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_45_1_3_Sample_9"] <- "S45_1min_3"

names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_45_5_1_Sample_10"] <- "S45_15min_1"
names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_45_5_2_Sample_11"] <- "S45_15min_2"
names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_45_5_3_Sample_12"] <- "S45_15min_3"

names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_47_1_1_Sample_1"] <- "S47_1min_1"
names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_47_1_2_Sample_2"] <- "S47_1min_2"
names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_47_1_3_Sample_3"] <- "S47_1min_3"

names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_47_5_1_Sample_4"] <- "S47_15min_1"
names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_47_5_2_Sample_5"] <- "S47_15min_2"
names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_47_5_3_Sample_6"] <- "S47_15min_3"

names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_ITS_1_1_Sample_13"] <- "ITS2_1min_1"
names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_ITS_1_2_Sample_14"] <- "ITS2_1min_2"
names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_ITS_1_3_Sample_15"] <- "ITS2_1min_3"

names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_ITS_5_1_Sample_16"] <- "ITS2_15min_1"
names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_ITS_5_2_Sample_17"] <- "ITS2_15min_2"
names(df)[names(df) == "LFQ.intensity.802_MadForAnja_OnBeadDigest_ITS_5_3_Sample_18"] <- "ITS2_15min_3"

#Removing bad runs
df <- df %>% dplyr::select(-ITS2_1min_3, -ITS2_5min_4, -S47_5min_4)# Bad runs

#grouping for easier handling

S45_1min <- names(df) %in% paste0("S45_1min_", c(1:4))
S45_3min <- names(df) %in% paste0("S45_3min_", c(1:4))
S45_5min <- names(df) %in% paste0("S45_5min_", c(1:4))
S45_15min <- names(df) %in% paste0("S45_15min_", c(1:4))

S47_1min <- names(df) %in% paste0("S47_1min_", c(1:4))
S47_5min <- names(df) %in% paste0("S47_5min_", c(1:4))
S47_15min <- names(df) %in% paste0("S47_15min_", c(1:4))

ITS2_1min <- names(df) %in% paste0("ITS2_1min_", c(1:4))
ITS2_5min <- names(df) %in% paste0("ITS2_5min_", c(1:4))
ITS2_15min <- names(df) %in% paste0("ITS2_15min_", c(1:4))

alls <- S45_1min | S45_3min | S45_5min | S45_15min |
  S47_1min | S47_5min | S47_15min |
  ITS2_1min | ITS2_5min | ITS2_15min
sum(alls)
names(df)[alls]

#log2 transformation
df[,alls] <- sapply(df[,alls], log2)
df[,alls][sapply(df[,alls], is.infinite)] <- NA_real_

#pre-filtering the data
df <- df %>% dplyr::filter(!Only.identified.by.site == "+" & !Reverse == "+" & !Potential.contaminant == "+")
df <- df %>% dplyr::filter(`Razor.+.unique.peptides`>=2)
#Filtering for detected at least 2 times in any group
nrow(df[rowSums(!sapply(df[alls], is.na))>0,])
df <- df[rowSums(!sapply(df[S45_1min], is.na))==sum(S45_1min) |
           rowSums(!sapply(df[S45_5min], is.na))==sum(S45_5min) |
           rowSums(!sapply(df[S45_15min], is.na))==sum(S45_15min) |
           rowSums(!sapply(df[S47_1min], is.na))==sum(S47_1min) |
           rowSums(!sapply(df[S47_5min], is.na))==sum(S47_5min) |
           rowSums(!sapply(df[S47_15min], is.na))==sum(S47_15min) |
           rowSums(!sapply(df[ITS2_1min], is.na))==sum(ITS2_1min) |
           rowSums(!sapply(df[ITS2_5min], is.na))==sum(ITS2_5min) |
           rowSums(!sapply(df[ITS2_15min], is.na))==sum(ITS2_15min),]

df.norm <- df

#imputation
width <- 0.3
downshift <- 1.8

u <- mean(unlist(df.norm[,alls]), na.rm = TRUE)
stdev <- sd(unlist(df.norm[,alls]), na.rm = TRUE)

df.imp <- df.norm
set.seed(1337)
df.imp[,alls][sapply(df.imp[,alls], is.na)] <- rnorm(sum(sapply(df.imp[,alls], is.na)), u-downshift*stdev, width*stdev)

rm(width, downshift, u, stdev)
#Saving interim results
fwrite(df, file=".\\res\\df.txt", quote = F, col.names = T, row.names = F, sep="\t")
fwrite(df.norm, file=".\\res\\df.norm.txt", quote = F, col.names = T, row.names = F, sep="\t")
fwrite(df.imp, file=".\\res\\df.imp.txt", quote = F, col.names = T, row.names = F, sep="\t")

#UMAP
library(umap)

set.seed(1337)
umap.sample <- umap(t(df.imp[,alls] - apply(df.imp[,alls], 1, median)))

df.umap.sample <- data.frame(umap_x = umap.sample$layout[,1],
                             umap_y = umap.sample$layout[,2],
                             Sample.name = unlist(dimnames(umap.sample$layout)[[1]])) %>%
  separate(Sample.name, c("bait","labtime","techr"), sep = "_", remove = F)

(umap.plot <- ggplot(df.umap.sample) +
  geom_text(aes(x = umap_x, y=umap_y, label = paste0(bait, "_", labtime), color = labtime) , show.legend = TRUE) +
  theme_classic()+theme(
    axis.line.x = element_line(colour = 'black', size=1.5, linetype='solid'), axis.ticks.length=unit(.5, "cm"), axis.ticks = element_line(size = 1.5),
    axis.line.y = element_line(colour = 'black', size=1.5, linetype='solid'), text = element_text(size=15, face="bold")))

ggsave(filename = ".\\res\\umap.plot_imp.pdf",umap.plot,
       device = "pdf", width=11, height=6, units = "in", useDingbats=F)
rm(umap.def2, umap.sample, umap.plot, df.umap.sample)

umap.def2 <- umap.defaults
umap.def2$n_neighbors <- 5
set.seed(1337)
umap.sample <- umap(t(df.imp[,S45_1min | S45_5min | S45_15min] - apply(df.imp[,S45_1min | S45_5min | S45_15min], 1, median)),umap.def2)

df.umap.sample <- data.frame(umap_x = umap.sample$layout[,1],
                             umap_y = umap.sample$layout[,2],
                             Sample.name = unlist(dimnames(umap.sample$layout)[[1]])) %>%
  separate(Sample.name, c("bait","labtime","techr"), sep = "_", remove = F)

(umap.plot <- ggplot(df.umap.sample) +
    geom_text(aes(x = umap_x, y=umap_y, label = paste0(bait, "_", labtime), color = labtime) , show.legend = TRUE) +
    theme_classic()+theme(
      axis.line.x = element_line(colour = 'black', size=1.5, linetype='solid'), axis.ticks.length=unit(.5, "cm"), axis.ticks = element_line(size = 1.5),
      axis.line.y = element_line(colour = 'black', size=1.5, linetype='solid'), text = element_text(size=15, face="bold")))

ggsave(filename = ".\\res\\umap.plot_imp_S45timecourse.pdf",umap.plot,
       device = "pdf", width=11, height=6, units = "in", useDingbats=F)
rm(umap.def2, umap.sample, umap.plot, df.umap.sample)

umap.def2 <- umap.defaults
umap.def2$n_neighbors <- 4
set.seed(1000)
umap.sample <- umap(t(df.imp[,S47_5min | S45_5min | ITS2_5min] - apply(df.imp[,S47_5min | S45_5min | ITS2_5min], 1, median)),umap.def2)

df.umap.sample <- data.frame(umap_x = umap.sample$layout[,1],
                             umap_y = umap.sample$layout[,2],
                             Sample.name = unlist(dimnames(umap.sample$layout)[[1]])) %>%
  separate(Sample.name, c("bait","labtime","techr"), sep = "_", remove = F)

(umap.plot <- ggplot(df.umap.sample) +
    geom_text(aes(x = umap_x, y=umap_y, label = paste0(bait, "_", labtime), color = labtime) , show.legend = TRUE) +
    theme_classic()+theme(
      axis.line.x = element_line(colour = 'black', size=1.5, linetype='solid'), axis.ticks.length=unit(.5, "cm"), axis.ticks = element_line(size = 1.5),
      axis.line.y = element_line(colour = 'black', size=1.5, linetype='solid'), text = element_text(size=15, face="bold")))

ggsave(filename = ".\\res\\umap.plot_imp_all_5min.pdf",umap.plot,
       device = "pdf", width=11, height=6, units = "in", useDingbats=F)
rm(umap.def2, umap.sample, umap.plot, df.umap.sample)


temp.long <- df %>% dplyr::select(Majority.protein.IDs, Gene.names, is.nol, S45_1min_1:S45_1min_3, S45_5min_1:S45_5min_4, S45_15min_1:S45_15min_3) %>%
  pivot_longer(S45_1min_1:S45_15min_3, names_to = "sp", values_to = "y") %>% separate(sp, into = c("target","timepoint","techr"), sep = "_")
temp.long$timepoint <- factor(temp.long$timepoint, levels = c("1min","5min","15min"))

temp.long.mean <- temp.long %>% group_by(Majority.protein.IDs, timepoint, is.nol) %>% do({
  temp <- as.data.frame(.)
  tibble(y = mean(temp$y, na.rm = T))
}) %>% ungroup() %>% as.data.frame()


#data preparation for LASSO
baits <- c("S47","S45", "ITS2")
labtimes <- c("1min","3min","5min","15min")
techrs <- c("1", "2", "3", "4")

df.long <- df.norm %>% dplyr::select(id, Majority.protein.IDs, Gene.names, S45_3min_1:ITS2_15min_3) %>%
  gather(Sample, Intensity, S45_3min_1:ITS2_15min_3) %>%
  separate(Sample, c("bait","labtime","techr"), sep = "_", remove = F)

df.imp.long <- df.imp %>% dplyr::select(id, Majority.protein.IDs, Gene.names, S45_3min_1:ITS2_15min_3) %>%
  gather(Sample, Intensity, S45_3min_1:ITS2_15min_3) %>%
  separate(Sample, c("bait","labtime","techr"), sep = "_", remove = F)

df.long$bait <- factor(df.long$bait, levels= baits)
df.long$labtime <- factor(df.long$labtime, levels= labtimes)
df.long$techr <- factor(df.long$techr, levels= techrs)

df.imp.long$bait <- factor(df.imp.long$bait, levels= baits)
df.imp.long$labtime <- factor(df.imp.long$labtime, levels= labtimes)
df.imp.long$techr <- factor(df.imp.long$techr, levels= techrs)

df.long <- df.long[order(df.long$techr),]
df.long <- df.long[order(df.long$labtime),]
df.long <- df.long[order(df.long$bait),]

df.imp.long <- df.imp.long[order(df.imp.long$techr),]
df.imp.long <- df.imp.long[order(df.imp.long$labtime),]
df.imp.long <- df.imp.long[order(df.imp.long$bait),]

# Experiment design specification
exp.des <- dplyr::select(df.imp.long,-Majority.protein.IDs,-Gene.names,-id,-Intensity) %>% distinct()
exp.des.m <- model.matrix(~ labtime + labtime:bait, exp.des) #,
exp.des.m <- exp.des.m[,c(1:5,7:9,11:12)]

rownames(exp.des.m) <- paste0(exp.des$bait,"_",exp.des$labtime)
colnames(exp.des.m) <- c("Intercept","lab3","lab5","lab15",
                         "lab1_S47", "lab5_S47", "lab15_S47",
                         "lab1_ITS2", "lab5_ITS2", "lab15_ITS2")


pdf(file = ".\\res\\ExpDesMat.pdf", width = 16, height = 6)
heatmap.2(exp.des.m, Rowv = F, Colv = F, margins = c(12,12), trace = "none",
          #cellnote = round(exp.des.effsc.m,2),
          notecol = "black")
dev.off()

exp.des.m.lasso <- exp.des.m[,-1]
exp.des.m.lasso <- scale(exp.des.m.lasso, T, F)

save.image(".\\backup\\CheckPoint2.RData")

require(glmnet)
require(selectiveInference)
require(gridExtra)
source(".\\V_common_fun.R")

#Running LASSO
sink('LassoSI_output.txt')

df.res <- df.imp %>% group_by(Majority.protein.IDs) %>% do({
  temp <- .
  temp_id <- temp$id
  print(temp$Gene.names)
  
  y <- dplyr::filter(df.imp.long, Majority.protein.IDs == temp$Majority.protein.IDs)
  y.unimp <- dplyr::filter(df.long, Majority.protein.IDs == temp$Majority.protein.IDs)
  y$is.valid <- !sapply(y.unimp$Intensity, is.na)
  
  #print(w.w)
  set.seed(1)
  LassoI_y <- fixedLassoInf_wrapper(mtx=exp.des.m.lasso,
                                    rhs=y$Intensity,
                                    alpha = 0.1,
                                    tailarea_rtol = 0.1,
                                    n_tries = 3,
                                    scales = NA_real_)
  
  pvs.min <- setNames(rep.int(NA_real_, ncol(exp.des.m.lasso)), colnames(exp.des.m.lasso))
  pvs.min[names(LassoI_y[["at_lambda_min"]][["res.min"]][["vars"]])] <- LassoI_y[["at_lambda_min"]][["res.min"]][["pv"]]
  temp[["lambda.min.converged"]] <- LassoI_y[["at_lambda_min"]][["res.min"]]$converged
  
  for (cname in colnames(exp.des.m.lasso)) {
    temp[[paste0(cname, ".min")]] <- LassoI_y[["at_lambda_min"]][["rescaled.glmnet.estimates.min"]][rownames(LassoI_y[["at_lambda_min"]][["rescaled.glmnet.estimates.min"]])==cname]
    temp[[paste0(cname, ".min.p")]] <- pvs.min[[cname]]
  }
  temp[["lambda.min"]] <- LassoI_y[["at_lambda_min"]][["lambda.min"]]
  
  #plotting
  ci0.1 <- data.frame(V1 = rep(0,ncol(exp.des.m.lasso)), V2 = rep(0,ncol(exp.des.m.lasso)))
  if(!is.null(LassoI_y[["at_lambda_min"]][["res.min"]][["ci"]])){
    ci0.1.temp <- as.data.frame(as.matrix(LassoI_y[["at_lambda_min"]][["res.min"]][["ci"]]))
    rownames(ci0.1.temp) <- LassoI_y[["at_lambda_min"]][["res.min"]][["vars"]]
    
    ci0.1$V1[sapply(rownames(ci0.1.temp),as.numeric)] <- ci0.1.temp$V1
    ci0.1$V2[sapply(rownames(ci0.1.temp),as.numeric)] <- ci0.1.temp$V2
  }
  rm(ci0.1.temp)
  
  resc_est <- c(LassoI_y[["at_lambda_min"]][["rescaled.glmnet.estimates.min"]][1],
                rep(0, length(LassoI_y[["at_lambda_min"]][["rescaled.glmnet.estimates.min"]])-1))
  resc_est[LassoI_y[["at_lambda_min"]][["rescaled.glmnet.estimates.min"]]@i+1] <- LassoI_y[["at_lambda_min"]][["rescaled.glmnet.estimates.min"]][LassoI_y[["at_lambda_min"]][["rescaled.glmnet.estimates.min"]]@i+1]
  
  ci0.1 <- rbind(data.frame(V1=resc_est[1], V2=resc_est[1]), ci0.1)
  
  y$resc_est_0.05 <- apply(as.data.frame(exp.des.m.lasso), 1, function(x){sum(ci0.1$V1*c(1,x))})
  y$resc_est_0.5 <- apply(as.data.frame(exp.des.m.lasso), 1, function(x){sum(resc_est*c(1,x))})
  y$resc_est_0.95 <- apply(as.data.frame(exp.des.m.lasso), 1, function(x){sum(ci0.1$V2*c(1,x))})
  
  y$sid <- factor(paste0(y$bait,"_", y$labtime),levels=unique(paste0(y$bait,"_", y$labtime)))
  
  # bplot.limits <- c(floor(min(y$Intensity))-2,
  #                   ceiling(max(y$Intensity))+4)
  # 
  # bx_data <- y %>% dplyr::select(sid, resc_est_0.05, resc_est_0.5, resc_est_0.95, bait, is.valid) %>% distinct()
  # bx_data.coef <- data.frame(coefx=rownames(LassoI_y[["at_lambda_min"]][["rescaled.glmnet.estimates.min"]]))
  # bx_data.coef$p <- rep(NA_real_,nrow(bx_data.coef))
  # bx_data.coef$p[LassoI_y[["at_lambda_min"]][["res.min"]][["vars"]]+1] <- LassoI_y[["at_lambda_min"]][["res.min"]][["pv"]]
  # bx_data.coef$c <- rep(NA_real_,nrow(bx_data.coef))
  # bx_data.coef$c[LassoI_y[["at_lambda_min"]][["rescaled.glmnet.estimates.min"]]@i+1] <- LassoI_y[["at_lambda_min"]][["rescaled.glmnet.estimates.min"]][LassoI_y[["at_lambda_min"]][["rescaled.glmnet.estimates.min"]]@i+1]
  # bx_data.coef$sig <- sapply(sapply(bx_data.coef$p,as.numeric),
  #                            function(x){ifelse(x<0.001,"***",ifelse(x<0.01,"**",ifelse(x<0.05,"*","n.s.")))})
  # 
  # bx_data.coef <- bx_data.coef %>% dplyr::filter(!is.na(c))
  # 
  # bplot <- ggplot() +
  #   # geom_line(data = y, aes(x=sid, y = Intensity, group = donor, color = donor)) +
  #   geom_errorbar(data = bx_data,
  #                 aes(x=sid, ymin = resc_est_0.05, ymax = resc_est_0.95, color = bait), width = 0.5) +
  #   geom_errorbar(data = bx_data,
  #                 aes(x=sid, ymin = resc_est_0.5, ymax = resc_est_0.5, color = bait), width = 0.75) +
  #   geom_point(data = y, aes(x=sid, y = Intensity, color = bait, fill=is.valid), shape = 21, size = 2) +
  #   geom_text(data = bx_data.coef, aes(x=1, y=bplot.limits[2],
  #                                      label = eval(paste0(bx_data.coef$coefx," = ", round(bx_data.coef$c,3),
  #                                                          ", p = ",
  #                                                          round(bx_data.coef$p,4), collapse = "\n"))),hjust = 0, vjust = 1) +
  #   scale_y_continuous(limits = bplot.limits,
  #                      breaks = seq(bplot.limits[1], bplot.limits[2],1)) +
  #   ylab(paste0("Intensity ",temp$Majority.protein.IDs)) +
  #   theme_bw()
  # 
  # # View(LassoI_y)
  # #stop("La")
  # LassoI_y_c.l <- as.data.frame(t(as.matrix(LassoI_y[["cvf"]][["glmnet.fit"]][["beta"]])))
  # LassoI_y_c.l$lnlambda <- log(LassoI_y[["cvf"]][["glmnet.fit"]][["lambda"]])
  # LassoI_y_c.l$cvm <- LassoI_y[["cvf"]][["cvm"]]
  # LassoI_y_c.l$cvsd <- LassoI_y[["cvf"]][["cvsd"]]
  # 
  # cvm.plot <- ggplot() +
  #   geom_vline(aes(xintercept=log(c(LassoI_y[["cvf"]][["lambda.min"]],LassoI_y[["cvf"]][["lambda.1se"]]))),
  #              color = "grey", linetype = "dashed")+
  #   geom_vline(aes(xintercept=log(LassoI_y[["at_lambda_min"]][["lambda.min"]])),
  #              color = "red", linetype = "dashed")+
  #   geom_segment(data = LassoI_y_c.l, aes(x=lnlambda, xend = lnlambda, y=cvm-cvsd, yend=cvm+cvsd)) +
  #   geom_segment(data = LassoI_y_c.l, aes(x=lnlambda+0.025, xend = lnlambda-0.025, y=cvm-cvsd, yend=cvm-cvsd)) +
  #   geom_segment(data = LassoI_y_c.l, aes(x=lnlambda+0.025, xend = lnlambda-0.025, y=cvm+cvsd, yend=cvm+cvsd)) +
  #   geom_point(data = LassoI_y_c.l, aes(x=lnlambda, y=cvm), size = 1.5, color = "red") +
  #   xlab("ln(lambda)") +
  #   ylab("Mean-Squared Error")+
  #   theme_bw()
  # 
  # LassoI_y_c.l <- LassoI_y_c.l %>% gather(coefx, value, 1:ncol(exp.des.m.lasso))
  # coef.plot <- ggplot() +
  #   geom_vline(aes(xintercept=log(c(LassoI_y[["cvf"]][["lambda.min"]],LassoI_y[["cvf"]][["lambda.1se"]]))),
  #              color = "grey", linetype = "dashed")+
  #   geom_vline(aes(xintercept=log(LassoI_y[["at_lambda_min"]][["lambda.min"]])),
  #              color = "red", linetype = "dashed")+
  #   geom_line(data = LassoI_y_c.l, aes(x=lnlambda, y=value, group=coefx, color = coefx), size = 1) +
  #   labs(color = "Coefficient") +
  #   xlab("ln(lambda)") +
  #   ylab("Coefficient value")+
  #   theme_bw()
  # 
  # ggsave(filename = paste0(".\\lasso_plots\\",temp$id,"_",ifelse(temp$Gene.names=="",temp$Majority.protein.IDs, temp$Gene.names),".pdf"),
  #        arrangeGrob(bplot, cvm.plot, coef.plot, nrow = 3, ncol = 1),
  #        device = "pdf", width=30, height=24, units = "in", useDingbats=F, limitsize = FALSE)
  
  rm(y, y.unimp, LassoI_y, pvs.min, bplot.limits, bplot, LassoI_y_c.l, cvm.plot, coef.plot, gna)
  #returning
  as_tibble(temp)
}) %>% ungroup() %>% as.data.frame()
sink()

save.image(".\\backup\\CheckPoint3.RData")

iis <- as.data.frame(fread("..\\ITS2 interactors.txt", header=TRUE, sep="\t"))
iis <- iis %>% separate_longer_delim(Gene.names, delim = ", ") %>% separate_longer_delim(Gene.names, delim = "/") %>%
  mutate(Gene.names = toupper(Gene.names))
#importing the data on nol proteins
nols <- as.data.frame(fread("..\\Protein categorization_Mus musculus_Trupej_20241206.txt", header=TRUE, sep="\t"))
nols <- nols %>% separate_longer_delim(Gene.names, delim = ", ") %>% separate_longer_delim(Gene.names, delim = "/") %>%
  mutate(Gene.names = toupper(Gene.names))

df.res <- df.res %>% rowwise() %>% do({
  temp <- .
  temp$nolloc = paste0(unique(nols$loc[nols$Gene.names %in% sapply(strsplit(temp$Gene.names, ";", T)[[1]], toupper)]), collapse = ";")
  as_tibble(temp)
}) %>% ungroup() %>% as.data.frame()

df.res <- df.res %>% rowwise() %>% do({
  temp <- .
  temp$its2i = paste0(unique(iis$loc[iis$Gene.names %in% sapply(strsplit(temp$Gene.names, ";", T)[[1]], toupper)]), collapse = ";")
  as_tibble(temp)
}) %>% ungroup() %>% as.data.frame()

table(df.res$nolloc)
table(df.res$its2i)

ptt <- function(x){
  pt <- 5
  return(sapply(x, function(x){ifelse(x>pt,pt,x)}))
}

S45.markers <- df.res$Gene.names[df.res$nolloc == "DFC"]
S47.markers <- df.res$Gene.names[df.res$nolloc == "FC"]
ITS2.markers <- df.res$Gene.names[df.res$its2i == "ITS2"]

(vp1 <- ggplot() +
  geom_hline(aes(yintercept = ptt(Inf)), linetype = "dashed") +
  geom_point(data=df.res, aes(x=lab5_S47.min, y=ptt(-log10(lab5_S47.min.p))), color = "grey") +
  geom_point(data=df.res %>% dplyr::filter(nolloc == "DFC"), aes(x=lab5_S47.min, y=ptt(-log10(lab5_S47.min.p))), size = 2, color = "darkred") +
  geom_point(data=df.res %>% dplyr::filter(nolloc == "FC"), aes(x=lab5_S47.min, y=ptt(-log10(lab5_S47.min.p))), size = 2, color = "black") +
  #geom_point(data=df.res %>% dplyr::filter(nolloc == "GC"), aes(x=lab5_S47.min, y=-log10(lab5_S47.min.p)), size = 2, color = "blue") +
  geom_text_repel(data=df.res %>% dplyr::filter(nolloc == "DFC" & Gene.names %in% S45.markers), aes(x=lab5_S47.min, y=ptt(-log10(lab5_S47.min.p)), label = Gene.names), color = "darkred") +
  geom_text_repel(data=df.res %>% dplyr::filter(nolloc == "FC" & Gene.names %in% S47.markers), aes(x=lab5_S47.min, y=ptt(-log10(lab5_S47.min.p)), label = Gene.names), color = "black") +
  xlab("45S - 47S, 5min, LASSO") + ylab("p (-log10)") +
  theme_classic())



