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


#importing the data
df.raw <- as.data.frame(fread(".\\proteinGroups.txt", header=TRUE, sep="\t"))
colnames(df.raw) <- gsub(" ",".",colnames(df.raw))

#analysis-relevant columns
df <- df.raw %>% dplyr::select(id, Majority.protein.IDs, Gene.names, `Razor.+.unique.peptides`, Unique.peptides,
                               LFQ.intensity.P_OnBeadDigest_45S_Rep1:LFQ.intensity.P_OnBeadDigest_45S_Rep4,
                               LFQ.intensity.P_OnBeadDigest_Efl_Rep1:LFQ.intensity.P_OnBeadDigest_Efl_Rep4,
                               LFQ.intensity.P_OnBeadDigest_Malat1_Rep1:LFQ.intensity.P_OnBeadDigest_Malat1_Rep4,
                               LFQ.intensity.P_OnBeadDigest_NORAD_Rep1:LFQ.intensity.P_OnBeadDigest_NORAD_Rep4,
                               Only.identified.by.site, Reverse, Potential.contaminant, Fasta.headers)

colnames(df)[6:21] <- sapply(colnames(df)[6:21], function(x){paste0(c(strsplit(x, "_", T)[[1]][3],
                                                                      substr(strsplit(x, "_", T)[[1]][4], start = 4, stop = 4)), collapse = "_")})

S45 <- names(df) %in% paste0("45S_",c(1:4))
EFL <- names(df) %in% paste0("Efl_",c(1:4))
MALAT1 <- names(df) %in% paste0("Malat1_",c(1:4))
NORAD <- names(df) %in% paste0("NORAD_",c(1:4))


alls <- S45 | EFL | MALAT1 | NORAD

sum(alls)
sum(S45)

df[,alls] <- sapply(df[,alls], log2)
df[,alls][sapply(df[,alls], is.infinite)] <- NA_real_

#annotating NOL proteins, alternative: https://www.biorxiv.org/content/10.1101/2024.07.03.601239v1
require(biomaRt)
mart <- useEnsembl("ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")

#View(listAttributes(mart))

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

#Quality check
#plotting number of quanted pgs per sample
(quant.pgs <- ggplot() +
  geom_col(aes(x=colnames(df[,alls]), y=apply(df[,alls],2,function(x){nrow(df)-sum(is.na(x))}))) +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("number of quantified PGs"))

ggsave(filename = ".\\qc\\n_quant_pgs.pdf",quant.pgs,
       device = "pdf", width=12, height=6, units = "in", useDingbats=F)
rm(quant.pgs)

#pre-filtering the data
df <- df %>% dplyr::filter(!Only.identified.by.site == "+" & !Reverse == "+" & !Potential.contaminant == "+")
df <- df %>% dplyr::filter(`Razor.+.unique.peptides`>=2)

#Filtering for detected at least 1 times in any group - just to get numbers
nrow(df[rowSums(!sapply(df[S45], is.na))>=1 |
           rowSums(!sapply(df[EFL], is.na))>=1 |
           rowSums(!sapply(df[MALAT1], is.na))>=1 |
           rowSums(!sapply(df[NORAD], is.na))>=1,])

#Filtering for detected at least 3 times in any group
df <- df[rowSums(!sapply(df[S45], is.na))>=3 |
           rowSums(!sapply(df[EFL], is.na))>=3 |
           rowSums(!sapply(df[MALAT1], is.na))>=3 |
           rowSums(!sapply(df[NORAD], is.na))>=3,]

#imputation
width <- 0.3
downshift <- 1.8

u <- mean(unlist(df[,alls]), na.rm = TRUE)
stdev <- sd(unlist(df[,alls]), na.rm = TRUE)

df.imp <- df
set.seed(1337)
df.imp[,alls][sapply(df.imp[,alls], is.na)] <- rnorm(sum(sapply(df.imp[,alls], is.na)), u-downshift*stdev, width*stdev)
rm(stdev, u, width, downshift)

fwrite(df, file=".\\res\\df.txt", quote = F, col.names = T, row.names = F, sep="\t")
fwrite(df.imp, file=".\\res\\df.imp.txt", quote = F, col.names = T, row.names = F, sep="\t")

#Plotting UMAP - all samples
library(umap)
set.seed(1337)
umap.sample <- umap(t(df.imp[,alls] - apply(df.imp[,alls], 1, median)))

df.umap.sample <- data.frame(umap_x = umap.sample$layout[,1],
                             umap_y = umap.sample$layout[,2],
                             Sample.name = unlist(dimnames(umap.sample$layout)[[1]])) %>%
  separate(Sample.name, c("bait","techr"), sep = "_", remove = F)

(umap.plot <- ggplot(df.umap.sample) +
    geom_point(aes(x = umap_x, y=umap_y, label = bait) , show.legend = TRUE) +
  geom_text(aes(x = umap_x, y=umap_y, label = bait) , show.legend = TRUE) +
  theme_classic()+theme(
    axis.line.x = element_line(colour = 'black', size=1.5, linetype='solid'), axis.ticks.length=unit(.5, "cm"), axis.ticks = element_line(size = 1.5),
    axis.line.y = element_line(colour = 'black', size=1.5, linetype='solid'), text = element_text(size=15, face="bold")))

ggsave(filename = ".\\res\\umap.plot_imp.pdf",umap.plot,
       device = "pdf", width=11, height=6, units = "in", useDingbats=F)
rm(umap.sample, df.umap.sample, umap.plot, umap.defaults2)

df.res <- df.imp %>% group_by(Majority.protein.IDs) %>% do({
  temp <- as.data.frame(.)
  temp$d = mean(unlist(temp %>% dplyr::select("45S_1":"45S_4"))) - mean(unlist(temp %>% dplyr::select("Efl_1":"NORAD_4")))
  temp$p = t.test(unlist(temp %>% dplyr::select("45S_1":"45S_4")), unlist(temp %>% dplyr::select("Efl_1":"NORAD_4")), var.equal = T)$p.value
  as_tibble(temp)
}) %>% ungroup() %>% as.data.frame() %>% dplyr::mutate(p.adj = p.adjust(p, method = "fdr"))

(tp <- ggplot() +
  geom_point(data=df.res, aes(x=d, y=-log10(p)), color = "grey") +
  geom_point(data=df.res %>% dplyr::filter(is.nol), aes(x=d, y=-log10(p)), color = "black") +
  geom_text_repel(data=df.res %>% dplyr::filter(is.nol & p.adj < 0.01 & d > 1), aes(x=d, y=-log10(p), label = Gene.names), color = "black") +
    geom_point(data=df.res %>% dplyr::filter(p.adj < 0.01 & d > 1), aes(x=d, y=-log10(p)), color = "red") +
  theme_classic())

ggsave(filename = ".\\res\\volcano_labels.pdf",tp,
       device = "pdf", width=12, height=6, units = "in", useDingbats=F)

df.res <- df.res %>% mutate(sig = (p.adj < 0.05 & d > 1))

table(df.res$sig)

fwrite(df.res, file=".\\res\\df.res.txt", quote = F, col.names = T, row.names = F, sep="\t")


