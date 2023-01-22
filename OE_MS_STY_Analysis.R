#-------Installing packages-------------------
library(ggplot2, datapasta, tidyverse, dplyr, ggrepel, cowplot, ggVennDiagram, readr,
        grid, ggpubr, DT, msigdbr, org.Hs.eg.db, clusterProfiler, enrichplot, 
        checkGeneSymbols, DOSE, ggnewscale, ggpattern, cowplot)

#-------Raw data analysis and wrangling for SYK-KO-MV411-Nuclear pYome----------------------------------
#-------Raw data analysis and wrangling for SYK-KO-MV411-Nuclear pYome

MV411KOYraw <- read_tsv(file.path("F:/(Frankfurt) Workspace/Analyses/LC-MSMS/Raw data/Phospho (STY)Sites_834.txt"))

#-------Raw data analysis and wrangling for SYK-KO-MV411-----------------------------------
#-------Raw data analysis and wrangling for SYK-KO-MV411

MV411KOraw <- read_tsv(file.path("F:/(Frankfurt) Workspace/Analyses/LC-MSMS/Raw data/proteinGroups_832.txt"))

#-------Raw data analysis and wrangling for SEM-Nuclear pYome----------------------------------
#-------Raw data analysis and wrangling for SEM-Nuclear pYome

SEMYraw <- read_tsv(file.path("F:/(Frankfurt) Workspace/Analyses/LC-MSMS/Raw data/Phospho (STY)Sites_1476.txt"))

SEMY0<- SEMYraw[,c("Ratio H/L normalized R1_1476___1",
                   "Ratio H/L normalized R1_1476___2",
                   "Ratio H/L normalized R1_1476___3",
                   "Ratio H/L normalized R2_1476___1",
                   "Ratio H/L normalized R2_1476___2",
                   "Ratio H/L normalized R2_1476___3",
                   "Gene names", "Protein names", "Intensity",
                   "Reverse", "Potential contaminant", "Localization prob",
                   "id", "Peptide IDs", "Protein")]

SEMY0<- SEMY0 %>% dplyr::filter(is.na(SEMY0$"Potential contaminant"))# remove contaminants
#SEMY0<- SEMY0 %>% dplyr::filter(is.na(SEMY0$"Only identified by site"))# remove identified only by sites
SEMY0<- SEMY0 %>% dplyr::filter(is.na(SEMY0$"Reverse"))# remove reverse
SEMY0<- subset(SEMY0, select= -c(Reverse,`Potential contaminant`))

#check the unique columns
cat(paste("duplicated Peptide IDs:", "", anyDuplicated(SEMY0$`Peptide IDs`), 
          "duplicated Protein names:", "",anyDuplicated(SEMY0$`Protein names`), 
          "duplicated Gene names:", "",anyDuplicated(SEMY0$`Gene names`),
          "duplicated ids:", "",anyDuplicated(SEMY0$"id"),
          sep = "\n")) 

#----combining the three multiplicities into one column----in the STY data, every read has 3 multiplicities (___1,2,3), 
#I wanted to pivot longer R1 separately from R2, so I separated them into two data frames then merged them again using a unique column id that I created, 
#this code looks very naive and I would be happy if you could shorten it

SEMYR1<- pivot_longer(data=SEMY0, cols=c(`Ratio H/L normalized R1_1476___1`,
                                       `Ratio H/L normalized R1_1476___2`,
                                       `Ratio H/L normalized R1_1476___3`),
                     values_to = "R1", names_to = "MultiplicityR1")
SEMYR1$"MultiplicityR1"<-gsub("Ratio H/L normalized R1_1476__","",SEMYR1$"MultiplicityR1")
SEMYR1$Gene_Multiplicity <- paste(SEMYR1$`Gene names`, SEMYR1$MultiplicityR1, SEMYR1$id)
cat(paste("duplicated R1 Multiplicities:", "", anyDuplicated(SEMYR1$Gene_Multiplicity)))


SEMYR2<- pivot_longer(data=SEMY0, cols=c(`Ratio H/L normalized R2_1476___1`,
                                         `Ratio H/L normalized R2_1476___2`,
                                         `Ratio H/L normalized R2_1476___3`),
                      values_to = "R2", names_to = "MultiplicityR2")
SEMYR2$"MultiplicityR2"<-gsub("Ratio H/L normalized R2_1476__","",SEMYR2$"MultiplicityR2")
SEMYR2$Gene_Multiplicity <- paste(SEMYR2$`Gene names`, SEMYR2$MultiplicityR2, SEMYR2$id)
cat(paste("duplicated R2 Multiplicities:", "", anyDuplicated(SEMYR2$Gene_Multiplicity)))

SEMY<- merge(SEMYR1, SEMYR2, on="Gene_Multiplicity")

cat(paste("duplicated Peptide IDs:", "", anyDuplicated(SEMY$`Peptide IDs`), 
          "duplicated Protein names:", "",anyDuplicated(SEMY$`Protein names`), 
          "duplicated Gene names:", "",anyDuplicated(SEMY$`Gene names`),
          "duplicated ids:", "",anyDuplicated(SEMY$"id"),
          sep = "\n")) 



SEMY[SEMY=="NaN"] <- NA#change NANs to Nas 
SEMY$"Gene names"<- sub(";.*", "\\1", SEMY$"Gene names")#take the first alias of gene name till semicolon
SEMY$"R1"<- as.numeric(SEMY$"R1")# as they were not numeric
SEMY$"R2"<- as.numeric(SEMY$"R2")# as they were not numeric
SEMY$"Intensity"<- as.numeric(SEMY$"Intensity")# as they were not numeric

SEMY<- SEMY %>% dplyr::filter(SEMY$"Localization prob">0.75)#filter for first class phosphosites
SEMY <- SEMY %>% mutate(R1=log2(R1), R2=log2(R2))
SEMY$"Intensity"<- log10(SEMY$"Intensity")# log10x

#SEMY$`Protein names` <- ifelse(is.na(SEMY$`Protein names`), SEMY$`Gene names`, SEMY$`Protein names`)#change all NA protein name to Gene name
SEMY$`Gene names` <- ifelse(is.na(SEMY$`Gene names`), "CALM3", SEMY$`Gene names`)#I looked for the name of this specific protein "P0DP25" manually,
#I would be happy if there is an automated/systemic way

SEMY$meansem <- rowMeans(subset(SEMY, select = c("R1", "R2")), na.rm = TRUE)
SEMY<- SEMY[!is.na(SEMY$meansem),]# remove invalid values

#-------Raw data analysis and wrangling for MV4-11-Nuclear pYome-------------------
#-------Raw data analysis and wrangling for MV4-11-Nuclear pYome
MV4Yraw <- read_tsv(file.path("F:/(Frankfurt) Workspace/Analyses/LC-MSMS/Raw data/Phospho (STY)Sites_19-653.txt"))


MV4Y0<- MV4Yraw[,c("Ratio H/L normalized R1_653___1",
                   "Ratio H/L normalized R1_653___2",
                   "Ratio H/L normalized R1_653___3",
                   "Ratio H/L normalized R2_653___1",
                   "Ratio H/L normalized R2_653___2",
                   "Ratio H/L normalized R2_653___3",
                   "Gene names", "Protein names", "Intensity",
                   "Reverse", "Potential contaminant", "Localization prob",
                   "id", "Peptide IDs", "Protein")]

MV4Y0<- MV4Y0 %>% dplyr::filter(is.na(MV4Y0$"Potential contaminant"))# remove contaminants
#MV4Y0<- MV4Y0 %>% dplyr::filter(is.na(MV4Y0$"Only identified by site"))# remove identified only by sites
MV4Y0<- MV4Y0 %>% dplyr::filter(is.na(MV4Y0$"Reverse"))# remove reverse
MV4Y0<- subset(MV4Y0, select= -c(Reverse,`Potential contaminant`))

#check the unique columns
cat(paste("duplicated Peptide IDs:", "", anyDuplicated(MV4Y0$`Peptide IDs`), 
          "duplicated Protein names:", "",anyDuplicated(MV4Y0$`Protein names`), 
          "duplicated Gene names:", "",anyDuplicated(MV4Y0$`Gene names`),
          "duplicated ids:", "",anyDuplicated(MV4Y0$"id"),
          sep = "\n")) 

#----combining the three multiplicities into one column----in the STY data, every read has 3 multiplicities (___1,2,3), 
#I wanted to pivot longer R1 separately from R2, so I separated them into two data frames then merged them again using a unique column id that I created, 
#this code looks very naive and I would be happy if you could shorten it

MV4YR1<- pivot_longer(data=MV4Y0, cols=c(`Ratio H/L normalized R1_653___1`,
                                         `Ratio H/L normalized R1_653___2`,
                                         `Ratio H/L normalized R1_653___3`),
                      values_to = "R1", names_to = "MultiplicityR1")
MV4YR1$"MultiplicityR1"<-gsub("Ratio H/L normalized R1_653__","",MV4YR1$"MultiplicityR1")
MV4YR1$Gene_Multiplicity <- paste(MV4YR1$`Gene names`, MV4YR1$MultiplicityR1, MV4YR1$id)
cat(paste("duplicated R1 Multiplicities:", "", anyDuplicated(MV4YR1$Gene_Multiplicity)))


MV4YR2<- pivot_longer(data=MV4Y0, cols=c(`Ratio H/L normalized R2_653___1`,
                                         `Ratio H/L normalized R2_653___2`,
                                         `Ratio H/L normalized R2_653___3`),
                      values_to = "R2", names_to = "MultiplicityR2")
MV4YR2$"MultiplicityR2"<-gsub("Ratio H/L normalized R2_653__","",MV4YR2$"MultiplicityR2")
MV4YR2$Gene_Multiplicity <- paste(MV4YR2$`Gene names`, MV4YR2$MultiplicityR2, MV4YR2$id)
cat(paste("duplicated R2 Multiplicities:", "", anyDuplicated(MV4YR2$Gene_Multiplicity)))

MV4Y<- merge(MV4YR1, MV4YR2, on="Gene_Multiplicity")


MV4Y[MV4Y=="NaN"] <- NA#change NANs to Nas 
MV4Y$"Gene names"<- sub(";.*", "\\1", MV4Y$"Gene names")#take the first alias of gene name till semicolon
MV4Y$"R1"<- as.numeric(MV4Y$"R1")# as they were not numeric
MV4Y$"R2"<- as.numeric(MV4Y$"R2")# as they were not numeric
MV4Y$"Intensity"<- as.numeric(MV4Y$"Intensity")# as they were not numeric

MV4Y<- MV4Y %>% dplyr::filter(MV4Y$"Localization prob">0.75)#filter for first class phosphosites
MV4Y <- MV4Y %>% mutate(R1=log2(R1), R2=log2(R2))
MV4Y$"Intensity"<- log10(MV4Y$"Intensity")# log10x

#SEMY$`Protein names` <- ifelse(is.na(SEMY$`Protein names`), SEMY$`Gene names`, SEMY$`Protein names`)#change all NA protein name to Gene name
MV4Y$`Gene names` <- ifelse(is.na(MV4Y$`Gene names`), "CALM3", MV4Y$`Gene names`)#I looked for the name of this specific protein "P0DP25" manually,
#I would be happy if there is an automated/systemic way

MV4Y$meanmv <- rowMeans(subset(MV4Y, select = c("R1", "R2")), na.rm = TRUE)
MV4Y<- MV4Y[!is.na(MV4Y$meanmv),]# remove invalid values



#-------Histograms------------
#-------Histograms


h1<- ggplot(MV4Y, aes(x=`meanmv`))  +
  geom_histogram(aes(y = ..density..), color = "grey30", fill = "white", bins=30) +
  geom_density(alpha = .2, fill = "cornflowerblue")+ theme_bw()+
  labs(y="Density",title="Histogram", x="Counts")

h2<- ggplot(SEMY, aes(x=`meansem`))  +
  geom_histogram(aes(y = ..density..), color = "grey30", fill = "white", bins=30) +
  geom_density(alpha = .2, fill = "cornflowerblue")+ theme_bw()+
  labs(y="Density",title="Histogram", x="Counts")

plot_grid(h1, h2)

#-------Main Plots-------------
#-------Main Plots

SEMY$Enrichment <- "No"
SEMY$Enrichment[SEMY$meansem> 0.5 ] <- "Higher"
SEMY$Enrichment[SEMY$meansem< -0.5 ] <- "Lower"
SEMY$Enrichment<- as.factor(SEMY$Enrichment)

MV4Y$Enrichment <- "No"
MV4Y$Enrichment[MV4Y$meanmv> 0.5 ] <- "Higher"
MV4Y$Enrichment[MV4Y$meanmv< -0.5 ] <- "Lower"
MV4Y$Enrichment<- as.factor(MV4Y$Enrichment)

sig_il_genes_SEM <- SEMY %>% dplyr::filter(`Gene names` %in% c("SYK", "MEF2C", "MEF2D", "MEIS1"))
sig_il_genes_SEM <- SEMY %>% dplyr::filter(`Gene names` %in% groups$targetsSYMBOL)
sig_il_genes_SEM <- dplyr::filter(dplyr::filter(SEMY, `Gene names` %in% groups$targetsSYMBOL), meansem< -0.5)


sig_il_genes_MV411 <- MV4Y %>% dplyr::filter(`Gene names` %in% c("SYK", "MEF2C", "MEF2D", "MEIS1"))
sig_il_genes_MV411 <- MV4Y %>% dplyr::filter(`Gene names` %in% groups$partnersSYMBOL)
sig_il_genes_MV411 <- dplyr::filter(dplyr::filter(MV4Y, `Gene names` %in% groups$partnersSYMBOL), meanmv< -0.5)

p1<- ggplot(SEMY, aes(x=`meansem`, y=Intensity, color=Enrichment)) + 
  geom_point(shape = 21, fill = "black", size=3, alpha = 9/10)+
  geom_point(data=sig_il_genes_SEM, aes(x=`meansem`, y=Intensity), color='red', size=3, alpha = 9/10)+
  theme_bw()+
  geom_vline(xintercept = c(-0.5,0.5), color="red")+
  scale_x_continuous(limits = c(-7.9, 7.9))+
  #scale_fill_manual(values = cols) + 
  #scale_size_manual(values = sizes) +
  labs(x="Log2 Normalized Ratio Ento/DMSO", y="-Log10Intensity", title="Nuclear-pYome in SEM")+
  geom_label_repel(data = sig_il_genes_SEM, aes(label=`Gene names`), box.padding = 0.5,
                   max.overlaps = Inf, color="black")+
  theme(legend.position="none") #+ theme(legend.position=c(0.2, 0.8))

p2<- ggplot(MV4Y, aes(x=`meanmv`, y=Intensity, color=Enrichment)) + 
  geom_point(shape = 21, fill = "black", size=3, alpha = 9/10)+
  geom_point(data=sig_il_genes_MV411,
             aes(x=`meanmv`, y=Intensity), color='red', 
             size=3, alpha = 9/10)+theme_bw()+
  geom_vline(xintercept = c(-0.5,0.5), color="red")+
  scale_x_continuous(limits = c(-4.5, 4.5))+
  #scale_fill_manual(values = cols) + 
  #scale_size_manual(values = sizes) +
  labs(x="Log2 Normalized Ratio Ento/DMSO", y="-Log10Intensity", title="Nuclear-pYome in MV4-11")+
  geom_label_repel(data = sig_il_genes_MV411, aes(label=`Gene names`), box.padding = 0.5,
                   max.overlaps = Inf, color="black")+
  theme(legend.position="none") #+ theme(legend.position=c(0.2, 0.8))


p3<- ggplot(SEMY, aes(x=`meansem`, y=`Intensity`)) + 
  geom_point(shape = 16, size=3)+
  geom_point(data=sig_il_genes_SEM, aes(x=`meansem`, y=Intensity), color='red', size=3, alpha = 9/10)+
  theme_bw()+ 
  geom_vline(xintercept = c(-1), color="red")+
  geom_vline(xintercept = c(0), color="black", size=1.5)+
  labs(x="Log2 Normalized Ratio Ento/DMSO", y="-Log10Intensity", title="Nuclear-pYome in SEM")+ 
  scale_x_continuous(breaks = seq(0, -8, by = -1), limits = c(-8, 0))+
  scale_y_continuous(breaks = seq(5.7, 10.5, by = 0.3), limits= c(5.7,10.3))+
  #scale_fill_manual(values = cols) + 
  #scale_size_manual(values = sizes) +
  geom_text_repel(segment.colour = NA, label=SEMY$`Gene names`)

p4<- ggplot(MV4Y, aes(x=`meanmv`, y=`Intensity`)) + 
  geom_point(shape = 16, size=3)+
  geom_point(data=sig_il_genes_MV411,
             aes(x=`meanmv`, y=Intensity), color='red',size=3, alpha = 9/10)+
  theme_bw()+ 
  geom_vline(xintercept = c(-1), color="red")+
  geom_vline(xintercept = c(0), color="black", size=1.5)+
  labs(x="Log2 Normalized Ratio Ento/DMSO", y="-Log10Intensity", title="Nuclear-pYome in MV4-11")+ 
  scale_x_continuous(breaks = seq(0, -4.1, by = -1), limits = c(-4.1, 0))+
  scale_y_continuous(breaks = seq(5.7, 10.5, by = 0.3), limits= c(5.7,10.3))+
  #scale_fill_manual(values = cols) + 
  #scale_size_manual(values = sizes) +
  geom_text_repel(segment.colour = NA, label=MV4Y$`Gene names`)


plot_grid(h1, p1, p3, h2, p2, p4, ncol=3, rel_widths=c(1,2,3))


#-------Correlating the two data sets----------
#-------Correlating the two data sets

Merged<- inner_join(SEMY, MV4Y, by= "Gene names")
Merged$Enrichment <- "No"
Merged$Enrichment[Merged$meanmv < -1 & Merged$meansem < -1] <- "High"
Significant<- dplyr::filter(Merged, Merged$meanmv< -1 & Merged$meansem < -1)
spearmanNucpYome<- cor.test(Merged$meansem, Merged$meanmv, method="spearman")
 
c1<- ggplot(Merged, aes(meansem, meanmv)) + 
  geom_point(shape = 21, fill = "black", color="gray", size=3, alpha = 3/10)+
  geom_vline(xintercept = c(0), color="black")+
  geom_hline(yintercept = c(0), color="black")+
  geom_smooth(method='lm', se=T)+
  theme_bw()+
  #stat_regline_equation(label.y = 6, aes(label = ..rr.label..))+
  labs(x="Log2 FC SEM", y="Log2 FC MV4-11", title="Nuclear-pYome")+
  annotate("text", x = -5, y = 1.3, col = "black", size=4.3,
           label = paste("Spearman r = ", 
                         format(round(spearmanNucpYome$estimate, 2), nsmall = 2)))+
  geom_text_repel(data=Significant, segment.colour = "black",
                  label=Significant$`Gene names`,
                  max.overlaps = 10,
                  force = 3,
                  arrow = arrow(length = unit(0.015, "npc")),
                  nudge_x= 0)
c1
#-------Venn Diagram----------
#-------Venn Diagram

x<- list(SEMY$`Gene names`, MV4Y$`Gene names`)
v1<- ggVennDiagram(x[1:2],
                   label_alpha = 0,
                   set_size = 5,
                   label_color = "black", label_size = 4,
                   stroke_size = 0.4,
                   edge_size = 2,
                   category.names = c("SEM", "MV4-11"),
                   show_intersect = F,
                   color = "black") +
  scale_fill_gradient(low = "white", high = "lightblue")+
  labs(title = "Total Quantified Phosphosites")+
  theme(legend.position ="none")


v1

#------ adjust SYMBOL and add ENTREZID------
#------ adjust SYMBOL and add ENTREZID
#SEMY
# first correct symbols and bind them to df
corrected_symbols_SEM<- HGNChelper::checkGeneSymbols(x = SEMY$`Gene names`)
SEMY$SYMBOL <- corrected_symbols_SEM$Suggested.Symbol

# second transform symbols to ENTREZID and bind them to df
ENTREZIDdfSEM<- clusterProfiler::bitr(geneID = SEMY$SYMBOL, 
                      fromType = "SYMBOL", 
                      toType = "ENTREZID", 
                      OrgDb = org.Hs.eg.db, 
                      drop = F)

SEMY<- SEMY %>%  left_join(ENTREZIDdfSEM, by=c("SYMBOL" = "SYMBOL"))

#MV4Y
corrected_symbols_MV4<- HGNChelper::checkGeneSymbols(x = MV4Y$`Gene names`)
MV4Y$SYMBOL <- corrected_symbols_MV4$Suggested.Symbol

# second transform symbols to ENTREZID and bind them to df
ENTREZIDdfMV4<- clusterProfiler::bitr(geneID = MV4Y$SYMBOL, 
                                   fromType = "SYMBOL", 
                                   toType = "ENTREZID", 
                                   OrgDb = org.Hs.eg.db, 
                                   drop = F)

MV4Y<- MV4Y %>%  left_join(ENTREZIDdfMV4, by=c("SYMBOL" = "SYMBOL"))


#----------Setting up a Universe to have suffecient genes for Over Representation analysis---------
#----------Setting up a Universe to have suffecient genes for Over Representation analysis

Universe<- as_tibble(full_join(SEMY0, MV4Y0, by='Gene names'))[,"Gene names"]
paste("Number of genes in the Universe joining SEM Nuc pYome and MV4 Nuc pYome = ", nrow(Universe))
Universe<- full_join(Universe, MV411KOYraw, by='Gene names')[,"Gene names"]
paste("Number of genes in the Universe joining SEM Nuc pYome, MV4 Nuc pYome and MV4 KO pYome = ", nrow(Universe))
Universe<- full_join(Universe, MV411KOraw, by='Gene names')[,"Gene names"]
paste("Number of genes in the Universe joining SEM Nuc pYome, MV4 Nuc pYome, MV4 KO pYome and MV4 KO = ", nrow(Universe))


Universe<- dplyr::rename(Universe, gene="Gene names")
Universe$gene<-gsub(";.*", "", Universe$gene)
corrected_symbols_Universe<- HGNChelper::checkGeneSymbols(x = Universe$gene)
Universe$SYMBOL <- corrected_symbols_Universe$Suggested.Symbol
ENTREZIDUniverse<- clusterProfiler::bitr(geneID = Universe$SYMBOL, 
                                      fromType = "SYMBOL", 
                                      toType = "ENTREZID", 
                                      OrgDb = org.Hs.eg.db, 
                                      drop = F)
Universe<- Universe %>%  left_join(ENTREZIDUniverse, by=c("SYMBOL" = "SYMBOL"))
Universe<- Universe[!duplicated(Universe$SYMBOL),]
paste("Final number of genes in the Universe = ", nrow(Universe))

#----------Nuclear data sets----------



#-------ORA_GO---------
#-------ORA_GO

SEMY_ORA_GO <- enrichGO(gene = filter(SEMY, meansem <= -1)$SYMBOL,
                    OrgDb = org.Hs.eg.db, 
                    keyType = "SYMBOL",
                    pvalueCutoff = 0.01, 
                    qvalueCutoff  = 0.01,
                    pAdjustMethod = "BH",
                    minGSSize = 5,
                    maxGSSize = Inf, 
                    universe=Universe$SYMBOL, 
                    ont = "ALL",
                    readable = T)

paste("Number of Enriched Pathways =", nrow(SEMY_ORA_GO@result), "Pathways")


SEMY_ORA_GO@result$Description<- substr(SEMY_ORA_GO@result$Description, 1, 37)


SEMY.dotplot_all<- enrichplot::dotplot(SEMY_ORA_GO, showCategory=20)
SEMY.dotplot_bp<- enrichplot::dotplot(filter(SEMY_ORA_GO, ONTOLOGY == "BP"), showCategory  = 20)
SEMY.dotplot_mf<- enrichplot::dotplot(filter(SEMY_ORA_GO, ONTOLOGY == "MF"), showCategory  = 20)
SEMY.dotplot_cc<- enrichplot::dotplot(filter(SEMY_ORA_GO, ONTOLOGY == "CC"), showCategory  = 20)


cowplot::plot_grid(plotlist = list(
  "BP" = SEMY.dotplot_bp,
  "MF" = SEMY.dotplot_mf),
  labels = c("BP", "MF"),
  ncol = 2) 


enrichplot::dotplot(SEMY_ORA_GO, showCategory=20)
enrichplot::upsetplot(SEMY_ORA_GO, 15) #number of genes overlapping between pathways
graphics::barplot(SEMY_ORA_GO, showCategory = 15, title = "ALL")


enrichplot::emapplot(pairwise_termsim(SEMY_ORA_GO),showCategory = 20,
                     color = "p.adjust", layout = "kk", legend_n = 5,
                     group_category = T, nCluster = 2, cex_label_category=0.7)+
  ggtitle("Enrichment of two sets of networks of RNA processing and Hematopoiesis ")



enrichplot::cnetplot(SEMY_ORA_GO, showCategory = 10, categorySize="pvalue",
         foldChange = NULL,#gene_list
         circular=F,
         node_label_size=2,
         cex_label_category=1,
         cex_label_gene=0.5,
         cex_gene=0.5,
         colorEdge=T)

enrichplot::cnetplot(SEMY_ORA_GO, showCategory = 5,
                 foldChange = list(filter(SEMY, meansem<= -0.5)$SYMBOL),
                 circular=T,
                 node_label_size=2,
                 cex_label_category=1,
                 cex_label_gene=0.5,
                 cex_gene=0.5,
                 colorEdge=T)




#-------ORA_KEGG -----
#-------ORA_KEGG 



SEMY_KEGG <- enrichKEGG(gene = filter(SEMY, meansem <= -1)$ENTREZID,
                        organism = "hsa",
                        universe = Universe$ENTREZID,
                        pvalueCutoff = 0.01)
                        
paste("Number of Enriched Pathways =", nrow(SEMY_KEGG@result), "Pathways")

SEMY_KEGG@result[1:20,] %>%
  mutate(Description=substr(Description, 1, 38)) %>%
  ggplot(
    aes(-log(p.adjust), reorder(Description, -log(p.adjust)),
        size=Count,
        color=qvalue))+
  geom_point(alpha=0.8)+
  scale_colour_gradient(low = "red2",  high = "mediumblue") +
  labs(title="KEGG enrichment analysis", y="Pathways")+
  theme_bw()



#-------ORA_Reactome------
#-------ORA_Reactome

SEMY_reactome<- ReactomePA::enrichPathway(gene = dplyr::filter(SEMY, meansem <= -1)$ENTREZID, 
                          pvalueCutoff = 0.1, 
                          pAdjustMethod = "BH",
                          minGSSize = 5,
                          maxGSSize = Inf,
                          universe = Universe$ENTREZID
                          )

paste("Number of Enriched Pathways =", nrow(SEMY_reactome@result), "Pathways")



SEMY_reactome@result$Description<- substr(SEMY_reactome@result$Description, 1, 50)

SEMY_reactome %>% 
  as_tibble() %>% 
  head(20) %>% 
  mutate(Description=substr(Description, 1, 38)) %>%
ggplot(
    aes(-log(p.adjust), reorder(Description, -log(p.adjust)),
        size=Count,
        color=p.adjust))+
  geom_point(alpha=0.8)+
  scale_colour_gradient(low = "red2",  high = "mediumblue") +
  labs(title="Reactome enrichment analysis", y="Pathways")+
  theme_bw()

enrichplot::dotplot(SEMY_reactome)
enrichplot::upsetplot(SEMY_reactome, 15)
graphics::barplot(SEMY_reactome, showCategory = 15, title = "ALL")
enrichplot::emapplot(pairwise_termsim(SEMY_reactome),showCategory = 15,color = "p.adjust", layout = "kk")
enrichplot::cnetplot(SEMY_reactome, showCategory = 5, categorySize="pvalue",
         foldChange = NULL,#gene_list
         circular=F,
         node_label_size=2,
         cex_label_category=0.7,
         cex_label_gene=0.5,
         cex_gene=0.5,
         colorEdge=T)
enrichplot::cnetplot(SEMY_reactome, showCategory = 5,
         foldChange = list(filter(MV4Y, meanmv<= -0.5)$SYMBOL),
         circular=T,
         node_label_size=2,
         cex_label_category=1,
         cex_label_gene=0.5,
         cex_gene=0.5,
         colorEdge=T)

#-------GSEA_GO--------
#-------GSEA_GO


SEMYList<- SEMY %>%
  dplyr::select(SYMBOL, meansem) %>% 
  mutate(meansem = meansem * -1) %>% #*-1 as we are interested in the down regulated sites
  deframe() 

SEMY_GSEA_GO <- clusterProfiler::gseGO(geneList = sort(SEMYList, decreasing = T), 
                       ont = "ALL",
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       minGSSize = 20, 
                       maxGSSize = 1000, 
                       pvalueCutoff = 0.1, 
                       pAdjustMethod = "BH", 
                       by = "fgsea")

paste("Number of Enriched Pathways =", nrow(SEMY_GSEA_GO@result), "Pathways")

SEMY_GSEA_GO@result$Description<- substr(SEMY_GSEA_GO@result$Description, 1, 50)




enrichplot::dotplot(SEMY_GSEA_GO, showCategory=20)
enrichplot::upsetplot(SEMY_GSEA_GO, 15) #number of genes overlapping between pathways
enrichplot::emapplot(pairwise_termsim(SEMY_GSEA_GO),showCategory = 15,color = "p.adjust", layout = "kk")
enrichplot::cnetplot(SEMY_GSEA_GO, showCategory = 2, categorySize="pvalue",
                     foldChange = NULL,#gene_list
                     circular=F,
                     node_label_size=2,
                     cex_label_category=1,
                     cex_label_gene=0.5,
                     cex_gene=0.5,
                     colorEdge=T)
enrichplot::cnetplot(SEMY_GSEA_GO, showCategory = 2,
                     foldChange =SEMYList,
                     circular=F,
                     node_label_size=2,
                     cex_label_category=1,
                     cex_label_gene=0.5,
                     cex_gene=0.5,
                     colorEdge=T)


SEMY_GSEA_GO %>%
  as_tibble() %>%
  filter(ONTOLOGY == c("CC", "BP", "MF")) %>% 
  arrange(desc(NES)) %>%
  head(20) %>%
  mutate(Description=substr(Description, 1, 50)) %>%
  ggplot(
    aes(NES, reorder(Description, NES),
        size=setSize,
        color=NES))+
  geom_point(alpha=0.8)+
  scale_colour_gradient(low = "red2",  high = "mediumblue") +
  labs(title="GSEA enrichment analysis-Upregulated", y="")+
  theme_bw()

SEMY_GSEA_GO %>%
  as_tibble() %>%
  filter(ONTOLOGY == c("CC", "BP", "MF")) %>% 
  arrange(NES) %>%
  head(20) %>%
  mutate(Description=substr(Description, 1, 50)) %>%
  ggplot(
    aes(NES, reorder(Description, NES),
        size=setSize,
        color=NES))+
  geom_point(alpha=0.8)+
  scale_colour_gradient(low = "red2",  high = "mediumblue") +
  labs(title="GSEA enrichment analysis-Downregulated", y="")+
  theme_bw()

#-------GSEA_Reactome_Pathway------------
#-------GSEA_Reactome_Pathway


SEMYListreactome<- SEMY %>%
  dplyr::select(ENTREZID, meansem) %>% 
  mutate(meansem = meansem * -1) %>% #*-1 as we are interested in the downregulated sites
  deframe() 

SEMY_GSEA_Reactome <- ReactomePA::gsePathway(geneList = sort(SEMYListreactome, decreasing = T), 
                       minGSSize = 20, 
                       maxGSSize = 1000, 
                       pvalueCutoff = 0.1, 
                       pAdjustMethod = "BH")

paste("Number of Enriched Pathways =", nrow(SEMY_GSEA_Reactome@result), "Pathways")



SEMY_GSEA_Reactome %>%
  as_tibble() %>%
  arrange(desc(NES)) %>%
  head(20) %>%
  mutate(Description=substr(Description, 1, 38)) %>%
  ggplot(
    aes(-log(p.adjust), reorder(Description, -log(p.adjust)),
        size=setSize,
        color=NES))+
  geom_point(alpha=0.8)+
  scale_colour_gradient(low = "red2",  high = "mediumblue") +
  labs(title="GSEA enrichment analysis-Upregulated", y="")+
  theme_bw()

SEMY_GSEA_Reactome %>%
  as_tibble() %>%
  arrange(NES) %>%
  head(20) %>%
  mutate(Description=substr(Description, 1, 38)) %>%
  ggplot(
    aes(-log(p.adjust), reorder(Description, -log(p.adjust)),
        size=setSize,
        color=NES))+
  geom_point(alpha=0.8)+
  scale_colour_gradient(low = "red2",  high = "mediumblue") +
  labs(title="GSEA enrichment analysis-Downregulated", y="")+
  theme_bw()

#-------GSEA_MITOCARTA--(I used it here as a negative control for nuclear proteins)----
#-------GSEA_MITOCARTA--(I used it here as a negative control for nuclear proteins)

mitocarta<- read_excel("C:/Downloads/Human.MitoCarta3.0.xls", sheet=4)%>% 
  #first we create a list where every element is a pathway with its genes in an character vector
  dplyr::select(MitoPathway, Genes) %>% 
  {
    split(.$Genes, .$MitoPathway)
  } %>% 
  #next we use purrrs map (similar to lapply()) to split every character vector into individual genes
  purrr::map(~ stringr::str_split(string = ., pattern = ",")[[1]]) %>% 
  # and remove the whitespaces 
  purrr::map(~ stringr::str_trim(string = ., side = "both")) %>%
  # unlist into vector (cave: numbers the names)
  unlist(use.names = T) %>% 
  # transform into tibble
  as_tibble(rownames = "pathway") %>%
  # and remove the numbering
  mutate(pathway = stringr::str_remove(string = pathway, "\\d{1,}$"))



SEMYMitocarta<- clusterProfiler::GSEA(geneList = sort(SEMYList, decreasing = T),
                                        pvalueCutoff = 0.1, 
                                        pAdjustMethod = "BH", 
                                        TERM2GENE = mitocarta)

paste("Number of Enriched Pathways =", nrow(SEMYMitocarta@result), "Pathways")

#-------Leading edges (Automated)--------
#-------Leading edges (Automated)


SEMY_ORA_GO@result$Description<- substr(SEMY_ORA_GO@result$Description, 1, 35)

#Step 1 (create termsim)
SEMYGSEAmitocartatermism<- SEMY_ORA_GO %>% 
  enrichplot::pairwise_termsim()



#Step 2 (plot tree to check clusters)
SEMYGSEAmitocartatermism@termsim %>% 
  t() %>% 
  as.dist() %>% 
  {
    hclust(1- .,method = "ward.D2") 
  } %>% 
  plot()


#Step 3 (Cut tree)
SEMYGSEAmitocartatermism_cluster<- SEMYGSEAmitocartatermism@termsim %>% 
  t() %>% 
  as.dist() %>% 
  {
    hclust(1- .,method = "ward.D2") 
  } %>% 
  cutree(k = 2)

#Step4 make a table with Pathways, Gene names and Clusters through combination
SEMYGSEAmitocartatermism_cluster <- data.frame(as.list(SEMYGSEAmitocartatermism_cluster))
SEMYGSEAmitocartatermism_cluster <- pivot_longer(SEMYGSEAmitocartatermism_cluster, cols = everything())
SEMYGSEAmitocartatermism_cluster<- dplyr::rename(SEMYGSEAmitocartatermism_cluster, "Description" = "name")
SEMYGSEAmitocartatermism_cluster<- dplyr::rename(SEMYGSEAmitocartatermism_cluster, "cluster" = "value")
SEMYGSEAmitocartatermism_cluster$Description<-gsub(".", " ", SEMYGSEAmitocartatermism_cluster$Description, fixed=TRUE)

SEMYGSEAmitocartatermism_cluster.nested <- SEMY_ORA_GO %>% 
  as_tibble() %>% 
  left_join(SEMYGSEAmitocartatermism_cluster, by="Description") %>% 
  dplyr::select(ID, Description, geneID, cluster) %>% 
  group_by(cluster) %>% 
  mutate(cluster_name = case_when(
    cluster == 1 ~ "Hematopoiesis",
    cluster == 2 ~ "RNAprocessing"
  )) %>%
  separate_rows(geneID, convert = TRUE) %>% 
  mutate(value= geneID) %>%
  tidyr::pivot_wider(names_from = "cluster_name", values_from = "geneID") %>% 
  rename(gene = value) %>%
  mutate(RNAprocessing= ifelse(is.na(RNAprocessing) == T, "No", "Yes"))%>%
  mutate(Hematopoiesis= ifelse(is.na(Hematopoiesis) == T, "No", "Yes"))

SEMYGSEAmitocartatermism_cluster.nested<- SEMYGSEAmitocartatermism_cluster.nested[!duplicated(SEMYGSEAmitocartatermism_cluster.nested$gene), ]

SEMYGSEAmitocartatermism_cluster.nested[,-c(1,3,7)] %>%
  rename_with(.cols = -gene, .fn = ~ stringr::str_trunc(., width = 30, side = "right")) %>%
  DT::datatable(filter = "top", 
                options = list(autoWidth = F, scrollX = T), 
                class = "compact hover row-border stripe dt-left cell-border nowrap")
