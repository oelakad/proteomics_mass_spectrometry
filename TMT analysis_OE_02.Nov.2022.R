#------Loading packages------------
#Loading packages
pacman::p_load(sva, edgeR, limma, ggplot2, factoextra, readr, rmarkdown)

#------group description----
#group description

groups<- as.vector("1\n127N\nKOPN-8\n2\n127C\nRCH-ACV\n3\n128N\nREH\n4\n128C\nRS4;11\n5\n129N\nSD-1\n6\n129C\nSEM\n7\n130N\n697\n8\n131\nInternal standard")
groups<- str_split(groups, "\n", n=24)
groups<- as.data.frame(groups)
colnames(groups)<- "V1"
groups<- as.data.frame(matrix(groups$V1, ncol = 3, byrow = TRUE))
colnames(groups)<- c("Number", "Labeling", "Cell line")
groups

#------Loading data----
#Loading data
TMTraw <- read_tsv(file.path("F:/(Frankfurt) Workspace/Analyses/LC-MSMS/Raw data/proteinGroups_TMT_Oct.22.txt")) 

TMT<- TMTraw %>% as_tibble() %>%
  filter(is.na(`Only identified by site`))%>% 
  filter(is.na(Reverse))%>% 
  filter(is.na(`Potential contaminant`)) %>% subset(select=c("Protein IDs",
                                                             "Reporter intensity 1 R1_996",
                                                             "Reporter intensity 2 R1_996",
                                                             "Reporter intensity 3 R1_996",
                                                             "Reporter intensity 4 R1_996",
                                                             "Reporter intensity 5 R1_996",
                                                             "Reporter intensity 6 R1_996",
                                                             "Reporter intensity 7 R1_996",
                                                             "Reporter intensity 8 R1_996",
                                                             "Reporter intensity 1 R2_996",
                                                             "Reporter intensity 2 R2_996",
                                                             "Reporter intensity 3 R2_996",
                                                             "Reporter intensity 4 R2_996",
                                                             "Reporter intensity 5 R2_996",
                                                             "Reporter intensity 6 R2_996",
                                                             "Reporter intensity 7 R2_996",
                                                             "Reporter intensity 8 R2_996"))


names(TMT) <- gsub(x=names(TMT), pattern = "Reporter intensity ", replacement = "")
names(TMT) <- gsub(x=names(TMT), pattern = "_996", replacement = "")
TMT<- TMT %>% column_to_rownames("Protein IDs")

names(TMT)
head(TMT,3)

boxplot(log2(TMT), col = rep(c('red', 'blue'), each = 8), 
        notch = TRUE, main = "RAW data: Exp1 (red), Exp2 (blue)",
        xlab = 'TMT Samples', ylab = 'log2 of Intensity')

plotDensities(log2(TMT), col = rep(c('red', 'blue'), 8), 
              main = 'Raw data')


plotMDS(log2(TMT), col = rep(c("red", "blue"), each = 8), 
        main = "RAW clusters group by TMT experiment")


fviz_dist(get_dist(t(na.omit(TMT)), method = "pearson"), lab_size = 13)

fviz_pca_ind(PCA(t(na.omit(TMT))))

#------SL----
#SL

exp1_raw <- TMT[1:8]
exp2_raw <- TMT[9:16]

# first basic normalization is to adjust each TMT experiment to equal signal per channel
# figure out the global scaling value
target <- mean(c(colSums(exp1_raw), colSums(exp2_raw)))

# do the sample loading normalization before the IRS normalization
# there is a different correction factor for each column
norm_facs1 <- target / colSums(exp1_raw)
exp1_sl <- sweep(exp1_raw, 2, norm_facs1, FUN = "*")
norm_facs2 <- target / colSums(exp2_raw)
exp2_sl <- sweep(exp2_raw, 2, norm_facs2, FUN = "*")


# make a pre-IRS data frame after sample loading normalizations
data_sl <- cbind(exp1_sl, exp2_sl)

# see what the SL normalized data look like
boxplot(log2(data_sl), col = rep(c("red", "blue"), each = 8), 
        notch = TRUE, main = "Sample Loading (SL) normalized data: \nExp1 (red), Exp2 (blue)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity')

plotDensities(log2(data_sl), col = rep(c("red", "blue"), 8), main = "SL normalization")

plotMDS(log2(data_sl), col = rep(c("red", "blue"), each = 8), 
        main = "SL clusters group by TMT experiment")

fviz_dist(get_dist(t(na.omit(data_sl)), method = "pearson"), lab_size = 13)

fviz_pca_ind(PCA(t(data_sl)))

head(data_sl, 3)


#------IRS----
#IRS

# make new data frame with row sums from each frame
irs <- tibble(rowSums(exp1_sl), rowSums(exp2_sl))
colnames(irs) <- c("sum1", "sum2")

# here I have used the internal control column instead of calculating average!!
irs$average1 <- exp1_raw$`8 R1`
irs$average2 <- exp2_raw$`8 R2`

# compute the scaling factor vectors
irs$fac1 <- irs$average1 / irs$sum1
irs$fac2 <- irs$average2 / irs$sum2


# make new data frame with IRS normalized data
data_irs <- exp1_sl * irs$fac1
data_irs <- cbind(data_irs, exp2_sl * irs$fac2)

# see what the IRS data look like
boxplot(log2(data_irs), col = rep(c("red", "blue"), each = 8), 
        main = "Internal Reference Scaling (IRS) normalized data: \nExp1 (red), Exp2 (green), Exp3 (blue)",
        xlab = 'IRS of SL Sample', ylab = 'log2 of Intensity', notch = TRUE)

# can also look at density plots (like a distribution histogram)    
plotDensities(log2(data_irs), col = rep(c("red", "blue"), ), main = "IRS of SL data")


plotMDS(log2(data_irs), col = rep(c("red", "blue"), each = 8), 
        main = "IRS of SL clusters group by TMT experiment")

fviz_dist(get_dist(t(na.omit(data_irs)), method = "pearson"), lab_size = 13)

fviz_pca_ind(PCA(t(na.omit(data_irs))))

head(data_irs, 3)


#------TMM----
#TMM

# perform TMM on the SL-normed data and visualize resulting distributions
data_irs<- na.omit(data_irs)
sl_tmm <- calcNormFactors(data_irs)
data_sl_tmm <- sweep(data_irs, 2, sl_tmm, FUN = "/") # data after SL and TMM on original scale

boxplot(log2(data_sl_tmm), col = rep(c("red", "blue"), each = 8), 
        notch = TRUE, main = "SL/TMM: \nExp1 (red), Exp2 (green), Exp3 (blue)",
        xlab = 'TMT Sample', ylab = 'log2 of Intensity')

plotDensities(log2(data_sl_tmm), col = rep(c("red", "blue"), 8), main = "SL/TMM normalization")


plotMDS(log2(data_sl_tmm), col = rep(c("red", "blue"), each = 8), 
        main = "SL/TMM clusters group by TMT experiment")

fviz_dist(get_dist(t(na.omit(data_sl_tmm)), method = "pearson"), lab_size = 13)

fviz_pca_ind(PCA(t(na.omit(data_sl_tmm))))

head(data_sl_tmm, 3)
