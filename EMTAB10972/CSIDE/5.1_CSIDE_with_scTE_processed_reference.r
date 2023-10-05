library(spacexr)
library(Matrix)
library(doParallel)
library(readr)
library(ggplot2)

setwd('/work/project/ext_014/EMTAB10972/')

#Customize dataframes and create reference
#Counts
Sys.setenv(VROOM_CONNECTION_SIZE=2000000)
counts <- as.data.frame(read_csv('CSIDE/Reference/sc_counts.csv.gz'))
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames

#cell_types
ct_source <- read_csv('rawdata/EMTAB10974_annotation.csv.gz')
cell_types <- setNames(ct_source[[6]], ct_source[[1]])
cell_types <- as.factor(cell_types) # convert to factor data type

#nUMI
col_sum <- colSums(counts)
meta_data <- read.delim("rawdata/meta_org1-12.tsv.gz")
meta_data <- meta_data[,11,drop = FALSE]
meta_data[,2] <- meta_data[,1]
meta_data[,1] <- rownames(meta_data)
rownames(meta_data) <- NULL
meta_data[,2] <- col_sum
colnames(meta_data) <- c('barcode','nUMI')
nUMI <- setNames(meta_data[[2]], meta_data[[1]])
reference <- Reference(counts, cell_types, nUMI)

setwd('/work/project/ext_014/EMTAB10972/CSIDE/')

#Generated from h5ad_conversion.ipynb, spatial data preparation
counts <- as.data.frame(read_csv('S1/counts.csv.gz'))
coords <- read.csv('S1/location.csv.gz')
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
rownames(coords) <- coords[,1]; coords[,1] <- NULL # Move first column to rownames
colnames(coords) <- c('x','y')
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
slice <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(slice@counts) # pixels to be used (a list of barcode names). 
plot_puck_continuous(slice, barcodes, slice@nUMI, ylimit = c(0,round(quantile(slice@nUMI,0.9))),title ='plot of nUMI') 

#RCTD creation
myRCTD1 <- create.RCTD(slice, reference, max_cores = 2,gene_cutoff = 0,fc_cutoff = 0, gene_cutoff_reg = 0,fc_cutoff_reg = 0, UMI_min = 1,counts_MIN = 1,CELL_MIN_INSTANCE = 1,keep_reference = T)
myRCTD1 <- run.RCTD(myRCTD1, doublet_mode = 'full')
weights <- as.data.frame(read_csv('S1/ct_proportion.csv.gz'))
rownames(weights) <- weights[,1]
weights[,1] <- NULL
myRCTD1 <- import_weights(myRCTD1, weights)

saveRDS(myRCTD1, file = "saved_rds/myRCTD1_scTEref_iw.rds")
#myRCTD1 <- readRDS('saved_rds/myRCTD1_scTEref_iw.rds')
# From h5ad_conversion.ipynb, extracting regions of the spots
S1 <- read_csv('S1/build_designmatrix_regions.csv.gz',show_col_types = FALSE)
Neural_crest <- S1[S1$regions=='Neural crest',][['barcodes']]
Diencephalon_mesencephalon <- S1[S1$regions=='Diencephalon–mesencephalon',][['barcodes']]
Telencephalon <- S1[S1$regions=='Telencephalon',][['barcodes']]
Rhombencephalon <- S1[S1$regions=='Rhombencephalon',][['barcodes']]
Mesenchymal <- S1[S1$regions=='Mesenchymal',][['barcodes']]
#Retina not included in the list as it does not exist in S1
#Retina <- S1[S1$regions=='Retina',][['barcodes']]
#names(Retina) <- 'Retina'
Epithelial <- S1[S1$regions=='Epithelial',][['barcodes']]
Choroid_plexus <- S1[S1$regions=='Choroid Plexus',][['barcodes']]
region_list <- list(Neural_crest,Diencephalon_mesencephalon,Telencephalon,Rhombencephalon,Mesenchymal,Epithelial,Choroid_plexus)
explanatory.variable1 <- build.designmatrix.regions(myRCTD1, region_list)

#Build.designmatrix.regions for S2
#Generated from h5ad_conversion.ipynb, spatial data preparation
counts <- as.data.frame(read_csv('S2/counts.csv.gz'))
coords <- read.csv('S2/location.csv.gz')
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
rownames(coords) <- coords[,1]; coords[,1] <- NULL # Move first column to rownames
colnames(coords) <- c('x','y')
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
slice <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(slice@counts) # pixels to be used (a list of barcode names). 
plot_puck_continuous(slice, barcodes, slice@nUMI, ylimit = c(0,round(quantile(slice@nUMI,0.9))),title ='plot of nUMI') 

#RCTD creation
myRCTD2 <- create.RCTD(slice, reference, max_cores = 2,gene_cutoff = 0,fc_cutoff = 0, gene_cutoff_reg = 0,fc_cutoff_reg = 0, UMI_min = 1,counts_MIN = 1,CELL_MIN_INSTANCE = 1,keep_reference = T)
myRCTD2 <- run.RCTD(myRCTD1, doublet_mode = 'full')
weights <- as.data.frame(read_csv('S2/ct_proportion.csv.gz'))
rownames(weights) <- weights[,1]
weights[,1] <- NULL
myRCTD1 <- import_weights(myRCTD1, weights)

saveRDS(myRCTD2, file = "saved_rds/myRCTD2_scTEref_iw.rds")
#myRCTD2 <- readRDS('saved_rds/myRCTD2_scTEref_iw.rds')
# From h5ad_conversion.ipynb, extracting regions of the spots
S2 <- read_csv('S2/build_designmatrix_regions.csv.gz',show_col_types = FALSE)
Neural_crest <- S2[S2$regions=='Neural crest',][['barcodes']]
Diencephalon_mesencephalon <- S2[S2$regions=='Diencephalon–mesencephalon',][['barcodes']]
Telencephalon <- S2[S2$regions=='Telencephalon',][['barcodes']]
Rhombencephalon <- S2[S2$regions=='Rhombencephalon',][['barcodes']]
Mesenchymal<- S2[S2$regions=='Mesenchymal',][['barcodes']]
#Retina not included in the list as it does not exist in S2
#Retina2 <- S2[S2$regions=='Retina',][['barcodes']]
Epithelial <- S2[S2$regions=='Epithelial',][['barcodes']]
Choroid_plexus <- S2[S2$regions=='Choroid Plexus',][['barcodes']]
region_list <- list(Neural_crest,Diencephalon_mesencephalon,Telencephalon,Rhombencephalon,Mesenchymal,Epithelial,Choroid_plexus)
explanatory.variable2 <- build.designmatrix.regions(myRCTD2, region_list)

#Build.designmatrix.regions for S3
#Generated from h5ad_conversion.ipynb, spatial data preparation
counts <- as.data.frame(read_csv('S3/counts.csv.gz'))
coords <- read.csv('S3/location.csv.gz')
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
rownames(coords) <- coords[,1]; coords[,1] <- NULL # Move first column to rownames
colnames(coords) <- c('x','y')
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
slice <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(slice@counts) # pixels to be used (a list of barcode names). 
plot_puck_continuous(slice, barcodes, slice@nUMI, ylimit = c(0,round(quantile(slice@nUMI,0.9))),title ='plot of nUMI') 

#RCTD creation
myRCTD3 <- create.RCTD(slice, reference, max_cores = 2,gene_cutoff = 0,fc_cutoff = 0, gene_cutoff_reg = 0,fc_cutoff_reg = 0, UMI_min = 1,counts_MIN = 1,CELL_MIN_INSTANCE = 1,keep_reference = T)
myRCTD3 <- run.RCTD(myRCTD3, doublet_mode = 'full')
weights <- as.data.frame(read_csv('S3/ct_proportion.csv.gz'))
rownames(weights) <- weights[,1]
weights[,1] <- NULL
myRCTD3 <- import_weights(myRCTD3, weights)

saveRDS(myRCTD3, file = "saved_rds/myRCTD3_scTEref_iw.rds")
#myRCTD3 <- readRDS('saved_rds/myRCTD3_scTEref_iw.rds')
# From h5ad_conversion.ipynb, extracting regions of the spots
S3 <- read_csv('S3/build_designmatrix_regions.csv.gz',show_col_types = FALSE)
Neural_crest <- S3[S3$regions=='Neural crest',][['barcodes']]
Diencephalon_mesencephalon <- S3[S3$regions=='Diencephalon–mesencephalon',][['barcodes']]
Telencephalon <- S3[S3$regions=='Telencephalon',][['barcodes']]
Rhombencephalon <- S3[S3$regions=='Rhombencephalon',][['barcodes']]
Mesenchymal <- S3[S3$regions=='Mesenchymal',][['barcodes']]
#Retina not applied due to only length of 1
#Retina <- S3[S3$regions=='Retina',][['barcodes']]
Epithelial <- S3[S3$regions=='Epithelial',][['barcodes']]
Choroid_plexus <- S3[S3$regions=='Choroid Plexus',][['barcodes']]
region_list <- list(Neural_crest,Diencephalon_mesencephalon,Telencephalon,Rhombencephalon,Mesenchymal,Epithelial,Choroid_plexus)
explanatory.variable3 <- build.designmatrix.regions(myRCTD3, region_list)

#Merge RCTD objects
RCTD_list <- list(myRCTD1,myRCTD2,myRCTD3)
replicate_names <- c('S1','S2','S3')
myRCTD.reps <- merge_RCTD_objects(RCTD_list,replicate_names)
saveRDS(myRCTD.reps, file = "saved_rds/myRCTD_reps_scTEref_iw.rds")
#myRCTD.reps <- readRDS('saved_rds/myRCTD_reps_scTEref_iw.rds')

#Explanatory variables preparation for each region
exvar_list <- list(explanatory.variable1, explanatory.variable2, explanatory.variable3)
exvar_list_nc <- c()
exvar_list_dm <- c()
exvar_list_t <- c()
exvar_list_r <- c()
exvar_list_m <- c()
exvar_list_e <- c()
exvar_list_cp <- c()
for (x in 1:3) {
    temp <- data.frame(exvar_list[x])
    vec <- temp[['X1']]
    names(vec) <- rownames(temp)
    exvar_list_nc[[x]] <- vec
}
for (x in 1:3) {
    temp <- data.frame(exvar_list[x])
    vec <- temp[['X2']]
    names(vec) <- rownames(temp)
    exvar_list_dm[[x]] <- vec
}
for (x in 1:3) {
    temp <- data.frame(exvar_list[x])
    vec <- temp[['X3']]
    names(vec) <- rownames(temp)
    exvar_list_t[[x]] <- vec
}
for (x in 1:3) {
    temp <- data.frame(exvar_list[x])
    vec <- temp[['X4']]
    names(vec) <- rownames(temp)
    exvar_list_r[[x]] <- vec
}
for (x in 1:3) {
    temp <- data.frame(exvar_list[x])
    vec <- temp[['X5']]
    names(vec) <- rownames(temp)
    exvar_list_m[[x]] <- vec
}
for (x in 1:3) {
    temp <- data.frame(exvar_list[x])
    vec <- temp[['X6']]
    names(vec) <- rownames(temp)
    exvar_list_e[[x]] <- vec
}
for (x in 1:3) {
    temp <- data.frame(exvar_list[x])
    vec <- temp[['X7']]
    names(vec) <- rownames(temp)
    exvar_list_cp[[x]] <- vec
}

#This session is to add one more item called doublet_mode, which is checked by aggregate_cell_types and is not included
#in the config (the one in config is called RCTDmode). Probably a bug in the original code and might be corrected in the future.
#At current version, run.CSIDE.replicates would not work if doublet_mode is specified to be FALSE.
doublet_mode <- myRCTD.reps@RCTD.reps[[1]]@config$RCTDmode
names(doublet_mode) <- 'doublet_mode'
myRCTD.reps@RCTD.reps[[1]]@config <- append(myRCTD.reps@RCTD.reps[[1]]@config, doublet_mode)
myRCTD.reps@RCTD.reps[[2]]@config <- append(myRCTD.reps@RCTD.reps[[2]]@config, doublet_mode)
myRCTD.reps@RCTD.reps[[3]]@config <- append(myRCTD.reps@RCTD.reps[[3]]@config, doublet_mode)

#For Neural crest region
#mg less than cell_type_threshold of 1
cell_types <- c('mesen','tel_NPC','tel_neuron','nc','pns','dien_NPC','dien_neuron','hind_neuron','hind_NPC','choroid_plexus','epi','retina')
myRCTD.reps_nc <- run.CSIDE.replicates(myRCTD.reps,cell_types, explanatory.variable.replicates = exvar_list_nc, cell_type_threshold = 1,doublet_mode = F, weight_threshold = 0.7,fdr = 0.1,population_de = T)
saveRDS(myRCTD.reps_nc, file = "saved_rds/myRCTD_reps_nc_scTEref_iw.rds")
#myRCTD.reps_nc <- readRDS('saved_rds/myRCTD_reps_nc_scTEref_iw.rds')

#For Dien-mesen region
#mg less than cell_type_threshold of 1
cell_types <- c('mesen','tel_NPC','tel_neuron','nc','pns','dien_NPC','dien_neuron','hind_neuron','hind_NPC','choroid_plexus','epi','retina')
myRCTD.reps_dm <- run.CSIDE.replicates(myRCTD.reps,cell_types, explanatory.variable.replicates = exvar_list_dm, cell_type_threshold = 1,doublet_mode = F, weight_threshold = 0.8,fdr = 0.1, population_de = T)
saveRDS(myRCTD.reps_dm, file = "saved_rds/myRCTD_reps_dm_scTEref_iw.rds")
#myRCTD.reps_dm <- readRDS('saved_rds/myRCTD_reps_dm_scTEref_iw.rds')

#For Rhombencephalon region
#mg less than minimum cell_type_threshold of 1 
cell_types <- c('mesen','tel_NPC','tel_neuron','nc','pns','dien_NPC','dien_neuron','hind_neuron','hind_NPC','choroid_plexus','epi','retina')
myRCTD.reps_r <- run.CSIDE.replicates(myRCTD.reps,cell_types, explanatory.variable.replicates = exvar_list_r, cell_type_threshold = 1,doublet_mode = F, weight_threshold = 0.8,fdr = 0.1, population_de = T)
saveRDS(myRCTD.reps_r, file = "saved_rds/myRCTD_reps_r_scTEref_iw.rds")
#myRCTD.reps_r <- readRDS('saved_rds/myRCTD_reps_r_scTEref_iw.rds')

#For Mesenchymal region
# mg not passing cell_type_threshold, hind_NPC,hind_neuron epi, dien_neuron, dien_NPC, pns out of bound
cell_types <- c('mesen','tel_NPC','tel_neuron','nc','choroid_plexus','retina')
myRCTD.reps_m <- run.CSIDE.replicates(myRCTD.reps,cell_types, explanatory.variable.replicates = exvar_list_m, cell_type_threshold = 1,doublet_mode = F, weight_threshold = 0.25,fdr = 0.1, population_de = T)
saveRDS(myRCTD.reps_m, file = "saved_rds/myRCTD_reps_m_scTEref_iw.rds")
#myRCTD.reps_m <- readRDS('saved_rds/myRCTD_reps_m_scTEref_iw.rds')

#For Epithelial region
#epi, cp, retina, hind_neuron,dien_neuron,dien_NPC OOB, mg not passing cell_type_threshold
cell_types <- c('mesen','tel_NPC','tel_neuron','nc','pns','hind_NPC')
myRCTD.reps_e <- run.CSIDE.replicates(myRCTD.reps,cell_types, explanatory.variable.replicates = exvar_list_e, cell_type_threshold = 1,doublet_mode = F, weight_threshold = 0.3,fdr = 0.1,  population_de = T)
saveRDS(myRCTD.reps_e, file = "saved_rds/myRCTD_reps_e_scTEref_iw.rds")
#myRCTD.reps_e <- readRDS('saved_rds/myRCTD_reps_e_scTEref_iw.rds')

#For Choroid Plexus region
cell_types <- c('mesen','tel_NPC','tel_neuron','nc','pns','dien_NPC','dien_neuron','hind_neuron','hind_NPC','choroid_plexus','epi','retina')
myRCTD.reps_cp <- run.CSIDE.replicates(myRCTD.reps,cell_types, explanatory.variable.replicates = exvar_list_cp, cell_type_threshold = 1,doublet_mode = F, weight_threshold = 0.6,fdr = 0.1, population_de = T)
saveRDS(myRCTD.reps_cp, file = "saved_rds/myRCTD_reps_cp_scTEref_iw.rds")
#myRCTD.reps_cp <- readRDS('saved_rds/myRCTD_reps_cp_scTEref_iw.rds')

#For Telecephalon region
# mg not passing cell_type_threshold, hind_NPC, dien_neuron , hind_neuron oob
cell_types <- c('mesen','tel_NPC','tel_neuron','nc','pns','dien_NPC','choroid_plexus','epi','retina')
myRCTD.reps_t <- run.CSIDE.replicates(myRCTD.reps,cell_types, explanatory.variable.replicates = exvar_list_t, cell_type_threshold = 1,doublet_mode = F, weight_threshold = 0.35,fdr = 0.1, population_de = T)
saveRDS(myRCTD.reps_t, file = "saved_rds/myRCTD_reps_t_scTEref_iw.rds")
#myRCTD.reps_t <- readRDS('saved_rds/myRCTD_reps_t_scTEref_iw.rds')


