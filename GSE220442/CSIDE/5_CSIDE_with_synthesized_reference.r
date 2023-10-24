library(spacexr)
library(Matrix)
library(doParallel)
library(readr)
library(ggplot2)

setwd('/work/project/ext_014/GSE220442/')

#Customize dataframes and create reference
#Sys.setenv(VROOM_CONNECTION_SIZE=1000000)
counts <- as.data.frame(read_csv('rawdata/counts.csv.gz'))
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
meta_data <- read_csv("rawdata/MTG_pure_annotation.csv.gz")
cell_types <- setNames(meta_data[[2]], meta_data[[1]])
cell_types <- as.factor(cell_types) # convert to factor data type
col_sum <- colSums(counts)
meta_data[,5] <- col_sum
colnames(meta_data) <- c('barcode','celltype','sample_id','donor','nUMI')
nUMI <- setNames(meta_data[[5]], meta_data[[1]])
reference <- Reference(counts, cell_types, nUMI)

setwd('/work/project/ext_014/GSE220442/CSIDE/')

#Build.designmatrix.regions for CT1
#Generated from h5ad_conversion.ipynb, spatial data preparation
counts <- as.data.frame(read_csv('CT1/counts.csv.gz'))
coords <- read.csv('CT1/location.csv.gz')
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
rownames(coords) <- coords[,1]; coords[,1] <- NULL # Move first column to rownames
colnames(coords) <- c('x','y')
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
slice <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(slice@counts) # pixels to be used (a list of barcode names). 
plot_puck_continuous(slice, barcodes, slice@nUMI, ylimit = c(0,round(quantile(slice@nUMI,0.9))),title ='plot of nUMI') 

#RCTD creation
cell_type_profiles <- as.data.frame(read_csv('Reference/cell_type_profiles.csv.gz'))
rownames(cell_type_profiles) <- cell_type_profiles[,1]
cell_type_profiles[,1] <- NULL
myRCTD1 <- create.RCTD(slice, reference, max_cores = 2,gene_cutoff = 0,fc_cutoff = 0, gene_cutoff_reg = 0,fc_cutoff_reg = 0, UMI_min = 1,counts_MIN = 1,CELL_MIN_INSTANCE = 1,keep_reference = T,cell_type_profiles = cell_type_profiles)
myRCTD1 <- run.RCTD(myRCTD1, doublet_mode = 'full')
weights <- as.data.frame(read_csv('CT1/ct_proportion.csv.gz'))
rownames(weights) <- weights[,1]
weights[,1] <- NULL
myRCTD1 <- import_weights(myRCTD1, weights)

saveRDS(myRCTD1, file = "saved_rds/myRCTD1.rds")
#myRCTD1 <- readRDS('saved_rds/myRCTD1.rds')
# From h5ad_conversion.ipynb, extracting regions of the spots
CT1 <- read_csv('CT1/build_designmatrix_regions.csv.gz',show_col_types = FALSE)
L1 <- CT1[CT1$regions=='Layer 1',][['barcodes']]
L2 <- CT1[CT1$regions=='Layer 2',][['barcodes']]
L3 <- CT1[CT1$regions=='Layer 3',][['barcodes']]
L4 <- CT1[CT1$regions=='Layer 4',][['barcodes']]
L5 <- CT1[CT1$regions=='Layer 5',][['barcodes']]
L6 <- CT1[CT1$regions=='Layer 6',][['barcodes']]
WM <- CT1[CT1$regions=='White Matter',][['barcodes']]
region_list <- list(L1,L2,L3,L4,L5,L6,WM)
explanatory.variable1 <- build.designmatrix.regions(myRCTD1, region_list)

#Build.designmatrix.regions for CT2
#Generated from h5ad_conversion.ipynb, spatial data preparation
counts <- as.data.frame(read_csv('CT2/counts.csv.gz'))
coords <- read.csv('CT2/location.csv.gz')
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
rownames(coords) <- coords[,1]; coords[,1] <- NULL # Move first column to rownames
colnames(coords) <- c('x','y')
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
slice <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(slice@counts) # pixels to be used (a list of barcode names). 
plot_puck_continuous(slice, barcodes, slice@nUMI, ylimit = c(0,round(quantile(slice@nUMI,0.9))),title ='plot of nUMI') 

#RCTD creation
myRCTD2 <- create.RCTD(slice, reference, max_cores = 2,gene_cutoff = 0,fc_cutoff = 0, gene_cutoff_reg = 0,fc_cutoff_reg = 0, UMI_min = 1,counts_MIN = 1,CELL_MIN_INSTANCE = 1,keep_reference = T,cell_type_profiles = cell_type_profiles)
myRCTD2 <- run.RCTD(myRCTD2, doublet_mode = 'full')
weights <- as.data.frame(read_csv('CT2/ct_proportion.csv.gz'))
rownames(weights) <- weights[,1]
weights[,1] <- NULL
myRCTD2 <- import_weights(myRCTD2, weights)

saveRDS(myRCTD2, file = "saved_rds/myRCTD2.rds")
#myRCTD2 <- readRDS('saved_rds/myRCTD2.rds')
# From h5ad_conversion.ipynb, extracting regions of the spots
CT2 <- read_csv('CT2/build_designmatrix_regions.csv.gz',show_col_types = FALSE)
L1 <- CT2[CT2$regions=='Layer 1',][['barcodes']]
L2 <- CT2[CT2$regions=='Layer 2',][['barcodes']]
L3 <- CT2[CT2$regions=='Layer 3',][['barcodes']]
L4 <- CT2[CT2$regions=='Layer 4',][['barcodes']]
L5 <- CT2[CT2$regions=='Layer 5',][['barcodes']]
L6 <- CT2[CT2$regions=='Layer 6',][['barcodes']]
WM <- CT2[CT2$regions=='White Matter',][['barcodes']]
region_list <- list(L1,L2,L3,L4,L5,L6,WM)
explanatory.variable2 <- build.designmatrix.regions(myRCTD2, region_list)

#Build.designmatrix.regions for CT3
#Generated from h5ad_conversion.ipynb, spatial data preparation
counts <- as.data.frame(read_csv('CT3/counts.csv.gz'))
coords <- read.csv('CT3/location.csv.gz')
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
rownames(coords) <- coords[,1]; coords[,1] <- NULL # Move first column to rownames
colnames(coords) <- c('x','y')
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
slice <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(slice@counts) # pixels to be used (a list of barcode names). 
plot_puck_continuous(slice, barcodes, slice@nUMI, ylimit = c(0,round(quantile(slice@nUMI,0.9))),title ='plot of nUMI') 

#RCTD creation
myRCTD3 <- create.RCTD(slice, reference, max_cores = 2,gene_cutoff = 0,fc_cutoff = 0, gene_cutoff_reg = 0,fc_cutoff_reg = 0, UMI_min = 1,counts_MIN = 1,CELL_MIN_INSTANCE = 1,keep_reference = T,cell_type_profiles = cell_type_profiles)
myRCTD3 <- run.RCTD(myRCTD3, doublet_mode = 'full')
weights <- as.data.frame(read_csv('CT3/ct_proportion.csv.gz'))
rownames(weights) <- weights[,1]
weights[,1] <- NULL
myRCTD3 <- import_weights(myRCTD3, weights)

saveRDS(myRCTD3, file = "saved_rds/myRCTD3.rds")
#myRCTD3 <- readRDS('saved_rds/myRCTD3.rds')
# From h5ad_conversion.ipynb, extracting regions of the spots
CT3 <- read_csv('CT3/build_designmatrix_regions.csv.gz',show_col_types = FALSE)
L1 <- CT3[CT3$regions=='Layer 1',][['barcodes']]
L2 <- CT3[CT3$regions=='Layer 2',][['barcodes']]
L3 <- CT3[CT3$regions=='Layer 3',][['barcodes']]
L4 <- CT3[CT3$regions=='Layer 4',][['barcodes']]
L5 <- CT3[CT3$regions=='Layer 5',][['barcodes']]
L6 <- CT3[CT3$regions=='Layer 6',][['barcodes']]
WM <- CT3[CT3$regions=='White Matter',][['barcodes']]
region_list <- list(L1,L2,L3,L4,L5,L6,WM)
explanatory.variable3 <- build.designmatrix.regions(myRCTD3, region_list)

#Build.designmatrix.regions for AD1
#Generated from h5ad_conversion.ipynb, spatial data preparation
counts <- as.data.frame(read_csv('AD1/counts.csv.gz'))
coords <- read.csv('AD1/location.csv.gz')
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
rownames(coords) <- coords[,1]; coords[,1] <- NULL # Move first column to rownames
colnames(coords) <- c('x','y')
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
slice <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(slice@counts) # pixels to be used (a list of barcode names). 
plot_puck_continuous(slice, barcodes, slice@nUMI, ylimit = c(0,round(quantile(slice@nUMI,0.9))),title ='plot of nUMI') 

#RCTD creation
myRCTD4 <- create.RCTD(slice, reference, max_cores = 2,gene_cutoff = 0,fc_cutoff = 0, gene_cutoff_reg = 0,fc_cutoff_reg = 0, UMI_min = 1,counts_MIN = 1,CELL_MIN_INSTANCE = 1,keep_reference = T,cell_type_profiles = cell_type_profiles)
myRCTD4 <- run.RCTD(myRCTD4, doublet_mode = 'full')
weights <- as.data.frame(read_csv('AD1/ct_proportion.csv.gz'))
rownames(weights) <- weights[,1]
weights[,1] <- NULL
myRCTD4 <- import_weights(myRCTD4, weights)

saveRDS(myRCTD4, file = "saved_rds/myRCTD4.rds")
#myRCTD4 <- readRDS('saved_rds/myRCTD4.rds')
# From h5ad_conversion.ipynb, extracting regions of the spots
AD1 <- read_csv('AD1/build_designmatrix_regions.csv.gz',show_col_types = FALSE)
L1 <- AD1[AD1$regions=='Layer 1',][['barcodes']]
L2 <- AD1[AD1$regions=='Layer 2',][['barcodes']]
L3 <- AD1[AD1$regions=='Layer 3',][['barcodes']]
L4 <- AD1[AD1$regions=='Layer 4',][['barcodes']]
L5 <- AD1[AD1$regions=='Layer 5',][['barcodes']]
L6 <- AD1[AD1$regions=='Layer 6',][['barcodes']]
WM <- AD1[AD1$regions=='White Matter',][['barcodes']]
region_list <- list(L1,L2,L3,L4,L5,L6,WM)
explanatory.variable4 <- build.designmatrix.regions(myRCTD4, region_list)

#Build.designmatrix.regions for AD2
#Generated from h5ad_conversion.ipynb, spatial data preparation
counts <- as.data.frame(read_csv('AD2/counts.csv.gz'))
coords <- read.csv('AD2/location.csv.gz')
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
rownames(coords) <- coords[,1]; coords[,1] <- NULL # Move first column to rownames
colnames(coords) <- c('x','y')
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
slice <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(slice@counts) # pixels to be used (a list of barcode names). 
plot_puck_continuous(slice, barcodes, slice@nUMI, ylimit = c(0,round(quantile(slice@nUMI,0.9))),title ='plot of nUMI') 

#RCTD creation
myRCTD5 <- create.RCTD(slice, reference, max_cores = 2,gene_cutoff = 0,fc_cutoff = 0, gene_cutoff_reg = 0,fc_cutoff_reg = 0, UMI_min = 1,counts_MIN = 1,CELL_MIN_INSTANCE = 1,keep_reference = T,cell_type_profiles = cell_type_profiles)
myRCTD5 <- run.RCTD(myRCTD5, doublet_mode = 'full')
weights <- as.data.frame(read_csv('AD2/ct_proportion.csv.gz'))
rownames(weights) <- weights[,1]
weights[,1] <- NULL
myRCTD5 <- import_weights(myRCTD5, weights)

saveRDS(myRCTD5, file = "saved_rds/myRCTD5.rds")
#myRCTD5 <- readRDS('saved_rds/myRCTD5.rds')
# From h5ad_conversion.ipynb, extracting regions of the spots
AD2 <- read_csv('AD2/build_designmatrix_regions.csv.gz',show_col_types = FALSE)
#No Layer 1 in AD2
#L1 <- AD2[AD2$regions=='Layer 1',][['barcodes']]
L2 <- AD2[AD2$regions=='Layer 2',][['barcodes']]
L3 <- AD2[AD2$regions=='Layer 3',][['barcodes']]
L4 <- AD2[AD2$regions=='Layer 4',][['barcodes']]
L5 <- AD2[AD2$regions=='Layer 5',][['barcodes']]
L6 <- AD2[AD2$regions=='Layer 6',][['barcodes']]
WM <- AD2[AD2$regions=='White Matter',][['barcodes']]
region_list <- list(L2,L3,L4,L5,L6,WM)
explanatory.variable5 <- build.designmatrix.regions(myRCTD5, region_list)
#to modify the matrix into 7 columns
L1col <- rep(c(0),each=3439)
explanatory.variable5 <- cbind(explanatory.variable5, L1col)
explanatory.variable5[,7] <- explanatory.variable5[,6]
explanatory.variable5[,6] <- explanatory.variable5[,5]
explanatory.variable5[,5] <- explanatory.variable5[,4]
explanatory.variable5[,4] <- explanatory.variable5[,3]
explanatory.variable5[,3] <- explanatory.variable5[,2]
explanatory.variable5[,2] <- explanatory.variable5[,1]
explanatory.variable5[,1] <- L1col
colnames(explanatory.variable5) <- NULL

#Build.designmatrix.regions for AD3
#Generated from h5ad_conversion.ipynb, spatial data preparation
counts <- as.data.frame(read_csv('AD3/counts.csv.gz'))
coords <- read.csv('AD3/location.csv.gz')
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
rownames(coords) <- coords[,1]; coords[,1] <- NULL # Move first column to rownames
colnames(coords) <- c('x','y')
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
slice <- SpatialRNA(coords, counts, nUMI)
barcodes <- colnames(slice@counts) # pixels to be used (a list of barcode names). 
plot_puck_continuous(slice, barcodes, slice@nUMI, ylimit = c(0,round(quantile(slice@nUMI,0.9))),title ='plot of nUMI') 

#RCTD creation
myRCTD6 <- create.RCTD(slice, reference, max_cores = 2,gene_cutoff = 0,fc_cutoff = 0, gene_cutoff_reg = 0,fc_cutoff_reg = 0, UMI_min = 1,counts_MIN = 1,CELL_MIN_INSTANCE = 1,keep_reference = T,cell_type_profiles = cell_type_profiles)
myRCTD6 <- run.RCTD(myRCTD6, doublet_mode = 'full')
weights <- as.data.frame(read_csv('AD3/ct_proportion.csv.gz'))
rownames(weights) <- weights[,1]
weights[,1] <- NULL
myRCTD6 <- import_weights(myRCTD6, weights)

saveRDS(myRCTD6, file = "saved_rds/myRCTD6.rds")
#myRCTD6 <- readRDS('saved_rds/myRCTD6.rds')
# From h5ad_conversion.ipynb, extracting regions of the spots
AD3 <- read_csv('AD3/build_designmatrix_regions.csv.gz',show_col_types = FALSE)
L1 <- AD3[AD3$regions=='Layer 1',][['barcodes']]
L2 <- AD3[AD3$regions=='Layer 2',][['barcodes']]
L3 <- AD3[AD3$regions=='Layer 3',][['barcodes']]
L4 <- AD3[AD3$regions=='Layer 4',][['barcodes']]
L5 <- AD3[AD3$regions=='Layer 5',][['barcodes']]
L6 <- AD3[AD3$regions=='Layer 6',][['barcodes']]
WM <- AD3[AD3$regions=='White Matter',][['barcodes']]
region_list <- list(L1,L2,L3,L4,L5,L6,WM)
explanatory.variable6 <- build.designmatrix.regions(myRCTD6, region_list)

#CREATE RCTD_replicates for all regions except L1 (L1 does not exist in AD2)
#RCTD_list <- list(myRCTD1,myRCTD2,myRCTD3,myRCTD4,myRCTD5,myRCTD6)
replicate_names <- c('CT1','CT2','CT3','AD1','AD2','AD3')
myRCTD.reps <- merge_RCTD_objects(RCTD_list,replicate_names)
saveRDS(myRCTD.reps, file = "saved_rds/myRCTD_reps_im.rds")
#myRCTD.reps <- readRDS('saved_rds/myRCTD_reps_im.rds')
#Create RCTD_replicates for L1
RCTD_list <- list(myRCTD1,myRCTD2,myRCTD3,myRCTD4,myRCTD6)
replicate_names <- c('CT1','CT2','CT3','AD1','AD3')
myRCTD.reps_for_L1 <- merge_RCTD_objects(RCTD_list,replicate_names)
saveRDS(myRCTD.reps_for_L1, file = "saved_rds/myRCTD_reps_im_for_L1.rds")
#myRCTD.reps_for_L1 <- readRDS('saved_rds/myRCTD_reps_im_for_L1.rds')

#Explanatory variables preparation for each region
exvar_list <- list(explanatory.variable1, explanatory.variable2, explanatory.variable3,explanatory.variable4,explanatory.variable5,explanatory.variable6)
exvar_list_L1 <- c()
exvar_list_L2 <- c()
exvar_list_L3 <- c()
exvar_list_L4 <- c()
exvar_list_L5 <- c()
exvar_list_L6 <- c()
exvar_list_WM <- c()
#L1 does not exist in AD2
for (x in 1:6) {
    temp <- data.frame(exvar_list[x])
    vec <- temp[['X1']]
    names(vec) <- rownames(temp)
    exvar_list_L1[[x]] <- vec
}
for (x in 1:6) {
    temp <- data.frame(exvar_list[x])
    vec <- temp[['X2']]
    names(vec) <- rownames(temp)
    exvar_list_L2[[x]] <- vec
}
for (x in 1:6) {
    temp <- data.frame(exvar_list[x])
    vec <- temp[['X3']]
    names(vec) <- rownames(temp)
    exvar_list_L3[[x]] <- vec
}
for (x in 1:6) {
    temp <- data.frame(exvar_list[x])
    vec <- temp[['X4']]
    names(vec) <- rownames(temp)
    exvar_list_L4[[x]] <- vec
}
for (x in 1:6) {
    temp <- data.frame(exvar_list[x])
    vec <- temp[['X5']]
    names(vec) <- rownames(temp)
    exvar_list_L5[[x]] <- vec
}
for (x in 1:6) {
    temp <- data.frame(exvar_list[x])
    vec <- temp[['X6']]
    names(vec) <- rownames(temp)
    exvar_list_L6[[x]] <- vec
}
for (x in 1:6) {
    temp <- data.frame(exvar_list[x])
    vec <- temp[['X7']]
    names(vec) <- rownames(temp)
    exvar_list_WM[[x]] <- vec
}

#This session is to add one more item called doublet_mode, which is checked by aggregate_cell_types and is not included
#in the config (the one in config is called RCTDmode). Probably a bug in the original code and might be corrected in the future.
#At current version, run.CSIDE.replicates would not work if doublet_mode is specified to be FALSE.
doublet_mode <- myRCTD.reps@RCTD.reps[[1]]@config$RCTDmode
names(doublet_mode) <- 'doublet_mode'
#for all regions except L1
myRCTD.reps@RCTD.reps[[1]]@config <- append(myRCTD.reps@RCTD.reps[[1]]@config, doublet_mode)
myRCTD.reps@RCTD.reps[[2]]@config <- append(myRCTD.reps@RCTD.reps[[2]]@config, doublet_mode)
myRCTD.reps@RCTD.reps[[3]]@config <- append(myRCTD.reps@RCTD.reps[[3]]@config, doublet_mode)
myRCTD.reps@RCTD.reps[[4]]@config <- append(myRCTD.reps@RCTD.reps[[4]]@config, doublet_mode)
myRCTD.reps@RCTD.reps[[5]]@config <- append(myRCTD.reps@RCTD.reps[[5]]@config, doublet_mode)
myRCTD.reps@RCTD.reps[[6]]@config <- append(myRCTD.reps@RCTD.reps[[6]]@config, doublet_mode)
#for L1
doublet_mode <- myRCTD.reps_for_L1@RCTD.reps[[1]]@config$RCTDmode
names(doublet_mode) <- 'doublet_mode'
myRCTD.reps_for_L1@RCTD.reps[[1]]@config <- append(myRCTD.reps_for_L1@RCTD.reps[[1]]@config, doublet_mode)
myRCTD.reps_for_L1@RCTD.reps[[2]]@config <- append(myRCTD.reps_for_L1@RCTD.reps[[2]]@config, doublet_mode)
myRCTD.reps_for_L1@RCTD.reps[[3]]@config <- append(myRCTD.reps_for_L1@RCTD.reps[[3]]@config, doublet_mode)
myRCTD.reps_for_L1@RCTD.reps[[4]]@config <- append(myRCTD.reps_for_L1@RCTD.reps[[4]]@config, doublet_mode)
myRCTD.reps_for_L1@RCTD.reps[[5]]@config <- append(myRCTD.reps_for_L1@RCTD.reps[[5]]@config, doublet_mode)

exvar_list_L1[5] <- exvar_list_L1[6]
exvar_list_L1[6] <- NULL

#For L1 region
cell_types <- c('Astro', 'Endo', 'Exc', 'Inh', 'Micro', 'OPC','Oligo')
myRCTD.reps_L1_CSIDE <- run.CSIDE.replicates(myRCTD.reps_for_L1,cell_types, explanatory.variable.replicates = exvar_list_L1, cell_type_threshold = 1,fdr = 0.1,doublet_mode = F, weight_threshold = 0.8)
saveRDS(myRCTD.reps_L1_CSIDE, file = "saved_rds/myRCTD_reps_L1_CSIDE_im_ft.rds")
#myRCTD.reps_L1_CSIDE <- readRDS('saved_rds/myRCTD_reps_L1_CSIDE_im_ft.rds')
myRCTD.reps_L1_popin <- CSIDE.population.inference(myRCTD.reps_L1_CSIDE,fdr = 0.1)
saveRDS(myRCTD.reps_L1_popin, file = "saved_rds/myRCTD_reps_L1_popin_im_ft.rds")
#myRCTD.reps_L1_popin <- readRDS('saved_rds/myRCTD_reps_L1_popin_im.rds')

#Oligo does not exist in CT2 after CSIDE run
meta.design.matrix <- data.frame('AD' = c(0,0,1,1))
myRCTD.reps_meta_L1 <- CSIDE.population.inference(myRCTD.reps_L1_CSIDE,fdr = 0.1,meta = TRUE,meta.design.matrix = meta.design.matrix,meta.test_var = 'AD')
saveRDS(myRCTD.reps_meta_L1, file = "saved_rds/myRCTD_reps_L1_meta_im_ft.rds")
#myRCTD.reps_meta_L1 <- readRDS('saved_rds/myRCTD_reps_L1_meta_im.rds')

#For L2 region
myRCTD.reps_L2 <- run.CSIDE.replicates(myRCTD.reps,cell_types, explanatory.variable.replicates = exvar_list_L2, cell_type_threshold = 1,fdr = 0.1,doublet_mode = F, weight_threshold = 0.8)
saveRDS(myRCTD.reps_L2, file = "saved_rds/myRCTD_reps_L2_im.rds")
#myRCTD.reps_L2 <- readRDS('saved_rds/myRCTD_reps_L2_im.rds')
myRCTD.reps_L2_popin <- CSIDE.population.inference(myRCTD.reps_L2,fdr = 0.1)
saveRDS(myRCTD.reps_L2_popin, file = "saved_rds/myRCTD_reps_L2_popin_im.rds")
#myRCTD.reps_L2_popin <- readRDS('saved_rds/myRCTD_reps_L2_popin_im.rds')

meta.design.matrix <- data.frame('AD' = c(0,0,0,1,1,1))
myRCTD.reps_meta_L2 <- CSIDE.population.inference(myRCTD.reps_L2,fdr = 0.1,meta = TRUE,meta.design.matrix = meta.design.matrix,meta.test_var = 'AD')
saveRDS(myRCTD.reps_meta_L2, file = "saved_rds/myRCTD_reps_L2_meta_im.rds")
#myRCTD.reps_meta_L2 <- readRDS('saved_rds/myRCTD_reps_L2_meta_im.rds')

#For L3 region
myRCTD.reps_L3 <- run.CSIDE.replicates(myRCTD.reps,cell_types, explanatory.variable.replicates = exvar_list_L3, cell_type_threshold = 1,fdr = 0.1,doublet_mode = F, weight_threshold = 0.8)
saveRDS(myRCTD.reps_L3, file = "saved_rds/myRCTD_reps_L3_im.rds")
#myRCTD.reps_L3 <- readRDS('saved_rds/myRCTD_reps_L3_im.rds')
myRCTD.reps_L3_popin <- CSIDE.population.inference(myRCTD.reps_L3,fdr = 0.1)
saveRDS(myRCTD.reps_L3_popin, file = "saved_rds/myRCTD_reps_L3_popin_im.rds")
#myRCTD.reps_L3_popin <- readRDS('saved_rds/myRCTD_reps_L3_popin_im.rds')

meta.design.matrix <- data.frame('AD' = c(0,0,0,1,1,1))
myRCTD.reps_meta_L3 <- CSIDE.population.inference(myRCTD.reps_L3,fdr = 0.1,meta = TRUE,meta.design.matrix = meta.design.matrix,meta.test_var = 'AD')
saveRDS(myRCTD.reps_meta_L3, file = "saved_rds/myRCTD_reps_L3_meta_im.rds")
#myRCTD.reps_meta_L3 <- readRDS('saved_rds/myRCTD_reps_L3_meta_im.rds')

#For L4 region
myRCTD.reps_L4 <- run.CSIDE.replicates(myRCTD.reps,cell_types, explanatory.variable.replicates = exvar_list_L4, cell_type_threshold = 1,fdr = 0.1,doublet_mode = F, weight_threshold = 0.8)
saveRDS(myRCTD.reps_L4, file = "saved_rds/myRCTD_reps_L4_im.rds")
#myRCTD.reps_L4 <- readRDS('saved_rds/myRCTD_reps_L4_im.rds')
myRCTD.reps_L4_popin <- CSIDE.population.inference(myRCTD.reps_L4,fdr = 0.1)
saveRDS(myRCTD.reps_L4_popin, file = "saved_rds/myRCTD_reps_L4_popin_im.rds")
#myRCTD.reps_L4_popin <- readRDS('saved_rds/myRCTD_reps_L4_popin_im.rds')

meta.design.matrix <- data.frame('AD' = c(0,0,0,1,1,1))
myRCTD.reps_meta_L4 <- CSIDE.population.inference(myRCTD.reps_L4,fdr = 0.1,meta = TRUE,meta.design.matrix = meta.design.matrix,meta.test_var = 'AD')
saveRDS(myRCTD.reps_meta_L4, file = "saved_rds/myRCTD_reps_L4_meta_im.rds")
#myRCTD.reps_meta_L4 <- readRDS('saved_rds/myRCTD_reps_L4_meta_im.rds')

#For L5 region
myRCTD.reps_L5 <- run.CSIDE.replicates(myRCTD.reps,cell_types, explanatory.variable.replicates = exvar_list_L5, cell_type_threshold = 1,fdr = 0.1,doublet_mode = F, weight_threshold = 0.8)
saveRDS(myRCTD.reps_L5, file = "saved_rds/myRCTD_reps_L5_im.rds")
#myRCTD.reps_L5 <- readRDS('saved_rds/myRCTD_reps_L5_im.rds')
myRCTD.reps_L5_popin <- CSIDE.population.inference(myRCTD.reps_L5,fdr = 0.1)
saveRDS(myRCTD.reps_L5_popin, file = "saved_rds/myRCTD_reps_L5_popin_im.rds")
#myRCTD.reps_L5_popin <- readRDS('saved_rds/myRCTD_reps_L5_popin_im.rds')

meta.design.matrix <- data.frame('AD' = c(0,0,0,1,1,1))
myRCTD.reps_meta_L5 <- CSIDE.population.inference(myRCTD.reps_L5,fdr = 0.1,meta = TRUE,meta.design.matrix = meta.design.matrix,meta.test_var = 'AD')
saveRDS(myRCTD.reps_meta_L5, file = "saved_rds/myRCTD_reps_L5_meta_im.rds")
#myRCTD.reps_meta_L5 <- readRDS('saved_rds/myRCTD_reps_L5_meta_im.rds')

#For L6 region
myRCTD.reps_L6 <- run.CSIDE.replicates(myRCTD.reps,cell_types, explanatory.variable.replicates = exvar_list_L6, cell_type_threshold = 1,fdr = 0.1,doublet_mode = F, weight_threshold = 0.8)
saveRDS(myRCTD.reps_L6, file = "saved_rds/myRCTD_reps_L6_im.rds")
#myRCTD.reps_L6 <- readRDS('saved_rds/myRCTD_reps_L6_im.rds')
myRCTD.reps_L6_popin <- CSIDE.population.inference(myRCTD.reps_L6,fdr = 0.1)
saveRDS(myRCTD.reps_L6_popin, file = "saved_rds/myRCTD_reps_L6_popin_im.rds")
#myRCTD.reps_L6_popin <- readRDS('saved_rds/myRCTD_reps_L6_popin_im.rds')

meta.design.matrix <- data.frame('AD' = c(0,0,0,1,1,1))
myRCTD.reps_meta_L6 <- CSIDE.population.inference(myRCTD.reps_L6,fdr = 0.1,meta = TRUE,meta.design.matrix = meta.design.matrix,meta.test_var = 'AD')
saveRDS(myRCTD.reps_meta_L6, file = "saved_rds/myRCTD_reps_L6_meta_im.rds")
#myRCTD.reps_meta_L6 <- readRDS('saved_rds/myRCTD_reps_L6_meta_im.rds')

#For WM region
myRCTD.reps_WM <- run.CSIDE.replicates(myRCTD.reps,cell_types, explanatory.variable.replicates = exvar_list_WM, cell_type_threshold = 1,fdr = 0.1,doublet_mode = F, weight_threshold = 0.8)
saveRDS(myRCTD.reps_WM, file = "saved_rds/myRCTD_reps_WM_im.rds")
#myRCTD.reps_WM <- readRDS('saved_rds/myRCTD_reps_WM_im.rds')
myRCTD.reps_WM_popin <- CSIDE.population.inference(myRCTD.reps_WM,fdr = 0.1)
saveRDS(myRCTD.reps_WM_popin, file = "saved_rds/myRCTD_reps_WM_popin_im.rds")
#myRCTD.reps_WM_popin <- readRDS('saved_rds/myRCTD_reps_WM_popin_im.rds')

meta.design.matrix <- data.frame('AD' = c(0,0,0,1,1,1))
myRCTD.reps_meta_WM <- CSIDE.population.inference(myRCTD.reps_WM,fdr = 0.1,meta = TRUE,meta.design.matrix = meta.design.matrix,meta.test_var = 'AD')
saveRDS(myRCTD.reps_meta_WM, file = "saved_rds/myRCTD_reps_WM_meta_im.rds")
#myRCTD.reps_meta_WM <- readRDS('saved_rds/myRCTD_reps_WM_meta_im.rds')
