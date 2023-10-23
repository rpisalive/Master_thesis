library(plyr)
library(dplyr)
library(tibble)
library(spacexr)

setwd('/work/project/ext_014/GSE220442/CSIDE')

#Load CISDE population inference results
myRCTD_reps_L1 <- readRDS('saved_rds/myRCTD_reps_L1_popin_im_ft.rds')
#myRCTD_reps_L2 <- readRDS('saved_rds/myRCTD_reps_L2_popin_im.rds')
#myRCTD_reps_L3 <- readRDS('saved_rds/myRCTD_reps_L3_popin_im.rds')
#myRCTD_reps_L4 <- readRDS('saved_rds/myRCTD_reps_L4_popin_im.rds')
#myRCTD_reps_L5 <- readRDS('saved_rds/myRCTD_reps_L5_popin_im.rds')
#myRCTD_reps_L6 <- readRDS('saved_rds/myRCTD_reps_L6_popin_im.rds')
#myRCTD_reps_WM <- readRDS('saved_rds/myRCTD_reps_WM_popin_im.rds')

#Add column to identify non-significant genes
for (cell_type in names(myRCTD_reps_L1@population_de_results)){
    myRCTD_reps_L1@population_de_results[[cell_type]]$cell_name <- cell_type
    myRCTD_reps_L1@population_de_results[[cell_type]]$sig_cat <- 'non_sig'
}
#Add column to identify significant genes
for (cell_type in names(myRCTD_reps_L1@population_sig_gene_df)){
    if (dim(myRCTD_reps_L1@population_sig_gene_df[[cell_type]])[1] == 0){
        myRCTD_reps_L1@population_sig_gene_df[[cell_type]] <- NULL
    } else {
        myRCTD_reps_L1@population_sig_gene_df[[cell_type]]$cell_name <- cell_type
        myRCTD_reps_L1@population_sig_gene_df[[cell_type]]$sig_cat <- 'sig'
        }
}
#Add column to identify non-significant genes
for (cell_type in names(myRCTD_reps_L2@population_de_results)){
    myRCTD_reps_L2@population_de_results[[cell_type]]$cell_name <- cell_type
    myRCTD_reps_L2@population_de_results[[cell_type]]$sig_cat <- 'non_sig'
}
#Add column to identify significant genes
for (cell_type in names(myRCTD_reps_L2@population_sig_gene_df)){
    if (dim(myRCTD_reps_L2@population_sig_gene_df[[cell_type]])[1] == 0){
        myRCTD_reps_L2@population_sig_gene_df[[cell_type]] <- NULL
    } else {
        myRCTD_reps_L2@population_sig_gene_df[[cell_type]]$cell_name <- cell_type
        myRCTD_reps_L2@population_sig_gene_df[[cell_type]]$sig_cat <- 'sig'
        }
}
#Add column to identify non-significant genes
for (cell_type in names(myRCTD_reps_L3@population_de_results)){
    myRCTD_reps_L3@population_de_results[[cell_type]]$cell_name <- cell_type
    myRCTD_reps_L3@population_de_results[[cell_type]]$sig_cat <- 'non_sig'
}
#Add column to identify significant genes
for (cell_type in names(myRCTD_reps_L3@population_sig_gene_df)){
    if (dim(myRCTD_reps_L3@population_sig_gene_df[[cell_type]])[1] == 0){
        myRCTD_reps_L3@population_sig_gene_df[[cell_type]] <- NULL
    } else {
        myRCTD_reps_L3@population_sig_gene_df[[cell_type]]$cell_name <- cell_type
        myRCTD_reps_L3@population_sig_gene_df[[cell_type]]$sig_cat <- 'sig'
        }
}
#Add column to identify non-significant genes
for (cell_type in names(myRCTD_reps_L4@population_de_results)){
    myRCTD_reps_L4@population_de_results[[cell_type]]$cell_name <- cell_type
    myRCTD_reps_L4@population_de_results[[cell_type]]$sig_cat <- 'non_sig'
}
#Add column to identify significant genes
for (cell_type in names(myRCTD_reps_L4@population_sig_gene_df)){
    if (dim(myRCTD_reps_L4@population_sig_gene_df[[cell_type]])[1] == 0){
        myRCTD_reps_L4@population_sig_gene_df[[cell_type]] <- NULL
    } else {
        myRCTD_reps_L4@population_sig_gene_df[[cell_type]]$cell_name <- cell_type
        myRCTD_reps_L4@population_sig_gene_df[[cell_type]]$sig_cat <- 'sig'
        }
}
#Add column to identify non-significant genes
for (cell_type in names(myRCTD_reps_L5@population_de_results)){
    myRCTD_reps_L5@population_de_results[[cell_type]]$cell_name <- cell_type
    myRCTD_reps_L5@population_de_results[[cell_type]]$sig_cat <- 'non_sig'
}
#Add column to identify significant genes
for (cell_type in names(myRCTD_reps_L5@population_sig_gene_df)){
    if (dim(myRCTD_reps_L5@population_sig_gene_df[[cell_type]])[1] == 0){
        myRCTD_reps_L5@population_sig_gene_df[[cell_type]] <- NULL
    } else {
        myRCTD_reps_L5@population_sig_gene_df[[cell_type]]$cell_name <- cell_type
        myRCTD_reps_L5@population_sig_gene_df[[cell_type]]$sig_cat <- 'sig'
        }
}
#Add column to identify non-significant genes
for (cell_type in names(myRCTD_reps_L6@population_de_results)){
    myRCTD_reps_L6@population_de_results[[cell_type]]$cell_name <- cell_type
    myRCTD_reps_L6@population_de_results[[cell_type]]$sig_cat <- 'non_sig'
}
#Add column to identify significant genes
for (cell_type in names(myRCTD_reps_L6@population_sig_gene_df)){
    if (dim(myRCTD_reps_L6@population_sig_gene_df[[cell_type]])[1] == 0){
        myRCTD_reps_L6@population_sig_gene_df[[cell_type]] <- NULL
    } else {
        myRCTD_reps_L6@population_sig_gene_df[[cell_type]]$cell_name <- cell_type
        myRCTD_reps_L6@population_sig_gene_df[[cell_type]]$sig_cat <- 'sig'
        }
}
#Add column to identify non-significant genes
for (cell_type in names(myRCTD_reps_WM@population_de_results)){
    myRCTD_reps_WM@population_de_results[[cell_type]]$cell_name <- cell_type
    myRCTD_reps_WM@population_de_results[[cell_type]]$sig_cat <- 'non_sig'
}
#Add column to identify significant genes
for (cell_type in names(myRCTD_reps_WM@population_sig_gene_df)){
    if (dim(myRCTD_reps_WM@population_sig_gene_df[[cell_type]])[1] == 0){
        myRCTD_reps_WM@population_sig_gene_df[[cell_type]] <- NULL
    } else {
        myRCTD_reps_WM@population_sig_gene_df[[cell_type]]$cell_name <- cell_type
        myRCTD_reps_WM@population_sig_gene_df[[cell_type]]$sig_cat <- 'sig'
        }
}

# for L1
combined <- data.frame()
for (cell_type in names(myRCTD_reps_L1@population_de_results)){
    if (cell_type %in% names(myRCTD_reps_L1@population_sig_gene_df)){
        assign(paste0(cell_type,'_non_sig'),tibble::rownames_to_column(myRCTD_reps_L1@population_de_results[[cell_type]],'genes'))
        assign(paste0(cell_type,'_sig'),tibble::rownames_to_column(myRCTD_reps_L1@population_sig_gene_df[[cell_type]],'genes'))
        assign(paste0(cell_type,'_non_sig'),anti_join(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig'))), by = "genes"))
        assign(cell_type, rbind.fill(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig')))))
    } else {
        assign(cell_type,tibble::rownames_to_column(myRCTD_reps_L1@population_de_results[[cell_type]],'genes'))
        }
    combined <- rbind.fill(combined, eval(as.name(cell_type)))
    write.csv(combined, 'vol_plot_data/L1/combined_across_samples.csv', row.names=FALSE)
}
# for L2
combined <- data.frame()
for (cell_type in names(myRCTD_reps_L2@population_de_results)){
    if (cell_type %in% names(myRCTD_reps_L2@population_sig_gene_df)){
        assign(paste0(cell_type,'_non_sig'),tibble::rownames_to_column(myRCTD_reps_L2@population_de_results[[cell_type]],'genes'))
        assign(paste0(cell_type,'_sig'),tibble::rownames_to_column(myRCTD_reps_L2@population_sig_gene_df[[cell_type]],'genes'))
        assign(paste0(cell_type,'_non_sig'),anti_join(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig'))), by = "genes"))
        assign(cell_type, rbind.fill(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig')))))
    } else {
        assign(cell_type,tibble::rownames_to_column(myRCTD_reps_L2@population_de_results[[cell_type]],'genes'))
        }
    combined <- rbind.fill(combined, eval(as.name(cell_type)))
    write.csv(combined, 'vol_plot_data/L2/combined_across_samples.csv', row.names=FALSE)
}
# for L3
combined <- data.frame()
for (cell_type in names(myRCTD_reps_L3@population_de_results)){
    if (cell_type %in% names(myRCTD_reps_L3@population_sig_gene_df)){
        assign(paste0(cell_type,'_non_sig'),tibble::rownames_to_column(myRCTD_reps_L3@population_de_results[[cell_type]],'genes'))
        assign(paste0(cell_type,'_sig'),tibble::rownames_to_column(myRCTD_reps_L3@population_sig_gene_df[[cell_type]],'genes'))
        assign(paste0(cell_type,'_non_sig'),anti_join(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig'))), by = "genes"))
        assign(cell_type, rbind.fill(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig')))))
    } else {
        assign(cell_type,tibble::rownames_to_column(myRCTD_reps_L3@population_de_results[[cell_type]],'genes'))
        }
    combined <- rbind.fill(combined, eval(as.name(cell_type)))
    write.csv(combined, 'vol_plot_data/L3/combined_across_samples.csv', row.names=FALSE)
}
# for L4
combined <- data.frame()
for (cell_type in names(myRCTD_reps_L4@population_de_results)){
    if (cell_type %in% names(myRCTD_reps_L4@population_sig_gene_df)){
        assign(paste0(cell_type,'_non_sig'),tibble::rownames_to_column(myRCTD_reps_L4@population_de_results[[cell_type]],'genes'))
        assign(paste0(cell_type,'_sig'),tibble::rownames_to_column(myRCTD_reps_L4@population_sig_gene_df[[cell_type]],'genes'))
        assign(paste0(cell_type,'_non_sig'),anti_join(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig'))), by = "genes"))
        assign(cell_type, rbind.fill(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig')))))
    } else {
        assign(cell_type,tibble::rownames_to_column(myRCTD_reps_L4@population_de_results[[cell_type]],'genes'))
        }
    combined <- rbind.fill(combined, eval(as.name(cell_type)))
    write.csv(combined, 'vol_plot_data/L4/combined_across_samples.csv', row.names=FALSE)
}
# for L5
combined <- data.frame()
for (cell_type in names(myRCTD_reps_L5@population_de_results)){
    if (cell_type %in% names(myRCTD_reps_L5@population_sig_gene_df)){
        assign(paste0(cell_type,'_non_sig'),tibble::rownames_to_column(myRCTD_reps_L5@population_de_results[[cell_type]],'genes'))
        assign(paste0(cell_type,'_sig'),tibble::rownames_to_column(myRCTD_reps_L5@population_sig_gene_df[[cell_type]],'genes'))
        assign(paste0(cell_type,'_non_sig'),anti_join(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig'))), by = "genes"))
        assign(cell_type, rbind.fill(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig')))))
    } else {
        assign(cell_type,tibble::rownames_to_column(myRCTD_reps_L5@population_de_results[[cell_type]],'genes'))
        }
    combined <- rbind.fill(combined, eval(as.name(cell_type)))
    write.csv(combined, 'vol_plot_data/L5/combined_across_samples.csv', row.names=FALSE)
}
# for L6
combined <- data.frame()
for (cell_type in names(myRCTD_reps_L6@population_de_results)){
    if (cell_type %in% names(myRCTD_reps_L6@population_sig_gene_df)){
        assign(paste0(cell_type,'_non_sig'),tibble::rownames_to_column(myRCTD_reps_L6@population_de_results[[cell_type]],'genes'))
        assign(paste0(cell_type,'_sig'),tibble::rownames_to_column(myRCTD_reps_L6@population_sig_gene_df[[cell_type]],'genes'))
        assign(paste0(cell_type,'_non_sig'),anti_join(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig'))), by = "genes"))
        assign(cell_type, rbind.fill(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig')))))
    } else {
        assign(cell_type,tibble::rownames_to_column(myRCTD_reps_L6@population_de_results[[cell_type]],'genes'))
        }
    combined <- rbind.fill(combined, eval(as.name(cell_type)))
    write.csv(combined, 'vol_plot_data/L6/combined_across_samples.csv', row.names=FALSE)
}
# for WM
combined <- data.frame()
for (cell_type in names(myRCTD_reps_WM@population_de_results)){
    if (cell_type %in% names(myRCTD_reps_WM@population_sig_gene_df)){
        assign(paste0(cell_type,'_non_sig'),tibble::rownames_to_column(myRCTD_reps_WM@population_de_results[[cell_type]],'genes'))
        assign(paste0(cell_type,'_sig'),tibble::rownames_to_column(myRCTD_reps_WM@population_sig_gene_df[[cell_type]],'genes'))
        assign(paste0(cell_type,'_non_sig'),anti_join(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig'))), by = "genes"))
        assign(cell_type, rbind.fill(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig')))))
    } else {
        assign(cell_type,tibble::rownames_to_column(myRCTD_reps_WM@population_de_results[[cell_type]],'genes'))
        }
    combined <- rbind.fill(combined, eval(as.name(cell_type)))
    write.csv(combined, 'vol_plot_data/WM/combined_across_samples.csv', row.names=FALSE)
}


