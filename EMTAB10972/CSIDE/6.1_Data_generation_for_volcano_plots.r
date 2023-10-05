library(plyr)
library(dplyr)
library(tibble)
library(spacexr)

setwd('/work/project/ext_014/EMTAB10972/CSIDE')

#Load CISDE results
myRCTD_reps_nc <- readRDS('saved_rds/myRCTD_reps_nc_im.rds')
myRCTD_reps_dm <- readRDS('saved_rds/myRCTD_reps_dm_im.rds')
myRCTD_reps_r <- readRDS('saved_rds/myRCTD_reps_r_im.rds')
myRCTD_reps_t <- readRDS('saved_rds/myRCTD_reps_t_im.rds')
myRCTD_reps_cp <- readRDS('saved_rds/myRCTD_reps_cp_im.rds')
myRCTD_reps_e <- readRDS('saved_rds/myRCTD_reps_e_im.rds')
myRCTD_reps_m <- readRDS('saved_rds/myRCTD_reps_m_im.rds')

#Add column to identify non-significant genes & cell type
for (cell_type in names(myRCTD_reps_nc@population_de_results)){
    myRCTD_reps_nc@population_de_results[[cell_type]]$cell_name <- cell_type
    myRCTD_reps_nc@population_de_results[[cell_type]]$sig_cat <- 'non_sig'
}
#Add column to identify significant genes & cell type
for (cell_type in names(myRCTD_reps_nc@population_sig_gene_df)){
    if (dim(myRCTD_reps_nc@population_sig_gene_df[[cell_type]])[1] == 0){
        myRCTD_reps_nc@population_sig_gene_df[[cell_type]] <- NULL
    } else {
        myRCTD_reps_nc@population_sig_gene_df[[cell_type]]$cell_name <- cell_type
        myRCTD_reps_nc@population_sig_gene_df[[cell_type]]$sig_cat <- 'sig'
        }
}
#Add column to identify non-significant genes & cell type
for (cell_type in names(myRCTD_reps_dm@population_de_results)){
    myRCTD_reps_dm@population_de_results[[cell_type]]$cell_name <- cell_type
    myRCTD_reps_dm@population_de_results[[cell_type]]$sig_cat <- 'non_sig'
}
#Add column to identify significant genes & cell type
for (cell_type in names(myRCTD_reps_dm@population_sig_gene_df)){
    if (dim(myRCTD_reps_dm@population_sig_gene_df[[cell_type]])[1] == 0){
        myRCTD_reps_dm@population_sig_gene_df[[cell_type]] <- NULL
    } else {
        myRCTD_reps_dm@population_sig_gene_df[[cell_type]]$cell_name <- cell_type
        myRCTD_reps_dm@population_sig_gene_df[[cell_type]]$sig_cat <- 'sig'
        }
}
#Add column to identify non-significant genes & cell type
for (cell_type in names(myRCTD_reps_r@population_de_results)){
    myRCTD_reps_r@population_de_results[[cell_type]]$cell_name <- cell_type
    myRCTD_reps_r@population_de_results[[cell_type]]$sig_cat <- 'non_sig'
}
#Add column to identify significant genes & cell type
for (cell_type in names(myRCTD_reps_r@population_sig_gene_df)){
    if (dim(myRCTD_reps_r@population_sig_gene_df[[cell_type]])[1] == 0){
        myRCTD_reps_r@population_sig_gene_df[[cell_type]] <- NULL
    } else {
        myRCTD_reps_r@population_sig_gene_df[[cell_type]]$cell_name <- cell_type
        myRCTD_reps_r@population_sig_gene_df[[cell_type]]$sig_cat <- 'sig'
        }
}
#Add column to identify non-significant genes & cell type
for (cell_type in names(myRCTD_reps_t@population_de_results)){
    myRCTD_reps_t@population_de_results[[cell_type]]$cell_name <- cell_type
    myRCTD_reps_t@population_de_results[[cell_type]]$sig_cat <- 'non_sig'
}
#Add column to identify significant genes & cell type
for (cell_type in names(myRCTD_reps_t@population_sig_gene_df)){
    if (dim(myRCTD_reps_t@population_sig_gene_df[[cell_type]])[1] == 0){
        myRCTD_reps_t@population_sig_gene_df[[cell_type]] <- NULL
    } else {
        myRCTD_reps_t@population_sig_gene_df[[cell_type]]$cell_name <- cell_type
        myRCTD_reps_t@population_sig_gene_df[[cell_type]]$sig_cat <- 'sig'
        }
}
#Add column to identify non-significant genes & cell type
for (cell_type in names(myRCTD_reps_cp@population_de_results)){
    myRCTD_reps_cp@population_de_results[[cell_type]]$cell_name <- cell_type
    myRCTD_reps_cp@population_de_results[[cell_type]]$sig_cat <- 'non_sig'
}
#Add column to identify significant genes & cell type
for (cell_type in names(myRCTD_reps_cp@population_sig_gene_df)){
    if (dim(myRCTD_reps_cp@population_sig_gene_df[[cell_type]])[1] == 0){
        myRCTD_reps_cp@population_sig_gene_df[[cell_type]] <- NULL
    } else {
        myRCTD_reps_cp@population_sig_gene_df[[cell_type]]$cell_name <- cell_type
        myRCTD_reps_cp@population_sig_gene_df[[cell_type]]$sig_cat <- 'sig'
        }
}
#Add column to identify non-significant genes & cell type
for (cell_type in names(myRCTD_reps_e@population_de_results)){
    myRCTD_reps_e@population_de_results[[cell_type]]$cell_name <- cell_type
    myRCTD_reps_e@population_de_results[[cell_type]]$sig_cat <- 'non_sig'
}
#Add column to identify significant genes & cell type
for (cell_type in names(myRCTD_reps_e@population_sig_gene_df)){
    if (dim(myRCTD_reps_e@population_sig_gene_df[[cell_type]])[1] == 0){
        myRCTD_reps_e@population_sig_gene_df[[cell_type]] <- NULL
    } else {
        myRCTD_reps_e@population_sig_gene_df[[cell_type]]$cell_name <- cell_type
        myRCTD_reps_e@population_sig_gene_df[[cell_type]]$sig_cat <- 'sig'
        }
}
#Add column to identify non-significant genes & cell type
for (cell_type in names(myRCTD_reps_m@population_de_results)){
    myRCTD_reps_m@population_de_results[[cell_type]]$cell_name <- cell_type
    myRCTD_reps_m@population_de_results[[cell_type]]$sig_cat <- 'non_sig'
}
#Add column to identify significant genes & cell type
for (cell_type in names(myRCTD_reps_m@population_sig_gene_df)){
    if (dim(myRCTD_reps_m@population_sig_gene_df[[cell_type]])[1] == 0){
        myRCTD_reps_m@population_sig_gene_df[[cell_type]] <- NULL
    } else {
        myRCTD_reps_m@population_sig_gene_df[[cell_type]]$cell_name <- cell_type
        myRCTD_reps_m@population_sig_gene_df[[cell_type]]$sig_cat <- 'sig'
        }
}

# for nc
combined <- data.frame()
for (cell_type in names(myRCTD_reps_nc@population_de_results)){
    if (cell_type %in% names(myRCTD_reps_nc@population_sig_gene_df)){
        assign(paste0(cell_type,'_non_sig'),tibble::rownames_to_column(myRCTD_reps_nc@population_de_results[[cell_type]],'genes'))
        assign(paste0(cell_type,'_sig'),tibble::rownames_to_column(myRCTD_reps_nc@population_sig_gene_df[[cell_type]],'genes'))
        assign(paste0(cell_type,'_non_sig'),anti_join(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig'))), by = "genes"))
        assign(cell_type, rbind.fill(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig')))))
    } else {
        assign(cell_type,tibble::rownames_to_column(myRCTD_reps_nc@population_de_results[[cell_type]],'genes'))
        }
    combined <- rbind.fill(combined, eval(as.name(cell_type)))
    write.csv(combined, 'vol_plot_data/nc/nc_combined_im.csv', row.names=FALSE)
}
#for dm
combined <- data.frame()
for (cell_type in names(myRCTD_reps_dm@population_de_results)){
    if (cell_type %in% names(myRCTD_reps_dm@population_sig_gene_df)){
        assign(paste0(cell_type,'_non_sig'),tibble::rownames_to_column(myRCTD_reps_dm@population_de_results[[cell_type]],'genes'))
        assign(paste0(cell_type,'_sig'),tibble::rownames_to_column(myRCTD_reps_dm@population_sig_gene_df[[cell_type]],'genes'))
        assign(paste0(cell_type,'_non_sig'),anti_join(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig'))), by = "genes"))
        assign(cell_type, rbind.fill(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig')))))
    } else {
        assign(cell_type,tibble::rownames_to_column(myRCTD_reps_dm@population_de_results[[cell_type]],'genes'))
        }
    combined <- rbind.fill(combined, eval(as.name(cell_type)))
    write.csv(combined, 'vol_plot_data/dm/dm_combined_im.csv', row.names=FALSE)
}
#for r
combined <- data.frame()
for (cell_type in names(myRCTD_reps_r@population_de_results)){
    if (cell_type %in% names(myRCTD_reps_r@population_sig_gene_df)){
        assign(paste0(cell_type,'_non_sig'),tibble::rownames_to_column(myRCTD_reps_r@population_de_results[[cell_type]],'genes'))
        assign(paste0(cell_type,'_sig'),tibble::rownames_to_column(myRCTD_reps_r@population_sig_gene_df[[cell_type]],'genes'))
        assign(paste0(cell_type,'_non_sig'),anti_join(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig'))), by = "genes"))
        assign(cell_type, rbind.fill(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig')))))
    } else {
        assign(cell_type,tibble::rownames_to_column(myRCTD_reps_r@population_de_results[[cell_type]],'genes'))
        }
    combined <- rbind.fill(combined, eval(as.name(cell_type)))
    write.csv(combined, 'vol_plot_data/r/r_combined_im.csv', row.names=FALSE)
}
#for t
combined <- data.frame()
for (cell_type in names(myRCTD_reps_t@population_de_results)){
    if (cell_type %in% names(myRCTD_reps_t@population_sig_gene_df)){
        assign(paste0(cell_type,'_non_sig'),tibble::rownames_to_column(myRCTD_reps_t@population_de_results[[cell_type]],'genes'))
        assign(paste0(cell_type,'_sig'),tibble::rownames_to_column(myRCTD_reps_t@population_sig_gene_df[[cell_type]],'genes'))
        assign(paste0(cell_type,'_non_sig'),anti_join(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig'))), by = "genes"))
        assign(cell_type, rbind.fill(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig')))))
    } else {
        assign(cell_type,tibble::rownames_to_column(myRCTD_reps_t@population_de_results[[cell_type]],'genes'))
        }
    combined <- rbind.fill(combined, eval(as.name(cell_type)))
    write.csv(combined, 'vol_plot_data/t/t_combined_im.csv', row.names=FALSE)
}
#for cp
combined <- data.frame()
for (cell_type in names(myRCTD_reps_cp@population_de_results)){
    if (cell_type %in% names(myRCTD_reps_cp@population_sig_gene_df)){
        assign(paste0(cell_type,'_non_sig'),tibble::rownames_to_column(myRCTD_reps_cp@population_de_results[[cell_type]],'genes'))
        assign(paste0(cell_type,'_sig'),tibble::rownames_to_column(myRCTD_reps_cp@population_sig_gene_df[[cell_type]],'genes'))
        assign(paste0(cell_type,'_non_sig'),anti_join(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig'))), by = "genes"))
        assign(cell_type, rbind.fill(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig')))))
    } else {
        assign(cell_type,tibble::rownames_to_column(myRCTD_reps_cp@population_de_results[[cell_type]],'genes'))
        }
    combined <- rbind.fill(combined, eval(as.name(cell_type)))
    write.csv(combined, 'vol_plot_data/cp/cp_combined_im.csv', row.names=FALSE)
}
#for m
combined <- data.frame()
for (cell_type in names(myRCTD_reps_m@population_de_results)){
    if (cell_type %in% names(myRCTD_reps_m@population_sig_gene_df)){
        assign(paste0(cell_type,'_non_sig'),tibble::rownames_to_column(myRCTD_reps_m@population_de_results[[cell_type]],'genes'))
        assign(paste0(cell_type,'_sig'),tibble::rownames_to_column(myRCTD_reps_m@population_sig_gene_df[[cell_type]],'genes'))
        assign(paste0(cell_type,'_non_sig'),anti_join(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig'))), by = "genes"))
        assign(cell_type, rbind.fill(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig')))))
    } else {
        assign(cell_type,tibble::rownames_to_column(myRCTD_reps_m@population_de_results[[cell_type]],'genes'))
        }
    combined <- rbind.fill(combined, eval(as.name(cell_type)))
    write.csv(combined, 'vol_plot_data/m/m_combined_im.csv', row.names=FALSE)
}
#for e
combined <- data.frame()
for (cell_type in names(myRCTD_reps_e@population_de_results)){
    if (cell_type %in% names(myRCTD_reps_e@population_sig_gene_df)){
        assign(paste0(cell_type,'_non_sig'),tibble::rownames_to_column(myRCTD_reps_e@population_de_results[[cell_type]],'genes'))
        assign(paste0(cell_type,'_sig'),tibble::rownames_to_column(myRCTD_reps_e@population_sig_gene_df[[cell_type]],'genes'))
        assign(paste0(cell_type,'_non_sig'),anti_join(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig'))), by = "genes"))
        assign(cell_type, rbind.fill(eval(as.name(paste0(cell_type,'_non_sig'))), eval(as.name(paste0(cell_type,'_sig')))))
    } else {
        assign(cell_type,tibble::rownames_to_column(myRCTD_reps_e@population_de_results[[cell_type]],'genes'))
        }
    combined <- rbind.fill(combined, eval(as.name(cell_type)))
    write.csv(combined, 'vol_plot_data/e/e_combined_im.csv', row.names=FALSE)
}


