### This is to genenerate defattr file to open in ChimeraX to visualise the functional scores in 3D structure
setwd("/PATH/to/base_dir")
load("pos_hn.RData")
load("pos_fg.RData")
load("pos_fp.RData") 

####
write_defattr <- function(df, column_name, file_name, attribute_name) { 
  dna1 <- c(
    paste0("\t", "/A:", df[!is.na(df[[column_name]]), "Pos"], "\t", df[!is.na(df[[column_name]]), column_name]),
    paste0("\t", "/B:", df[!is.na(df[[column_name]]), "Pos"], "\t", df[!is.na(df[[column_name]]), column_name])
  )
  
  # Combine the text components with the user-defined attribute name
  write_text <- c(paste0("attribute: ", attribute_name, "\t"), "recipient: residues\t", dna1)
  
  # Write to file
  write.table(write_text, file = file_name, quote = FALSE, row.names = FALSE, col.names = FALSE)
}

# Writing three files with different columns and user-defined attribute names
write_defattr(pos_hn, "median_mean_scaled_ss_fitness", "hnf1a_dbd_scaled_ss_fitness.defattr", "scaled_ss_fitness")
write_defattr(pos_hn, "median_scaled_fitness_b1h", "hnf1a_dbd_scaled_fitness_b1h.defattr", "scaled_fitness_b1h")
write_defattr(pos_hn, "median_residual_scaled_BA", "hnf1a_dbd_scaled_residual.defattr", "scaled_residual_BA")
write_defattr(pos_hn, "abs_median_residual_scaled_BA", "hnf1a_abs_dbd_scaled_residual.defattr", "abs_median_residual_scaled_BA")
write_defattr(pos_hn, "median_scal_resi_Decayresidual", "hnf1a_abs_dbd_distance_corrected.defattr", "median_scal_resi_Decayresidual")

### for foxg1, foxp1, it is Chain C.  
# for foxg1 
write_defattr <- function(df, column_name, file_name, attribute_name) {
  # Create the dna1 vector
  dna1 <- c(
    paste0("\t", "/C:", df[!is.na(df[[column_name]]), "Pos"], "\t", df[!is.na(df[[column_name]]), column_name]),
  )
  
  # Combine the text components with the user-defined attribute name
  write_text <- c(paste0("attribute: ", attribute_name, "\t"), "recipient: residues\t", dna1)
  
  # Write to file
  write.table(write_text, file = file_name, quote = FALSE, row.names = FALSE, col.names = FALSE)
}
# Writing three files with different columns and user-defined attribute names
write_defattr(pos_fg, "median_abundance", "fg_dbd_scaled_ss_fitness.defattr", "scaled_ss_fitness")
write_defattr(pos_fg, "median_binding", "fg_dbd_scaled_fitness_b1h.defattr", "scaled_fitness_b1h")
write_defattr(pos_fg, "median_residual", "fg_dbd_scaled_residual.defattr", "scaled_residual_BA")
write_defattr(pos_fg, "median_abs_residual", "fg_abs_dbd_scaled_residual.defattr", "abs_median_residual_scaled_BA")

write_defattr(pos_fp, "median_abundance", "fp_dbd_scaled_ss_fitness.defattr", "scaled_ss_fitness")
write_defattr(pos_fp, "median_binding", "fp_dbd_scaled_fitness_b1h.defattr", "scaled_fitness_b1h")
write_defattr(pos_fp, "median_residual", "fp_dbd_scaled_residual.defattr", "scaled_residual_BA")
write_defattr(pos_fp, "median_abs_residual", "fp_abs_dbd_scaled_residual.defattr", "abs_median_residual_scaled_BA")
