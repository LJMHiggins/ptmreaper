# Peptide ANOVA

peptide_anova <- function(){
  #loop over cells
  column_names <- c("Df1", "Df2", "Sum Sq1", "Sum Sq2", "Mean Sq1", "Mean Sq2",
                    "F value1", "F value2", "Pr(>F)1", "Pr(>F)2", "Cell.id", "Prot.id")

  # Create an empty dataframe with specified column names
  anova_results <- data.frame(matrix(ncol = length(column_names), nrow = 0))
  colnames(anova_results) <- column_names

  prot_anova_sub <- subset(prot_anova, treatment.1 != "Veh_HCL")

  for (cell in unique(prot_anova_sub$cell.id)){
    df_cell_sub <- subset(prot_anova_sub, cell.id == cell)
    for (prot in unique(df_cell_sub$protein.group)){
      df_sub_x <- subset(df_cell_sub, protein.group == prot)

      aov_result <- aov(Prot_quantity ~ condition, data = df_sub_x)
      aov_result_sum <- unlist(summary(aov_result))
      aov_result_sum["Cell.id"] <- cell
      aov_result_sum["prot.id"] <- prot

      aov_result_sum <- as.data.frame(t(aov_result_sum))

      anova_results <- rbind(anova_results, aov_result_sum)

    }
  }

  anova_results$qvalue <- p.adjust(anova_results$`Pr(>F)1`, method = "fdr")
}
