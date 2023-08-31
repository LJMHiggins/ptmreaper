
#' Perform protools2 limma and label output
#'
#' @param df.design
#' @param df.comparisons
#' @param df.normalized.areas
#' @param combipept_ptm
#'
#' @return
#' @export
#'
#' @examples
labelled_protools2_limma <- function(df.design,
                                     df.comparisons,
                                     df.normalized.areas,
                                     combipept_ptm){

  comparisons <- protools2::compare.conditions.by.limma(df.design = df.design,
                                                        df.comparisons = df.comparisons,
                                                        df.normalized.areas = df.normalized.areas)
  # Create modification site label, with prot, location, and type.
  combipept_ptm$site_label <- paste0(combipept_ptm$...25, " - ", gsub(" ", "", combipept_ptm$pep_mod))
  for (i in 1:length(comparisons)) {
    comparisons[[i]] %>%
      rename(db_id = protein) -> comparisons[[i]]
    comparisons[[i]] <- merge(comparisons[[i]], combipept_ptm, by.x = "db_id")
  }
  return(comparisons)
}

#' Title
#'
#' @param comparisons List of limma comparisons. Generated from "labelled_protools2_limma".
#'
#' @return
#' @export
#'
#' @examples
plot_all_volcano <- function(comparisons){

  vol_list <- vector(mode = "list", length(comparisons))

  for (i in 1:length(comparisons)) {

    comparisons[[i]]$diffexpressed <- "NO"
    comparisons[[i]]$diffexpressed[comparisons[[i]]$FDR < 0.05] <- "YES"

    comparisons[[i]]$delabel <- NA
    comparisons[[i]]$delabel[comparisons[[i]]$diffexpressed != "NO"] <- comparisons[[i]]$site_label[comparisons[[i]]$diffexpressed != "NO"]
    plo <- ggplot(data=comparisons[[i]], aes(x=difference.test.vs.control, y=-log10(FDR), label=delabel)) +
      geom_point(alpha=0.5) +
      theme_minimal() +
      geom_text_repel( max.overlaps=10)+
      geom_vline(xintercept=c(-0.5, 0.5), col="red") +
      geom_hline(yintercept=-log10(0.25), col="Blue")+
      geom_hline(yintercept=-log10(0.05), col="red") +
      ggtitle(names(comparisons[i]))+
      theme(plot.title = element_text(hjust = 0.5))
    print(plo)

    vol_list[[i]] <- protools2::identify_differences_in_comparison_plus_volcano(
      comparisons[[i]],fold.cutoff = 0, qval.cutoff = 0.25)
  }
  return(vol_list)
}

