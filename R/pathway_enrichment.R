
#' Pathway enrichment of modified peptides
#'
#' @description
#' Pathway enrichment of modified peptides. Note, this approach uses protein id
#' taht modified peptide is assignmed to. It is therefore a proxy of the type of
#' proteins being modified, and not necessarily representative of the functional
#' consequences of the modification.
#'
#'
#' @param comparisons List of comparison data.frames from protools2 limma comparison.
#' @param combipept Combipept data from pescal output. Requires whole sheet.
#' @param databases Databases to use in pathway enrichment search.
#' @param graph.heading String. Title for graphs.
#'
#' @return List of pathway enrichment data from protools2::pathway_enrichment function.
#' @export
#'
#' @examples
ptmprot_pathway_enrichment <- function(comparisons,
                                       combipept,
                                       databases,
                                       graph.heading){
  ## Read in proteins from combipept - to match peptide ids to proteins ##
  combipept_filt <- combipept[29:31]
  colnames(combipept_filt) <- c("protein", "protein.protein", "db_id")

  .modpep_enrichment <- function(comp_df, combipept_filt, databases) {
    hits_df <- comp_df[1:5]
    hits_df %>%
      rename(protein = acc_no) -> hits_df

    hits_df2 <- merge(hits_df, combipept_filt, by="db_id")

    hits_df2$db_id <- hits_df2$protein.y

    hits_df2 %>%
      rename(protein = db_id) %>%
      select(1:4)-> hits_df3

    diff_peps <- protools2::identify_differences_in_comparison_plus_volcano(hits_df3,fold.cutoff = 0, qval.cutoff = 0.25)

    increased_aps <- diff_peps$df.increased$sites
    decreased_aps <- diff_peps$df.decreased$sites

    background_aps <- diff_peps

    pe <- protools2::pathway_enrichment(increased.peptides = increased_aps,
                                        decreased.peptides = decreased_aps,
                                        background.data = hits_df3,
                                        graph.heading = graph.heading,
                                        prot_dbs = databases)
    return(pe)
  }
  ###### process search #####

  all_pe_results <- list()
  for (n in names(comparisons)) {
    tryCatch({
      pe <- .modpep_enrichment(comp_df = comparisons[[n]],
                               combipept_filt = combipept_filt,
                               databases = databases)

      all_pe_results[[n]] <- pe
    }, error=function(e) all_pe_results[[n]] <- NA
    )
  }
  return(all_pe_results)
}
