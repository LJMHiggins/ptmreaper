#' Calculate peptide ID False Discover Rate
#'
#' @param main.pescal.dir
#' @param decoy.pescal.dir
#' @param peptides
#'
#' @return
#' @export
#'
#' @examples
compute_FDR <- function(main.pescal.dir,
                        decoy.pescal.dir,
                        peptides){
  .calculate_FDR <- function(pep_id,
                             combipept.main,
                             combipept.decoy) {
    ## Word of warning - function variable cannot have same name as a column of the dataframe.
    peptide_max.scr <- combipept.main %>%
      filter(db_id == pep_id) %>%
      pull(max_scr)
    peptide_max.ppm <- combipept.main %>%
      filter(db_id == pep_id) %>%
      pull(max_ppm)

    spm.target.hits.scr <- nrow(combipept.main %>%
                                  filter(max_scr > peptide_max.scr))
    spm.decoy.hits.scr <- nrow(combipept.decoy %>%
                                 filter(max_scr > peptide_max.scr))
    fdr_max.scr <- (spm.decoy.hits.scr / spm.target.hits.scr) * 100

    spm.target.hits.ppm <- nrow(combipept.main %>%
                                  filter(max_ppm > peptide_max.ppm))
    spm.decoy.hits.ppm <- nrow(combipept.decoy %>%
                                 filter(max_ppm > peptide_max.ppm))
    fdr_max.ppm <- (spm.decoy.hits.ppm / spm.target.hits.ppm) * 100

    fdrs <- list(fdr.msc.scr = fdr_max.scr,
                 fdr.ppm = fdr_max.ppm)
    return(fdrs)
  }
  combipept.main <- ptmreaper::read_combipept(main.pescal.dir)
  combipept.decoy <- ptmreaper::read_combipept(decoy.pescal.dir)

  fdr.results <- foreach::foreach(i=1:nrow(peptides), .combine = "rbind") %do% {
    pep.for.calc <- peptides[i,1]
    site <- peptides[i,"site_label"]
    db_id <- peptides[i, "db_id"]
    fdrs <- .calculate_FDR(pep.for.calc,
                           combipept.main = combipept.main,
                           combipept.decoy = combipept.decoy)
    out <-data.frame(
      "db_id" = db_id,
      "Site"= site,
      "FDR max_scr" = fdrs[["fdr.msc.scr"]],
      "FDR max_ppm" = fdrs[["fdr.ppm"]])

    return(out)
  }
}


#' Plot mascot score ranks
#'
#' @param main.pescal.dir
#' @param decoy.pescal.dir
#'
#' @return
#' @export
#'
#' @examples
plot_mscr_ranks <- function(main.pescal.dir,
                            decoy.pescal.dir){
  combipept.main <- ptmreaper::read_combipept(main.pescal.dir)
  combipept.decoy <- ptmreaper::read_combipept(decoy.pescal.dir)

  plot.combi.df <- rbind(
    combipept.main %>%
      arrange(max_scr) %>%
      mutate(Database = rep("SwissProt", nrow(combipept.main))) %>%
      mutate(rank = scale(seq(1:nrow(combipept.main)), center=FALSE)),
    combipept.decoy %>%
      arrange(max_scr) %>%
      mutate(Database = rep("Decoy_SwissProt", nrow(combipept.decoy)))%>%
      mutate(rank = scale(seq(1:nrow(combipept.decoy)),center=FALSE))
  )

  ggplot(data = plot.combi.df, aes(x=rank, y=max_scr, color = Database)) +
    geom_point()+
    theme_bw()+
    theme(axis.text.x = element_blank(),
          text = element_text(size = 20))
}
