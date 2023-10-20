
#' PTM reap
#'
#' @description
#' Wrapper function for reading combipept data, output areas, filtering modified
#' peptides and return normalised results (if selected)
#'
#'
#' @param ptm_pattern String. Ptm type to retrieve from combipept data. Must be
#' exact match.
#' @param pescal.path Path to pescal output file.
#' @param normalise Logical. Whether to normalise or return raw values.
#' @param delta_score_cutoff Numeric. Max_delta_score cut off value for ptm peptides.
#'
#' @return
#' @export
#'
#' @examples
ptm_reap <- function(ptm_pattern,
                     pescal.path,
                     normalise = TRUE,
                     delta_score_cutoff = 10){
  ptm_pep <- ptmreaper::read_combipept(pescal.path = pescal.path) %>%
    subset(grepl(x = pep_mod, pattern = ptm_pattern, fixed = T)) %>% # subset on acetyl
    subset(max_delta_score > delta_score_cutoff)

  output_areas <- ptmreaper::read_output_areas(pescal.path = pescal.path)

  output_areas_ptm <- output_areas %>%
    subset(output_areas$db_id %in% ptm_pep$db_id)

  if (normalise == TRUE){
    results <- ptmreaper::ptm_norm(output_areas_ptm)
    return(results)
  } else {
    return(output_areas_ptm)
  }
}

#' Read combipept data
#'
#' @description
#' Function to read in combipept sheet from pescal output.
#'
#' @param pescal.path Path with location of of pescal output file.
#'
#' @return Dataframe with combipept data with "acc_no", "protein", "peptide",
#' "pep_mod", "...25", "max_scr", "mean_scr", "max_delta_score", "xmod_pos" and
#' "db_id" columns.
#'
#' @export
#'
#' @examples
read_combipept <- function(pescal.path){
  # Function for reading combipept data from pescal output
  # Combipept sheet is for extracting ids for modified peptides
  require(dplyr)
  df <-
    readxl::read_excel(pescal.path, sheet = "combiPeptData") %>% subset(
      select = c(
        "acc_no",
        "protein",
        "peptide",
        "pep_mod",
        "...25",
        "max_ppm",
        "max_scr",
        "mean_scr",
        "max_delta_score",
        "xmod_pos",
        "db_id"
      ))
  return(df)
}

#' Read output areas data
#'
#' @param pescal.path Path with location of of pescal output file.
#'
#' @return Dataframe with output areas sheet from pescal output file.
#' @export
#'
#' @examples
read_output_areas <- function(pescal.path){
  # Function for reading output areasd data from pescal output
  require(dplyr)
  df <-
    readxl::read_excel(pescal.path, sheet = "output_areas")
  return(df)
}

#' Normalise ptm peptide output areas
#'
#' @param output_areas_df Dataframe. Output areas data from pescal output.
#'
#' @return List of dataframes with normalised data. Note, recent focus
#' has shifted to df.norm.log2.centered.na.impute.v2.
#'
#' @export
#'
#' @examples
ptm_norm <- function(output_areas_df){
  # Function includes scaling so enter unscaled data
  cols <- colnames(dplyr::select_if(output_areas_df, is.numeric))
  df.areas.n <- data.frame(ids=output_areas_df$db_id,
                           scale(output_areas_df[,cols],center = F,
                                 scale =  colSums(output_areas_df[,cols])))

  cols <- colnames(dplyr::select_if(df.areas.n, is.numeric))
  df.areas.n[df.areas.n == 0] <- NA
  df.areas.n[,cols] <- df.areas.n[,cols]*1000000
  df.norm <- df.areas.n
  ### Copied - ids changed ###

  df.norm.log2.centered <- data.frame(ids=df.norm$ids,
                                      scale(log2(df.norm[,cols]),scale = F))

  df.norm.log2.centered.scaled <- data.frame(ids=df.norm$ids,
                                             scale(log2(df.norm[,cols])))

  df.norm.log2.centered.scaled.na.imputed <- df.norm.log2.centered.scaled
  df.norm.log2.centered.scaled.na.imputed[is.na(df.norm.log2.centered.scaled.na.imputed)] <- min(df.norm.log2.centered.scaled.na.imputed[,cols], na.rm = T)/5

  ## Second imputation method - due to error of log2 norm / 5
  df.norm.log2.centered.na.impute.v2 <- df.norm.log2.centered
  df.norm.log2.centered.na.impute.v2[cols] <- lapply(
    df.norm.log2.centered.na.impute.v2[cols], function(col){
      replace(col, is.na(col), min(col, na.rm = TRUE) - 1)
    }
  )
  return(list(normalized.data=df.norm,
              normalized.plus.log2.cent.data = df.norm.log2.centered,
              df.norm.log2.centered.na.impute.v2 = df.norm.log2.centered.na.impute.v2,
              normalized.plus.log2.cent.scaled.data = df.norm.log2.centered.scaled,
              df.norm.log2.centered.scaled.na.imputed = df.norm.log2.centered.scaled.na.imputed))
}
