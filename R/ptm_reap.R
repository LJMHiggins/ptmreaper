
ptm_reap <- function(ptm_pattern, pescal.path, normalise = TRUE){
  ptm_pep <- read_combipept(x = pescal.path) %>%
    subset(grepl(x = pep_mod, pattern = ptm_pattern, fixed = T)) %>% # subset on acetyl
    subset(max_delta_score > 10)

  output_areas <- read_output_areas(x = pescal.path)

  output_areas_ptm <- output_areas %>%
    subset(output_areas$db_id %in% acetyl_pep$db_id)

  if (normalise == TRUE){
    results <- ptm_norm(output_areas_ptm)
    return(results)
  } else {
    return(output_areas_ptm)
  }


}

read_combipept <- function(x){
  # Function for reading combipept data from pescal output
  # Combipept sheet is for extracting ids for modified peptides
  require(dplyr)
  df <-
    readxl::read_excel(x, sheet = "combiPeptData") %>% subset(
      select = c(
        "acc_no",
        "protein",
        "peptide",
        "pep_mod",
        "...25",
        "max_scr",
        "mean_scr",
        "max_delta_score",
        "xmod_pos",
        "db_id"
      ))
  return(df)
}

read_output_areas <- function(x){
  # Function for reading output areasd data from pescal output
  require(dplyr)
  df <-
    readxl::read_excel(x, sheet = "output_areas")
  return(df)
}

ptm_norm <- function(df){
  # Function includes scaling so enter unscaled data
  cols <- colnames(dplyr::select_if(df, is.numeric))
  df.areas.n <- data.frame(ids=df$db_id,
                           scale(df[,cols],center = F,
                                 scale =  colSums(df[,cols])))

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

  return(list(normalized.data=df.norm,
              normalized.plus.log2.cent.data=df.norm.log2.centered,
              normalized.plus.log2.cent.scaled.data=df.norm.log2.centered.scaled,
              df.norm.log2.centered.scaled.na.imputed=df.norm.log2.centered.scaled.na.imputed))
}
