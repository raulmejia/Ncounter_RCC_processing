relocate_row_names <- function(somedf,somecolname){
  rownames(somedf) <- somedf[, somecolname]
  resultdf <- somedf[ , - which( colnames(somedf) %in% somecolname) ]
  return(resultdf)
}

relocate_count_norm_samplesID <- function(anotherdf){
  resultdf <- as.data.frame(anotherdf[,"count_norm_vec"] )
  rownames( resultdf) <- rownames( anotherdf)
  colnames( resultdf ) <- as.character(anotherdf[1, "count_norm_samplesID"])
  return(resultdf)
}
relocate_count_samplesID <- function(anotherdf){
  resultdf <- as.data.frame(anotherdf[,"count_vec"] )
  rownames( resultdf) <- rownames( anotherdf)
  colnames( resultdf ) <- as.character(anotherdf[1, "count_samplesID"])
  return(resultdf)
}

NachoNorm2matrix <- function(NachoObj){
  GeneNames <- NachoObj$nacho$Name
  count_norm_vec <- NachoObj$nacho$Count_Norm
  count_norm_samplesID <- NachoObj$nacho[ , myIDcolname]
  Long_df_norm <- data.frame( GeneNames ,
                              count_norm_samplesID , 
                              count_norm_vec )
  Long_df_norm$count_norm_samplesID <- as.factor( Long_df_norm$count_norm_samplesID )
  List_dfs_norm <- split(Long_df_norm, Long_df_norm$count_norm_samplesID)
  List_dfs_norm_rownames <- lapply( List_dfs_norm, relocate_row_names, "GeneNames" )
  
  List_dfs_norm_rownanmes_colnames <- lapply(List_dfs_norm_rownames , relocate_count_norm_samplesID )
  
  result_mat <- do.call( cbind, List_dfs_norm_rownanmes_colnames)
  result_mat <- as.matrix( result_mat )
  return( result_mat)
}

Nacho_Orig_count_2matrix <- function(NachoObj){
  GeneNames <- NachoObj$nacho$Name
  count_vec <- NachoObj$nacho$Count
  count_samplesID <- NachoObj$nacho[ , myIDcolname]
  Long_df <- data.frame( GeneNames ,
                              count_samplesID , 
                              count_vec )
  Long_df$count_samplesID <- as.factor( Long_df$count_samplesID )
  List_dfs <- split(Long_df, Long_df$count_samplesID)
  List_dfs_rownames <- lapply( List_dfs, relocate_row_names, "GeneNames" )
  
  List_dfs_rownanmes_colnames <- lapply(List_dfs_rownames , relocate_count_samplesID )
  
  result_mat <- do.call( cbind, List_dfs_rownanmes_colnames)
  result_mat <- as.matrix( result_mat )
  return( result_mat)
}
