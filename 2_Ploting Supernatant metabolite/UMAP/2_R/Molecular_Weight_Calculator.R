# Molecular Weight Calculator (unfinished)

get_mol_weight <- function(element){
  # define molecular weight
  weight <- c(C = 12, H = 1, O = 16, N = 14, S = 32, P = 31)
  return(weight[element])
}
formula_count <- function(formula){
  elements <- str_extract_all(formula, "[A-Z][a-z]*\\d*")[[1]]
  elements <- str_extract_all(elements, "[A-Z][a-z]*|\\d+")
  elements <- unlist(elements)
  
  max_length <- max(sum(!grepl("\\d+", elements)), sum(grepl("\\d+", elements)))
  df <- data.frame(
    Element = rep(NA, max_length),
    Count = rep(NA, max_length)
  )
  
  
  df$Element[1:sum(!grepl("\\d+", elements))] <- elements[!grepl("\\d+", elements)]
  df$Count[1:sum(grepl("\\d+", elements))] <- as.numeric(elements[grepl("\\d+", elements)])
  df$Count[is.na(df$Count | df$Count == "")] <- 1
  return(df)
  df$Element <- sapply(df$Element, get_mol_weight)
  df$Element <- as.numeric(df$Element)
  df$Count <- as.numeric(df$Count)
  df$weight <- df$Element * df$Count
  return(sum(df$weight))
}