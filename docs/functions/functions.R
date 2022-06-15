#!/usr/bin/Rscript

#extract_samples <- function(x, regex_samples){
  # in our case, regex is "\\d{4}"
 # out <- str_extract(x, regex_samples)
  # returns list of sample names
  #return(out)
#}

extract_sample_names <- function(x, regex_samples = "\\d{4}"){
  require(stringr)
  out <- str_extract(x, regex_samples)
  # returns list of sample names
  return(out)
}

read_in_files_metadata_initial_tally <- function(file_list, metadata, regex_samples = "\\d{4}"){
  require(DT)
  require(stringr)
  require(dplyr)
  require(tidyr)
  #assuming files in table format and have a header - should be the case for neofuse output anyway
  data_lst = lapply(X = file_list, FUN = read.table, header = T)
  #print(data_lst)
  # in our case, regex is "\\d{4}", extract sample names from file list
  Sample <- unlist(lapply(X = file_list, FUN = extract_sample_names))
  #print(Sample)
  # create a new list with the sample ID
  new_list <- list()
  for (i in 1:length(data_lst)){
    new_list[[i]] <- cbind(data_lst[[i]], Sample=rep(Sample[i], nrow(data_lst[[i]])))
    #print(paste("length of list:", length(new_list)))
    #print(new_list)
  }
  # make this all into one big data frame
  df_complete <- data.frame()
  for (i in 1:length(new_list)){
    df_complete <- rbind.data.frame(df_complete, new_list[[i]])
  }
  # get rid of numbers in brackets
  df_complete$Gene1 <- str_replace_all(df_complete$Gene1, pattern = "\\(\\d++\\)", replacement = "")
  df_complete$Gene2 <- str_replace_all(df_complete$Gene2, pattern = "\\(\\d++\\)", replacement = "")
  # ensure sample column is integer type
  #df_complete$Sample <- as.numeric(levels(df_complete$Sample))[df_complete$Sample]
  df_complete$Sample <- as.integer(df_complete$Sample)
  # In metadata file, replace Tumour with CAF and Normal with TAN
  metadata$Condition <- str_replace(metadata$Condition, "Tumour", "CAF")
  metadata$Condition <- str_replace(metadata$Condition, "Normal", "TAN")
  # ensure Sample column is equal for joining 
  colnames(metadata)[1] <- "Sample"
  # complete join with metadata
  neofuse_results_with_metadata <- df_complete %>% full_join(metadata, by = "Sample")
  output_datatable <- DT::datatable(neofuse_results_with_metadata, rownames = F, width = 1500, height = 1000, options=list(scrollX=T, scrollY = 800))
  print(paste("Total number of neoantigens identified without filtering: ", sum(!is.na(neofuse_results_with_metadata$Fusion)), sep = ""))
  # drop any rows with NA under Fusion column
  neofuse_results_with_metadata_without_na <- neofuse_results_with_metadata %>% drop_na("Fusion")
  print(paste("Total number of neoantigens identified in CAF samples without filtering: ", 
              nrow(neofuse_results_with_metadata_without_na %>% filter(Condition == "CAF")), sep = ""))
  print(paste("Total number of neoantigens identified in TAN samples without filtering: ", 
              nrow(neofuse_results_with_metadata_without_na %>% filter(Condition == "TAN")), sep = ""))
  outputs_list <- list("DataTable" = output_datatable, "Initial_Results" = neofuse_results_with_metadata_without_na)
  return(outputs_list)
}



neoantigens_unique_CAF <- function(neofuse_results_one_dataframe){
#take in full dataframe of results across samples, return dataframe with Fusion peptides that are unique to the CAFs
require(tidyr)
require(dplyr)
#produce one tibble per patient for filtering
    neoantigens_per_patient_split <- neofuse_results_one_dataframe %>% 
                                  group_by(Patient) %>% 
                                  group_split()
    #print(neoantigens_per_patient_split)
    output_dataframe <- data.frame()
    for (i in 1:length(neoantigens_per_patient_split)){
      group_dataframe <- neoantigens_per_patient_split[[i]]
      # group by condition
      neoantigens_unique_caf_df <- group_dataframe[group_dataframe$Condition == "CAF" & 
                                                     !(group_dataframe$Fusion_Peptide %in% group_dataframe$Fusion_Peptide[group_dataframe$Condition == "TAN"]),]
                                                     #group_dataframe$Fusion_Peptide %in% group_dataframe$Fusion_Peptide[group_dataframe$Condition == "TAN"],]
      #print(neoantigens_unique_caf_df)
      output_dataframe <- rbind.data.frame(output_dataframe, neoantigens_unique_caf_df)
    }
    output_dataframe_without_na <- output_dataframe %>% drop_na("Fusion")
    print(paste("Total number of neoantigens identified in CAF samples only: ", 
                nrow(output_dataframe_without_na), sep = ""))
    return(output_dataframe_without_na)
}

neoantigens_unique_CAF_gene_level <- function(neofuse_results_one_dataframe){
  #take in full dataframe of results across samples, return dataframe with Fusion peptides that are unique to the CAFs
  require(tidyr)
  require(dplyr)
  #produce one tibble per patient for filtering
  neoantigens_per_patient_split <- neofuse_results_one_dataframe %>% 
    group_by(Patient) %>% 
    group_split()
  #print(neoantigens_per_patient_split)
  output_dataframe <- data.frame()
  for (i in 1:length(neoantigens_per_patient_split)){
    group_dataframe <- neoantigens_per_patient_split[[i]]
    # group by condition
    neoantigens_unique_caf_df <- group_dataframe[group_dataframe$Condition == "CAF" & 
                                                   !(group_dataframe$Fusion %in% group_dataframe$Fusion[group_dataframe$Condition == "TAN"]),]
    #group_dataframe$Fusion_Peptide %in% group_dataframe$Fusion_Peptide[group_dataframe$Condition == "TAN"],]
    #print(neoantigens_unique_caf_df)
    output_dataframe <- rbind.data.frame(output_dataframe, neoantigens_unique_caf_df)
  }
  output_dataframe_without_na <- output_dataframe %>% drop_na("Fusion")
  print(paste("Total number of neoantigens identified in CAF samples only: ", 
              nrow(output_dataframe_without_na), sep = ""))
  return(output_dataframe_without_na)
}

neoantigens_unique_CAF_gene_level_fusion_level <- function(neofuse_results_one_dataframe, fusions_called){
  #take in full dataframe of results across samples, return dataframe with Fusion peptides that are unique to the CAFs
  require(tidyr)
  require(dplyr)
  #produce one tibble per patient for filtering
  fusions_per_patient_split <- fusions_called %>% 
    group_by(Patient) %>% 
    group_split()
  #print(neoantigens_per_patient_split)
  fusions_unique_caf <- data.frame()
 # output_dataframe <- data.frame()
  for (i in 1:length(fusions_per_patient_split)){
    group_dataframe <- fusions_per_patient_split[[i]]
    # group by condition
    #filter out repeated gene1
    fusions_unique_caf_df <- group_dataframe[group_dataframe$Condition == "CAF" & 
                                                   !(group_dataframe$gene1 %in% group_dataframe$gene1[group_dataframe$Condition == "TAN"]),]
    #filter out repeated gene2
    fusions_unique_caf_df <- group_dataframe[group_dataframe$Condition == "CAF" & 
                                               !(group_dataframe$gene2 %in% group_dataframe$gene2[group_dataframe$Condition == "TAN"]),]
    #group_dataframe$Fusion_Peptide %in% group_dataframe$Fusion_Peptide[group_dataframe$Condition == "TAN"],]
    #print(neoantigens_unique_caf_df)
    fusions_unique_caf <- rbind.data.frame(fusions_unique_caf, fusions_unique_caf_df)
    #print(fusions_unique_caf)
  }
  fusions_unique_without_na <- fusions_unique_caf %>% drop_na(gene1, gene2)
  # rename gene1, gene2 cols
  x <- which(colnames(fusions_unique_without_na) == "gene1")
  y <- which(colnames(fusions_unique_without_na) == "gene2")
  colnames(fusions_unique_without_na)[x] <- "Gene1"
  colnames(fusions_unique_without_na)[y] <- "Gene2"
  fusions_unique_caf <- semi_join(neofuse_results_one_dataframe, fusions_unique_without_na, by = c("Gene1", "Gene2", "Patient", "Sample"))
  #print(paste("Total number of neoantigens identified in CAF samples only: ", 
   #           nrow(fusions_unique_without_na), sep = ""))
  #only keep fusions found in the unique df
  
  return(fusions_unique_caf)
}

# MAIN FUNCTION
read_in_fusions_combine_filter_tan <- function(file_list, metadata, regex_samples = "\\d{4}"){
  require(DT)
  require(stringr)
  require(dplyr)
  require(tidyr)
  #assuming files in table format and have a header - should be the case for neofuse output anyway
  data_lst = lapply(X = file_list, FUN = read.table, header = T, comment.char="", check.names = F)
  # in our case, regex is "\\d{4}", extract sample names from file list
  Sample <- unlist(lapply(X = file_list, FUN = extract_sample_names))
  # create a new list with the sample ID
  new_list <- list()
  for (i in 1:length(data_lst)){
    colnames(data_lst[[i]])[1] <- "gene1"
    new_list[[i]] <- cbind(data_lst[[i]], Sample=rep(Sample[i], nrow(data_lst[[i]])))
  }
  # make this all into one big data frame
  df_complete <- data.frame()
  for (i in 1:length(new_list)){
    df_complete <- rbind.data.frame(df_complete, new_list[[i]])
  }
  # get rid of numbers in brackets
  df_complete$gene1 <- str_replace_all(df_complete$gene1, pattern = "\\(\\d++\\)", replacement = "")
  df_complete$gene2 <- str_replace_all(df_complete$gene2, pattern = "\\(\\d++\\)", replacement = "")
  
  # ensure sample column is integer type
  df_complete$Sample <- as.integer(df_complete$Sample)
  
  # ensure Sample column is equal for joining 
  colnames(metadata)[1] <- "Sample"
  # complete join with metadata
  fusions_with_metadata <- df_complete %>% full_join(metadata, by = "Sample")
  fusions_per_patient_split <- fusions_with_metadata %>% 
    group_by(Patient) %>% 
    group_split()
  #print(neoantigens_per_patient_split)
  fusions_unique_caf <- data.frame()
  # output_dataframe <- data.frame()
  for (i in 1:length(fusions_per_patient_split)){
    group_dataframe <- fusions_per_patient_split[[i]]
    # group by condition
    #filter out repeated gene1
    fusions_unique_caf_df <- group_dataframe[group_dataframe$Condition == "CAF" & 
                                               !(group_dataframe$gene1 %in% group_dataframe$gene1[group_dataframe$Condition == "TAN"]),]
    #filter out repeated gene2
    fusions_unique_caf_df <- group_dataframe[group_dataframe$Condition == "CAF" & 
                                               !(group_dataframe$gene2 %in% group_dataframe$gene2[group_dataframe$Condition == "TAN"]),]
    #group_dataframe$Fusion_Peptide %in% group_dataframe$Fusion_Peptide[group_dataframe$Condition == "TAN"],]
    #print(neoantigens_unique_caf_df)
    fusions_unique_caf <- rbind.data.frame(fusions_unique_caf, fusions_unique_caf_df)
    #print(fusions_unique_caf)
  }
  fusions_unique_without_na <- fusions_unique_caf %>% drop_na(gene1, gene2)
  # drop any rows with NA under Fusion column
  x <- which(colnames(fusions_unique_without_na) == "gene1")
  y <- which(colnames(fusions_unique_without_na) == "gene2")
  colnames(fusions_unique_without_na)[x] <- "Gene1"
  colnames(fusions_unique_without_na)[y] <- "Gene2"
  #neofuse_results_with_metadata_without_na <- neofuse_results_with_metadata %>% drop_na("Fusion")
  return(fusions_unique_without_na)
}

# SAME FUNCTION AS ABOVE BUT LOOK FOR FUSIONS IN TANS NOT CAF

read_in_fusions_combine_filter_only_keep_tan_fusions <- function(file_list, metadata, regex_samples = "\\d{4}"){
  file_list <- files_to_read
  metadata <- metadata
  library(DT)
  library(stringr)
  library(dplyr)
  library(tidyr)
  #assuming files in table format and have a header - should be the case for neofuse output anyway
  data_lst = lapply(X = file_list, FUN = read.table, header = T, comment.char="", check.names = F)
  #print(data_lst)
  # in our case, regex is "\\d{4}", extract sample names from file list
  Sample <- unlist(lapply(X = file_list, FUN = extract_sample_names))
  # create a new list with the sample ID
  new_list <- list()
  for (i in 1:length(data_lst)){
    colnames(data_lst[[i]])[1] <- "Fusion"
    new_list[[i]] <- cbind(data_lst[[i]], Sample=rep(Sample[i], nrow(data_lst[[i]])))
    #print(paste("length of list:", length(new_list)))
    #print(new_list)
  }
  # make this all into one big data frame
  df_complete <- data.frame()
  for (i in 1:length(new_list)){
    df_complete <- rbind.data.frame(df_complete, new_list[[i]])
  }
  # get rid of numbers in brackets
  df_complete$Fusion <- str_replace_all(df_complete$Fusion, pattern = "\\(\\d++\\)", replacement = "")
  df_complete$Gene1 <- str_replace_all(df_complete$Gene1, pattern = "\\(\\d++\\)", replacement = "")
  df_complete$Gene2 <- str_replace_all(df_complete$Gene2, pattern = "\\(\\d++\\)", replacement = "")
  
  # ensure sample column is integer type
  #df_complete$Sample <- as.numeric(levels(df_complete$Sample))[df_complete$Sample]
  df_complete$Sample <- as.integer(df_complete$Sample)
  # In metadata file, replace Tumour with CAF and Normal with TAN
  metadata$Condition <- str_replace(metadata$Condition, "Tumour", "CAF")
  metadata$Condition <- str_replace(metadata$Condition, "Normal", "TAN")
  # ensure Sample column is equal for joining 
  colnames(metadata)[1] <- "Sample"
  # complete join with metadata
  fusions_with_metadata <- df_complete %>% full_join(metadata, by = "Sample")
  fusions_per_patient_split <- fusions_with_metadata %>% 
    group_by(Patient) %>% 
    group_split()
  #print(neoantigens_per_patient_split)
  fusions_unique_caf <- data.frame()
  # output_dataframe <- data.frame()
  for (i in 1:length(fusions_per_patient_split)){
    group_dataframe <- fusions_per_patient_split[[i]]
    # group by condition
    #filter out repeated gene1
    fusions_unique_caf_df <- group_dataframe[group_dataframe$Condition == "TAN" & 
                                               !(group_dataframe$Gene1 %in% group_dataframe$Gene1[group_dataframe$Condition == "CAF"]),]
    #filter out repeated gene2
    fusions_unique_caf_df <- group_dataframe[group_dataframe$Condition == "TAN" & 
                                               !(group_dataframe$Gene2 %in% group_dataframe$Gene2[group_dataframe$Condition == "CAF"]),]
    #group_dataframe$Fusion_Peptide %in% group_dataframe$Fusion_Peptide[group_dataframe$Condition == "TAN"],]
    #print(neoantigens_unique_caf_df)
    fusions_unique_caf <- rbind.data.frame(fusions_unique_caf, fusions_unique_caf_df)
    #print(fusions_unique_caf)
  }
  fusions_unique_without_na <- fusions_unique_caf %>% drop_na(Gene1, Gene2)
  # drop any rows with NA under Fusion column
  x <- which(colnames(fusions_unique_without_na) == "gene1")
  y <- which(colnames(fusions_unique_without_na) == "gene2")
  colnames(fusions_unique_without_na)[x] <- "Gene1"
  colnames(fusions_unique_without_na)[y] <- "Gene2"
  #neofuse_results_with_metadata_without_na <- neofuse_results_with_metadata %>% drop_na("Fusion")
  return(fusions_unique_without_na)
}

# Make function to take in CAF secreted factors table from Wu et al paper and return factors specific to a particular cancer
filter_by_cell_type <- function(x, Cell_type_selection, return_df = F){
  require(stringr)
  require(dplyr)
  split_object <- str_split(x$`Cell type`, pattern = ",")
  df <- data.frame()
  df_small <- data.frame()
  num <- which(colnames(x) == "Cell type")
  for (i in 1:length(split_object)){
    Cell_type <- split_object[[i]]
    if (length(split_object[[i]]) > 1){
      z <- x[i,-c(num)]
      reps <- length(Cell_type)
      z_repeated <- z %>% slice(rep(1:n(), each = reps))
      df_small <- cbind.data.frame(z_repeated, Cell_type)
    } else {
      z <- x[i,-c(num)]
      reps <- 1
      df_small <- cbind.data.frame(z, Cell_type)
    }
    df <- rbind.data.frame(df,df_small)
  }
  out <- df %>% filter(Cell_type_selection == Cell_type)
  rownames(out) <- NULL
  if (return_df == T){
    out <- list(out, df)
  }
  return(out)
}

# MAIN FUNCTION FINAL
read_in_fusions_combine_filter_out <- function(file_list, metadata, regex_samples = "\\d{4}", keep = "CAF", filter_out = "TAN", Confidence = c("high", "medium")){
  require(stringr)
  require(dplyr)
  require(tidyr)
  #assuming files in table format and have a header - should be the case for neofuse output anyway
  data_lst = lapply(X = file_list, FUN = read.table, header = T, comment.char="", check.names = F)
  # in our case, regex is "\\d{4}", extract sample names from file list
  Sample <- unlist(lapply(X = file_list, FUN = extract_sample_names))
  new_list <- list()
  for (i in 1:length(data_lst)){
    #colnames(data_lst[[i]])[1] <- "gene1"
    new_list[[i]] <- cbind(data_lst[[i]], Sample=rep(Sample[i], nrow(data_lst[[i]])))
  }
  
  # make this all into one big data frame
  df_complete <- data.frame()
  for (i in 1:length(new_list)){
    df_complete <- rbind.data.frame(df_complete, new_list[[i]])
  }
  print(paste("total number of fusion neoantigens before filtering: ", nrow(df_complete), sep = " "))
  # get rid of numbers in brackets
  df_complete$Fusion <- str_replace_all(df_complete$Fusion, pattern = "\\(\\d++\\)", replacement = "")
  df_complete$Gene1 <- str_replace_all(df_complete$Gene1, pattern = "\\(\\d++\\)", replacement = "")
  df_complete$Gene2 <- str_replace_all(df_complete$Gene2, pattern = "\\(\\d++\\)", replacement = "")
  # ensure sample column is integer type
  df_complete$Sample <- as.integer(df_complete$Sample)
  # ensure Sample column is equal for joining 
  colnames(metadata)[1] <- "Sample"
  # complete join with metadata
  fusions_with_metadata <- df_complete %>% full_join(metadata, by = "Sample")
  fusions_per_patient_split <- fusions_with_metadata %>% 
    group_by(Patient) %>% 
    group_split()
  fusions_unique_keep <- data.frame()
  # output_dataframe <- data.frame()
  for (i in 1:length(fusions_per_patient_split)){
    group_dataframe <- fusions_per_patient_split[[i]]
    # filter out repeated fusion
    fusions_unique_keep_df <- group_dataframe[group_dataframe$Condition == keep & 
                                                !(group_dataframe$Fusion %in% group_dataframe$Fusion[group_dataframe$Condition == filter_out]),]
    #filter out repeated gene1
    #fusions_unique_keep_df <- group_dataframe[group_dataframe$Condition == keep & 
     #                                           !(group_dataframe$Gene1 %in% group_dataframe$Gene1[group_dataframe$Condition == filter_out]),]
    #filter out repeated gene2
    #fusions_unique_keep_df <- group_dataframe[group_dataframe$Condition == keep & 
     #                                           !(group_dataframe$Gene2 %in% group_dataframe$Gene2[group_dataframe$Condition == filter_out]),]
    fusions_unique_keep <- rbind.data.frame(fusions_unique_keep, fusions_unique_keep_df)
    #print(fusions_unique_keep)
  }
  fusions_unique_desired_confidence <- fusions_unique_keep[(fusions_unique_keep$Confidence %in% Confidence),]
  outputs <- list(fusions_with_metadata = fusions_with_metadata, fusions_unique_desired_confidence = fusions_unique_desired_confidence)
  return(outputs)
}


read_in_fusions_combine <- function(file_list, metadata, regex_samples = "\\d{4}"){
  require(DT)
  require(stringr)
  require(dplyr)
  require(tidyr)
  #assuming files in table format and have a header - should be the case for neofuse output anyway
  data_lst = lapply(X = file_list, FUN = read.table, header = T, comment.char= "", check.names = F)
  #print(data_lst)
  # in our case, regex is "\\d{4}", extract sample names from file list
  Sample <- unlist(lapply(X = file_list, FUN = extract_sample_names))
  # create a new list with the sample ID
  new_list <- list()
  for (i in 1:length(data_lst)){
    #colnames(data_lst[[i]])[1] <- "gene1"
    new_list[[i]] <- cbind(data_lst[[i]], Sample=rep(Sample[i], nrow(data_lst[[i]])))
    #print(paste("length of list:", length(new_list)))
    #print(new_list)
  }
  # make this all into one big data frame
  df_complete <- data.frame()
  for (i in 1:length(new_list)){
    df_complete <- rbind.data.frame(df_complete, new_list[[i]])
  }
  # ensure sample column is integer type
  df_complete$Sample <- as.integer(df_complete$Sample)
  # ensure Sample column is equal for joining 
  colnames(metadata)[1] <- "Sample"
  # ensure all colnames are correct for downstream processing, capitalise them
  gene1_colnum <- which(colnames(df_complete) == "#gene1")
  colnames(df_complete)[gene1_colnum] <- "Gene1"
  gene2_colnum <- which(colnames(df_complete) == "gene2")
  colnames(df_complete)[gene2_colnum] <- "Gene2"
  bpoint1_colnum <- which(colnames(df_complete) == "breakpoint1")
  colnames(df_complete)[bpoint1_colnum] <- "Breakpoint1"
  bpoint2_colnum <- which(colnames(df_complete) == "breakpoint2")
  colnames(df_complete)[bpoint2_colnum] <- "Breakpoint2"
  # complete join with metadata
  fusions_with_metadata <- df_complete %>% full_join(metadata, by = "Sample")
  return(fusions_with_metadata)
}

