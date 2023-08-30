summarise_iqtree <- function(file_path) {
  # open iqtree file
  iqtree_text <- readLines(file_path)
  # abstract iqtree info
  iqtree_summarise <- list()

  # the best model and its criteria
  line_best_model <- grep("Best-fit model according to", iqtree_text)
  iqtree_summarise$model_crireia <- sub(".*to\\s+(.*):.*", "\\1", iqtree_text[line_best_model])
  iqtree_summarise$best_model <- sub(".*:\\s*(.*)", "\\1", iqtree_text[line_best_model])
  line_outgroup <- grep("Tree is UNROOTED", iqtree_text)
  iqtree_summarise$outgroup <- sub(".*'\\s*(.*)'.*", "\\1", iqtree_text[line_outgroup])
  rm(line_best_model, line_outgroup)

  # Find the line numbers for the substitution model and likelihood tree
  line_sub <- grep("SUBSTITUTION PROCESS", iqtree_text)
  line_likelihood <- grep("MAXIMUM LIKELIHOOD TREE", iqtree_text)
  
  # Extract the substitution model information
  model_info <- iqtree_text[(line_sub+3):(line_likelihood-1)]
  model_substitution <- sub(".*: (.*)", "\\1", model_info[1])
  model_components <- unlist(strsplit(model_substitution, "\\+"))
  
  para_num <- length(model_components)
  iqtree_summarise$rate_type <- "None"
  if (grepl("R|G", model_components[para_num])){
    if ((grepl("I", model_components[para_num - 1]))){
      iqtree_summarise$rate_type <- paste(model_components[(para_num - 1):para_num], collapse = "+")
    } else {
      iqtree_summarise$rate_type <- model_components[para_num]
    }
  } else if ((grepl("I", model_components[para_num]))){
    iqtree_summarise$rate_type <- model_components[para_num]
  }
  
  # Check if the model includes a rate heterogeneity mixture
  if (grepl("(m|M)ixture", model_info[1])){
    # Find the line numbers for the rate heterogeneity and site proportion and rates sections
    line_RHAS <-  grep("Model of rate heterogeneity:", model_info)
    
    # Extract the rate heterogeneity and site proportion and rates information
    iqtree_summarise$model_info <- read.table(text = model_info[3:(line_RHAS-1)], header = TRUE)
    if (grepl("R", model_components[length(model_components)])){
      line_table <- grep("Category", model_info)
      block_SPR <- model_info[line_table:(length(model_info) - 1)]
      iqtree_summarise$RHAS <- read.table(text = block_SPR, header = TRUE)
      rm(line_table, block_SPR)
    } else if (grepl("G", model_components[length(model_components)])){
      line_table <- grep("Category", model_info)
      block_SPR <- model_info[line_table:(length(model_info) - 2)]
      iqtree_summarise$RHAS <- read.table(text = block_SPR, header = TRUE)
      iqtree_summarise$G_alpha <-  as.numeric(sub(".*: ([0-9.]+)", "\\1", model_info[line_table - 2]))
      rm(line_table, block_SPR)
    } else { 
      iqtree_summarise$RHAS <- NULL
    }
    
  } else {
    # Find the line numbers for the rate parameter R, state frequencies, rate matrix Q, rate heterogeneity, and site proportion and rates sections
    line_R <- grep("Rate parameter R", model_info)
    line_F <- grep("State frequencies", model_info)
    line_Q <- grep("Rate matrix Q", model_info)
    line_RHAS <- grep("Model of rate heterogeneity:", model_info)
    
    # Extract the rate parameter R, state frequencies, and rate matrix Q information
    block_R <- model_info[(line_R+2):(line_F-2)]
    block_Q  <- model_info[(line_Q+2):(line_RHAS-1)]
    
    iqtree_summarise$R <- as.numeric(sub(".*: ([0-9.]+)", "\\1", block_R))
    names(iqtree_summarise$R) <- sub("(.*): [0-9.]+", "\\1", block_R)
    
    if (grepl("equal frequencies", model_info[line_F])){
      iqtree_summarise$F <- rep(1/4, 4)
    } else {
      block_F  <- model_info[(line_F+2):(line_Q-2)]
      iqtree_summarise$F <- as.numeric(sub(".*= ([0-9.]+)", "\\1", block_F))
      rm(block_F)
    }
    names(iqtree_summarise$F) <- c("A", "C", "G", "T")
    
    iqtree_summarise$Q <- read.table(text = block_Q, header = FALSE, row.names = 1)
    colnames(iqtree_summarise$Q) <- rownames(iqtree_summarise$Q)
    
    # Create a parameter string for the model
    parameter <- paste0(model_components[1], "{", paste(as.character(iqtree_summarise$R), collapse = ","),  "}+",
                        model_components[2], "{", paste(as.character(iqtree_summarise$F), collapse = ","),  "}")
    
    # Store the model information in a data frame
    iqtree_summarise$model_info <- data.frame("No" = 1,
                                              "Component" = model_components[1],
                                              "Rate" = 1.0000,
                                              "Weight" = 1.0000,
                                              "Parameters" = parameter)
    
    
    # Extract the rate heterogeneity and site proportion and rates information
    if (grepl("R", model_components[length(model_components)])){
      line_table <- grep("Category", model_info)
      block_SPR <- model_info[line_table:(length(model_info) - 1)]
      iqtree_summarise$RHAS <- read.table(text = block_SPR, header = TRUE)
      rm(line_table, block_SPR)
    } else if (grepl("G", model_components[length(model_components)])){
      line_table <- grep("Category", model_info)
      block_SPR <- model_info[line_table:(length(model_info) - 2)]
      iqtree_summarise$RHAS <- read.table(text = block_SPR, header = TRUE)
      iqtree_summarise$G_alpha <-  as.numeric(sub(".*: ([0-9.]+)", "\\1", model_info[line_table - 2]))
      rm(line_table, block_SPR)
    } else { 
      iqtree_summarise$RHAS <- NULL
    }
    
    # Remove unnecessary variables
    rm(line_sub, line_R, line_F, line_Q, line_RHAS)
    rm(block_R, block_Q)
  }
  # tree length
  line_length <- grep("Total tree length", iqtree_text)
  iqtree_summarise$total_tree_length <- as.numeric(sub(".*: ([0-9.]+)", "\\1", iqtree_text[line_length]))
  iqtree_summarise$internal_tree_length <- as.numeric(sub(".*: ([0-9.]+).*", "\\1", iqtree_text[line_length + 1]))
  rm(line_length)

  # statistics

  iqtree_summarise$log_likelihood <- as.numeric(sub(".*: ([0-9.-]+) .*", "\\1", iqtree_text[line_likelihood + 3]))
  iqtree_summarise$log_likelihood_sd <- as.numeric(sub(".*: ([0-9.-]+) \\(s.e. ([0-9.]+)\\)", "\\2", iqtree_text[line_likelihood + 3]))
  iqtree_summarise$log_likelihood_Unconstrained <- as.numeric(sub(".*: ([0-9.-]+)", "\\1", iqtree_text[line_likelihood + 4]))
  iqtree_summarise$num_free_params <- as.numeric(sub(".*: ([0-9]+)", "\\1", iqtree_text[line_likelihood + 5]))
  iqtree_summarise$AIC <- as.numeric(sub(".*: ([0-9.-]+)", "\\1", iqtree_text[line_likelihood + 6]))
  iqtree_summarise$AICc <- as.numeric(sub(".*: ([0-9.-]+)", "\\1", iqtree_text[line_likelihood + 7]))
  iqtree_summarise$BIC <- as.numeric(sub(".*: ([0-9.-]+)", "\\1", iqtree_text[line_likelihood + 8]))
  # iqtree_summarise$num_params <- sub(".*: ([0-9]+)", "\\1", iqtree_text[line_likelihood + 11])
  # iqtree_summarise$samp_size <- sub(".*: ([0-9]+)", "\\1", iqtree_text[line_likelihood + 12])
  rm(line_likelihood)

  return(iqtree_summarise)
}

short_description <- function(summarise_list) {
  cat("Best model: ", summarise_list$best_model, "\n")
  cat("Tree length: ", summarise_list$total_tree_length, "(inner: ", summarise_list$internal_tree_length, ")\n", sep = "")
  cat("Log-likelihood: ", summarise_list$log_likelihood, "\n")
  cat("Number of free parameters: ", summarise_list$num_free_params, "\n")
  cat("AIC, AICc, BIC: ", summarise_list$AIC, summarise_list$AICc, summarise_list$BIC, "\n")
  cat("Model parameters: \n")
  print(summarise_list$model_info$Parameters)
  cat("Rate heterogeneity: \n")
  if (nrow(summarise_list$RHAS) > 0) {
    cat("    Rate:", summarise_list$RHAS$Relative_rate, "\n")
    cat("    Prop:", summarise_list$RHAS$Proportion, "\n")
  } else {
    cat("    None\n")
  }
}

short_description_output <- function(sum_one, sum_mix, filename) {
  # Open a connection to the output file
  output_file <- file(filename, open = "wt")
  
  # Define a helper function to write to both the console and the output file
  write_output <- function(...) {
    cat(..., "\n")
    cat(..., "\n", file = output_file)
  }
  
  # Write the output for the first summary
  write_output("The phylogenic analysis result for one model:")
  write_output("Best model: ", sum_one$best_model)
  write_output("Tree length: ", sum_one$total_tree_length, "(inner: ", sum_one$internal_tree_length, ")", sep = "")
  write_output("Log-likelihood: ", sum_one$log_likelihood)
  write_output("Number of free parameters: ", sum_one$num_free_params)
  write_output("AIC, AICc, BIC: ", sum_one$AIC, sum_one$AICc, sum_one$BIC)
  write_output("Model parameters: ")
  print(sum_one$model_info$Parameters)
  cat("Model parameters: \n", file = output_file)
  write.table(sum_one$model_info$Parameters, file = output_file, row.names = FALSE)
  
  write_output("Rate heterogeneity: ")
  if (nrow(sum_one$RHAS) > 0) {
    write_output("    Rate:", sum_one$RHAS$Relative_rate)
    write_output("    Prop:", sum_one$RHAS$Proportion)
  } else {
    write_output("    None")
  }
  
  # Write a separator
  write_output("---------------------------------------------")
  
  # Write the output for the second summary
  write_output("The phylogenic analysis result for mix model:")
  write_output("Best model: ", sum_mix$best_model)
  write_output("Tree length: ", sum_mix$total_tree_length, "(inner: ", sum_mix$internal_tree_length, ")", sep = "")
  write_output("Log-likelihood: ", sum_mix$log_likelihood)
  write_output("Number of free parameters: ", sum_mix$num_free_params)
  write_output("AIC, AICc, BIC: ", sum_mix$AIC, sum_mix$AICc, sum_mix$BIC)
  write_output("Model parameters: ")
  print(sum_mix$model_info$Parameters)
  cat("Model parameters: \n", file = output_file)
  write.table(sum_mix$model_info$Parameters, file = output_file, row.names = FALSE)
  
  write_output("Rate heterogeneity: ")
  if (nrow(sum_mix$RHAS) > 0) {
    write_output("    Rate:", sum_mix$RHAS$Relative_rate)
    write_output("    Prop:", sum_mix$RHAS$Proportion)
  } else {
    write_output("    None")
  }
  
  # Close the connection to the output file
  close(output_file)
}

summary_table <- function(sum_model){
  data.frame("Model_string" = sum_model$best_model,
             "Rates" = sum_model$rate_type,
             "Likelihood" = sum_model$log_likelihood,
             "Unconstrained_likelihood" = sum_model$log_likelihood_Unconstrained,
             "parameters" = paste(sum_model$model_info$Parameters, collapse = ";"),
             "AIC" = sum_model$AIC,
             "AICc" = sum_model$AICc,
             "BIC" = sum_model$BIC)
}

