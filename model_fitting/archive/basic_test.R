# Load necessary libraries
library(lattice)
library(MASS)
require(pscl) # alternatively can use package ZIM for zero-inflated models
library(lmtest)
library(marginaleffects)
library(ggcorrplot)
library(performance)
library(see)
library(topmodels)
library(statmod)
library(assessor)
library(dplyr)
library(optparse)  # For command-line argument parsing

# Model functions
poisson_model <- function(data) {
  glm(formula_nonzi, family = 'poisson', data = data)
}

zero_inflated_poisson_model <- function(data) {
  zeroinfl(formula, dist = 'poisson', data = data)
}

negative_binomial_model <- function(data) {
  glm.nb(formula_nonzi, data = data)
}

zero_inflated_negative_binomial_model <- function(data) {
  zeroinfl(formula, dist = 'negbin', data = data)
}

read_data <- function(dataset) {
  if (dataset == 'cohen_multi') {
    data <- read.csv("/home/mcn26/palmer_scratch/tabula_data/formatted/COHEN_MIXED_CELL.tsv", sep='\t', header = T)[-1]
  } else if (dataset == 'cohen_retina') {
  data <- read.csv("/home/mcn26/palmer_scratch/tabula_data/formatted/COHEN_RETINA.tsv", sep='\t', header = T)[-1]
  data <- data[data$reads_transfection_BC>0, ]
  }
  return(data)
}

# Parse arguments
option_list <- list(
  make_option(c("--dataset"), type = "character"),
  make_option(c("--model_choice"), type = "character"),
  make_option(c("--out_file"), type = "character")
)

args <- parse_args(OptionParser(option_list = option_list))

# Main function
main <- function(args) {

  data <- read_data(args$dataset)



  model_functions <- list(
    poisson = poisson_model,
    zi_poisson = zero_inflated_poisson_model,
    negative_binomial = negative_binomial_model,
    zi_negative_binomial = zero_inflated_negative_binomial_model
  )

  if (!args$model_choice %in% names(model_functions)) {
    stop("Invalid model choice")
  }

  model_function <- model_functions[[args$model_choice]]
  fit_model <- model_function(data)

  # write out model file
  saveRDS(fit_model, file = args$out_file)

}

# Run the main function
formula <- UMIs_MPRA_BC ~ cell_type_annotation + rep_id + CRE_id | rep_id
formula_nonzi <- UMIs_MPRA_BC ~ cell_type_annotation + rep_id + CRE_id
main(args)
