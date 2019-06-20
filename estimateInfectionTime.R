#!/usr/bin/Rscript

# This script produces estimates of infection time based on
#  - identify_founders.tab
#  - optionally: bounds.csv
#  - optionally: viral_loads.csv

# The user needs to specify which models to use
# - Two different structures are available:
# -- The full model using viral loads, bounds and the relevant estimate
# -- The mutation rate calibration model which only scales the estimate
# - Three different timepoint/region combinations can be considered:
# -- 1 month nflg
# -- 1 month v3
# -- 6 month with pooled v3 and nflg
# - The user needs to specify which estimate to use:
# -- PFitter
# -- Synonymous PFitter
# -- Multifounder PFitter
# -- Multifounder Synonymous PFitter
# -- The other estimators are not supported by this script, but the user can modify this script to use them.

# An important note on file formats:
# Information from the three different files must be linked together.
# The identify_founders.tab is indexed by the name of the input fasta file.
# The viral load and bounds files MUST be indexed by THE SAME unique patient identifier that occurs in the name of the input file to the pipeline that will end up as the index for the identify_founders.tab file. The names of the files are scanned with regular expressions for the patient identifier, so care must be taken to prevent spurious matches. For example, if a patient identifier is just the number 1, then you can't have an identifier 10 since the 1 will match 10 when using regular expression. This can be prevented by left padding the ids with 0's or prepending a character string to them.
# For the Full model, the identifiers will be read from the bounds file. If the identifier cannot be found in the index columns of both the viral load and identify_founders.tab files, no timing estimate will be produced for it.
# If the bounds file is specified, then an additional estimate will be produced in a column named time_since_infection_bounded that is a variation of the time_since_infection column that is restricted so that it lies within the bounds.

cat('\n\n')
cat('Running estimateInfectionTime.R\n-------------------------------\n')

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option("--model_structure", 
            default = "full",
            help = "Which model structure do you want to use? (full or slope). full: The full model with viral loads, bounds and interactions. slope: The model that only modifies the estimate by multiplying it with a constant."),
make_option("--identify_founders.tab",
            help = "Path (Relative or Full) to and name of the identify_founders.tab file. Always required."),
make_option("--vl_file",
            help = "Path (Relative or Full) to and name of the file with the viral loads. Required for the full model only."),
make_option("--bounds_file",
            default = NULL,
            help = "Path (Relative or Full) to and name of the file with the bounds. Required for the full model and optional for the slope model."),
make_option("--estimator",
            help = "Which estimator from the identify_founders.tab file should be used? pfitter, synonymous_pfitter, multifounder_pfitter, or multifounder_synonymous_pfitter. The match is not case sensitive."),
make_option("--time_and_region",
            default = '1m_nflg',
            help = "Three different datasets were used to construct calibrated models. These varied based on the region of the genome and the time since infection. Choose which model you want to use by specifying the time and region: '1m_v3', '1m_nflg', '6m_pooled'. Since many of the assumptions break before 6 months, it is recommended to use either the 1m_v3 or 1m_nflg options."),
make_option("--path_to_model_coefficients_file",
            default = '.',
            help = "Path to the folder containing the model coefficients file. Because of how the rest of the pipeline works, you usually have to run it from within the hiv-founder-id folder, thus this defaults to '.'. If you are running this script from a different folder, specify the path to the model coefficients file with this variable. You can include the name of the file into this variable, if it is not included, then it is assumed that the file is named 'model_coefficients.csv'"),
make_option("--verbose",
            action = "store_true",
            default = FALSE,
            help = "Produces verbose output"),
make_option("--debug",
            action = "store_true",
            default = FALSE,
            help = "Prints out the structure detected for the vl_file and bounds_file as well as all the options that were selected.")
           )

opt <- parse_args(OptionParser(option_list = option_list,
                  description = "Predicts times of HIV infection from sequence data from samples very soon after infection. Comments at the start of this script provides some guidelines about the required file formats."))

# Processing settings relevant to ALL model structures:
identify_founders_tab <- read.delim(opt$identify_founders.tab, sep = '\t', stringsAsFactors = FALSE)

estimator_column_name <- switch(tolower(opt$estimator),
                                "pfitter" = "PFitter.time.est",
                                "synonymous_pfitter" = "Synonymous.PFitter.time.est",
                                "multifounder_pfitter" = "multifounder.PFitter.time.est",
                                "multifounder_synonymous_pfitter" = "multifounder.Synonymous.PFitter.time.est",
                                "synonymous_multifounder_pfitter" = "multifounder.Synonymous.PFitter.time.est",
                                "error")

if (estimator_column_name == 'error'){
  stop('invalid estimator specified: must be one of pfitter, synonymous_pfitter, multifounder_pfitter, or multifounder_synonymous_pfitter')
}

time_and_region <- opt$time_and_region
stopifnot(!is.null(time_and_region))
stopifnot(!is.na(time_and_region))
stopifnot(time_and_region %in% c('1m_v3', '1m_nflg', '6m_pooled'))

# remove trailing / from path:
mod_co_file <- gsub('\\/$', '', opt$path_to_model_coefficients_file)
if (file.exists(file.path(mod_co_file, 'model_coefficients.csv'))){
  all_coeffs <- read.csv(file.path(mod_co_file, 'model_coefficients.csv'), stringsAsFactors = FALSE)
} else if (file.exists(mod_co_file)) {
  all_coeffs <- read.csv(mod_co_file, stringsAsFactors = FALSE)
} else {
  stop(paste('Neither ', file.path(mod_co_file, 'model_coefficients.csv') , ' nor ', file.exists(mod_co_file) , ' exists.', sep = ''))
}

if (opt$debug){
  cat('\n\n=====================\nDEBUG\n=====================\n\n')
  cat("\nopt:\n----\n\n")
  print(opt)
  cat("\n\nBounds file:\n------------\n\n")
  if (is.null(opt$bounds_file)){
    cat('\nBounds file option not specified.\n')
  } else if (file.exists(opt$bounds_file)){
    print(str(read.csv(opt$bounds_file, stringsAsFactors = FALSE)))
  } else {
    cat("\nStrucutre of Bounds File:\n")
    cat("\nValue specified for bounds file, ", opt$bounds_file, " does not point to a file that exists.\n\n")
  }
  if (is.null(opt$vl_file)){
    cat('\nViral loads file option not specified.\n')
  } else if (file.exists(opt$vl_file)){
    cat("\nStrucutre of Viral Loads file:\n")
    print(str(read.csv(opt$vl_file, stringsAsFactors = FALSE)))
  } else {
    cat("\nValue specified for viral loads file, ", opt$vl_file, " does not point to a file that exists.\n\n")
  }
  cat('\n\nEnd of Debugging output\n\n')
}

# enforces bounds after estimation has already been performed
bounds_enforcer <- function(bounds_file, results){
  bounds_dat <- read.csv(bounds_file, stringsAsFactors = FALSE)
  results$time_since_infection_bounded <- NA_real_
  for (i in 1:nrow(bounds_dat)){
    
    bounds_indx <- i
    bounds_identifier <- bounds_dat[i, 1]
    lower_bound <- bounds_dat[i, 'lower']
    upper_bound <- bounds_dat[i, 'upper']
  
    results_indx <- which(grepl(as.character(bounds_identifier), as.character(results$identify_founders_tab_indx)))
    results_identifier <- as.character(results$identify_founders_tab_indx)[ results_indx ]

    time_since_infection <- results$time_since_infection[ results_indx ]
    time_since_infection_bounded <- time_since_infection
    if (time_since_infection < lower_bound) {
      time_since_infection_bounded <- lower_bound
    }
    if (time_since_infection > upper_bound) {
      time_since_infection_bounded <- upper_bound
    }

    results$time_since_infection_bounded[results_indx] <- time_since_infection_bounded
  }
  results
}

# predicts infection time using the full model
predict_with_full_model <- function(opt, identify_founders_tab, estimator_column_name, time_and_region, all_coeffs){

  coeffs <- subset(all_coeffs, y_type == 'time_since_infection' & model_structure == 'full')
  coeffs <- coeffs[,c(-1, -2)]

  which_row <- (coeffs$time_and_region == time_and_region) & (coeffs$estimator_column_name == estimator_column_name)
  rel_coeff <- coeffs[which_row, c('term', 'coeff')]

  coeff_est       <- rel_coeff[rel_coeff$term == 'est',       'coeff']
  coeff_Intercept <- rel_coeff[rel_coeff$term == 'Intercept', 'coeff']
  coeff_lPVL      <- rel_coeff[rel_coeff$term == 'lPVL',      'coeff']
  coeff_lPVL_est  <- rel_coeff[rel_coeff$term == 'lPVL:est',  'coeff']
  coeff_bound     <- rel_coeff[rel_coeff$term == 'bound',     'coeff']
  coeff_bound_est <- rel_coeff[rel_coeff$term == 'bound:est', 'coeff']

  vl_dat <- read.csv(opt$vl_file, stringsAsFactors = FALSE)
  bounds_dat <- read.csv(opt$bounds_file, stringsAsFactors = FALSE)

  results <- NULL
  for (i in 1:nrow(bounds_dat)){
    
    bounds_indx <- i
    bounds_identifier <- bounds_dat[i, 1]
    lower_bound <- bounds_dat[i, 'lower']
    upper_bound <- bounds_dat[i, 'upper']

    vl_indx <- which(grepl(as.character(bounds_identifier), as.character(vl_dat[,1])))
    vl_identifier <- as.character(vl_dat[,1])[ vl_indx ]
    
    ift_indx <- which(grepl(as.character(bounds_identifier), as.character(identify_founders_tab[,1])))
    ift_identifier <- as.character(identify_founders_tab[,1])[ ift_indx ]

    est <- identify_founders_tab[ift_indx, estimator_column_name]
    vl <- vl_dat[vl_indx, 'vl']
    bound <- bounds_dat[bounds_indx, 'upper']

    time_since_infection <- coeff_Intercept + 
                              coeff_est * est + 
                              vl * coeff_lPVL + 
                              vl * est * coeff_lPVL_est + 
                              bound * coeff_bound + 
                              bound * est * coeff_bound_est

    time_since_infection_bounded <- time_since_infection
    if (time_since_infection < lower_bound) {
      time_since_infection_bounded <- lower_bound
    }
    if (time_since_infection > upper_bound) {
      time_since_infection_bounded <- upper_bound
    }

    results <- rbind(results, 
      data.frame(identify_founders_tab_indx = ift_identifier,
                 bounds_indx = bounds_identifier,
                 vl_indx = vl_identifier,
                 time_since_infection = time_since_infection,
                 time_since_infection_bounded = time_since_infection_bounded,
                 stringsAsFactors = FALSE),
      stringsAsFactors = FALSE)
  }

  results <- bounds_enforcer(opt$bounds_file, results)
  results
}

# predicts infection time using the slope only model
predict_with_slope_model <- function(opt, identify_founders_tab, estimator_column_name, time_and_region, all_coeffs){
  
  coeffs <- subset(all_coeffs, y_type == 'time_since_infection' & model_structure == 'slope')
  coeffs <- coeffs[,c(-1, -2)]

  which_row <- (coeffs$time_and_region == time_and_region) & (coeffs$estimator_column_name == estimator_column_name)
  coeff_est <- coeffs[which_row, 'coeff']

  results <- NULL

  for (i in 1:nrow(identify_founders_tab)){
    est <- identify_founders_tab[i, estimator_column_name]
    indx <- identify_founders_tab[i, 1]
    results <- rbind(results, 
      data.frame(identify_founders_tab_indx = indx,
                 time_since_infection = est * coeff_est,
                 stringsAsFactors = FALSE),
      stringsAsFactors = FALSE)
  }
  if (!is.null(opt$bounds_file)){
    results <- bounds_enforcer(opt$bounds_file, results)
  }
  results
}

# main loop
if (opt$model_structure == 'full'){
  results <- predict_with_full_model(opt = opt, 
                          identify_founders_tab = identify_founders_tab, 
                          estimator_column_name = estimator_column_name,
                          time_and_region = time_and_region,
                          all_coeffs = all_coeffs)
} else {
  results <- predict_with_slope_model(opt = opt, 
                                      identify_founders_tab = identify_founders_tab, 
                                      estimator_column_name = estimator_column_name,
                                      time_and_region = time_and_region,
                                      all_coeffs = all_coeffs)
}

bound_label <- ifelse(is.null(opt$bounds_file), 'unbounded', 'bounded')
bound_status <- ifelse(is.null(opt$bounds_file), 'not bounded since no bounds were supplied', 'bounded based on the specified bounds')
out_name <- file.path(dirname(opt$identify_founders.tab), 
                             paste('infection_time_estimates_', 
                                   opt$model_structure,
                                   '_',
                                   opt$time_and_region,
                                   '_',
                                   opt$estimator,
                                   '_',
                                   bound_label,
                                   '.csv', sep = ''))

cat(paste('Infection time estimation completed for ', nrow(results), ' samples based on the metrics in ', opt$identify_founders.tab, '\n',
          'The ', opt$model_structure, ' model calibrated on the ', opt$time_and_region, ' samples with the "', opt$estimator, '" input variable was used. Estimates were ', bound_status, '.\n\n', sep = ''))

if (is.null(opt$bounds_file)){
  print(results)
} else {
  print(results[, c('identify_founders_tab_indx', 'time_since_infection', 'time_since_infection_bounded')])
}


cat(paste("\n\nResults written to ", out_name , '\n\n', sep = ''))
write.csv(results, out_name, 
          row.names = FALSE)
































































