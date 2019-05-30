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
# -- The other estimators are not supported by this script, but the user can modify this script to use them by just extending the approach.

# An important note on formats:
# Information from the three different files must be linked together.
# The identify_founders.tab is indexed by the name of the input fasta file
# The viral load and bounds files MUST be indexed by THE SAME unique patient identifier that occurs in the name of the input file to the pipeline that will end up as the index for the identify_founders.tab file.
# The identifiers will be read from the bounds file. If the identifier cannot be found in the index columns of both the viral load and identify_founders.tab files, no timing estimate will be produced for it.

suppressPackageStartupMessages(library("optparse"))

option_list <- list(
make_option("--model_structure", 
            default = "full",
            help = "Which model structure do you want to use? full: The full model with viral loads, bounds and interactions. slope: The model that only modifies the estimate by multiplying it with a constant."),
make_option("--identify_founders.tab",
            help = "Path (Relative or Full) to and name of the identify_founders.tab file."),
make_option("--vl_file",
            help = "Path (Relative or Full) to and name of the file with the viral loads."),
make_option("--bounds_file",
            help = "Path (Relative or Full) to and name of the file with the bounds."),
make_option("--estimator",
            help = "Which estimator from the identify_founders.tab file should be used? PFitter, ..."),
make_option("--time_and_region",
            default = '1m_nflg',
            help = "Choose a time and region from which to choose parameters: '1m_v3', '1m_nflg', '6m_pooled'"),
make_option("--verbose",
            action = "store_true",
            default = FALSE,
            help = "Produces verbose output")
           )

opt <- parse_args(OptionParser(option_list = option_list,
                  description = "Predicts times of HIV infection from sequence data during acute infection",
                  epilogue = "I want donuts"))

# Processing settings relevant to ALL model structures:
identify_founders_tab <- read.delim(opt$identify_founders.tab, sep = '\t', stringsAsFactors = FALSE)

estimator_column_name <- switch(tolower(opt$estimator),
                                pfitter = "PFitter.time.est",
                                "synonymous pfitter" = "Synonymous.PFitter.time.est",
                                "error")
if (estimator_column_name == 'error'){
  stop('invalid estimator specified')
}

time_and_region <- opt$time_and_region
stopifnot(!is.null(time_and_region))
stopifnot(!is.na(time_and_region))
stopifnot(time_and_region %in% c('1m_v3', '1m_nflg', '6m_pooled'))



predict_with_full_model <- function(opt, identify_founders_tab, estimator_column_name, time_and_region){

  coeffs <- data.frame(time_and_region = rep("", 72),
                       estimator_column_name = rep("", 72),
                       term = rep("", 72),
                       coeff = rep(NA_real_, 72),
                       stringsAsFactors = FALSE)

  coeffs[1 , ] <- list('1m_v3', "PFitter.time.est",                             'est',       0.02)
  coeffs[2 , ] <- list('1m_v3', "PFitter.time.est",                             'Intercept', 38.60)
  coeffs[3 , ] <- list('1m_v3', "PFitter.time.est",                             'lPVL',      2.31)
  coeffs[4 , ] <- list('1m_v3', "PFitter.time.est",                             'lPVL:est',  -0.00)
  coeffs[5 , ] <- list('1m_v3', "PFitter.time.est",                             'bound',     -0.02)
  coeffs[6 , ] <- list('1m_v3', "PFitter.time.est",                             'bound:est', 0.00)
  coeffs[7 , ] <- list('1m_v3', "Synonymous.PFitter.time.est",                  'est',       0.14)
  coeffs[8 , ] <- list('1m_v3', "Synonymous.PFitter.time.est",                  'Intercept', 36.11)
  coeffs[9 , ] <- list('1m_v3', "Synonymous.PFitter.time.est",                  'lPVL',      2.97)
  coeffs[10, ] <- list('1m_v3', "Synonymous.PFitter.time.est",                  'lPVL:est',  -0.02)
  coeffs[11, ] <- list('1m_v3', "Synonymous.PFitter.time.est",                  'bound',     -0.07)
  coeffs[12, ] <- list('1m_v3', "Synonymous.PFitter.time.est",                  'bound:est', 0.00)
  coeffs[13, ] <- list('1m_v3', "multifounder.PFitter.time.est",                'est',       0.03)
  coeffs[14, ] <- list('1m_v3', "multifounder.PFitter.time.est",                'Intercept', 46.87)
  coeffs[15, ] <- list('1m_v3', "multifounder.PFitter.time.est",                'lPVL',      1.07)
  coeffs[16, ] <- list('1m_v3', "multifounder.PFitter.time.est",                'lPVL:est',  0.00)
  coeffs[17, ] <- list('1m_v3', "multifounder.PFitter.time.est",                'bound',     0.03)
  coeffs[18, ] <- list('1m_v3', "multifounder.PFitter.time.est",                'bound:est', 0.00)
  coeffs[19, ] <- list('1m_v3', "multifounder.Synonymous.PFitter.time.est",     'est',       0.95)
  coeffs[20, ] <- list('1m_v3', "multifounder.Synonymous.PFitter.time.est",     'Intercept', 31.20)
  coeffs[21, ] <- list('1m_v3', "multifounder.Synonymous.PFitter.time.est",     'lPVL',      2.07)
  coeffs[22, ] <- list('1m_v3', "multifounder.Synonymous.PFitter.time.est",     'lPVL:est',  -0.06)
  coeffs[23, ] <- list('1m_v3', "multifounder.Synonymous.PFitter.time.est",     'bound',     0.08)
  coeffs[24, ] <- list('1m_v3', "multifounder.Synonymous.PFitter.time.est",     'bound:est', 0.00)
  coeffs[25, ] <- list('1m_nflg', "PFitter.time.est",                           'est',       -0.00)
  coeffs[26, ] <- list('1m_nflg', "PFitter.time.est",                           'Intercept', 41.60)
  coeffs[27, ] <- list('1m_nflg', "PFitter.time.est",                           'lPVL',      0.14)
  coeffs[28, ] <- list('1m_nflg', "PFitter.time.est",                           'lPVL:est',  0.00)
  coeffs[29, ] <- list('1m_nflg', "PFitter.time.est",                           'bound',     0.04)
  coeffs[30, ] <- list('1m_nflg', "PFitter.time.est",                           'bound:est', -0.00)
  coeffs[31, ] <- list('1m_nflg', "Synonymous.PFitter.time.est",                'est',       -0.14)
  coeffs[32, ] <- list('1m_nflg', "Synonymous.PFitter.time.est",                'Intercept', 44.61)
  coeffs[33, ] <- list('1m_nflg', "Synonymous.PFitter.time.est",                'lPVL',      0.01)
  coeffs[34, ] <- list('1m_nflg', "Synonymous.PFitter.time.est",                'lPVL:est',  0.00)
  coeffs[35, ] <- list('1m_nflg', "Synonymous.PFitter.time.est",                'bound',     0.02)
  coeffs[36, ] <- list('1m_nflg', "Synonymous.PFitter.time.est",                'bound:est', 0.00)
  coeffs[37, ] <- list('1m_nflg', "multifounder.PFitter.time.est",              'est',       -0.03)
  coeffs[38, ] <- list('1m_nflg', "multifounder.PFitter.time.est",              'Intercept', 42.51)
  coeffs[39, ] <- list('1m_nflg', "multifounder.PFitter.time.est",              'lPVL',      0.03)
  coeffs[40, ] <- list('1m_nflg', "multifounder.PFitter.time.est",              'lPVL:est',  0.00)
  coeffs[41, ] <- list('1m_nflg', "multifounder.PFitter.time.est",              'bound',     0.04)
  coeffs[42, ] <- list('1m_nflg', "multifounder.PFitter.time.est",              'bound:est', -0.00)
  coeffs[43, ] <- list('1m_nflg', "multifounder.Synonymous.PFitter.time.est",   'est',       -0.46)
  coeffs[44, ] <- list('1m_nflg', "multifounder.Synonymous.PFitter.time.est",   'Intercept', 48.86)
  coeffs[45, ] <- list('1m_nflg', "multifounder.Synonymous.PFitter.time.est",   'lPVL',      -0.27)
  coeffs[46, ] <- list('1m_nflg', "multifounder.Synonymous.PFitter.time.est",   'lPVL:est',  0.02)
  coeffs[47, ] <- list('1m_nflg', "multifounder.Synonymous.PFitter.time.est",   'bound',     0.00)
  coeffs[48, ] <- list('1m_nflg', "multifounder.Synonymous.PFitter.time.est",   'bound:est', 0.00)
  coeffs[49, ] <- list('6m_pooled', "PFitter.time.est",                         'est',       0.01)
  coeffs[50, ] <- list('6m_pooled', "PFitter.time.est",                         'Intercept', 188.42)
  coeffs[51, ] <- list('6m_pooled', "PFitter.time.est",                         'lPVL',      -1.16)
  coeffs[52, ] <- list('6m_pooled', "PFitter.time.est",                         'lPVL:est',  -0.00)
  coeffs[53, ] <- list('6m_pooled', "PFitter.time.est",                         'bound',     0.02)
  coeffs[54, ] <- list('6m_pooled', "PFitter.time.est",                         'bound:est', 0.00)
  coeffs[55, ] <- list('6m_pooled', "Synonymous.PFitter.time.est",              'est',       0.24)
  coeffs[56, ] <- list('6m_pooled', "Synonymous.PFitter.time.est",              'Intercept', 178.07)
  coeffs[57, ] <- list('6m_pooled', "Synonymous.PFitter.time.est",              'lPVL',      -0.80)
  coeffs[58, ] <- list('6m_pooled', "Synonymous.PFitter.time.est",              'lPVL:est',  -0.01)
  coeffs[59, ] <- list('6m_pooled', "Synonymous.PFitter.time.est",              'bound',     0.06)
  coeffs[60, ] <- list('6m_pooled', "Synonymous.PFitter.time.est",              'bound:est', -0.00)
  coeffs[61, ] <- list('6m_pooled', "multifounder.PFitter.time.est",            'est',       -0.03)
  coeffs[62, ] <- list('6m_pooled', "multifounder.PFitter.time.est",            'Intercept', 192.38)
  coeffs[63, ] <- list('6m_pooled', "multifounder.PFitter.time.est",            'lPVL',      -1.33)
  coeffs[64, ] <- list('6m_pooled', "multifounder.PFitter.time.est",            'lPVL:est',  0.00)
  coeffs[65, ] <- list('6m_pooled', "multifounder.PFitter.time.est",            'bound',     0.02)
  coeffs[66, ] <- list('6m_pooled', "multifounder.PFitter.time.est",            'bound:est', 0.00)
  coeffs[67, ] <- list('6m_pooled', "multifounder.Synonymous.PFitter.time.est", 'est',       0.30)
  coeffs[68, ] <- list('6m_pooled', "multifounder.Synonymous.PFitter.time.est", 'Intercept', 178.16)
  coeffs[69, ] <- list('6m_pooled', "multifounder.Synonymous.PFitter.time.est", 'lPVL',      -0.84)
  coeffs[70, ] <- list('6m_pooled', "multifounder.Synonymous.PFitter.time.est", 'lPVL:est',  -0.01)
  coeffs[71, ] <- list('6m_pooled', "multifounder.Synonymous.PFitter.time.est", 'bound',     0.06)
  coeffs[72, ] <- list('6m_pooled', "multifounder.Synonymous.PFitter.time.est", 'bound:est', -0.00)

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

    vl_indx <- which(grepl(as.character(bounds_identifier), as.character(vl_dat[,1])))
    vl_identifier <- as.character(vl_dat[,1])[ vl_indx ]
    
    ift_indx <- which(grepl(as.character(bounds_identifier), as.character(identify_founders_tab[,1])))
    ift_identifier <- as.character(identify_founders_tab[,1])[ ift_indx ]

    est <- identify_founders_tab[ift_indx, estimator_column_name]
    vl <- vl_dat[vl_indx, 'vl']
    bound <- bounds_dat[bounds_indx, 'upper']

    results <- rbind(results, 
      data.frame(bounds_indx = bounds_identifier,
                 vl_indx = vl_identifier,
                 identify_founders_tab_indx = ift_identifier,
                 time_since_infection = coeff_Intercept + 
                   coeff_est * est + 
                   vl * coeff_lPVL + 
                   vl * est * coeff_lPVL_est + 
                   bound * coeff_bound + 
                   bound * est * coeff_bound_est,
                 stringsAsFactors = FALSE),
      stringsAsFactors = FALSE)
  }

  results
}

predict_with_slope_model <- function(opt, identify_founders_tab, estimator_column_name, time_and_region){
  coeffs <- data.frame(time_and_region = rep("", 12),
                       estimator_column_name = rep("", 12),
                       coeff = rep(NA_real_, 12),
                       stringsAsFactors = FALSE)

  coeffs[ 1, ] <- list("1m_v3",     "PFitter.time.est",                         0.08)
  coeffs[ 2, ] <- list("1m_v3",     "Synonymous.PFitter.time.est",              0.48)
  coeffs[ 3, ] <- list("1m_v3",     "multifounder.PFitter.time.est",            0.45)
  coeffs[ 4, ] <- list("1m_v3",     "multifounder.Synonymous.PFitter.time.est", 0.94)
  coeffs[ 5, ] <- list("1m_nflg",   "PFitter.time.est",                         0.23)
  coeffs[ 6, ] <- list("1m_nflg",   "Synonymous.PFitter.time.est",              1.04)
  coeffs[ 7, ] <- list("1m_nflg",   "multifounder.PFitter.time.est",            0.40)
  coeffs[ 8, ] <- list("1m_nflg",   "multifounder.Synonymous.PFitter.time.est", 1.70)
  coeffs[ 9, ] <- list("6m_pooled", "PFitter.time.est",                         0.54)
  coeffs[10, ] <- list("6m_pooled", "Synonymous.PFitter.time.est",              2.03)
  coeffs[11, ] <- list("6m_pooled", "multifounder.PFitter.time.est",            0.82)
  coeffs[12, ] <- list("6m_pooled", "multifounder.Synonymous.PFitter.time.est", 2.65)

  which_row <- (coeffs$time_and_region == time_and_region) & (coeffs$estimator_column_name == estimator_column_name)
  coeff_est <- coeffs[which_row, 'coeff']

  results <- NULL

  for (i in 1:nrow(identify_founders_tab)){
    est <- identify_founders_tab[i, estimator_column_name]
    indx <- identify_founders_tab[i, 1]
    results <- rbind(results, 
      data.frame(indx = indx,
                 time_since_infection = est * coeff_est,
                 stringsAsFactors = FALSE),
      stringsAsFactors = FALSE)
  }
  results
}

if (opt$model_structure == 'full'){
  results <- predict_with_full_model(opt = opt, 
                          identify_founders_tab = identify_founders_tab, 
                          estimator_column_name = estimator_column_name,
                          time_and_region = time_and_region)
} else {
  results <- predict_with_slope_model(opt = opt, 
                                      identify_founders_tab = identify_founders_tab, 
                                      estimator_column_name = estimator_column_name,
                                      time_and_region = time_and_region)
}

out_name <- file.path(dirname(opt$identify_founders.tab), 
                             paste('infection_time_estimates_', 
                                   opt$model_structure,
                                   '_',
                                   opt$time_and_region,
                                   '.csv', sep = ''))
print(results)
print(paste("Results written to ", out_name, sep = ''))
write.csv(results, out_name, 
          row.names = FALSE)


