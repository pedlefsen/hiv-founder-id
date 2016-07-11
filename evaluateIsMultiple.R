## First do all the stuff in README.postprocessing.txt.
library( "ROCR" ) # for "prediction" and "performance"
library( "glmnet" ) # for cv.glmnet
library( "parallel" ) # for mcapply

source( "readIdentifyFounders_safetosource.R" );
source( "getArtificialBounds_safetosource.R" );
source( "getResultsByRegionAndTime_safetosource.R" );
source( "writeResultsTables_safetosource.R" );
source( "summarizeCovariatesOnePerParticipant_safetosource.R" );

GOLD.STANDARD.DIR <- "/fh/fast/edlefsen_p/bakeoff/gold_standard";
RESULTS.DIR <- "/fh/fast/edlefsen_p/bakeoff_analysis_results/";

#' Evaluate isMultiple estimates and produce results tables.
#'
#' This function runs the BakeOff results analysis for the isMultiple results.
#' This will also look for plasma viral load measurements in files
#' called
#' /fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/caprisa_002_viralloads.csv
#' and
#' /fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/rv217_viralloads.csv. These
#' files have three important columns: ptid,viralload,timepoint.  For
#' each ptid the viral loads at timepoints 1,2,3 correspond to the
#' gold-standard, 1m, and 6m time points.  Viral loads are not logged
#' in the input file.
#'
#' @param use.glm.validate evaluate predicted values from leave-one-out cross-validation.
#' @param include.bounds.in.glm include the corresponding prior bounds in the regression equations used for cross-validation.
#' @param include.bounds.in.lasso include the corresponding prior bounds in the regression equations used for cross-validation.
#' @param include.helpful.additional.cols.in.glm include helpful.additional.cols in the glm (but don't force them to stay in the lasso)? By default this is the opposite of include.bounds.in.glm).
#' @param helpful.additional.cols [TODO: extra cols to be included in the glm. By default, "pred.days", which is a special code meaning that the results should be computed for _each_ of the sets of preditions in the set: { same-time-and-region, same-time, same-region, all-times-and-regions }.]
#' @param results.dirname the subdirectory of "/fh/fast/edlefsen_p/bakeoff/analysis_sequences" and also of "/fh/fast/edlefsen_p/bakeoff_analysis_results"
#' @param force.recomputation if FALSE (default) and if there is a saved version called isMultiple.results.by.region.and.time.Rda (under bakeoff_analysis_results/results.dirname), then that file will be loaded; otherwise the results will be recomputed and saved in that location.
#' @param partition.bootstrap.seed the random seed to use when bootstrapping samples by selecting one partition number per ptid, repeatedly; we do it this way because there are an unequal number of partitions, depending on sampling depth.
#' @param partition.bootstrap.samples the number of bootstrap replicates to conduct; the idea is to get an estimate of the variation in estimates and results (errors) across these samples.
#' @param partition.bootstrap.num.cores the number of cores to run the boostrap replicates on (defaults to all of the cores returned by parallel::detectCores()).
#' @return NULL
#' @export

evaluateIsMultiple <- function (
                             use.bounds = TRUE,
                             use.glm.validate = TRUE,
                             use.lasso.validate = TRUE,
                             include.bounds.in.glm = TRUE,
                             include.bounds.in.lasso = TRUE,
                             include.helpful.additional.cols.in.glm = !include.bounds.in.glm,
                             helpful.additional.cols = c(),
                             results.dirname = "raw_edited_20160216",
                             force.recomputation = FALSE,
                             partition.bootstrap.seed = 98103,
                             partition.bootstrap.samples = 100,
                             partition.bootstrap.num.cores = detectCores(),
                             regions = c( "nflg", "v3" ),
                             times = c( "1m", "6m" )
                             # regions = c( "nflg", "v3", "rv217_v3" ),
                             # times = c( "1m", "6m", "1m6m" )
                            )
{
    is.multiple.results.by.region.and.time.Rda.filename <-
        paste( RESULTS.DIR, results.dirname, "/isMultiple.results.by.region.and.time.Rda", sep = "" );

    ## Read in the gold standards.
    # From this file we read in indicators of whether to use the multiple- or single-founder true profile.
    rv217.gold.standards.in <- read.csv( paste( GOLD.STANDARD.DIR, "/rv217/RV217_gold_standards.csv", sep = "" ) );
    rv217.gold.is.multiple <- rv217.gold.standards.in[ , "gold.is.multiple" ];
    names( rv217.gold.is.multiple ) <- rv217.gold.standards.in[ , "ptid" ];
    
    caprisa002.gold.standards.in <- read.csv( paste( GOLD.STANDARD.DIR, "/caprisa_002/caprisa_002_gold_standards.csv", sep = "" ) );
    caprisa002.gold.is.multiple <- caprisa002.gold.standards.in[ , "gold.is.multiple" ];
    names( caprisa002.gold.is.multiple ) <- caprisa002.gold.standards.in[ , "ptid" ];

    rv217.pvl.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/rv217_viralloads.csv" );
    caprisa002.pvl.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/caprisa_002_viralloads.csv" );
    
    ## Sometimes there are multiple entries for one ptid/sample, eg
    ## for the NFLGs there are often right half and left half (RH and
    ## LH) and sometimes additionally NFLG results.  If so, for
    ## something to be called single founder, all of the estimates
    ## (across regions, for a given test) must agree that it is one
    ## founder.
    compute.estimates.is.one.founder.per.person <- function ( estimates.is.one.founder ) {
        identify.founders.ptids <- rownames( estimates.is.one.founder );
            
        estimates.is.one.founder.per.person <- 
            t( sapply( unique( identify.founders.ptids ), function ( .ptid ) {
            .ptid.subtable <- estimates.is.one.founder[ identify.founders.ptids == .ptid, , drop = FALSE ];
            # Use the AND condition, meaning it's only called single-founder if all rows call it single-founder.
            if( nrow( .ptid.subtable ) == 1 ) {
              return( .ptid.subtable );
            }
            apply( .ptid.subtable, 2, function( .column ) { as.numeric( all( as.logical( .column ) ) ) } );
        } ) );
            
        colnames( estimates.is.one.founder.per.person ) <- colnames( estimates.is.one.founder );
        return( estimates.is.one.founder.per.person );
    } # compute.estimates.is.one.founder.per.person (..)

    bound.and.evaluate.is.multiple.results.per.ppt <-
        function ( estimates.is.one.founder.per.person, gold.is.one.founder.per.person, results.covars.per.person.with.extra.cols, the.time, the.artificial.bounds = NA, ppt.suffix.pattern = "\\..+", return.lasso.coefs = TRUE ) {
          ## TODO: REMOVE.  DEBUGGING.
          #estimates.is.one.founder.per.person.in <- estimates.is.one.founder.per.person;
          
       ## Special: the ppt names might have suffices in results.per.person; if so, strip off the suffix for purposes of matching ppts to the covars, etc.
       ppt.names <- rownames( estimates.is.one.founder.per.person );
       ppt.suffices <- NULL;
       if( !is.na( ppt.suffix.pattern ) && ( length( grep( ppt.suffix.pattern, ppt.names ) ) > 0 ) ) {
         # if one has it, they should all have it.
         stopifnot( length( grep( ppt.suffix.pattern, ppt.names ) ) == length( ppt.names ) );
           .ppt.names <- ppt.names;
           ppt.names <- gsub( ppt.suffix.pattern, "", .ppt.names );
           names( ppt.names ) <- .ppt.names;
           ppt.suffices <- gsub( paste( "^.+?(", ppt.suffix.pattern, ")", sep = "" ), "\\1", .ppt.names, perl = TRUE );
           names( ppt.suffices ) <- .ppt.names;
       }
       all.ptids <- unique( ppt.names );

       mode( estimates.is.one.founder.per.person ) <- "numeric";
       
       if( use.glm.validate || use.lasso.validate ) {
           results.covars.per.person.with.extra.cols <-
               cbind( estimates.is.one.founder.per.person, results.covars.per.person.with.extra.cols[ , setdiff( colnames( results.covars.per.person.with.extra.cols ), colnames( estimates.is.one.founder.per.person ) ) , drop = FALSE ] );
        
        .keep.cols <-
           grep( "num.*\\.seqs|totalbases|upper|lower", colnames( results.covars.per.person.with.extra.cols ), value = TRUE, perl = TRUE, invert = TRUE );
           # There are redundancies because the mut.rate.coef for DS is identical to PFitter's and for Bayesian it is very similar.
        .keep.cols <-
            grep( "Star[pP]hy\\.mut\\.rate\\.coef", .keep.cols, value = TRUE, invert = TRUE );
        # Also exclude this strange test.
        .keep.cols <-
            grep( "DS\\.Star[pP]hy\\.is\\.starlike", .keep.cols, value = TRUE, invert = TRUE );
        # Also exclude this, which is based on the strange test.
        .keep.cols <-
            grep( "StarPhy\\.is\\.one\\.founder", .keep.cols, value = TRUE, invert = TRUE );

        # For COB and infer, use only the real-data sources (mtn003 or hvtn502). So exclude the "mtn003" and "hvtn502" ones.
        .keep.cols <-
            grep( "\\.(one|six)months?\\.", .keep.cols, value = TRUE, invert = TRUE );
           
        ## Try removing some variables that are rarely selected or are too correlated (eg diversity is highly correlated with sd.entropy, max.hd, insites.is.one.founder, insites.founders) [SEE BELOW WHERE WE ADD TO THIS PROGRAMMATICALLY]
           .donotkeep.cols <- c( "multifounder.DS.Starphy.R", "PFitter.chi.sq.stat", "Synonymous.DS.StarPhy.R", "StarPhy.is.one.founder", "DS.Starphy.fits", "DS.Starphy.is.starlike", "insites.is.one.founder", "inf.to.priv.ratio" );
           # NOTE ALSO MORE (formerly excluded below because cor > 0.9):
           .donotkeep.cols <- c( .donotkeep.cols, "StarPhy.founders", "InSites.founders", "PFitter.max.hd",   "PFitter.mean.hd", "sd.entropy",       "mean.entropy",     "inf.sites" );

           .keep.cols <- setdiff( .keep.cols, .donotkeep.cols );
           
           single.cols <- grep( "\\.is\\.|fits", .keep.cols, perl = TRUE, value = TRUE );
           mut.rate.coef.cols <- grep( "mut\\.rate\\.coef", .keep.cols, value = TRUE );
           
           all.additional.cols <- setdiff( .keep.cols, c( single.cols, mut.rate.coef.cols ) );
          
        if( use.lasso.validate ) {
            keep.cols <- c( all.additional.cols, single.cols, mut.rate.coef.cols );
            estimate.cols <- setdiff( keep.cols, all.additional.cols );
        } else {
            keep.cols <- c( helpful.additional.cols, single.cols, mut.rate.coef.cols );
            estimate.cols <- setdiff( keep.cols, helpful.additional.cols );
        }

        # Don't evaluate estimators that have no variation at
        # all. Note that we also redo this after holding
        # out each person, in case a value has no variation among the
        # subset excluding that person.  But we do it here first.
        estimators.to.exclude <-
          apply( results.covars.per.person.with.extra.cols[ , estimate.cols, drop = FALSE ], 2, function ( .col ) {
            return( ( var( .col ) == 0 ) || ( sum( !is.na( .col ) ) <= 1 ) );
          } );
        estimate.cols <- setdiff( estimate.cols, names( which( estimators.to.exclude ) ) );
           
        # Also always evaluate no-estimate: "none".
        estimate.cols <- c( "none", estimate.cols );
          
        results.covars.per.person <-
            results.covars.per.person.with.extra.cols[ , keep.cols, drop = FALSE ];
        results.covars.per.person.df <-
            data.frame( results.covars.per.person );

        regression.df <- cbind( data.frame( is.one.founder = gold.is.one.founder.per.person[ rownames( results.covars.per.person.df ) ] ), results.covars.per.person.df, lapply( the.artificial.bounds[ grep( "(one|six)month", names( the.artificial.bounds ), invert = TRUE, value = TRUE ) ], function( .mat ) { .mat[ rownames( results.covars.per.person.df ), , drop = FALSE ] } ) );
        
        ## Ok build a regression model with no intercept, including only the helpful.additional.cols, and also the lower and upper bounds associated with either 5 weeks or 30 weeks, depending on the.time (if there's a 1m sample, uses "mtn003").
        if( the.time == "6m" ) {
            .lower.bound.colname <- "sampledwidth_uniform_hvtn502.lower";
            .upper.bound.colname <- "sampledwidth_uniform_hvtn502.upper";
        } else if( the.time == "1m.6m" ) {
            .lower.bound.colname <- "sampledwidth_uniform_1mmtn003_6mhvtn502.lower";
            .upper.bound.colname <- "sampledwidth_uniform_1mmtn003_6mhvtn502.upper";
        } else {
            .lower.bound.colname <- "sampledwidth_uniform_mtn003.lower";
            .upper.bound.colname <- "sampledwidth_uniform_mtn003.upper";
        }
        
        if( use.glm.validate ) {
            ## Note that there might be multiple rows per ppt in the regression.df and in this prediction output matrix; the values will be filled in using leave-one-ptid-out xv, ie in each iteration there might be multiple rows filled in, since multiple rows correspond to one held-out ptid.
            glm.validation.estimates.is.one.founder.per.person <- matrix( NA, nrow = nrow( regression.df ), ncol = length( estimate.cols ) );
        }
        if( use.lasso.validate ) {
            ## See note above (for use.glm.validate)
            lasso.validation.estimates.is.one.founder.per.person <-
                matrix( NA, nrow = nrow( regression.df ), ncol = length( estimate.cols ) );
            if( return.lasso.coefs ) {
                ## This is really a 3D array, but I'm just lazily representing it directly this way.  Note this is by removed ppt, not by regression row (there might be multiple rows per ppt in the regression.df and the predictio noutput matrices).
                lasso.validation.estimates.is.one.founder.per.person.coefs <-
                    as.list( rep( NA, length( all.ptids ) ) );
                names( lasso.validation.estimates.is.one.founder.per.person.coefs ) <-
                    all.ptids;
            }
        }
        for( .ptid.i in 1:length( all.ptids ) ) {
            the.ptid <- all.ptids[ .ptid.i ];
            the.rows.for.ptid <- which( ppt.names == the.ptid );
            the.rows.excluding.ptid <- which( ppt.names != the.ptid );
            ## TODO: REMOVE
            print( paste( "PTID", .ptid.i, "removed:", the.ptid, "rows:(", paste( the.rows.for.ptid, collapse = ", " ), ")" ) );
            if( use.lasso.validate && return.lasso.coefs ) {
                .lasso.validation.estimates.is.one.founder.per.person.coefs.row <-
                    as.list( rep( NA, length( estimate.cols ) ) );
                names( .lasso.validation.estimates.is.one.founder.per.person.coefs.row ) <-
                    estimate.cols;
            }
            regression.df.without.ptid.i <-
                regression.df[ the.rows.excluding.ptid, , drop = FALSE ];

            # Don't evaluate estimators that have no variation at
            # all. Note that we also do this above, before holding
            # out each person.  But we do it here too.
            .estimators.to.exclude <- # exclude the first ("none") when doing this:
              apply( regression.df.without.ptid.i[ , estimate.cols[ -1 ], drop = FALSE ], 2, function ( .col ) {
                return( ( var( .col ) == 0 ) || ( sum( !is.na( .col ) ) <= 1 ) );
              } );
            .estimate.cols <- setdiff( estimate.cols, names( which( .estimators.to.exclude ) ) );
           
            for( .col.i in 1:length( .estimate.cols ) ) {
                .estimate.colname <- .estimate.cols[ .col.i ];
                ##print( .estimate.colname );
                if( use.glm.validate ) {
                  # covariates for glm
                  .covariates.glm <- c();
                  if( include.helpful.additional.cols.in.glm ) {
                      .covariates.glm <-
                          c( .covariates.glm, helpful.additional.cols );
                  }
                  if( include.bounds.in.glm ) {
                      .covariates.glm <-
                        #c( .covariates.glm, .lower.bound.colname, .upper.bound.colname );
                        c( .covariates.glm, .upper.bound.colname );
                  }
                  # glm:v
                  .covars.to.exclude <- apply( regression.df.without.ptid.i, 2, function ( .col ) {
                      return( ( var( .col ) == 0 ) || ( sum( !is.na( .col ) ) <= 1 ) );
                  } );
                  .retained.covars <-
                    setdiff( colnames( regression.df.without.ptid.i ), names( which( .covars.to.exclude ) ) );
                  if( ( .estimate.colname == "none" ) || ( .estimate.colname %in% .retained.covars ) ) {
                    if( .estimate.colname == "none" ) {
                          .cv.glm <- intersect( .retained.covars, .covariates.glm );
                          if( length( .cv.glm ) == 0 ) {
                              .cv.glm <- "1";
                          }
                          .formula <- as.formula( paste( "is.one.founder ~ ", paste( .cv.glm, collapse = "+" ) ) );
                          .ptids.to.exclude <- c();
                    } else {
                        .formula <- as.formula( paste( "is.one.founder ~ ", paste( intersect( .retained.covars, c( .covariates.glm, .estimate.colname ) ), collapse = "+" ) ) );

                      .rows.to.exclude.helper <-
                        is.na( regression.df.without.ptid.i[ , .estimate.colname ] );
                      .rows.to.exclude <- which( .rows.to.exclude.helper );
                      if( length( .rows.to.exclude ) == 0 ) {
                        .ptids.to.exclude <- c();
                      } else {
                        .ptids.to.exclude <-
                          rownames( regression.df.without.ptid.i )[ .rows.to.exclude ];
                      }                
                        
                    }
                    
                    .df <- regression.df.without.ptid.i[ , .retained.covars, drop = FALSE ];

                    if( length( .ptids.to.exclude ) > 0 ) {
                      .df <- .df[ -.rows.to.exclude, , drop = FALSE ];
                    }
                    .lm <-
                        glm( .formula, family = "binomial", data = .df );
                    ## TODO: REMOVE
                    ##print( summary( .lm ) );
                    # Ok so now put the predicted values at their appropriate places in the matrix.
                    for( .row.i in the.rows.for.ptid ) {
                        .pred.value.glm <- predict( .lm, regression.df[ .row.i, , drop = FALSE ], type = "response" );
                        glm.validation.estimates.is.one.founder.per.person[ .row.i, .col.i ] <- 
                            .pred.value.glm;
                    }
                  } # If the estimate can be used
    
                } # End if use.glm.validate
    
                if( use.lasso.validate ) {
                  # covariates for lasso
                  .covariates.lasso <- c( all.additional.cols );
                  if( include.bounds.in.lasso ) {
                      .covariates.lasso <-
                        c( .covariates.lasso, .lower.bound.colname, .upper.bound.colname );
                        #c( .covariates.lasso, .upper.bound.colname );
                  }
                  # lasso:
                  if( .estimate.colname == "none" ) {
                      .mat1 <- as.matrix( regression.df.without.ptid.i[ , .covariates.lasso, drop = FALSE ] );
                        .ptids.to.exclude <- c();
                  } else {
                      .mat1 <- as.matrix( regression.df.without.ptid.i[ , c( .covariates.lasso, .estimate.colname ) ] );
                      .rows.to.exclude.helper <-
                        is.na( regression.df.without.ptid.i[ , .estimate.colname ] );
                      .rows.to.exclude <- which( .rows.to.exclude.helper );
                      if( length( .rows.to.exclude ) == 0 ) {
                        .ptids.to.exclude <- c();
                      } else {
                        .ptids.to.exclude <-
                          rownames( regression.df.without.ptid.i )[ .rows.to.exclude ];
                      }                
                  }
                  .out <- regression.df.without.ptid.i[[ "is.one.founder" ]];
                  
                  .covars.to.exclude <- apply( .mat1[ , setdiff( colnames( .mat1 ), .estimate.colname ) ], 2, function ( .col ) {
                      return( ( var( .col ) == 0 ) || ( sum( !is.na( .col ) ) <= 1 ) );
                  } );
                  .mat1 <- .mat1[ , setdiff( colnames( .mat1 ), names( which( .covars.to.exclude ) ) ), drop = FALSE ];
                  .covars.to.exclude <- names( which( .covars.to.exclude ) );
                  if( length( .ptids.to.exclude ) > 0 ) {
                      ## TODO: REMOVE
                      #print( paste( "excluding samples (", paste( .ptids.to.exclude, collapse = ", " ), ") due to NAs in", .estimate.colname ) );
                      .mat1 <- .mat1[ -.rows.to.exclude, , drop = FALSE ];
                      .out <- .out[ -.rows.to.exclude ];
                   }
                   # Exclude covars that are too highly correlated with the estimate.
                    COR.THRESHOLD <- 0.9;
                    if( .estimate.colname != "none" ) {
                      .covars.to.consider <-
                        setdiff( colnames( .mat1 ), c( .covars.to.exclude, .estimate.colname ) );
                      .covar.is.too.highly.correlated.with.upper.bound <-
                        sapply( .covars.to.consider, function( .covar.colname ) {
                            #print( .covar.colname );
                          .cor <- 
                              cor( .mat1[ , .estimate.colname ], .mat1[ , .covar.colname ], use = "pairwise" );
                            #print( .cor );
                          if( is.na( .cor ) || ( abs( .cor ) >= COR.THRESHOLD ) ) {
                            TRUE
                          } else {
                            FALSE
                          }
                        } );
                      .covars.to.exclude <- c( .covars.to.exclude,
                          names( which( .covar.is.too.highly.correlated.with.upper.bound ) ) );
                    }
                   # Exclude covars that are too highly correlated with the upper bound.
                   if( include.bounds.in.lasso ) {
                     .covars.to.consider <-
                       setdiff( colnames( .mat1 ), c( .covars.to.exclude, .estimate.colname, .upper.bound.colname ) );
                     if( length( .covars.to.consider ) > 0 ) {
                       .covar.is.too.highly.correlated.with.upper.bound <-
                         sapply( .covars.to.consider, function( .covar.colname ) {
                              #print( .covar.colname );
                            .cor <- 
                                cor( .mat1[ , .upper.bound.colname ], .mat1[ , .covar.colname ], use = "pairwise" );
                              #print( .cor );
                            if( is.na( .cor ) || ( abs( .cor ) >= COR.THRESHOLD ) ) {
                              TRUE
                            } else {
                              FALSE
                            }
                          } );
                        .covars.to.exclude <- c( .covars.to.exclude,
                            names( which( .covar.is.too.highly.correlated.with.upper.bound ) ) );
                     } # End if there are any remaining covars to consider removing.
                    } # End if include.bounds.in.lasso
                                              
                  # Exclude covars that are too highly correlated with each other.
                  .covars.to.consider <-
                    setdiff( colnames( .mat1 ), c( .covars.to.exclude, .estimate.colname, .lower.bound.colname, .upper.bound.colname ) );
                  ## Process these in reverse order to ensure that we prioritize keeping those towards the top.
                  .covars.to.consider <- rev( .covars.to.consider );
                  .new.covars.to.exclude <- rep( FALSE, length( .covars.to.consider ) );
                  names( .new.covars.to.exclude ) <- .covars.to.consider;
                  for( .c.i in 1:length( .covars.to.consider ) ) {
                    .covar.colname <- .covars.to.consider[ .c.i ];
                    #print( .covar.colname );
                    # Only consider those not already excluded.
                    .cor <- 
                      cor( .mat1[ , .covar.colname ], .mat1[ , names( which( !.new.covars.to.exclude ) ) ], use = "pairwise" );
                    #print( .cor );
                        if( ( length( .cor ) > 0 ) && any( .cor[ !is.na( .cor ) ] < 1 & .cor[ !is.na( .cor ) ] >= COR.THRESHOLD ) ) {
                          .new.covars.to.exclude[ .c.i ] <- TRUE;
                        } else {
                          .new.covars.to.exclude[ .c.i ] <- FALSE;
                        }
                  } # End foreach of the .covars.to.consider

                  .covars.to.exclude <- c( .covars.to.exclude, names( which( .new.covars.to.exclude ) ) );
                  .retained.covars <- setdiff( colnames( .mat1 ), .covars.to.exclude );
                  if( ( .estimate.colname == "none" ) || ( .estimate.colname %in% .retained.covars ) ) {
                    # penalty.factor = 0 to force the .estimate.colname variable.
                    
                    tryCatch( {
                      cv.glmnet.fit <- cv.glmnet( .mat1, .out, family = "binomial",
                                                 penalty.factor = as.numeric( colnames( .mat1 ) != .estimate.colname ) );
                    ## TODO: REMOVE
                    ##print( coef( cv.glmnet.fit, s = "lambda.min" ) );
                      if( return.lasso.coefs ) {
                        .lasso.validation.estimates.is.one.founder.per.person.coefs.cell <-
                            coef( cv.glmnet.fit, s = "lambda.min" );
                      
                        .lasso.validation.estimates.is.one.founder.per.person.coefs.row[[ .col.i ]] <-
                            .lasso.validation.estimates.is.one.founder.per.person.coefs.cell;
                      }
                      # Ok so now put the predicted values at their appropriate places in the matrix.
                      for( .row.i in the.rows.for.ptid ) {
                        .pred.value.lasso <- predict( cv.glmnet.fit, newx = as.matrix( regression.df[ .row.i, colnames( .mat1 ), drop = FALSE ] ), s = "lambda.min", type = "response");
                                lasso.validation.estimates.is.one.founder.per.person[ .row.i, .col.i ] <- 
                                    .pred.value.lasso;
                      }
                     },
                     error = function( e ) {
                        if( .estimate.colname == "none" ) {
                            warning( paste( "lasso failed with error", e, "\nReverting to simple regression with only an intrecept." ) );
                            .formula <- as.formula( paste( "is.one.founder ~ 1" ) );
                        } else {
                            warning( paste( "lasso failed with error", e, "\nReverting to simple regression vs", .estimate.colname ) );
                        }
                        .formula <- as.formula( paste( "is.one.founder ~", .estimate.colname ) );
                                        # Ok so now put the predicted values at their appropriate places in the matrix.
                        .lm <- glm( .formula, family = "binomial", data = regression.df.without.ptid.i );
                        ## TODO: REMOVE
                        ##print( "fallback from lasso to glm\n", summary( .lm ) ); 
                        for( .row.i in the.rows.for.ptid ) {
                            .pred.value.lasso <-
                                predict( .lm, regression.df[ .row.i, , drop = FALSE ], type = "response" );
                            lasso.validation.estimates.is.one.founder.per.person[ .row.i, .col.i ] <-
                                .pred.value.lasso;
                        }
                      },
                      finally = {
                      }
                    );
                  } # End if the estimate variable is usable
                } # End if use.lasso.validate
    
            } # End foreach .col.i
            if( use.lasso.validate && return.lasso.coefs ) {
                lasso.validation.estimates.is.one.founder.per.person.coefs[[ .ptid.i ]] <- 
                    .lasso.validation.estimates.is.one.founder.per.person.coefs.row;
            } # End if use.lasso.validate
        } # End foreach .ptid.i
        if( use.glm.validate ) {
            colnames( glm.validation.estimates.is.one.founder.per.person ) <-
                paste( "glm.validation.results", estimate.cols, sep = "." );
            rownames( glm.validation.estimates.is.one.founder.per.person ) <-
                rownames( regression.df );
            estimates.is.one.founder.per.person <-
                cbind( estimates.is.one.founder.per.person,
                      glm.validation.estimates.is.one.founder.per.person );
        }
        if( use.lasso.validate ) {
            colnames( lasso.validation.estimates.is.one.founder.per.person ) <-
                paste( "lasso.validation.results", estimate.cols, sep = "." );
            rownames( lasso.validation.estimates.is.one.founder.per.person ) <-
                rownames( regression.df );
            estimates.is.one.founder.per.person <-
                cbind( estimates.is.one.founder.per.person,
                      lasso.validation.estimates.is.one.founder.per.person );
        }
      } # End if use.glm.validate || use.lasso.validate

      # If any column is all NAs, remove it.
      estimates.is.one.founder.per.person <-
          estimates.is.one.founder.per.person[ , apply( estimates.is.one.founder.per.person, 2, function ( .col ) { !all( is.na( .col ) ) } ), drop = FALSE ];
            
      ## For fairness in evaluating when some methods
      ## completely fail to give a result, (so it's NA
      ## presently), we change all of these estimates
      ## from NA to 1 (meaning yes, it's single-founder).
                                        #stopifnot( sum( is.na( estimates.is.one.founder.per.person ) ) == 0 );
      if( sum( is.na( estimates.is.one.founder.per.person ) ) > 0 ) {
        warning( paste( "NAs!", sum( is.na( estimates.is.one.founder.per.person ) ) ) );
        estimates.is.one.founder.per.person.oneNAs <- apply( estimates.is.one.founder.per.person, 1:2, function( .value ) {
            if( is.na( .value ) ) {
                1
            } else {
                .value
            }
        } );
        estimates.is.one.founder.per.person <- estimates.is.one.founder.per.person.oneNAs;
      } # End if any are NA, replace with 1 (meaning single-founder).
       
      ## EVALUATION
        
        isMultiple.aucs <- 
            sapply( 1:ncol( estimates.is.one.founder.per.person ), function( .col.i ) {
              #print( .col.i );
              #print( as.numeric( estimates.is.one.founder.per.person[ , .col.i ] ) );
                if( sum( sapply( as.numeric( estimates.is.one.founder.per.person[ , .col.i ] ), function ( .x ) { !is.null( .x ) && !is.na( .x ) } ) ) > 1 ) {
                  
                  performance( prediction( as.numeric( estimates.is.one.founder.per.person[ , .col.i ] ), gold.is.one.founder.per.person ), measure = "auc" )@y.values[[ 1 ]];
                } else {
                  print( paste( "COL", .col.i, "--> WARNING: sum( sapply( as.numeric( estimates.is.one.founder.per.person[ , .col.i ] ), function ( .x ) { !is.null( .x ) && !is.na( .x ) } ) ) is", sum( sapply( as.numeric( estimates.is.one.founder.per.person[ , .col.i ] ), function ( .x ) { !is.null( .x ) && !is.na( .x ) } ) ) ) );
                  0;
                }
            } );
       names( isMultiple.aucs ) <- colnames( estimates.is.one.founder.per.person );

       isMultiple.aucs.list <- list( AUC = isMultiple.aucs );
       if( use.lasso.validate && return.lasso.coefs ) {
           isMultiple.aucs.list <- c( isMultiple.aucs.list, list( lasso.coefs = lasso.validation.estimates.is.one.founder.per.person.coefs ) );
       }
       
       results.list <- list( unbounded = isMultiple.aucs.list );

       ## If there are multi-timepoint and multi-region predictors, also include per-region, per-timepoint diffs.by.stat results.
       if( !is.null( ppt.suffices ) ) {
           unique.ppt.suffices <- unique( ppt.suffices );
           .results.by.suffix <- lapply( unique.ppt.suffices, function ( .ppt.suffix ) {
               ## TODO: REMOVE
               #print( paste( "suffix:", .ppt.suffix ) );
               .isMultiple.aucs <- 
                   sapply( 1:ncol( estimates.is.one.founder.per.person ), function( .col.i ) {
                                        #print( .col.i );
                     if( sum( sapply( as.numeric( estimates.is.one.founder.per.person[ ppt.suffices == .ppt.suffix, .col.i ] ), function( .x ) { !is.null( .x ) && !is.na( .x ) } ) ) > 1 ) {
                       performance( prediction( as.numeric( estimates.is.one.founder.per.person[ ppt.suffices == .ppt.suffix, .col.i ] ), gold.is.one.founder.per.person[ ppt.suffices == .ppt.suffix ] ), measure = "auc" )@y.values[[ 1 ]];
                     } else {
                       print( paste( "COL", .col.i, "--> WARNING: sum( sapply( as.numeric( estimates.is.one.founder.per.person[ ppt.suffices == .ppt.suffix, .col.i ] ), function( .x ) { !is.null( .x ) && !is.na( .x ) } ) ) is", sum( sapply( as.numeric( estimates.is.one.founder.per.person[ ppt.suffices == .ppt.suffix, .col.i ] ), function( .x ) { !is.null( .x ) && !is.na( .x ) } ) ) ) );
                       0;
                     }
                   } );
               names( .isMultiple.aucs ) <- colnames( estimates.is.one.founder.per.person );

               .isMultiple.aucs.list <- list( AUC = .isMultiple.aucs );
       
               return( .isMultiple.aucs.list );
           } );
           names( .results.by.suffix ) <- paste( "unbounded", unique.ppt.suffices, sep = "" );
           results.list <- c( results.list, .results.by.suffix );
       } # End if there are suffices, also include results by suffix.

        return( results.list );
    } # bound.and.evaluate.is.multiple.results.per.ppt (..)
    
    get.is.multiple.results.for.region.and.time <- function ( the.region, the.time, partition.size ) {
        identify.founders.tab.file <- paste( RESULTS.DIR, results.dirname, the.region, the.time, "identify_founders.tab", sep = "/" );
        stopifnot( file.exists( identify.founders.tab.file ) );

        cat( identify.founders.tab.file, fill = T );
        
        if( length( grep( "^(.*?)\\/[^\\/]+$", identify.founders.tab.file ) ) == 0 ) {
            identify.founders.tab.file.path <- ".";
        } else {
            identify.founders.tab.file.path <-
                gsub( "^(.*?)\\/[^\\/]+$", "\\1", identify.founders.tab.file );
        }
        identify.founders.tab.file.short <-
            gsub( "^.*?\\/?([^\\/]+?)$", "\\1", identify.founders.tab.file, perl = TRUE );
        identify.founders.tab.file.short.nosuffix <-
            gsub( "^([^\\.]+)(\\..+)?$", "\\1", identify.founders.tab.file.short, perl = TRUE );
        identify.founders.tab.file.suffix <-
            gsub( "^([^\\.]+)(\\..+)?$", "\\2", identify.founders.tab.file.short, perl = TRUE );
    
        ## identify-founders results
        identify.founders.study <- readIdentifyFounders( identify.founders.tab.file );
        single.colnames <- grep( "\\.is\\.|fits", colnames( identify.founders.study ), perl = TRUE, value = TRUE );

        ## _All_ of these are "is.one.founder", so note that in the comparison they need to be reversed.
        estimates.is.one.founder <-
            identify.founders.study[ , single.colnames, drop = FALSE ];

        if( is.na( partition.size ) ) {
            ## Now the issue is that there are multiple input files per ppt, eg for the NFLGs ther are often "LH" and "RH" files.  What to do?  The number of sequences varies.  Do a weighted average.
            ## Sometimes there are multiple entries for one ptid/sample, eg for the NFLGs there are often right half and left half (RH and LH) and sometimes additionally NFLG results.  If so, for something to be called single founder, all of the estimates (across regions, for a given test) must agree that it is one founder.
            estimates.is.one.founder.per.person <-
                compute.estimates.is.one.founder.per.person( estimates.is.one.founder );
            
            if( the.region == "v3" ) {
                gold.is.one.founder.per.person <-
                    1-caprisa002.gold.is.multiple[ rownames( estimates.is.one.founder.per.person ) ];
            } else {
                stopifnot( ( the.region == "nflg" ) || ( length( grep( "rv217", the.region ) ) > 0 ) );
                gold.is.one.founder.per.person <-
                    1-rv217.gold.is.multiple[ rownames( estimates.is.one.founder.per.person ) ];
            }
            
            # Note that we use the pvl at the earliest time ie for "1m6m" we use timepoint 2.
            if( the.region == "nflg" || ( length( grep( "rv217", the.region ) ) > 0 ) ) {
                pvl.at.the.time <- sapply( rownames( identify.founders.study ), function( .ptid ) { as.numeric( as.character( rv217.pvl.in[ ( rv217.pvl.in[ , "ptid" ] == .ptid ) & ( rv217.pvl.in[ , "timepoint" ] == ifelse( the.time == "6m", 3, 2 ) ), "viralload" ] ) ) } );
            } else {
                stopifnot( the.region == "v3" );
                pvl.at.the.time <- sapply( rownames( identify.founders.study ), function( .ptid ) { as.numeric( as.character( caprisa002.pvl.in[ ( caprisa002.pvl.in[ , "ptid" ] == .ptid ) & ( caprisa002.pvl.in[ , "timepoint" ] == ifelse( the.time == "6m", 3, 2 ) ), "viralload" ] ) ) } );
            }
            ## Add log plasma viral load (lPVL).
            identify.founders.study.with.lPVL <- cbind( identify.founders.study, log( pvl.at.the.time ) );
            colnames( identify.founders.study.with.lPVL )[ ncol( identify.founders.study.with.lPVL ) ] <- "lPVL";
            
            results.covars.per.person.with.extra.cols <-
              summarizeCovariatesOnePerParticipant( identify.founders.study.with.lPVL );
           
          if( use.bounds ) {
              the.artificial.bounds <- getArtificialBounds( the.region, the.time, results.dirname );
              # Only keep the "sampledwidth" bounds.
              the.artificial.bounds <-
                the.artificial.bounds[ grep( "sampledwidth", names( the.artificial.bounds ), value = TRUE ) ];

              return( list( results.per.person = estimates.is.one.founder.per.person, gold.is.one.founder.per.person = gold.is.one.founder.per.person, results.covars.per.person.with.extra.cols = results.covars.per.person.with.extra.cols, bounds = the.artificial.bounds, evaluated.results = bound.and.evaluate.is.multiple.results.per.ppt( estimates.is.one.founder.per.person, gold.is.one.founder.per.person, results.covars.per.person.with.extra.cols, the.time, the.artificial.bounds ) ) );
          } else {
              return( list( estimates.is.one.founder.per.person = estimates.is.one.founder.per.person, gold.is.one.founder.per.person = gold.is.one.founder.per.person, results.covars.per.person.with.extra.cols = results.covars.per.person.with.extra.cols, evaluated.results = bound.and.evaluate.is.multiple.results.per.ppt( estimates.is.one.founder.per.person, gold.is.one.founder.per.person, results.covars.per.person.with.extra.cols, the.time ) ) );
          }
        } else { # else !is.na( partition.size )
            ## ERE I AM.
            stop( "NOT IMPLEMENTED: EVALUATING PARTITIONS FOR ISMULTIPLE" );
        }
    } # get.is.multiple.results.for.region.and.time ( the.region, the.time, partition.size )

    getIsMultipleResultsByRegionAndTime <- function ( partition.size = NA ) {
        getResultsByRegionAndTime( gold.standard.varname = "gold.is.one.founder.per.person", get.results.for.region.and.time.fn = get.is.multiple.results.for.region.and.time, evaluate.results.per.person.fn = bound.and.evaluate.is.multiple.results.per.ppt, partition.size = partition.size, regions = regions, times = times )
    } # getIsMultipleResultsByRegionAndTime (..)
    
    if( force.recomputation || !file.exists( is.multiple.results.by.region.and.time.Rda.filename ) ) {
        results.by.region.and.time <- getIsMultipleResultsByRegionAndTime();
        save( results.by.region.and.time, file = is.multiple.results.by.region.and.time.Rda.filename );
    } else {
        # loads results.by.region.and.time
        load( file = is.multiple.results.by.region.and.time.Rda.filename );
    }

    writeResultsTables( results.by.region.and.time, "_evaluateIsMultiple.tab", regions = regions, results.are.bounded = TRUE );

    if( FALSE ) {
      # AUC is 1.0 for 1m, nflg, "Synonymous.PFitter.is.starlike"-- but even for "none" it's good (AUC 0.965!)
      .uses <- sapply( 1: length( results.by.region.and.time[[1]][[2]][[ 5 ]][[1]][[2]] ), function( .i ) { .dgCMatrix <- ( results.by.region.and.time[[1]][[2]][[ 5 ]][[1]][[2]][[.i]][[ "Synonymous.PFitter.is.starlike" ]] ); .rv <- as.logical( .dgCMatrix ); names( .rv ) <- rownames( .dgCMatrix ); return( .rv ) } ); .use.threshold.for.inclusion <- 0.2; .vars.to.include <- names( which( apply( .uses, 1, mean ) > .use.threshold.for.inclusion ) );
      .vars.to.include
# [1] "(Intercept)"                       "lPVL"                             
# [3] "diversity"                         "priv.sites"                       
# [5] "inf.sites.clusters"                "DS.Starphy.R"                     
# [7] "sampledwidth_uniform_mtn003.upper" "Synonymous.PFitter.is.starlike"
      .uses <- sapply( 1: length( results.by.region.and.time[[1]][[2]][[ 5 ]][[1]][[2]] ), function( .i ) { .dgCMatrix <- ( results.by.region.and.time[[1]][[2]][[ 5 ]][[1]][[2]][[.i]][[ "none" ]] ); .rv <- as.logical( .dgCMatrix ); names( .rv ) <- rownames( .dgCMatrix ); return( .rv ) } ); .use.threshold.for.inclusion <- 0.2; .vars.to.include <- names( which( apply( .uses, 1, mean ) > .use.threshold.for.inclusion ) );
      .vars.to.include
# [1] "(Intercept)"                         
# [2] "lPVL"                                
# [3] "diversity"                           
# [4] "priv.sites"                          
# [5] "DS.Starphy.R"                        
# [6] "multifounder.Synonymous.DS.StarPhy.R"
# [7] "sampledwidth_uniform_mtn003.upper"
## Now for the multi-region, multi-time results, where we can get an AUC of 0.933 using lasso with "multifounder.PFitter.mut.rate.coef", but even "none" gets 0.927.
      #sort( ( results.by.region.and.time[[3]][[1]][[1]][[1]][[5]][[1]][[1]] ))
      .uses <- sapply( 1: length( results.by.region.and.time[[3]][[1]][[1]][[1]][[ 5 ]][[1]][[2]] ), function( .i ) { .dgCMatrix <- ( results.by.region.and.time[[3]][[1]][[1]][[1]][[ 5 ]][[1]][[2]][[.i]][[ "multifounder.PFitter.mut.rate.coef" ]] ); .rv <- as.logical( .dgCMatrix ); names( .rv ) <- rownames( .dgCMatrix ); return( .rv ) } ); .use.threshold.for.inclusion <- 0.2; .vars.to.include <- names( which( apply( .uses, 1, mean ) > .use.threshold.for.inclusion ) );
      .vars.to.include
# [1] "(Intercept)"                                  
# [2] "lPVL"                                         
# [3] "diversity"                                    
# [4] "priv.sites"                                   
# [5] "inf.sites.clusters"                           
# [6] "DS.Starphy.R"                                 
# [7] "multifounder.Synonymous.DS.StarPhy.R"         
# [8] "sampledwidth_uniform_1mmtn003_6mhvtn502.upper"
# [9] "multifounder.PFitter.mut.rate.coef"                 
      .uses <- sapply( 1: length( results.by.region.and.time[[3]][[1]][[1]][[1]][[ 5 ]][[1]][[2]] ), function( .i ) { .dgCMatrix <- ( results.by.region.and.time[[3]][[1]][[1]][[1]][[ 5 ]][[1]][[2]][[.i]][[ "none" ]] ); .rv <- as.logical( .dgCMatrix ); names( .rv ) <- rownames( .dgCMatrix ); return( .rv ) } ); .use.threshold.for.inclusion <- 0.2; .vars.to.include <- names( which( apply( .uses, 1, mean ) > .use.threshold.for.inclusion ) );
      .vars.to.include
# [1] "(Intercept)"                                  
# [2] "lPVL"                                         
# [3] "diversity"                                    
# [4] "priv.sites"                                   
# [5] "inf.sites.clusters"                           
# [6] "DS.Starphy.R"                                 
# [7] "multifounder.Synonymous.DS.StarPhy.R"         
# [8] "sampledwidth_uniform_1mmtn003_6mhvtn502.upper"
}
    
    # Return the file name.
    return( output.table.path );
} # evaluateIsMultiple (... )

