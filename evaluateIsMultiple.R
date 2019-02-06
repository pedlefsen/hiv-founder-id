## First do all the stuff in README.postprocessing.txt.
library( "ROCR" ) # for "prediction" and "performance"
library( "glmnet" ) # for cv.glmnet
library( "parallel" ) # for mcapply

source( "readIdentifyFounders_safetosource.R" );
source( "getArtificialBounds_safetosource.R" );
source( "getResultsByRegionAndTime_safetosource.R" );
source( "writeResultsTables_safetosource.R" );
source( "summarizeCovariatesOnePerParticipant_safetosource.R" );

GOLD.STANDARD.DIR <- "/fh/fast/edlefsen_p/bakeoff/gold_standard/";
#RESULTS.DIR <- "/fh/fast/edlefsen_p/bakeoff_analysis_results/";
RESULTS.DIR <- "/fh/fast/edlefsen_p/bakeoff_merged_analysis_sequences_results/";

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
#' @param results.dirname the subdirectory of RESULTS.DIR
#' @param force.recomputation if FALSE (default) and if there is a saved version called isMultiple.results.by.region.and.time.Rda (under bakeoff_analysis_results/results.dirname), then that file will be loaded; otherwise the results will be recomputed and saved in that location.
#' @param partition.bootstrap.seed the random seed to use when bootstrapping samples by selecting one partition number per ptid, repeatedly; we do it this way because there are an unequal number of partitions, depending on sampling depth.
#' @param partition.bootstrap.samples the number of bootstrap replicates to conduct; the idea is to get an estimate of the variation in estimates and results (errors) across these samples.
#' @param partition.bootstrap.num.cores the number of cores to run the boostrap replicates on (defaults to all of the cores returned by parallel::detectCores()).
#' @return the filename of the Rda output. If you load( filename ), it will add "results.by.region.and.time" to your environment.
#' @export

evaluateIsMultiple <- function (
                             use.bounds = TRUE,
                             use.glm.validate = TRUE,
                             use.lasso.validate = TRUE,
                             include.bounds.in.glm = TRUE,
                             include.bounds.in.lasso = TRUE,
                             include.helpful.additional.cols.in.glm = !include.bounds.in.glm,
                             helpful.additional.cols = c( "diversity", "priv.sites", "multifounder.Synonymous.DS.Starphy.R", "DS.Starphy.R", "inf.sites.clusters", "lPVL" ),
                             #results.dirname = "raw_edited_20160216",
                             results.dirname = "raw_fixed",
                             force.recomputation = TRUE,
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
    rv217.gold.standards.in <- read.csv( paste( GOLD.STANDARD.DIR, "rv217/RV217_gold_standards.csv", sep = "" ) );
    rv217.gold.is.multiple <- rv217.gold.standards.in[ , "gold.is.multiple" ];
    names( rv217.gold.is.multiple ) <- rv217.gold.standards.in[ , "ptid" ];
    
    caprisa002.gold.standards.in <- read.csv( paste( GOLD.STANDARD.DIR, "caprisa_002/caprisa_002_gold_standards.csv", sep = "" ) );
    caprisa002.gold.is.multiple <- caprisa002.gold.standards.in[ , "gold.is.multiple" ];
    names( caprisa002.gold.is.multiple ) <- caprisa002.gold.standards.in[ , "ptid" ];

    rv217.pvl.in <- read.csv( paste( GOLD.STANDARD.DIR, "rv217/rv217_viralloads.csv", sep = "" ) );
    caprisa002.pvl.in <- read.csv( paste( GOLD.STANDARD.DIR, "caprisa_002/caprisa_002_viralloads.csv", sep = "" ) );
    
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

       # Do not allow use of information about the time of sampling.
       if( "6m.not.1m" %in% colnames( results.covars.per.person.with.extra.cols ) ) {
         .forbidden.column.i <- which( "6m.not.1m" == colnames( results.covars.per.person.with.extra.cols ) );
         results.covars.per.person.with.extra.cols <-
           results.covars.per.person.with.extra.cols[ , -.forbidden.column.i, drop = FALSE ];
       }
       
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
        .keep.cols <-
            grep( "Star[pP]hy\\.fits", .keep.cols, value = TRUE, invert = TRUE );
        .keep.cols <-
            grep( "Star[pP]hy\\.founders", .keep.cols, value = TRUE, invert = TRUE );

        # For COB and infer, use only the real-data sources (mtn003 or hvtn502). So exclude the "mtn003" and "sixmonths" ones.
        .keep.cols <-
            grep( "\\.(one|six)months?\\.", .keep.cols, value = TRUE, invert = TRUE );
        ## Try removing some variables that are rarely selected or are too correlated (eg diversity is highly correlated with sd.entropy, max.hd, insites.is.one.founder, insites.founders) [SEE BELOW WHERE WE ADD TO THIS PROGRAMMATICALLY]
           .donotkeep.cols <- c( "multifounder.DS.Starphy.R", "PFitter.chi.sq.stat", "Synonymous.DS.StarPhy.R", "StarPhy.is.one.founder", "DS.Starphy.fits", "DS.Starphy.is.starlike", "insites.is.one.founder", "inf.to.priv.ratio" );
           # NOTE ALSO MORE (formerly excluded below because cor > 0.9):
           .donotkeep.cols <- c( .donotkeep.cols, "StarPhy.founders", "InSites.founders", "PFitter.max.hd",   "PFitter.mean.hd", "sd.entropy",       "mean.entropy",     "inf.sites" );
           # This one is problematic, because it is almost always 1; this is a problem because it breaks cv.glmnet, and also because it effectively acts as an intercept.
           .donotkeep.cols <- c( .donotkeep.cols, "StarPhy.founders", "InSites.founders" );
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
        
        ## Ok build a regression model with an intercept, including only the helpful.additional.cols, and also the lower and upper bounds associated with either 5 weeks or 30 weeks, depending on the.time (if there's a 1m sample, uses "mtn003").
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
                      return( ( sum( !is.na( .col ) ) <= 1 ) || ( var( .col, na.rm = TRUE ) == 0 ) );
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
                        suppressWarnings( glm( .formula, family = "binomial", data = .df ) );
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
                      return( ( sum( !is.na( .col ) ) <= 1 ) || ( var( .col, na.rm = TRUE ) == 0 ) );
                  } );
                  .mat1 <- .mat1[ , setdiff( colnames( .mat1 ), names( which( .covars.to.exclude ) ) ), drop = FALSE ];
                  .covars.to.exclude <- names( which( .covars.to.exclude ) );
                  if( length( .ptids.to.exclude ) > 0 ) {
                      ## TODO: REMOVE
                      #print( paste( "excluding samples (", paste( .ptids.to.exclude, collapse = ", " ), ") due to NAs in", .estimate.colname ) );
                      .mat1 <- .mat1[ -.rows.to.exclude, , drop = FALSE ];
                      .out <- .out[ -.rows.to.exclude ];
                   }
                  
                  # At least one DF is needed for the cv.glmnet, and one for the intercept
                  MINIMUM.DF <- 2; # how much more should nrow( .lasso.mat ) be than ncol( .lasso.mat ) at minimum?
                  
                  MINIMUM.CORRELATION.WITH.OUTCOME <- 0.1;
                  cors.with.the.outcome <-
                    sapply( setdiff( colnames( .mat1 ), c( .covars.to.exclude, .estimate.colname ) ), function( .covar.colname ) {
                          .cor <- 
                              cor( .out, .mat1[ , .covar.colname ], use = "pairwise" );
                            return( .cor );
                        } );
                  sorted.cors.with.the.outcome <-
                    cors.with.the.outcome[ order( abs( cors.with.the.outcome ) ) ];
                  # Sort the columns of .mat1 by their correlation with the outcome. This is to ensure that more-relevant columns get selected when removing columns due to pairwise correlation among them.  We want to keep the best one.
                  if( .estimate.colname == "none" ) {
                    .mat1 <- .mat1[ , rev( names( sorted.cors.with.the.outcome ) ), drop = FALSE ];
                  } else {
                    .mat1 <- .mat1[ , c( .estimate.colname, rev( names( sorted.cors.with.the.outcome ) ) ), drop = FALSE ];
                  }
                  
                  # Exclude covars that are not sufficiently correlated with the outcome.
                      .covars.to.exclude <- c( .covars.to.exclude,
                          names( which( sapply( setdiff( colnames( .mat1 ), c( .covars.to.exclude, .estimate.colname ) ), function( .covar.colname ) {
                            #print( .covar.colname );
                          .cor <- 
                              cor( .out, .mat1[ , .covar.colname ], use = "pairwise" );
                            #print( .cor );
                          if( is.na( .cor ) || ( abs( .cor ) <= MINIMUM.CORRELATION.WITH.OUTCOME ) ) {
                            TRUE
                          } else {
                            FALSE
                          }
                        } ) ) ) );
                  
                    COR.THRESHOLD <- 0.8;
                   # Exclude covars that are too highly correlated with the estimate.
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
                  .needed.df <-
                    ( length( .retained.covars ) - ( nrow( .mat1 ) - MINIMUM.DF ) );
                  if( .needed.df > 0 ) {
                    # Then remove some covars until there is at least MINIMUM.DF degrees of freedom.
                    # They are in order, so just chop them off the end.
                    .mat1 <- .mat1[ , 1:( ncol( .mat1 ) - .needed.df ), drop = FALSE ];
                  }
                  if( ( .estimate.colname == "none" ) || ( .estimate.colname %in% .retained.covars ) ) {
                    # penalty.factor = 0 to force the .estimate.colname variable.
                    
                    tryCatch( {
                      tryCatch( {
                        cv.glmnet.fit <- cv.glmnet( .mat1, .out, family = "binomial", penalty.factor = as.numeric( colnames( .mat1 ) != .estimate.colname ), grouped = FALSE, nfold = length( .out ) );
                      }, error = function( e ) {
                        # First try removing the penalty factor, which for some reason can help even when it doesn't seem to matter.
                        cv.glmnet.fit <- cv.glmnet( .mat1, .out, family = "binomial", grouped = FALSE, nfold = length( .out ) );
                      } );
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
                            warning( paste( "ptid", .ptid.i, "col", .col.i, "lasso failed with error", e, "\nReverting to simple regression with only an intrecept." ) );
                            .formula <- as.formula( paste( "is.one.founder ~ 1" ) );
                        } else {
                            warning( paste( "ptid", .ptid.i, "col", .col.i, "lasso failed with error", e, "\nReverting to simple regression vs", .estimate.colname ) );
                        }
                        .formula <- as.formula( paste( "is.one.founder ~", .estimate.colname ) );
                                        # Ok so now put the predicted values at their appropriate places in the matrix.
                        .lm <- suppressWarnings( glm( .formula, family = "binomial", data = regression.df.without.ptid.i ) );
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
        #warning( paste( "NAs!", sum( is.na( estimates.is.one.founder.per.person ) ) ) );
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
              the.artificial.bounds <- getArtificialBounds( the.region, the.time, RESULTS.DIR, results.dirname );
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
        ## TODO: REMOVE
        #.m <- sapply( 1:length( results.by.region.and.time[[3]][[1]][[1]][["1m.6m"]][[5]][["unbounded"]][[2]] ), function( .i ) { .rv <- as.numeric( results.by.region.and.time[[3]][[1]][[1]][["1m.6m"]][[5]][["unbounded"]][[2]][[.i]][["multifounder.Synonymous.PFitter.is.poisson" ]] ) != 0; names( .rv ) <- rownames( results.by.region.and.time[[3]][[1]][[1]][["1m.6m"]][[5]][["unbounded"]][[2]][[.i]][["multifounder.Synonymous.PFitter.is.poisson" ]] ); return( .rv ); } );apply( .m, 1, mean )
    } else {
        # loads results.by.region.and.time
        load( file = is.multiple.results.by.region.and.time.Rda.filename );
    }

    writeResultsTables( results.by.region.and.time, "_evaluateIsMultiple.tab", regions = regions, results.are.bounded = TRUE );

    if( FALSE ) {
        get.uses <- function ( .varname = "none", regions = c( "nflg", "v3" ), times = c( "1m", "6m" ) ) {
            if( length( regions ) == 2 ) {
                .results.for.region <- results.by.region.and.time[[3]][[1]][[1]];
            } else {
                .results.for.region <- results.by.region.and.time[[ regions ]];
            }
            if( length( times ) == 2 ) {
                the.time <- "1m.6m";
            } else {
                the.time <- times;
            }
            .results.by.removed.ptid <-
                .results.for.region[[ the.time ]][[ "evaluated.results" ]][["unbounded"]][[ "lasso.coefs" ]];
            .uses <- sapply( 1:length( .results.by.removed.ptid ), function( .i ) { .dgCMatrix <- .results.by.removed.ptid[[.i]][[.varname]]; .rv <- as.logical( .dgCMatrix ); names( .rv ) <- rownames( .dgCMatrix ); return( .rv ); } );
            if( is.null( dim( .uses ) ) ) {
              table( names( which( unlist( .uses ) ) ) );
            } else {
              apply( .uses, 1, sum );
            }
        } # get.uses (..)

        ## Overall performance, when training everything together (but note never use of 6m.not.1m)
#        sort( unlist( results.by.region.and.time[[3]][[1]][[1]][[1]][[5]][["unbounded"]][[1]] ) )
#                                    lasso.validation.results.PFitter.is.poisson 
#                                                              0.9095098 
#                         lasso.validation.results.PFitter.mut.rate.coef 
#                                                              0.9111856 
#                           lasso.validation.results.PFitter.is.starlike 
#                                                              0.9237537 
#               lasso.validation.results.multifounder.PFitter.is.poisson 
#                                                              0.9291998 
#              lasso.validation.results.Synonymous.PFitter.mut.rate.coef 
#                                                              0.9308756 
#                 lasso.validation.results.Synonymous.PFitter.is.poisson 
#                                                              0.9359028 
#                lasso.validation.results.Synonymous.PFitter.is.starlike 
#                                                              0.9363217 
#    lasso.validation.results.multifounder.Synonymous.PFitter.is.poisson 
#                                                              0.9388354 
#            lasso.validation.results.multifounder.PFitter.mut.rate.coef 
#                                                              0.9400922 
#                        lasso.validation.results.InSites.is.one.founder 
#                                                              0.9409300 
#                                          lasso.validation.results.none 
#                                                              0.9426058 
# lasso.validation.results.multifounder.Synonymous.PFitter.mut.rate.coef 
#                                                              0.9455383 
        ## 1m evaluated on v3 data from caprisa002, when training everything together (but note never use of 6m.not.1m)
#        sort( unlist( results.by.region.and.time[[3]][[1]][[1]][[1]][[5]][["unbounded.1m.v3"]][[1]] ) )
#              glm.validation.results.multifounder.PFitter.mut.rate.coef 
#                                                              0.9333333 
#                           glm.validation.results.PFitter.mut.rate.coef 
#                                                              0.9733333 
#                         lasso.validation.results.PFitter.mut.rate.coef 
#                                                              0.9733333 
#                                          lasso.validation.results.none 
#                                                              0.9866667 
#                            lasso.validation.results.PFitter.is.poisson 
#                                                              0.9866667 
#                           lasso.validation.results.PFitter.is.starlike 
#                                                              0.9866667 
#                 lasso.validation.results.Synonymous.PFitter.is.poisson 
#                                                              0.9866667 
#                lasso.validation.results.Synonymous.PFitter.is.starlike 
#                                                              0.9866667 
#               lasso.validation.results.multifounder.PFitter.is.poisson 
#                                                              0.9866667 
#    lasso.validation.results.multifounder.Synonymous.PFitter.is.poisson 
#                                                              0.9866667 
#              lasso.validation.results.Synonymous.PFitter.mut.rate.coef 
#                                                              0.9866667 
#                        lasso.validation.results.InSites.is.one.founder 
#                                                              1.0000000 
#            lasso.validation.results.multifounder.PFitter.mut.rate.coef 
#                                                              1.0000000 
# lasso.validation.results.multifounder.Synonymous.PFitter.mut.rate.coef 
#                                                              1.0000000 
#        sort( unlist( results.by.region.and.time[[3]][[1]][[1]][[1]][[5]][["unbounded.1m.nflg"]][[1]] ) )
#   glm.validation.results.multifounder.Synonymous.PFitter.mut.rate.coef 
#                                                              0.9115385 
#                           lasso.validation.results.PFitter.is.starlike 
#                                                              0.9115385 
#               lasso.validation.results.multifounder.PFitter.is.poisson 
#                                                              0.9269231 
#              glm.validation.results.multifounder.PFitter.mut.rate.coef 
#                                                              0.9307692 
#    lasso.validation.results.multifounder.Synonymous.PFitter.is.poisson 
#                                                              0.9538462 
#                lasso.validation.results.Synonymous.PFitter.is.starlike 
#                                                              0.9576923 
#            lasso.validation.results.multifounder.PFitter.mut.rate.coef 
#                                                              0.9653846 
#                                          lasso.validation.results.none 
#                                                              0.9692308 
#                        lasso.validation.results.InSites.is.one.founder 
#                                                              0.9692308 
#              lasso.validation.results.Synonymous.PFitter.mut.rate.coef 
#                                                              0.9692308 
#                 lasso.validation.results.Synonymous.PFitter.is.poisson 
#                                                              0.9730769 
# lasso.validation.results.multifounder.Synonymous.PFitter.mut.rate.coef 
#                                                              0.9730769 
#        sort( unlist( results.by.region.and.time[[3]][[1]][[1]][[1]][[5]][["unbounded.6m.v3"]][[1]] ) )
#              lasso.validation.results.Synonymous.PFitter.mut.rate.coef 
#                                                             0.93055556 
#            lasso.validation.results.multifounder.PFitter.mut.rate.coef 
#                                                             0.93055556 
#                                          lasso.validation.results.none 
#                                                             0.95833333 
#                        lasso.validation.results.InSites.is.one.founder 
#                                                             0.95833333 
#                            lasso.validation.results.PFitter.is.poisson 
#                                                             0.95833333 
#                           lasso.validation.results.PFitter.is.starlike 
#                                                             0.95833333 
#               lasso.validation.results.multifounder.PFitter.is.poisson 
#                                                             0.95833333 
#                         lasso.validation.results.PFitter.mut.rate.coef 
#                                                             0.95833333 
#                           glm.validation.results.PFitter.mut.rate.coef 
#                                                             0.97222222 
#                 lasso.validation.results.Synonymous.PFitter.is.poisson 
#                                                             0.97222222 
#                lasso.validation.results.Synonymous.PFitter.is.starlike 
#                                                             0.97222222 
#    lasso.validation.results.multifounder.Synonymous.PFitter.is.poisson 
#                                                             0.97222222 
# lasso.validation.results.multifounder.Synonymous.PFitter.mut.rate.coef 
#                                                             0.98611111
#         sort( unlist( results.by.region.and.time[[3]][[1]][[1]][[1]][[5]][["unbounded.6m.nflg"]][[1]] ) )        
#                glm.validation.results.Synonymous.PFitter.mut.rate.coef 
#                                                              0.9375000 
# OF NOTE ALSO:
# lasso.validation.results.multifounder.Synonymous.PFitter.mut.rate.coef 
#                                                              0.8666667 
#        get.uses( "multifounder.Synonymous.PFitter.mut.rate.coef" )
#                                           (Intercept) 
#                                            57 
# multifounder.Synonymous.PFitter.mut.rate.coef 
#                                            57 
#                                     diversity 
#                                            57 
#                                    priv.sites 
#                                            57 
#          multifounder.Synonymous.DS.StarPhy.R 
#                                            47 
#                                  DS.Starphy.R 
#                                            57 
# sampledwidth_uniform_1mmtn003_6mhvtn502.lower 
#                                            57 
#                            inf.sites.clusters 
#                                            57 
# sampledwidth_uniform_1mmtn003_6mhvtn502.upper 
#                                            57 
#                                   v3_not_nflg 
#                                            57 
#                                          lPVL 
#                                            56 

### TODO: Describe that model, what does it do to make such good decisions, and why does it go wrong for the nflg.6m results, whereas this one does well: (note I've verified that the named coefs are the same in these two lists (above and below here):
#        get.uses( "Synonymous.PFitter.mut.rate.coef" )
#                                   (Intercept) 
#                                            57 
#              Synonymous.PFitter.mut.rate.coef 
#                                            57 
#                                     diversity 
#                                            57 
#          multifounder.Synonymous.DS.StarPhy.R 
#                                            57 
#                                    priv.sites 
#                                            56 
#                                  DS.Starphy.R 
#                                            57 
# sampledwidth_uniform_1mmtn003_6mhvtn502.lower 
#                                            57 
#                            inf.sites.clusters 
#                                            57 
# sampledwidth_uniform_1mmtn003_6mhvtn502.upper 
#                                            57 
#                                   v3_not_nflg 
#                                            56 
#                                          lPVL 
#                                            56 
#        setdiff( names( get.uses( "multifounder.Synonymous.PFitter.mut.rate.coef" ) ), "multifounder.Synonymous.PFitter.mut.rate.coef" )
#  [1] "(Intercept)"                                  
#  [2] "diversity"                                    
#  [3] "priv.sites"                                   
#  [4] "multifounder.Synonymous.DS.StarPhy.R"         
#  [5] "DS.Starphy.R"                                 
#  [6] "sampledwidth_uniform_1mmtn003_6mhvtn502.lower"
#  [7] "inf.sites.clusters"                           
#  [8] "sampledwidth_uniform_1mmtn003_6mhvtn502.upper"
#  [9] "v3_not_nflg"                                  
# [10] "lPVL"
#        evaluate.specific.isMultiple.model( model.vars = c( "multifounder.Synonymous.PFitter.mut.rate.coef", "diversity", "priv.sites", "multifounder.Synonymous.DS.StarPhy.R", "DS.Starphy.R", "sampledwidth_uniform_1mmtn003_6mhvtn502.lower", "inf.sites.clusters", "sampledwidth_uniform_1mmtn003_6mhvtn502.upper", "v3_not_nflg", "lPVL" ) );
#         Call:
# glm(formula = .formula, family = "binomial", data = regression.df)
# 
# Deviance Residuals: 
#      Min        1Q    Median        3Q       Max  
# -3.16057  -0.00125   0.04350   0.15376   1.45252  
# 
# Coefficients:
#                                                 Estimate Std. Error z value
# (Intercept)                                   -8.206e+00  4.896e+00  -1.676
# multifounder.Synonymous.PFitter.mut.rate.coef -2.214e+03  1.676e+03  -1.321
# diversity                                     -1.705e+03  6.415e+02  -2.658
# priv.sites                                    -2.100e-02  2.228e-02  -0.943
# multifounder.Synonymous.DS.StarPhy.R           7.537e-01  6.115e+00   0.123
# DS.Starphy.R                                   3.757e+00  3.259e+00   1.153
# sampledwidth_uniform_1mmtn003_6mhvtn502.lower  2.387e-02  1.375e-02   1.736
# inf.sites.clusters                             5.742e-01  5.078e-01   1.131
# sampledwidth_uniform_1mmtn003_6mhvtn502.upper  7.701e-03  7.843e-03   0.982
# v3_not_nflg                                    6.369e+00  3.279e+00   1.943
# lPVL                                           1.024e+00  5.327e-01   1.923
#                                               Pr(>|z|)   
# (Intercept)                                    0.09368 . 
# multifounder.Synonymous.PFitter.mut.rate.coef  0.18657   
# diversity                                      0.00787 **
# priv.sites                                     0.34591   
# multifounder.Synonymous.DS.StarPhy.R           0.90191   
# DS.Starphy.R                                   0.24892   
# sampledwidth_uniform_1mmtn003_6mhvtn502.lower  0.08260 . 
# inf.sites.clusters                             0.25816   
# sampledwidth_uniform_1mmtn003_6mhvtn502.upper  0.32619   
# v3_not_nflg                                    0.05207 . 
# lPVL                                           0.05448 . 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 128.807  on 106  degrees of freedom
# Residual deviance:  26.948  on  96  degrees of freedom
#   (1 observation deleted due to missingness)
# AIC: 48.948
# 
# Number of Fisher Scoring iterations: 9
#        evaluate.specific.isMultiple.model( model.vars, step = TRUE )
# Step:  AIC=44.79
# is.one.founder ~ multifounder.Synonymous.PFitter.mut.rate.coef + 
#     diversity + DS.Starphy.R + sampledwidth_uniform_1mmtn003_6mhvtn502.lower + 
#     inf.sites.clusters + v3_not_nflg + lPVL
# 
#                                                 Df Deviance     AIC
# <none>                                               28.788  44.788
# - inf.sites.clusters                             1   30.975  44.975
# - multifounder.Synonymous.PFitter.mut.rate.coef  1   31.098  45.098
# - lPVL                                           1   32.392  46.392
# - sampledwidth_uniform_1mmtn003_6mhvtn502.lower  1   32.849  46.849
# - DS.Starphy.R                                   1   35.362  49.362
# - v3_not_nflg                                    1   48.157  62.157
# - diversity                                      1   92.318 106.318
# 
# Call:
# glm(formula = is.one.founder ~ multifounder.Synonymous.PFitter.mut.rate.coef + 
#     diversity + DS.Starphy.R + sampledwidth_uniform_1mmtn003_6mhvtn502.lower + 
#     inf.sites.clusters + v3_not_nflg + lPVL, family = "binomial", 
#     data = regression.df)
# 
# Deviance Residuals: 
#      Min        1Q    Median        3Q       Max  
# -2.58921  -0.00305   0.05761   0.16963   1.74489  
# 
# Coefficients:
#                                                 Estimate Std. Error z value
# (Intercept)                                   -5.699e+00  3.424e+00  -1.665
# multifounder.Synonymous.PFitter.mut.rate.coef -2.273e+03  1.759e+03  -1.292
# diversity                                     -1.586e+03  5.312e+02  -2.986
# DS.Starphy.R                                   5.626e+00  2.605e+00   2.160
# sampledwidth_uniform_1mmtn003_6mhvtn502.lower  2.008e-02  1.070e-02   1.876
# inf.sites.clusters                             6.346e-01  4.787e-01   1.326
# v3_not_nflg                                    6.830e+00  2.771e+00   2.465
# lPVL                                           6.694e-01  3.938e-01   1.700
#                                               Pr(>|z|)   
# (Intercept)                                    0.09601 . 
# multifounder.Synonymous.PFitter.mut.rate.coef  0.19624   
# diversity                                      0.00282 **
# DS.Starphy.R                                   0.03080 * 
# sampledwidth_uniform_1mmtn003_6mhvtn502.lower  0.06071 . 
# inf.sites.clusters                             0.18496   
# v3_not_nflg                                    0.01370 * 
# lPVL                                           0.08915 . 
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 128.807  on 106  degrees of freedom
# Residual deviance:  28.788  on  99  degrees of freedom
#   (1 observation deleted due to missingness)
# AIC: 44.788
# 
# Number of Fisher Scoring iterations: 8
        
        ##
# evaluate.specific.model( c( "diversity", "v3_not_nflg", "sampledwidth_uniform_1mmtn003_6mhvtn502.lower", "lPVL" ) )
# 
# Call:
# glm(formula = .formula, family = "binomial", data = regression.df)
# 
# Deviance Residuals: 
#      Min        1Q    Median        3Q       Max  
# -2.26627  -0.00246   0.17435   0.38541   2.18787  
# 
# Coefficients:
#                                                 Estimate Std. Error z value
# (Intercept)                                   -9.624e-01  2.336e+00  -0.412
# diversity                                     -1.348e+03  2.952e+02  -4.566
# v3_not_nflg                                    4.873e+00  1.453e+00   3.355
# sampledwidth_uniform_1mmtn003_6mhvtn502.lower  1.553e-02  6.834e-03   2.272
# lPVL                                           4.605e-01  2.492e-01   1.848
#                                               Pr(>|z|)    
# (Intercept)                                   0.680363    
# diversity                                     4.97e-06 ***
# v3_not_nflg                                   0.000793 ***
# sampledwidth_uniform_1mmtn003_6mhvtn502.lower 0.023096 *  
# lPVL                                          0.064554 .  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 129.487  on 107  degrees of freedom
# Residual deviance:  47.513  on 103  degrees of freedom
# AIC: 57.513
# 
# Number of Fisher Scoring iterations: 7
#evaluate.specific.model( c( "multifounder.Synonymous.PFitter.mut.rate.coef", "diversity", "v3_not_nflg", "sampledwidth_uniform_1mmtn003_6mhvtn502.lower", "lPVL" ) )
# 
# Call:
# glm(formula = .formula, family = "binomial", data = regression.df)
# 
# Deviance Residuals: 
#      Min        1Q    Median        3Q       Max  
# -2.20021  -0.01191   0.09426   0.28203   1.78930  
# 
# Coefficients:
#                                                 Estimate Std. Error z value
# (Intercept)                                   -3.223e+00  2.690e+00  -1.198
# multifounder.Synonymous.PFitter.mut.rate.coef -2.803e+03  1.221e+03  -2.296
# diversity                                     -1.208e+03  3.007e+02  -4.016
# v3_not_nflg                                    5.246e+00  1.658e+00   3.164
# sampledwidth_uniform_1mmtn003_6mhvtn502.lower  2.134e-02  8.190e-03   2.606
# lPVL                                           7.908e-01  3.265e-01   2.422
#                                               Pr(>|z|)    
# (Intercept)                                    0.23081    
# multifounder.Synonymous.PFitter.mut.rate.coef  0.02170 *  
# diversity                                     5.92e-05 ***
# v3_not_nflg                                    0.00156 ** 
# sampledwidth_uniform_1mmtn003_6mhvtn502.lower  0.00917 ** 
# lPVL                                           0.01544 *  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 128.807  on 106  degrees of freedom
# Residual deviance:  40.271  on 101  degrees of freedom
#   (1 observation deleted due to missingness)
# AIC: 52.271
# 
# Number of Fisher Scoring iterations: 8
#evaluate.specific.model( c( "multifounder.Synonymous.PFitter.mut.rate.coef", "diversity", "v3_not_nflg", "sampledwidth_uniform_1mmtn003_6mhvtn502.upper", "lPVL" ) )
#### Conclusion: upper is better than lower!!
# Call:
# glm(formula = .formula, family = "binomial", data = regression.df)
# 
# Deviance Residuals: 
#      Min        1Q    Median        3Q       Max  
# -3.03133  -0.01292   0.09300   0.26800   2.35990  
# 
# Coefficients:
#                                                 Estimate Std. Error z value
# (Intercept)                                   -3.780e+00  2.630e+00  -1.437
# multifounder.Synonymous.PFitter.mut.rate.coef -2.794e+03  1.223e+03  -2.284
# diversity                                     -1.246e+03  3.271e+02  -3.807
# v3_not_nflg                                    5.726e+00  1.849e+00   3.098
# sampledwidth_uniform_1mmtn003_6mhvtn502.upper  1.700e-02  5.920e-03   2.872
# lPVL                                           6.845e-01  2.855e-01   2.398
#                                               Pr(>|z|)    
# (Intercept)                                    0.15072    
# multifounder.Synonymous.PFitter.mut.rate.coef  0.02235 *  
# diversity                                      0.00014 ***
# v3_not_nflg                                    0.00195 ** 
# sampledwidth_uniform_1mmtn003_6mhvtn502.upper  0.00408 ** 
# lPVL                                           0.01650 *  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# 
# (Dispersion parameter for binomial family taken to be 1)
# 
#     Null deviance: 128.807  on 106  degrees of freedom
# Residual deviance:  35.984  on 101  degrees of freedom
#   (1 observation deleted due to missingness)
# AIC: 47.984
# 
# Number of Fisher Scoring iterations: 8

    } # END IF FALSE
    
    # Return the file name.
    return( is.multiple.results.by.region.and.time.Rda.filename );
} # evaluateIsMultiple (... )

