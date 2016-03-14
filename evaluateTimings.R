## First do all the stuff in README.postprocessing.txt.

library( "parallel" ); # for mclapply
library( "glmnet" ); # for cv.glmnet

source( "readIdentifyFounders_safetosource.R" );
source( "getDaysSinceInfection_safetosource.R" );
source( "summarizeCovariatesOnePerParticipant_safetosource.R" );

#' Evaluate timings estimates and produce results tables.
#'
#' This function runs the BakeOff results analysis for the timings results.
#'
#' The "center of bounds" (COB) approach is the way that we do it at the
#' VTN, and the way it was done in RV144, etc: use the midpoint
#' between the bounds on the actual infection time computed from the
#' dates and results of the HIV positivity tests (antibody or PCR).
#' The typical approach is to perform antibody testing every X days
#' (historically this is 6 months in most HIV vaccine trials, except
#' during the vaccination phase there are more frequent visits and on
#' every visit HIV testing is conducted).  The (fake) bounds used here
#' are calculated in the createArtificialBoundsOnInfectionDate.R file.
#' The actual bounds would be too tight, since the participants were
#' detected HIV+ earlier in these people than what we expect to see in
#' a trial in which testing is conducted every X days.  For the center
#' of bounds approach we load the bounds files in subdirs of the
#' "bounds" subdirectory eg at
#' /fh/fast/edlefsen_p/bakeoff/analysis_sequences/bounds/nflg/1m/.
#' These files have names beginning with "artificialBounds_" and
#' ending with ".tab".
#'
#' @param use.bounds compute results for the COB approach, and also return evaluations of bounded versions of the other results.
#' @param use.infer compute results for the PREAST approach.
#' @param use.anchre compute results for the anchre approach (presently disabled).
#' @param use.glm.validate evaluate predicted values from leave-one-out cross-validation.
#' @param include.bounds.in.glm include the corresponding prior bounds in the regression equations used for cross-validation.
#' @param include.helpful.additional.cols.in.glm include "priv.sites" and "multifounder.Synonymous.PFitter.is.poisson" in the regression equations used for cross-validation.
#' @param results.dirname the subdirectory of "/fh/fast/edlefsen_p/bakeoff/analysis_sequences" and also of "/fh/fast/edlefsen_p/bakeoff_analysis_results"
#' @param force.recomputation if FALSE (default) and if there is a saved version called timings.results.by.region.and.time.Rda (under bakeoff_analysis_results/results.dirname), then that file will be loaded; otherwise the results will be recomputed and saved in that location.
#'
#' @return NULL
#' @export

evaluateTimings <- function (
                             use.bounds = TRUE,
                             use.infer = TRUE,
                             use.anchre = TRUE,
                             use.glm.validate = TRUE,
                             use.lasso.validate = TRUE,
                             include.bounds.in.glm = TRUE,
                             include.bounds.in.lasso = TRUE,
                             include.helpful.additional.cols.in.glm = !include.bounds.in.glm,
                             results.dirname = "raw_edited_20160216",
                             force.recomputation = FALSE,
                             partition.bootstrap.seed = 98103,
                             partition.bootstrap.samples = 100,
                             partition.bootstrap.num.cores = detectCores()
                            )
{
    timings.results.by.region.and.time.Rda.filename <-
        paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/timings.results.by.region.and.time.Rda", sep = "" );
    
    identify.founders.date.estimates <-
        c( "PFitter.time.est", "Synonymous.PFitter.time.est", "multifounder.PFitter.time.est", "multifounder.Synonymous.PFitter.time.est" ); 
    
    rv217.gold.standard.infection.dates.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/rv217_gold_standard_timings.csv" );
    rv217.gold.standard.infection.dates <- as.Date( as.character( rv217.gold.standard.infection.dates.in[,2] ), "%m/%d/%y" );
    names( rv217.gold.standard.infection.dates ) <- as.character( rv217.gold.standard.infection.dates.in[,1] );
    
    caprisa002.gold.standard.infection.dates.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/caprisa_002_gold_standard_timings.csv" );
    caprisa002.gold.standard.infection.dates <- as.Date( as.character( caprisa002.gold.standard.infection.dates.in[,2] ), "%Y/%m/%d" );
    names( caprisa002.gold.standard.infection.dates ) <- as.character( caprisa002.gold.standard.infection.dates.in[,1] );
    
    regions <- c( "nflg", "v3", "rv217_v3" );
    times <- c( "1m", "6m", "1m6m" );
    #times <- c( "1m6m" );
    
    rmse <- function( x, na.rm = FALSE ) {
        if( na.rm ) {
            x <- x[ !is.na( x ) ];
        }
        return( sqrt( mean( x ** 2 ) ) );
    }
    
    compute.results.one.per.ppt <- function ( results, weights ) {
        apply( results, 2, function ( .column ) {
            .rv <- 
            sapply( unique( rownames( results ) ), function( .ppt ) {
                .ppt.cells <- .column[ rownames( results ) == .ppt ];
                if( all( is.na( .ppt.cells ) ) ) {
                  return( NA );
                }
                .ppt.weights <- weights[ rownames( results ) == .ppt ];
                .ppt.weights <- .ppt.weights / sum( .ppt.weights, na.rm = TRUE );
                sum( .ppt.cells * .ppt.weights, na.rm = TRUE );
            } );
            names( .rv ) <- unique( rownames( results ) );
            return( .rv );
        } );
    } # compute.results.one.per.ppt

    get.infer.results.columns <- function ( the.region, the.time, the.ptids, partition.size = NA ) {
        ## Add to results: "infer" results.
        if( ( the.region == "v3" ) || ( the.region == "rv217_v3" ) ) {
            the.region.dir <- "v3_edited_20160216";
        } else {
            the.region.dir <- "nflg_copy_20160222";
        }
        if( is.na( partition.size ) ) {
            infer.results.directories <- dir( paste( "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/", results.dirname, "/", the.region.dir, "/", the.time, sep = "" ), "founder-inference-bakeoff_", full.name = TRUE );
        } else {
            #infer.results.directories <- dir( paste( "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/", results.dirname, "/", the.region.dir, "/", the.time, "/partitions", sep = "" ), "founder-inference-bakeoff_", full.name = TRUE );
            stop( "TODO: When there are some Infer results run on partitions, evaluate the results here" );
        }
        # Special: for v3, separate out the caprisa seqs from the rv217 seqs
        if( the.region == "v3" ) {
             infer.results.directories <- grep( "_100\\d\\d\\d$", infer.results.directories, value = TRUE );
        } else if( the.region == "rv217_v3" ) {
             infer.results.directories <- grep( "_100\\d\\d\\d$", infer.results.directories, value = TRUE, invert = TRUE );
        }
          
        infer.results.files <- sapply( infer.results.directories, dir, "outtoi.csv", full.name = TRUE );
        infer.results.list <-
            lapply( unlist( infer.results.files ), function( .file ) {
                .rv <- as.matrix( read.csv( .file, header = FALSE ), nrow = 1 );
                stopifnot( ncol( .rv ) == 3 );
                return( .rv );
            } );
        names( infer.results.list ) <- unlist( infer.results.files );
        
        if( length( infer.results.list ) == 0 ) {
            return( NULL );
        }
        infer.results <- do.call( rbind, infer.results.list );
        colnames( infer.results ) <- c( "Infer", "Infer.CI.high", "Infer.CI.low" );
        ## reorder them
        infer.results <- infer.results[ , c( "Infer", "Infer.CI.low", "Infer.CI.high" ), drop = FALSE ];
        infer.results.bounds.ptid <- gsub( "^.+_(\\d+)/.+$", "\\1", names( infer.results.list ) );
        rownames( infer.results ) <- infer.results.bounds.ptid;
        
        ## Separate it into separate tables by bounds type.
        infer.results.bounds.type <-
          gsub( "^.+_artificialBounds_(.+)_\\d+/.+$", "\\1", names( infer.results.list ) );
        infer.results.bounds.type[ grep( "csv$", infer.results.bounds.type ) ] <- NA;
        infer.results.nobounds.table <- infer.results[ is.na( infer.results.bounds.type ), , drop = FALSE ];
        infer.results.bounds.types <- setdiff( unique( infer.results.bounds.type ), NA );
        infer.results.bounds.tables <- lapply( infer.results.bounds.types, function ( .bounds.type ) {
          return( infer.results[ !is.na( infer.results.bounds.type ) & ( infer.results.bounds.type == .bounds.type ), , drop = FALSE ] );
        } );
        names( infer.results.bounds.tables ) <- infer.results.bounds.types;
        
        # Add just the estimates from infer (not the CIs) to the results table.
        new.results.columns <-
          matrix( NA, nrow = length( the.ptids ), ncol = 1 + length( infer.results.bounds.tables ) );
        rownames( new.results.columns ) <- the.ptids;
        colnames( new.results.columns ) <- c( "Infer.time.est", infer.results.bounds.types );
        
        .shared.ptids.nobounds <-
            intersect( rownames( new.results.columns ), rownames( infer.results.nobounds.table ) );
        .result.ignored <- sapply( .shared.ptids.nobounds, function( .ptid ) {
            .infer.subtable <-
                infer.results.nobounds.table[ rownames( infer.results.nobounds.table ) == .ptid, 1, drop = FALSE ];
            stopifnot( sum( rownames( infer.results.nobounds.table ) == .ptid ) == nrow( .infer.subtable ) );
            # But there might be fewer of these than there are rows in the results.table (if eg there are 3 input fasta files and infer results for only 2 of them).
            stopifnot( nrow( .infer.subtable ) <= sum( rownames( new.results.columns ) == .ptid ) );
            if( nrow( .infer.subtable ) < sum( rownames( new.results.columns ) == .ptid ) ) {
              .infer.subtable <- rbind( .infer.subtable, matrix( NA, nrow = ( sum( rownames( new.results.columns ) == .ptid ) - nrow( .infer.subtable ) ), ncol = ncol( .infer.subtable ) ) );
            }
            new.results.columns[ rownames( new.results.columns ) == .ptid, 1 ] <<-
                .infer.subtable;
            return( NULL );
        } );
        .result.ignored <- sapply( names( infer.results.bounds.tables ), function ( .bounds.type ) {
          .shared.ptids <-
              intersect( rownames( new.results.columns ), rownames( infer.results.bounds.tables[[ .bounds.type ]] ) );
          ..result.ignored <- sapply( .shared.ptids, function( .ptid ) {
              .infer.subtable <-
                infer.results.bounds.tables[[ .bounds.type ]][ rownames( infer.results.bounds.tables[[ .bounds.type ]] ) == .ptid, 1, drop = FALSE ];
              stopifnot( sum( rownames( infer.results.bounds.tables[[ .bounds.type ]] ) == .ptid ) == nrow( .infer.subtable ) );
              # But there might be fewer of these than there are rows in the results.table (if eg there are 3 input fasta files and infer results for only 2 of them).
              stopifnot( nrow( .infer.subtable ) <= sum( rownames( new.results.columns ) == .ptid ) );
              if( nrow( .infer.subtable ) < sum( rownames( new.results.columns ) == .ptid ) ) {
                .infer.subtable <- rbind( .infer.subtable, matrix( NA, nrow = ( sum( rownames( new.results.columns ) == .ptid ) - nrow( .infer.subtable ) ), ncol = ncol( .infer.subtable ) ) );
              }
              new.results.columns[ rownames( new.results.columns ) == .ptid, .bounds.type ] <<-
                  .infer.subtable;
              return( NULL );
          } );                       
          return( NULL );
        } );
          
        colnames( new.results.columns ) <- c( "Infer.time.est", paste( "Infer", gsub( "_", ".", infer.results.bounds.types ), "time.est", sep = "." ) );
        return( new.results.columns );
    } # get.infer.results.columns (..)

    get.anchre.results.columns <- function ( the.region, the.time, sample.dates.in, partition.size = NA ) {
        stopifnot( is.na( partition.size ) ); # TODO: Implement support for anchre on partitions.
        stopifnot( the.time == "1m6m" ); # There's only anchre results for longitudinal data.
        ## Add to results: "anchre" results. (only at 1m6m)
        anchre.results.directories <- dir( paste( "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/", results.dirname, "/", the.region.dir, "/1m6m", sep = "" ), "anchre", full.name = TRUE );
        if( length( anchre.results.directories ) == 0 ) {
            return( NULL );
        }
        anchre.results.files <-
            sapply( anchre.results.directories, dir, "mrca.csv", full.name = TRUE );
        anchre.results <- do.call( rbind,
            lapply( unlist( anchre.results.files ), function( .file ) {
                stopifnot( file.exists( .file ) );
                .file.short <-
                    gsub( "^.*?\\/?([^\\/]+?)$", "\\1", .file, perl = TRUE );
                .file.short.nosuffix <-
                    gsub( "^([^\\.]+)(\\..+)?$", "\\1", .file.short, perl = TRUE );
                .file.converted <-
                    paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region.dir, "/1m6m/", .file.short.nosuffix, ".anc2tsv.tab", sep = "" );
                # convert it.
                system( paste( "./anc2tsv.sh", .file, ">", .file.converted ) );
                stopifnot( file.exists( .file.converted ) );
                .rv <- as.matrix( read.delim( .file.converted, header = TRUE, sep = "\t" ), nrow = 1 );
                ## No negative dates!  Just call it NA.
                .rv <- apply( .rv, 1:2, function( .str ) { if( length( grep( "^-", .str ) ) > 0 ) { NA } else { .str } } );
                stopifnot( ncol( .rv ) == 4 );
                return( .rv );
            } ) );
        colnames( anchre.results ) <- c( "Anchre.r2t.est", "Anchre.est", "Anchre.CI.low", "Anchre.CI.high" );
        rownames( anchre.results ) <-
            gsub( "^.+_(\\d+)$", "\\1", names( unlist( anchre.results.files ) ) );
        # Special: for v3, only use caprisa seqs (not rv217, for now).
        if( the.region == "v3" ) {
            anchre.results <-
                anchre.results[ grep( "^100\\d\\d\\d", rownames( anchre.results ) ), , drop = FALSE ];
        } else if( the.region == "rv217_v3" ) {
            anchre.results <-
                anchre.results[ grep( "^100\\d\\d\\d", rownames( anchre.results ), invert = TRUE ), , drop = FALSE ];
        }
        # Add just the estimate from anchre.
        sample.dates <- as.Date( as.character( sample.dates.in[ , 2 ] ) );
        names( sample.dates ) <- sample.dates.in[ , 1 ];
        anchre.r2t.days.before.sample <- sapply( 1:nrow( anchre.results ), function( .i ) { 0 - as.numeric( as.Date( anchre.results[ .i, 1 ] ) - sample.dates[ rownames( anchre.results )[ .i ] ] ) } );
        names( anchre.r2t.days.before.sample ) <- rownames( anchre.results );
        anchre.days.before.sample <- sapply( 1:nrow( anchre.results ), function( .i ) { 0 - as.numeric( as.Date( anchre.results[ .i, 2 ] ) - sample.dates[ rownames( anchre.results )[ .i ] ] ) } );
        names( anchre.days.before.sample ) <- rownames( anchre.results );
        
        anchre.columns <- cbind( anchre.r2t.days.before.sample, anchre.days.before.sample );
        colnames( anchre.columns ) <- c( "Anchre.r2t.time.est", "Anchre.bst.time.est" );
        return( anchre.columns );
    } # get.anchre.results.columns (..)

    
    compute.diffs.by.stat <- function ( results.one.per.ppt, days.since.infection ) {
        diffs.by.stat <- 
            lapply( colnames( results.one.per.ppt ), function( .stat ) {
                .rv <- ( as.numeric( results.one.per.ppt[ , .stat ] ) - as.numeric( days.since.infection[ rownames( results.one.per.ppt ) ] ) );
                names( .rv ) <- rownames( results.one.per.ppt );
                return( .rv );
            } );
        names( diffs.by.stat ) <- colnames( results.one.per.ppt );
        return( diffs.by.stat );
    } # compute.diffs.by.stat ( results.one.per.ppt, days.since.infection )

    bound.and.evaluate.results.per.ppt <-
        function ( results.one.per.ppt, days.since.infection, results.covars.one.per.ppt.with.extra.cols, the.artificial.bounds = NA ) {
        ## Also include all of the date estimates, which is
        ## everything in "results.one.per.ppt" so
        ## far. (However we will exclude some bounded
        ## results with deterministic bounds from the
        ## optimization, see below). It's also redundant to
        ## use PFitter-based days estimates, since we are
        ## using the mutation rate coefs, which are a
        ## linear fn of the corresponding days ests (I have
        ## confirmed that the predictions are the same
        ## (within +- 0.6 days, due to the rounding that
        ## PFitter does to its days estimates).  So
        ## basically unless we added non-PFitter results
        ## (PREAST/infer, anchre, or center-of-bounds), there won't be
        ## anything to do here.
        days.est.cols <- colnames( results.one.per.ppt );
        days.est.cols <- grep( "deterministic", days.est.cols, invert = TRUE, value = TRUE );
        days.est.cols <- grep( "PFitter|Star[Pp]y", days.est.cols, invert = TRUE, perl = TRUE, value = TRUE );
    
      if( use.glm.validate || use.lasso.validate ) {
        results.covars.one.per.ppt.with.extra.cols <-
            cbind( results.one.per.ppt[ , days.est.cols, drop = FALSE ], results.covars.one.per.ppt.with.extra.cols );
        
        .keep.cols <-
            grep( "num.*\\.seqs|totalbases", colnames( results.covars.one.per.ppt.with.extra.cols ), value = TRUE, perl = TRUE, invert = TRUE );
        ### TODO: Something else.  Just trying to get down to a reasonable set; basically there are very highly clustered covariates here and it screws up the inference.
        ## Also remove all of the mut.rate.coef except for multifounder.Synonymous.PFitter.mut.rate.coef.
        
        # mut.rate.coef.keep.cols <- c( "multifounder.Synonymous.PFitter.mut.rate.coef",
        #                 grep( "mut\\.rate\\.coef", mut.rate.coef.keep.cols, invert = TRUE, value = TRUE ) )
        #mut.rate.coef.keep.cols <- c( "multifounder.Synonymous.PFitter.mut.rate.coef", "inf.to.priv.ratio", "priv.sites", "inf.sites.clusters", "InSites.founders", "multifounder.Synonymous.PFitter.is.poisson" );
        ## Keep only the mut.rate.coef cols and priv.sites and multifounder.Synonymous.PFitter.is.poisson.
        mut.rate.coef.cols <- grep( "mut\\.rate\\.coef", .keep.cols, value = TRUE );
        all.additional.cols <- setdiff( .keep.cols, mut.rate.coef.cols );
        helpful.additional.cols <- c( "priv.sites","multifounder.Synonymous.PFitter.is.poisson" );
          
        if( use.lasso.validate ) {
            keep.cols <- c( all.additional.cols, mut.rate.coef.cols, days.est.cols );
            estimate.cols <- setdiff( keep.cols, all.additional.cols );
        } else {
            keep.cols <- c( helpful.additional.cols, mut.rate.coef.cols, days.est.cols );
            estimate.cols <- setdiff( keep.cols, helpful.additional.cols );
        }
          
        results.covars.one.per.ppt <-
            results.covars.one.per.ppt.with.extra.cols[ , keep.cols, drop = FALSE ];
        results.covars.one.per.ppt.df <-
            data.frame( results.covars.one.per.ppt );
        
        regression.df <- cbind( data.frame( days.since.infection = days.since.infection[ rownames( results.covars.one.per.ppt.df ) ] ), results.covars.one.per.ppt.df, lapply( the.artificial.bounds, function( .mat ) { .mat[ rownames( results.covars.one.per.ppt.df ), , drop = FALSE ] } ) );
        
        ## Ok build a regression model with no intercept, including only the helpful.additional.cols, and also the lower and upper bounds associated with either 5 weeks or 30 weeks, depending on the.time (if there's a 1m sample, uses "5weeks").
        if( the.time == "6m" ) {
            .lower.bound.colname <- "uniform_30weeks.lower";
            .upper.bound.colname <- "uniform_30weeks.upper";
        } else {
            .lower.bound.colname <- "uniform_5weeks.lower";
            .upper.bound.colname <- "uniform_5weeks.upper";
        }
        
        if( use.glm.validate ) {
            glm.validation.results.one.per.ppt <- matrix( NA, nrow = nrow( results.covars.one.per.ppt.df ), ncol = length( estimate.cols ) );
        }
        if( use.lasso.validate ) {
            lasso.validation.results.one.per.ppt <- matrix( NA, nrow = nrow( results.covars.one.per.ppt.df ), ncol = length( estimate.cols ) );
        }
        for( .row.i in 1:nrow( regression.df ) ) {
            regression.df.without.row.i <-
                regression.df[ -.row.i, , drop = FALSE ];
            for( .col.i in 1:length( estimate.cols ) ) {
                .estimate.colname <- estimate.cols[ .col.i ];
                if( use.glm.validate ) {
                  # covariates for glm
                  .covariates.glm <- c();
                  if( include.helpful.additional.cols.in.glm ) {
                      .covariates.glm <-
                          c( .covariates.glm, helpful.additional.cols );
                  }
                  if( include.bounds.in.glm ) {
                      .covariates.glm <-
                          c( .covariates.glm, .lower.bound.colname, .upper.bound.colname );
                  }
                  # glm:
                  .covars.to.exclude <- apply( regression.df.without.row.i, 2, function ( .col ) {
                      return( ( var( .col ) == 0 ) || ( sum( !is.na( .col ) ) <= 1 ) );
                  } );
                  .retained.covars <- setdiff( colnames( regression.df.without.row.i ), names( which( .covars.to.exclude ) ) );
                  if( .estimate.colname %in% .retained.covars ) {
                    .df <- regression.df.without.row.i[ , .retained.covars, drop = FALSE ];
                    .formula <- as.formula( paste( "days.since.infection ~ 0 + ", paste( intersect( .retained.covars, c( .covariates.glm, .estimate.colname ) ), collapse = "+" ) ) );
                    .pred.value.glm <- predict( lm( .formula, data = regression.df.without.row.i ), regression.df[ .row.i, , drop = FALSE ] );
                    glm.validation.results.one.per.ppt[ .row.i, .col.i ] <- 
                        .pred.value.glm;
                  } # If the estimate can be used
    
                } # End if use.glm.validate
    
                if( use.lasso.validate ) {
                  # covariates for lasso
                  .covariates.lasso <- c( all.additional.cols );
                  if( include.bounds.in.lasso ) {
                      .covariates.lasso <-
                          c( .covariates.lasso, .lower.bound.colname, .upper.bound.colname );
                  }
                  # lasso:
                  .mat1 <- as.matrix( regression.df.without.row.i[ , c( .covariates.lasso, .estimate.colname ) ] );
                  .out <- regression.df.without.row.i[[ "days.since.infection" ]];
    
                  .covars.to.exclude <- apply( .mat1, 2, function ( .col ) {
                      return( ( var( .col ) == 0 ) || ( sum( !is.na( .col ) ) <= 1 ) );
                  } );
                  .retained.covars <- setdiff( colnames( .mat1 ), names( which( .covars.to.exclude ) ) );
                  if( .estimate.colname %in% .retained.covars ) {
                    .mat1 <- .mat1[ , setdiff( colnames( .mat1 ), names( which( .covars.to.exclude ) ) ), drop = FALSE ];
                    # penalty.factor = 0 to force the .estimate.colname variable.
    
                    tryCatch( {
                    cv.glmnet.fit <- cv.glmnet( .mat1, .out, intercept = FALSE,
                                               penalty.factor = as.numeric( colnames( .mat1 ) != .estimate.colname ) );
                    .pred.value.lasso <- predict( cv.glmnet.fit, newx = as.matrix( regression.df[ .row.i, colnames( .mat1 ), drop = FALSE ] ), s = "lambda.min" );
                            lasso.validation.results.one.per.ppt[ .row.i, .col.i ] <<- 
                        .pred.value.lasso;
                     },
                     error = function( e ) {
                         warning( paste( "lasso failed with error", e, "\nReverting to simple regression vs", .estimate.colname ) );
                         .formula <- as.formula( paste( "days.since.infection ~ 0 + ", .estimate.colname ) );
                         .pred.value.lasso <- predict( lm( .formula, data = regression.df.without.row.i ), regression.df[ .row.i, , drop = FALSE ] );
                            lasso.validation.results.one.per.ppt[ .row.i, .col.i ] <<- 
                        .pred.value.lasso;
                     },
                        finally = {
                        }
                        );
                  } # End if the estimate variable is usable
                } # End if use.lasso.validate
    
            } # End foreach .col.i
        } # End foreach .row.i
        if( use.glm.validate ) {
            colnames( glm.validation.results.one.per.ppt ) <-
                paste( "glm.validation.results", estimate.cols, sep = "." );
            rownames( glm.validation.results.one.per.ppt ) <-
                rownames( regression.df );
            results.one.per.ppt <-
                cbind( results.one.per.ppt,
                      glm.validation.results.one.per.ppt );
        }
        if( use.lasso.validate ) {
            colnames( lasso.validation.results.one.per.ppt ) <-
                paste( "lasso.validation.results", estimate.cols, sep = "." );
            rownames( lasso.validation.results.one.per.ppt ) <-
                rownames( regression.df );
            results.one.per.ppt <-
                cbind( results.one.per.ppt,
                      lasso.validation.results.one.per.ppt );
        }
      } # End if use.glm.validate || use.lasso.validate
      
      ## For fairness in evaluating when some methods
      ## completely fail to give a result, (so it's NA
      ## presently), we change all of these estimates
      ## from NA to 0.  When we put bounds on the
      ## results, below, they will be changed from 0 to
      ## a boundary endpoint if 0 is outside of the
      ## bounds.
      results.one.per.ppt.zeroNAs <- apply( results.one.per.ppt, 1:2, function( .value ) {
          if( is.na( .value ) ) {
              0
          } else {
              .value
          }
      } );
      
      ## unbounded results:
      diffs.by.stat.zeroNAs <- compute.diffs.by.stat( results.one.per.ppt.zeroNAs, days.since.infection );
      diffs.by.stat <- compute.diffs.by.stat( results.one.per.ppt, days.since.infection );
    
      unbounded.results <- 
        list( bias = lapply( diffs.by.stat, mean, na.rm = T ), se = lapply( diffs.by.stat, sd, na.rm = T ), rmse = lapply( diffs.by.stat, rmse, na.rm = T ), n = lapply( diffs.by.stat, function( .vec ) { sum( !is.na( .vec ) ) } ), bias.zeroNAs = lapply( diffs.by.stat.zeroNAs, mean, na.rm = T ), se.zeroNAs = lapply( diffs.by.stat.zeroNAs, sd, na.rm = T ), rmse.zeroNAs = lapply( diffs.by.stat.zeroNAs, rmse, na.rm = T ), n.zeroNAs = lapply( diffs.by.stat.zeroNAs, function( .vec ) { sum( !is.na( .vec ) ) } ) );
      
      ## bounded results:
      if( use.bounds ) {
    
        ## Ok, well, now we can also evaluate versions
        ## of variants of each method, but bounding the
        ## results.  Ie if the estimated time is within
        ## the bounds, that time is used, otherwise,
        ## the boundary.  Note we don't do this with
        ## the deterministic bounds, and we only do it
        ## for the time corresponding to the sample (
        ## 5weeks for "1m" and 30weeks for "6m" ) [each
        ## gets an additional week for the difference
        ## between the 2 weeks added at the beginning
        ## and 1 week subtracted at the end, for
        ## eclipse phase; also this accounts for a bit
        ## (~10%) additional variation in the time
        ## between visits at 6 months than at 1-2
        ## months -- a totally made up number as at
        ## this time I have no idea what the right
        ## number is, but this seems reasonable.]
        .artificial.bounds.to.use <-
            grep( ifelse( the.time == "6m", "30weeks", "5weeks" ), grep( "deterministic", names( the.artificial.bounds ), invert = TRUE, value = TRUE ), value = TRUE );
        results.one.per.ppt.bounded <-
          lapply( .artificial.bounds.to.use, function ( .artificial.bounds.name ) {
            .mat <-
            apply( results.one.per.ppt, 2, function ( .results.column ) {
            sapply( names( .results.column ), function ( .ppt ) {
              .value <- .results.column[ .ppt ];
              if( !( .ppt %in% rownames( the.artificial.bounds[[ .artificial.bounds.name ]] ) ) ) {
                  stop( paste( .ppt, "has no bounds!" ) );
              }
              if( !is.na( .value ) ) {
                .value.is.below.lb <- ( .value < the.artificial.bounds[[ .artificial.bounds.name ]][ .ppt, "lower" ] );
                if( .value.is.below.lb ) {
                  return( the.artificial.bounds[[ .artificial.bounds.name ]][ .ppt, "lower" ] );
                }
                .value.is.above.ub <- ( .value > the.artificial.bounds[[ .artificial.bounds.name ]][ .ppt, "upper" ] );
                if( .value.is.above.ub ) {
                  return( the.artificial.bounds[[ .artificial.bounds.name ]][ .ppt, "upper" ] );
                }
              }
              return( .value );
            } );
          } );
            rownames( .mat ) <- rownames( results.one.per.ppt );
            return( .mat );
        } );
        names( results.one.per.ppt.bounded ) <- .artificial.bounds.to.use;
    
        bounded.results.by.bound.type <- lapply( results.one.per.ppt.bounded, function ( .results.one.per.ppt ) {
          .diffs.by.stat <- compute.diffs.by.stat( .results.one.per.ppt, days.since.infection );
    
          return( list( bias = lapply( .diffs.by.stat, mean, na.rm = T ), se = lapply( .diffs.by.stat, sd, na.rm = T ), rmse = lapply( .diffs.by.stat, rmse, na.rm = T ), n = lapply( .diffs.by.stat, function( .vec ) { sum( !is.na( .vec ) ) } ) ) );
        } );
        return( c( list( unbounded = unbounded.results ), bounded.results.by.bound.type ) );
      } else {
        return( list( unbounded = unbounded.results ) );
      }
    } # bound.and.evaluate.results.per.ppt (..)

    get.timings.results.for.region.and.time <- function ( the.region, the.time, partition.size ) {
        .days.since.infection.filename <- paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/sampleDates.tbl", sep = "" );
        if( the.region == "v3" ) {
            days.since.infection <-
                getDaysSinceInfection(
                    .days.since.infection.filename,
                    caprisa002.gold.standard.infection.dates
                );
        } else {
            stopifnot( ( the.region == "nflg" ) || ( length( grep( "rv217", the.region ) ) > 0 ) );
            days.since.infection <-
                getDaysSinceInfection(
                    .days.since.infection.filename,
                    rv217.gold.standard.infection.dates
                );
        }
            
        ## identify-founders results; we always get and use these.
        if( is.na( partition.size ) ) {
            results.in <- readIdentifyFounders( paste( paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/identify_founders.tab", sep = "" ) ) );
        } else {
            results.in <- readIdentifyFounders( paste( paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/partitions/identify_founders.tab", sep = "" ) ), partition.size = partition.size );
        }
        
        results.covars.one.per.ppt.with.extra.cols <-
          summarizeCovariatesOnePerParticipant( results.in );
        ## TODO: Add to that the loading of the viralloads.csv files (in the gold_standards dirs).
        
        results <- results.in;
        
        days.colnames <- c( grep( "time", colnames( results ), value = T ), grep( "days", colnames( results ), value = T ) );
        
        days.est.colnames <- grep( "est", days.colnames, value = TRUE );
        days.est <- results[ , days.est.colnames, drop = FALSE ];
        lambda.est.colnames <-
            gsub( "PFitter\\.lambda\\.est", "PFitter.lambda", gsub( "(?:days|time|fits)", "lambda", days.est.colnames, perl = TRUE ) );
        stopifnot( all( lambda.est.colnames %in% colnames( results ) ) );
        days.est.colnames.nb <- gsub( "[^\\.]+\\.Star[Pp]hy", "PFitter", gsub( "(?:days|time|fits).*$", "nbases", days.est.colnames, perl = TRUE ) );
        days.est.nb <- results[ , days.est.colnames.nb, drop = FALSE ];
        days.est.colnames.nseq <- gsub( "[^\\.]+\\.Star[Pp]hy", "PFitter", gsub( "(?:days|time|fits).*$", "nseq", days.est.colnames, perl = TRUE ) );
        days.est.nseq <- results[ , days.est.colnames.nseq, drop = FALSE ];
        
        results <- results[ , days.est.colnames, drop = FALSE ];
        
        if( use.infer && is.na( partition.size ) ) { ## TODO: process infer results on partitions.
          infer.results.columns <- get.infer.results.columns( the.region, the.time, rownames( results ), partition.size );
          results <- cbind( results, infer.results.columns );
        } # End if use.infer
        
        if( use.anchre && ( the.time == "1m6m" ) && is.na( partition.size ) ) {  ## TODO: Handle the anchre results for the partitions
          if( the.region == "v3" ) {
              sample.dates.in <-
                  getDaysSinceInfection(
                      .days.since.infection.filename,
                      caprisa002.gold.standard.infection.dates,
                      return.sample.dates.in = TRUE
                  );
          } else {
              stopifnot( ( the.region == "nflg" ) || ( length( grep( "rv217", the.region ) ) > 0 ) );
              sample.dates.in <-
                  getDaysSinceInfection(
                      .days.since.infection.filename,
                      rv217.gold.standard.infection.dates,
                      return.sample.dates.in = TRUE
                  );
          }
          anchre.results.columns <-
              get.anchre.results.columns( the.region, the.time, sample.dates.in, partition.size );
            
          results <- cbind( results, anchre.results.columns );
        } # End if use.ancher and the.time is 1m6m, add anchre results too.
       
        if( is.na( partition.size ) ) {
          ## Now the issue is that there are multiple input files per ppt, eg for the NFLGs ther are often "LH" and "RH" files.  What to do?  The number of sequences varies.  Do a weighted average.
          .weights <- days.est.nseq*days.est.nb;
          results.one.per.ppt <-
              compute.results.one.per.ppt( results, .weights );
        
          if( use.bounds ) {
            ## The "center of bounds" approach is the way that
            ## we do it at the VTN, and the way it was done in
            ## RV144, etc: use the midpoint between the bounds
            ## on the actual infection time computed from the
            ## dates and results of the HIV positivity tests
            ## (antibody or PCR).  The typical approach is to
            ## perform antibody testing every X days
            ## (historically this is 6 months in most HIV
            ## vaccine trials, except during the vaccination
            ## phase there are more frequent visits and on
            ## every visit HIV testing is conducted).  The
            ## (fake) bounds used here are calculated in the
            ## createArtificialBoundsOnInfectionDate.R file.
            ## The actual bounds would be too tight, since the
            ## participants were detected HIV+ earlier in these
            ## people than what we expect to see in a trial in
            ## which testing is conducted every X days.  For
            ## the center of bounds approach we load the bounds
            ## files in subdirs of the "bounds" subdirectory eg
            ## at
            ## /fh/fast/edlefsen_p/bakeoff/analysis_sequences/bounds/nflg/1m/.
            ## These files have names beginning with
            ## "artificialBounds_" and ending with ".tab".
            bounds.subdirname <- "bounds";
            .artificial.bounds.dirname <-
                paste( "/fh/fast/edlefsen_p/bakeoff/analysis_sequences/", results.dirname, "/", bounds.subdirname, "/", the.region, "/", the.time, "/", sep = "" );
            artificial.bounds.filenames <-
                dir( .artificial.bounds.dirname, pattern = "artificialBounds_.*.tab", recursive = FALSE, full.names = TRUE );
            names( artificial.bounds.filenames ) <- gsub( "^.*artificialBounds_(.*).tab$", "\\1", artificial.bounds.filenames );
            
            the.artificial.bounds <- lapply( names( artificial.bounds.filenames ), function ( .artificial.bounds.name ) {
                .tbl <- read.table( artificial.bounds.filenames[[ .artificial.bounds.name ]], header = TRUE, sep = "\t" );
                # Special: for v3, only use caprisa seqs (not rv217, for now).
                if( the.region == "v3" ) {
                    .tbl <-
                        .tbl[ grep( "^100\\d\\d\\d", rownames( .tbl ) ), , drop = FALSE ];
                } else if( the.region == "rv217_v3" ) {
                    .tbl <-
                        .tbl[ grep( "^100\\d\\d\\d", rownames( .tbl ), invert = TRUE ), , drop = FALSE ];
                }
                
              return( .tbl );
            } );
            names( the.artificial.bounds ) <- names( artificial.bounds.filenames );
            
            center.of.bounds.table <- sapply( names( the.artificial.bounds ), function ( .artificial.bounds.name ) {
                round( apply( the.artificial.bounds[[ .artificial.bounds.name ]], 1, mean ) )
            } );
            colnames( center.of.bounds.table ) <-
              paste( "COB", gsub( "_", ".", colnames( center.of.bounds.table ) ), "time.est", sep = "." );
              
            results.one.per.ppt <-
              cbind( results.one.per.ppt, center.of.bounds.table[ rownames( results.one.per.ppt ), , drop = FALSE ] );
            
              return( list( results.one.per.ppt = results.one.per.ppt, days.since.infection = days.since.infection, results.covars.one.per.ppt.with.extra.cols = results.covars.one.per.ppt.with.extra.cols, bounds = the.artificial.bounds, evaluated.results = bound.and.evaluate.results.per.ppt( results.one.per.ppt, days.since.infection, results.covars.one.per.ppt.with.extra.cols, the.artificial.bounds ) ) );
          } else {
              return( list( results.one.per.ppt = results.one.per.ppt, days.since.infection = days.since.infection, results.covars.one.per.ppt.with.extra.cols = results.covars.one.per.ppt.with.extra.cols, evaluated.results = bound.and.evaluate.results.per.ppt( results.one.per.ppt, days.since.infection, results.covars.one.per.ppt.with.extra.cols ) ) );
          }
        } else { # else !is.na( partition.size )
            ## Here the multiple results per participant come from the partitions.  We want to evaluate each one, and summarize them afterwards.
            partition.id <- results.in[ , "partition.id" ];
            
            diffs.per.ppt.by.id <- apply( results, 2, function ( .column ) {
              .rv <- 
                lapply( unique( rownames( results ) ), function( .ppt ) {
                    .values <- .column[ rownames( results ) == .ppt ];
                    .partition.ids <- partition.id[ rownames( results ) == .ppt ];
                    sapply( unique( .partition.ids ), function( the.partition.id ) {
                        ..values <- .values[ .partition.ids == the.partition.id ];
                        return( as.numeric( ..values ) - as.numeric( days.since.infection[ .ppt ] ) );
                    } );
              } );
              names( .rv ) <- unique( rownames( results ) );
              return( .rv );
            } );
            
            mean.bias.per.ppt <- lapply( diffs.per.ppt.by.id, function( .lst ) { sapply( .lst, mean, na.rm = T ) } );
            sd.bias.per.ppt <- lapply( diffs.per.ppt.by.id, function( .lst ) { sapply( .lst, sd, na.rm = T ) } );
            median.mean.bias <- lapply( mean.bias.per.ppt, function( .lst ) { median( unlist( .lst ), na.rm = T ) } );
            min.mean.bias <- lapply( mean.bias.per.ppt, function( .lst ) { min( unlist( .lst ), na.rm = T ) } );
            max.mean.bias <- lapply( mean.bias.per.ppt, function( .lst ) { max( unlist( .lst ), na.rm = T ) } );
            median.sd.bias <- lapply( sd.bias.per.ppt, function( .lst ) { median( unlist( .lst ), na.rm = T ) } );
            min.sd.bias <- lapply( sd.bias.per.ppt, function( .lst ) { min( unlist( .lst ), na.rm = T ) } );
            max.sd.bias <- lapply( sd.bias.per.ppt, function( .lst ) { max( unlist( .lst ), na.rm = T ) } );
        
            num.partitions.per.ppt <- sapply( diffs.per.ppt.by.id[[1]], length );
        
            the.summary <- list( median.mean.bias = median.mean.bias, min.mean.bias = min.mean.bias, max.mean.bias = max.mean.bias,  median.sd.bias = median.sd.bias, min.sd.bias = min.sd.bias, max.sd.bias = max.sd.bias );
            
            results.per.ppt.by.id <- apply( results, 2, function ( .column ) {
              .rv <- 
                lapply( unique( rownames( results ) ), function( .ppt ) {
                    .values <- .column[ rownames( results ) == .ppt ];
                    .partition.ids <- partition.id[ rownames( results ) == .ppt ];
                    sapply( unique( .partition.ids ), function( the.partition.id ) {
                        return( .values[ .partition.ids == the.partition.id ] );
                    } );
                } );
              names( .rv ) <- unique( rownames( results ) );
              return( .rv );
            } );
        
            .the.partition.ids.by.sample <-
                sapply( 1:partition.bootstrap.samples, function( .sample.id ) {
                    sapply( num.partitions.per.ppt, sample, size = 1 );
                } );
            do.one.sample <- function ( .sample.id ) {
                 print( .sample.id );
        
                 .thissample.the.partition.ids <-
                     .the.partition.ids.by.sample[ , .sample.id ];
                .thissample.results.one.per.ppt <-
                    sapply( colnames( results ), function( est.name ) {
                        sapply( names( .thissample.the.partition.ids ), function( .ppt ) {
                            unname( results.per.ppt.by.id[[ est.name ]][[ .ppt ]][ .thissample.the.partition.ids[ .ppt ] ] )
                        } ) } );
                return( list( results.one.per.ppt = .thissample.results.one.per.ppt, evaluated.results = bound.and.evaluate.results.per.ppt( .thissample.results.one.per.ppt, days.since.infection, results.covars.one.per.ppt.with.extra.cols, the.artificial.bounds ) ) );
          } # do.one.sample (..)
            
          set.seed( partition.bootstrap.seed );
          bootstrap.results <- 
              mclapply( 1:partition.bootstrap.samples, do.one.sample, mc.cores = partition.bootstrap.num.cores );
        
            matrix.of.unbounded.results.rmses <- sapply( bootstrap.results, function( .results.for.bootstrap ) { .results.for.bootstrap[[ "evaluated.results" ]][[ "unbounded" ]]$rmse } );
            mode( matrix.of.unbounded.results.rmses ) <- "numeric";
            ## This uses the second position, which is the first of the unbounded ones, and for now the only one.  It's "5weeks" unless the.time is "1m" in which case it is "30weeks".
            matrix.of.bounded.results.rmses <- sapply( bootstrap.results, function( .results.for.bootstrap ) { .results.for.bootstrap[[ "evaluated.results" ]][[ 2 ]]$rmse } );
            mode( matrix.of.bounded.results.rmses ) <- "numeric";
            
            #hist( apply( matrix.of.unbounded.results.rmses, 1, diff ) )
            
          return( list( summary = the.summary, mean.bias.per.ppt.by.est = mean.bias.per.ppt, sd.bias.per.ppt.by.est = sd.bias.per.ppt, num.partitions.per.ppt = num.partitions.per.ppt, days.since.infection = days.since.infection, results.covars.one.per.ppt.with.extra.cols = results.covars.one.per.ppt.with.extra.cols, bounds = the.artificial.bounds, bootstrap.results = bootstrap.results, bootstrap.unbounded.rmse = matrix.of.unbounded.results.rmses, bootstrap.bounded.rmse = matrix.of.bounded.results.rmses ) );
        } # End if is.na( partition.size ) .. else ..
        
    } # get.timings.results.for.region.and.time (..)
    
    getTimingsResultsByRegionAndTime <- function ( partition.size = NA ) {
        if( !is.na( partition.size ) ) {
            regions <- "v3"; # Only v3 has partition results at this time.
        }
        timings.results.by.region.and.time <-
            lapply( regions, function( the.region ) {
                ## TODO: REMOVE
                cat( the.region, fill = T );
           timings.results.by.time <- 
               lapply( times, function( the.time ) {
                   ## TODO: REMOVE
                   cat( the.time, fill = T );
                   get.timings.results.for.region.and.time( the.region, the.time, partition.size );               
               } );
           names( timings.results.by.time ) <- times;

           timings.results.1m.6m <- lapply( setdiff( names( timings.results.by.time[[1]] ), "evaluated.results" ), function ( .varname ) {
               #print( .varname );
               if( .varname == "bounds" ) {
                   .rv <- 
                   lapply( names( timings.results.by.time[[ "1m" ]][[ .varname ]] ), function( .bounds.type ) {
                       #print( .bounds.type );
                     ..rv <- 
                         rbind(
                             timings.results.by.time[[ "1m" ]][[ .varname ]][[ .bounds.type ]],
                             timings.results.by.time[[ "6m" ]][[ .varname ]][[ .bounds.type ]]
                     );
                     rownames( ..rv ) <-
                         c( paste( rownames( timings.results.by.time[[ "1m" ]][[ .varname ]][[ .bounds.type ]] ), "1m", sep = "." ),
                           paste( rownames( timings.results.by.time[[ "6m" ]][[ .varname ]][[ .bounds.type ]] ), "6m", sep = "." ) );
                     return( ..rv );
                   } );
                   names( .rv ) <-
                       names( timings.results.by.time[[ "1m" ]][[ .varname ]] );
                   return( .rv );
               } else if( .varname == "days.since.infection" ) {
                   .rv <- c( 
                             timings.results.by.time[[ "1m" ]][[ .varname ]],
                             timings.results.by.time[[ "6m" ]][[ .varname ]]
                       );
                     names( .rv ) <-
                         c( paste( names( timings.results.by.time[[ "1m" ]][[ .varname ]] ), "1m", sep = "." ),
                           paste( names( timings.results.by.time[[ "6m" ]][[ .varname ]] ), "6m", sep = "." ) );
                   return( .rv );
               } else {
                     .rv <- 
                         rbind(
                             timings.results.by.time[[ "1m" ]][[ .varname ]],
                             timings.results.by.time[[ "6m" ]][[ .varname ]]
                     );
                     rownames( .rv ) <-
                         c( paste( rownames( timings.results.by.time[[ "1m" ]][[ .varname ]] ), "1m", sep = "." ),
                           paste( rownames( timings.results.by.time[[ "6m" ]][[ .varname ]] ), "6m", sep = "." ) );
                     return( .rv );
                 }
           } );
           names( timings.results.1m.6m ) <-
               setdiff( names( timings.results.by.time[[1]] ), "evaluated.results" );
           timings.results.1m.6m <-
               c( timings.results.1m.6m,
                 list( evaluated.results = bound.and.evaluate.results.per.ppt( timings.results.1m.6m[[ "results.one.per.ppt" ]], timings.results.1m.6m[[ "days.since.infection" ]], timings.results.1m.6m[[ "results.covars.one.per.ppt.with.extra.cols" ]], timings.results.1m.6m[[ "bounds" ]] ) ) );
           return( c( list( "1m.6m" = timings.results.1m.6m ), timings.results.by.time ) );
       } ); # End foreach the.region
      names( timings.results.by.region.and.time ) <- regions;
      return( timings.results.by.region.and.time );
    } # getTimingsResultsByRegionAndTime ( partition.size )
    
    if( force.recomputation || !file.exists( timings.results.by.region.and.time.Rda.filename ) ) {
        timings.results.by.region.and.time <- getTimingsResultsByRegionAndTime();
        save( timings.results.by.region.and.time, file = timings.results.by.region.and.time.Rda.filename );
    } else {
        # loads timings.results.by.region.and.time
        load( file = timings.results.by.region.and.time.Rda.filename );
    }
    
    # Make a table out of it. (one per study).
    results.table.by.region.and.time.and.bounds.type <-
        lapply( names( timings.results.by.region.and.time ), function( the.region ) {
            .rv <- 
                lapply( names( timings.results.by.region.and.time[[ the.region ]] ), function( the.time ) {
                  ..rv <- 
                    lapply( timings.results.by.region.and.time[[ the.region ]][[ the.time ]], function( results.by.bounds.type ) {
                      sapply( results.by.bounds.type, function( results.list ) { results.list } );
                    } );
                  names( ..rv ) <- names( timings.results.by.region.and.time[[ the.region ]][[ the.time ]] )
                  return( ..rv );
                } );
            names( .rv ) <- times;
            return( .rv );
        } );
    names( results.table.by.region.and.time.and.bounds.type ) <- regions;
    
    ## Write these out.
    bounds.types <- names( results.table.by.region.and.time.and.bounds.type[[1]][[1]] );
    .result.ignored <- sapply( regions, function ( the.region ) {
        ..result.ignored <- 
        sapply( times, function ( the.time ) {
          ...result.ignored <- 
            sapply( bounds.types, function ( the.bounds.type ) {
              out.file <- paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/", the.bounds.type, "_evaluateTimings.tab", sep = "" );
              ## TODO: REMOVE
              print( the.bounds.type );
              .tbl <-
                apply( results.table.by.region.and.time.and.bounds.type[[ the.region ]][[ the.time ]][[ "evaluated.results" ]][[ the.bounds.type ]], 1:2, function( .x ) { sprintf( "%0.2f", .x ) } );
              #print( .tbl );
              write.table( .tbl, quote = FALSE, file = out.file, sep = "\t" );
              return( NULL );
            } );
          return( NULL );
        } );
        return( NULL );
    } );
    
    ## For partition size == 10
    timings.results.by.region.and.time.p10 <- getTimingsResultsByRegionAndTime( partition.size = 10 );
    
    # Make a table out of it. (one per study).
    results.table.by.region.and.time.p10 <-
        lapply( names( timings.results.by.region.and.time.p10 ), function( the.region ) {
            .rv <- 
                lapply( names( timings.results.by.region.and.time.p10[[ the.region ]] ), function( the.time ) {
                    sapply( timings.results.by.region.and.time.p10[[ the.region ]][[ the.time ]], function( results.list ) { results.list } ) } );
            names( .rv ) <- times;
            return( .rv );
        } );
    names( results.table.by.region.and.time.p10 ) <- "v3";
    
    ## Write these out.
    .result.ignored <- sapply( "v3", function ( the.region ) {
        ..result.ignored <- 
        sapply( times, function ( the.time ) {
            out.file <- paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/", the.region, "/", the.time, "/evaluateTimings_p10.tab", sep = "" );
            .tbl <- apply( results.table.by.region.and.time.p10[[ the.region ]][[ the.time ]], 1:2, function( .x ) { sprintf( "%0.2f", .x ) } );
            write.table( .tbl, quote = FALSE, file = out.file, sep = "\t" );
            return( NULL );
        } );
        return( NULL );
    } );
    
    return( invisible( NULL ) );
} # evaluateTimings (..)

## Here is where the action is.
evaluateTimings();
