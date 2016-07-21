## First do all the stuff in README.postprocessing.txt.

library( "parallel" ); # for mclapply
library( "glmnet" ); # for cv.glmnet
library( "glmnetUtils" ); # for formula interface (cv.glmnet.formula): see https://github.com/Hong-Revo/glmnetUtils

source( "readIdentifyFounders_safetosource.R" );
source( "getDaysSinceInfection_safetosource.R" );
source( "getArtificialBounds_safetosource.R" );
source( "getResultsByRegionAndTime_safetosource.R" );
source( "writeResultsTables_safetosource.R" );
source( "summarizeCovariatesOnePerParticipant_safetosource.R" );

GOLD.STANDARD.DIR <- "/fh/fast/edlefsen_p/bakeoff/gold_standard";
RESULTS.DIR <- "/fh/fast/edlefsen_p/bakeoff_analysis_results/";

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
#' /fh/fast/edlefsen_p/bakeoff/analysis_sequences/raw_edited_20160216/bounds/nflg/1m/.
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
#' These files have names beginning with "artificialBounds_" and
#' ending with ".tab".
#'
#' @param use.bounds compute results for the COB approach, and also return evaluations of bounded versions of the other results.
#' @param use.infer compute results for the PREAST approach.
#' @param use.anchre compute results for the anchre approach.
#' @param use.glm.validate evaluate predicted values from leave-one-out cross-validation, using a model with one predictor, maybe with helpful.additional.cols or with the bounds.
#' @param use.lasso.validate evaluate predicted values from leave-one-out cross-validation, using a model with one predictor and a lasso-selected subset of other predictors, maybe with the bounds.
#' @param include.intercept if TRUE, include an intercept term, and for time-pooled analyses also include a term to shift the intercept for 6m.not.1m samples. Note that this analysis is affected by the low variance in the true dates in the training data, which for 1m samples has SD around 5 and for 6m samples has an SD around 10, so even the "none" results do pretty well, and it is difficult to improve the estimators beyond this.
#' @param helpful.additional.cols extra cols to be included in the glm and lasso: Note that interactions will be added only with members of the other set (helpful.additional.cols.with.interactions), not, with each other in this set.
#' @param helpful.additional.cols.with.interactions extra cols to be included in the glm both as-is and interacting with each other and with members of the helpful.additional.cols set.
#' @param results.dirname the subdirectory of "/fh/fast/edlefsen_p/bakeoff/analysis_sequences" and also of "/fh/fast/edlefsen_p/bakeoff_analysis_results"
#' @param force.recomputation if FALSE (default) and if there is a saved version called timings.results.by.region.and.time.Rda (under bakeoff_analysis_results/results.dirname), then that file will be loaded; otherwise the results will be recomputed and saved in that location.
#' @param partition.bootstrap.seed the random seed to use when bootstrapping samples by selecting one partition number per ptid, repeatedly; we do it this way because there are an unequal number of partitions, depending on sampling depth.
#' @param partition.bootstrap.samples the number of bootstrap replicates to conduct; the idea is to get an estimate of the variation in estimates and results (errors) across these samples.
#' @param partition.bootstrap.num.cores the number of cores to run the boostrap replicates on (defaults to all of the cores returned by parallel::detectCores()).
#' @return NULL
#' @export

evaluateTimings <- function (
                             use.bounds = TRUE,
                             use.infer = TRUE,
                             use.anchre = TRUE,
                             use.glm.validate = TRUE,
                             use.lasso.validate = TRUE,
                             include.intercept = FALSE,
                             helpful.additional.cols = c( "lPVL" ),
                             helpful.additional.cols.with.interactions = c( "v3_not_nflg", "X6m.not.1m" ),
                             results.dirname = "raw_edited_20160216",
                             force.recomputation = FALSE,
                             partition.bootstrap.seed = 98103,
                             partition.bootstrap.samples = 100,
                             partition.bootstrap.num.cores = detectCores(),
                             regions = c( "nflg", "v3" ),
                             times = c( "1m", "6m" )
                             #regions = c( "nflg", "v3", "rv217_v3" ),
                             #times = c( "1m", "6m", "1m6m" )
                            )
{
    if( include.intercept ) {
        results.by.region.and.time.Rda.filename <-
            paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/Timings.results.by.region.and.time.include.intercept.Rda", sep = "" );
    } else {
        results.by.region.and.time.Rda.filename <-
            paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/Timings.results.by.region.and.time.Rda", sep = "" );
    }
    
    rv217.gold.standard.infection.dates.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/rv217_gold_standard_timings.csv" );
    rv217.gold.standard.infection.dates <- as.Date( as.character( rv217.gold.standard.infection.dates.in[,2] ), "%m/%d/%y" );
    names( rv217.gold.standard.infection.dates ) <- as.character( rv217.gold.standard.infection.dates.in[,1] );
    
    caprisa002.gold.standard.infection.dates.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/caprisa_002_gold_standard_timings.csv" );
    caprisa002.gold.standard.infection.dates <- as.Date( as.character( caprisa002.gold.standard.infection.dates.in[,2] ), "%Y/%m/%d" );
    names( caprisa002.gold.standard.infection.dates ) <- as.character( caprisa002.gold.standard.infection.dates.in[,1] );

    ## TODO: REMOVE. Experimenting to evaluate the value of knowing whether there are multiple founders.
    # rv217.gold.standards.in <- read.csv( paste( GOLD.STANDARD.DIR, "/rv217/RV217_gold_standards.csv", sep = "" ) );
    # rv217.gold.is.multiple <- rv217.gold.standards.in[ , "gold.is.multiple" ];
    # names( rv217.gold.is.multiple ) <- rv217.gold.standards.in[ , "ptid" ];
    # 
    # caprisa002.gold.standards.in <- read.csv( paste( GOLD.STANDARD.DIR, "/caprisa_002/caprisa_002_gold_standards.csv", sep = "" ) );
    # caprisa002.gold.is.multiple <- caprisa002.gold.standards.in[ , "gold.is.multiple" ];
    # names( caprisa002.gold.is.multiple ) <- caprisa002.gold.standards.in[ , "ptid" ];

    rv217.pvl.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/rv217/rv217_viralloads.csv" );
    caprisa002.pvl.in <- read.csv( "/fh/fast/edlefsen_p/bakeoff/gold_standard/caprisa_002/caprisa_002_viralloads.csv" );
    
    rmse <- function( x, na.rm = FALSE ) {
        if( na.rm ) {
            x <- x[ !is.na( x ) ];
        }
        return( sqrt( mean( x ** 2 ) ) );
    }
    
    compute.results.per.person <- function ( results, weights ) {
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
    } # compute.results.per.person

    get.infer.results.columns <-
        function ( the.region, the.time, the.ptids, partition.size = NA ) {
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

        # For now only use the sampledwidth ones, so exclude the other artificialBounds results.
        infer.results.directories.sampledwidth <-
          grep( "sampledwidth", infer.results.directories, value = TRUE );
        infer.results.directories.unbounded <-
          grep( "artificialBounds", infer.results.directories, value = TRUE, invert = TRUE );
        infer.results.directories <- c( infer.results.directories.unbounded, infer.results.directories.sampledwidth );
        
        infer.results.files <- sapply( infer.results.directories, dir, "toi.csv", full.name = TRUE );
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

        ## Exclude old/outdated bounds
        infer.results.bounds.types <- 
          grep( "(one|six)month", infer.results.bounds.types, invert = TRUE );

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
        ## Add to results: "infer" results.
        if( ( the.region == "v3" ) || ( the.region == "rv217_v3" ) ) {
            the.region.dir <- "v3_edited_20160216";
        } else {
            the.region.dir <- "nflg_copy_20160222";
        }
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

    compute.diffs.by.stat <- function ( results.per.person, days.since.infection ) {
        diffs.by.stat <- 
            lapply( colnames( results.per.person ), function( .stat ) {
                .rv <- ( as.numeric( results.per.person[ , .stat ] ) - as.numeric( days.since.infection[ rownames( results.per.person ) ] ) );
                names( .rv ) <- rownames( results.per.person );
                return( .rv );
            } );
        names( diffs.by.stat ) <- colnames( results.per.person );
        return( diffs.by.stat );
    } # compute.diffs.by.stat ( results.per.person, days.since.infection )

    bound.and.evaluate.results.per.ppt <-
        function ( results.per.person, days.since.infection, results.covars.per.person.with.extra.cols, the.time, the.artificial.bounds = NA, ppt.suffix.pattern = "\\..+", return.lasso.coefs = TRUE ) {

       ## Special: the ppt names might have suffices in results.per.person; if so, strip off the suffix for purposes of matching ppts to the covars, etc.
       ppt.names <- rownames( results.per.person );
       ppt.suffices <- NULL;
       if( !is.na( ppt.suffix.pattern ) && ( length( grep( ppt.suffix.pattern, ppt.names ) ) > 0 ) ) {
         # if one has it, they should all have it.
         stopifnot( length( grep( ppt.suffix.pattern, ppt.names ) ) == length( ppt.names ) );
           .ppt.names <- ppt.names;
           ppt.names <- gsub( ppt.suffix.pattern, "", .ppt.names );
           names( ppt.names ) <- .ppt.names;
           ppt.suffices <- gsub( paste( "^.+?(", ppt.suffix.pattern, ")$", sep = "" ), "\\1", .ppt.names, perl = TRUE );
           names( ppt.suffices ) <- .ppt.names;
       }
       all.ptids <- unique( ppt.names );
        
        ## Also include all of the date estimates, which is
        ## everything in "results.per.person" so
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
        days.est.cols <- colnames( results.per.person );
        days.est.cols <- grep( "deterministic", days.est.cols, invert = TRUE, value = TRUE );
        days.est.cols <- grep( "PFitter|Star[Pp]hy", days.est.cols, invert = TRUE, perl = TRUE, value = TRUE );
        # Also exclude anything time-dependent with times we don't use anymore.
        days.est.cols <- grep( "(one|six)month", days.est.cols, invert = TRUE, perl = TRUE, value = TRUE );
       if( the.time == "1m.6m" ) {
         # Then keep the 1mmtn003.6mhvtn502 results, but not each separately.
        days.est.cols <- grep( "\\.mtn|\\.hvtn", days.est.cols, invert = TRUE, perl = TRUE, value = TRUE );
       }

### TODO: PUT THAT BACK.  TESTING.
       # Do not allow use of information about the time of sampling unless include.intercept is TRUE
#        if( !include.intercept && ( "6m.not.1m" %in% colnames( results.covars.per.person.with.extra.cols ) ) ) {
#          .forbidden.column.i <- which( "6m.not.1m" == colnames( results.covars.per.person.with.extra.cols ) );
#          results.covars.per.person.with.extra.cols <-
#            results.covars.per.person.with.extra.cols[ , -.forbidden.column.i, drop = FALSE ];
#        }
       
      if( use.glm.validate || use.lasso.validate ) {
        results.covars.per.person.with.extra.cols <-
            cbind( results.per.person[ , days.est.cols, drop = FALSE ], results.covars.per.person.with.extra.cols );
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
        #.donotkeep.cols <- c( "multifounder.Synonymous.PFitter.is.poisson", "Synonymous.PFitter.is.poisson", "multifounder.PFitter.is.poisson", "PFitter.is.poisson", "multifounder.Synonymous.PFitter.is.starlike", "Synonymous.PFitter.is.starlike", "multifounder.PFitter.is.starlike", "PFitter.is.starlike", "inf.sites", "mean.entropy", "PFitter.mean.hd", "inf.to.priv.ratio", "StarPhy.founders", "multifounder.DS.Starphy.R", "PFitter.chi.sq.stat", "Synonymous.DS.StarPhy.R", "StarPhy.is.one.founder", "multifounder.DS.Starphy.fits", "multifounder.Synonymous.DS.StarPhy.fits", "Synonymous.DS.StarPhy.fits", "DS.Starphy.fits", "DS.Starphy.is.starlike", "sd.entropy", "PFitter.max.hd", "InSites.is.one.founder", "InSites.founders" );
        ## TODO: REMOVE . TESTING
        #.donotkeep.cols <- c( .donotkeep.cols, "gold.is.multiple" );
        .donotkeep.cols <- c( grep( "DS\\.Star[pP]y\\.(fits|founders|is.starlike)", .keep.cols, value = TRUE ) );
        .donotkeep.cols <- c( .donotkeep.cols, "StarPhy.founders", "InSites.founders", "multifounder.DS.Starphy.fits" );
        .keep.cols <- setdiff( .keep.cols, .donotkeep.cols );

        ## Keep only the mut.rate.coef cols and priv.sites and multifounder.Synonymous.PFitter.is.poisson, and Infer and anchre cols.
        Infer.cols <- grep( "Infer", .keep.cols, value = TRUE );
        ## Exclude Infer cols for obsolete bounds.
        Infer.cols <- grep( "(one|six)month", Infer.cols, value = TRUE, invert = TRUE );
        ## Actually also exclude the ones that don't match the time
        if( the.time == "1m.6m" ) {
          Infer.cols <- grep( "\\.mtn|\\.hvtn", Infer.cols, value = TRUE, invert = TRUE );
        }
        anchre.cols <- grep( "anchre", .keep.cols, value = TRUE );
        mut.rate.coef.cols <- grep( "mut\\.rate\\.coef", .keep.cols, value = TRUE );
        COB.cols <- grep( "^COB", .keep.cols, value = TRUE );
        ## Exclude COB cols for obsolete bounds.
        COB.cols <- grep( "(one|six)month", COB.cols, value = TRUE, invert = TRUE );
        ## Actually also exclude the ones that don't match the time
        if( the.time == "1m.6m" ) {
          COB.cols <- grep( "\\.mtn|\\.hvtn", COB.cols, value = TRUE, invert = TRUE );
        }
        estimate.cols <- c( COB.cols, mut.rate.coef.cols, Infer.cols, anchre.cols );
        all.additional.cols <- setdiff( .keep.cols, estimate.cols );
        
        if( use.lasso.validate ) {
          keep.cols <- unique( c( helpful.additional.cols, helpful.additional.cols.with.interactions, all.additional.cols, estimate.cols ) );
        } else {
          keep.cols <- unique( c( helpful.additional.cols, helpful.additional.cols.with.interactions, estimate.cols ) );
      }
        
        # Don't evaluate estimators that have no variation at
        # all. Note that technically we could/should do this after holding
        # out each person, in case a value has no variation among the
        # subset excluding that person.  But we do it here only.
        estimators.to.exclude <-
          apply( results.covars.per.person.with.extra.cols[ , estimate.cols, drop = FALSE ], 2, function ( .col ) {
            return( ( sum( !is.na( .col ) ) <= 1 ) || ( var( .col, na.rm = TRUE ) == 0 ) );
          } );
        estimate.cols <- setdiff( estimate.cols, names( which( estimators.to.exclude ) ) );
           
        # Also always evaluate no-estimate: "none".
        estimate.cols <- c( "none", estimate.cols );

        # To make things consistent we actually change this _to_ X6m.not.1m.
        colnames( results.covars.per.person.with.extra.cols )[ colnames( results.covars.per.person.with.extra.cols ) == "6m.not.1m" ] <- "X6m.not.1m";
        # Some cols from "helpful.additional.cols" or "helpful.additional.cols.with.interactions" might be missing
        keep.cols <-
            intersect( keep.cols, colnames( results.covars.per.person.with.extra.cols ) );
        results.covars.per.person <-
            results.covars.per.person.with.extra.cols[ , keep.cols, drop = FALSE ];
        results.covars.per.person.df <-
            data.frame( results.covars.per.person );
        ## DO NOT Undo conversion of the colnames (X is added before "6m.not.1m").  We want it to be called "X6m.not.1m" so it can work in the regression formulas.
        #colnames( results.covars.per.person.df ) <- colnames( results.covars.per.person.with.extra.cols );

        ## TODO: REMOVE. TESTING inclusion of interactions with gold.is.multiple.
        # ..df <- sapply( 1:nrow( results.covars.per.person.df ), function( .ppt.i ) { ( results.covars.per.person.df[ .ppt.i, "gold.is.multiple" ] * results.covars.per.person.df[ .ppt.i, setdiff( all.additional.cols, "gold.is.multiple" ) ] ) } );
        # gold.is.multiple.interactions.df <- t( ..df );
        # colnames( gold.is.multiple.interactions.df ) <- paste( "gold.is.multiple", colnames( gold.is.multiple.interactions.df ), sep = ".." );
        # results.covars.per.person.df <- cbind( results.covars.per.person.df, gold.is.multiple.interactions.df );
        # keep.cols <- c( keep.cols, colnames( gold.is.multiple.interactions.df ) );
        # all.additional.cols <- c( all.additional.cols, colnames( gold.is.multiple.interactions.df ) );
        
        regression.df <- cbind( data.frame( days.since.infection = days.since.infection[ rownames( results.covars.per.person.df ) ] ), results.covars.per.person.df, lapply( the.artificial.bounds[ grep( "(one|six)month", names( the.artificial.bounds ), invert = TRUE, value = TRUE ) ], function( .mat ) { .mat[ rownames( results.covars.per.person.df ), , drop = FALSE ] } ) );
        
        ## Ok build a regression model, including only the helpful.additional.cols and helpful.additional.cols.with.interactions (and interactions among them), and also the lower and upper bounds associated with either 5 weeks or 30 weeks, depending on the.time (if there's a 1m sample, uses "mtn003").
        if( ( the.time == "6m" ) || ( the.time == "1m6m" ) ) {
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
          #glm.fit.statistics <- matrix( 0, nrow = length( all.ptids ), nrow( results.covars.per.person.df ), ncol = 0 );
          
            ## Note that there might be multiple rows per ppt in the regression.df and in this prediction output matrix; the values will be filled in using leave-one-ptid-out xv, ie in each iteration there might be multiple rows filled in, since multiple rows correspond to one held-out ptid.
          glm.validation.results.per.person <- matrix( NA, nrow = nrow( results.covars.per.person.df ), ncol = length( estimate.cols ) );
          glm.withbounds.validation.results.per.person <- matrix( NA, nrow = nrow( results.covars.per.person.df ), ncol = length( estimate.cols ) );
        }
        if( use.lasso.validate ) {
            ## These are again per-sample, see above Note.
            lasso.validation.results.per.person <- matrix( NA, nrow = nrow( results.covars.per.person.df ), ncol = length( estimate.cols ) );
            lasso.withbounds.validation.results.per.person <- matrix( NA, nrow = nrow( results.covars.per.person.df ), ncol = length( estimate.cols ) );
            # lasso.nointercept.validation.results.per.person <- matrix( NA, nrow = nrow( results.covars.per.person.df ), ncol = length( estimate.cols ) );
            # lasso.withbounds.nointercept.validation.results.per.person <- matrix( NA, nrow = nrow( results.covars.per.person.df ), ncol = length( estimate.cols ) );
            if( return.lasso.coefs ) {
                ## This is really a 3D array, but I'm just lazily representing it directly this way.  Note this is by removed ppt, not by regression row (there might be multiple rows per ppt in the regression.df and the prediction output matrices).
                # and these are per-ptid!
                lasso.validation.results.per.person.coefs <-
                    as.list( rep( NA, length( all.ptids ) ) );
                names( lasso.validation.results.per.person.coefs ) <-
                    all.ptids;
                lasso.withbounds.validation.results.per.person.coefs <-
                    as.list( rep( NA, length( all.ptids ) ) );
                names( lasso.withbounds.validation.results.per.person.coefs ) <-
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
                .lasso.validation.results.per.person.coefs.row <-
                    as.list( rep( NA, length( estimate.cols ) ) );
                names( .lasso.validation.results.per.person.coefs.row ) <-
                    estimate.cols;
                .lasso.withbounds.validation.results.per.person.coefs.row <-
                    as.list( rep( NA, length( estimate.cols ) ) );
                names( .lasso.withbounds.validation.results.per.person.coefs.row ) <-
                    estimate.cols;
            }
            regression.df.without.ptid.i <-
                regression.df[ the.rows.excluding.ptid, , drop = FALSE ];
            .out <- regression.df.without.ptid.i[[ "days.since.infection" ]];
            for( .col.i in 1:length( estimate.cols ) ) {
                .estimate.colname <- estimate.cols[ .col.i ];
                ## TODO: REMOVE
                #print( .estimate.colname );
                if( use.glm.validate ) {
                  # covariates for glm
                  .covariates.glm <-
                    c( helpful.additional.cols, helpful.additional.cols.with.interactions );
                  # covariates for glm.withbounds
                  .covariates.glm.withbounds <-
                    c( helpful.additional.cols, helpful.additional.cols.with.interactions, .upper.bound.colname );

                  .covars.to.exclude <- apply( regression.df.without.ptid.i, 2, function ( .col ) {
                      return( ( sum( !is.na( .col ) ) <= 1 ) || ( var( as.numeric( .col ), na.rm = TRUE ) == 0 ) );
                  } );
                  # Also exclude covars that are too highly correlated.
                  .retained.covars <-
                      setdiff( colnames( regression.df.without.ptid.i ), names( which( .covars.to.exclude ) ) );
                  if( ( .estimate.colname == "none" ) || ( .estimate.colname %in% .retained.covars ) ) {
                     ## ERE I AM.
                     .interactors <-
                         intersect( .retained.covars, helpful.additional.cols.with.interactions );
                     # Add interactions among the interactors.
                     if( length( .interactors ) > 1 ) {
                         ## NOTE that for now it is up to the caller to ensure that the df is not too high in this case.
                         .interactions.among.interactors <- c();
                         for( .interactor.i in 1:( length( .interactors ) - 1 ) ) {
                             .interactor <- .interactors[ .interactor.i ];
                             .remaining.interactors <-
                                 .interactors[ ( .interactor.i + 1 ):length( .interactors ) ];
                             .interactions.among.interactors <-
                                 c( .interactions.among.interactors,
                                   paste( .interactor, .remaining.interactors, collapse = "+", sep = ":" ) );
                         } # End foreach .interactor.i
                         .interactors <- c( .interactors, .interactions.among.interactors );
                     } # End if there's more than one "interactor"
                     
                     if( .estimate.colname == "none" ) {
                        .cv.glm <- intersect( .retained.covars, .covariates.glm );
                        if( length( .cv.glm ) == 0 ) {
                          .cv.glm <- 1;
                      }
                      ## SEE NOTE BELOW ABOUT USE OF AN INTERCEPT. WHEN we're using "none", then we always do use an intercept (only for the non-bounded result).
                        if( include.intercept ) {
                            .formula <- paste( "days.since.infection ~", paste( .cv.glm, collapse = "+" ) );
                            .formula.withbounds <- paste( "days.since.infection ~", paste( intersect( .retained.covars, .covariates.glm.withbounds ), collapse = "+" ) );
                        } else {
                            .cv.glm.nointercept <- intersect( .retained.covars, .covariates.glm );
                            if( length( .cv.glm.nointercept ) == 0 ) {
                                .formula.nointercept <- "days.since.infection ~ 1";  # _only_ intercept!
                            } else {
                                .formula.nointercept <- paste( "days.since.infection ~ 0 + ", paste( .cv.glm.nointercept, collapse = "+" ) );
                            }
                            .covariates.glm.withbounds.nointercept <- .covariates.glm.withbounds;
                            .formula.withbounds.nointercept <- paste( "days.since.infection ~ 0 + ", paste( intersect( .retained.covars, .covariates.glm.withbounds.nointercept ), collapse = "+" ) );
                            .formula <- .formula.nointercept;
                            .formula.withbounds <- .formula.withbounds.nointercept;
                        }
                    } else {
                        ## NOTE WE MUST BE CAREFUL ABOUT USING AN INTERCEPT, BECAUSE THEN WE JUST CHEAT BY PREDICTING EVERYTHING IS CENTERED AT the tight true bounds on the days.since.infection in our training data (when looking at one time at a time, and now with 6m.not.1m this should be true for both times as well.).
                      if( include.intercept ) {
                          .formula <- paste( "days.since.infection ~", paste( intersect( .retained.covars, c( .covariates.glm, .estimate.colname ) ), collapse = "+" ) );
                          .formula.withbounds <- paste( "days.since.infection ~", paste( intersect( .retained.covars, c( .covariates.glm.withbounds, .estimate.colname ) ), collapse = "+" ) );
                      } else {
                          .covariates.glm.nointercept <- .covariates.glm;
                          .formula.nointercept <- paste( "days.since.infection ~ 0 + ", paste( intersect( .retained.covars, c( .covariates.glm.nointercept, .estimate.colname ) ), collapse = "+" ) );
                          .formula <- .formula.nointercept;
                          .covariates.glm.withbounds.nointercept <- .covariates.glm.withbounds;
                          .formula.withbounds.nointercept <-
                              paste( "days.since.infection ~ 0 + ", paste( intersect( .retained.covars, c( .covariates.glm.withbounds.nointercept, .estimate.colname ) ), collapse = "+" ) );
                          .formula.withbounds <- .formula.withbounds.nointercept;
                      }
                    }
                     if( length( .interactors ) > 0 ) {
                         .interactees <-
                             setdiff( intersect( .retained.covars, .covariates.glm ), .interactors );
                         if( length( .interactees ) > 0 ) {
                             ## NOTE that for now it is up to the caller to ensure that the df is not too high in this case.
                             .interactions.formula.part.components <- c();
                             for( .interactor.i in 1:length( .interactors ) ) {
                                 .interactor <- .interactors[ .interactor.i ];
                                 .interactions.formula.part.components <-
                                     c( .interactions.formula.part.components,
                                       paste( .interactor, .interactees, collapse = "+", sep = ":" ) );
                             } # End foreach .interactor.i
                                                                      
                             .interactions.formula.part <-
                                 paste( .interactions.formula.part.components, collapse =" + " );
                             .formula <-
                                 paste( .formula, .interactions.formula.part, sep = " + " );
                         } # End if there's any "interactees"
                         ## Also add interactions with the estimate colname and everything included so far.
                         if( .estimate.colname != "none" ) {
                             .everything.included.so.far.in.formula <- 
                                 strsplit( .formula, split = "\\s*\\+\\s*" )[[1]];
                             ## TODO: Move this part down later
                             .everything.included.so.far.in.formula.withbounds <- 
                                 strsplit( .formula.withbounds, split = "\\s*\\+\\s*" )[[1]];

                             # Add interactions.
                             .new.interactions.with.estimator.in.formula <-
                                 sapply( .everything.included.so.far.in.formula[ -1 ], function ( .existing.part ) {
                                     paste( .existing.part, .estimate.colname, sep = ":" )
                                 } );
                             ## TODO: Move this part down later
                             .new.interactions.with.estimator.in.formula.withbounds <-
                                 sapply( .everything.included.so.far.in.formula.withbounds[ -1 ], function ( .existing.part ) {
                                     paste( .existing.part, .estimate.colname, sep = ":" )
                                 } );

                             .formula <- paste( c( .everything.included.so.far.in.formula, .new.interactions.with.estimator.in.formula ), collapse = " + " );
                         } # End if .estimate.colname != "none"

                         ## withbounds
                         .interactees.withbounds <-
                             setdiff( intersect( .retained.covars, .covariates.glm.withbounds ), .interactors );
                         if( length( .interactees.withbounds ) > 0 ) {
                             ## NOTE that for now it is up to the caller to ensure that the df is not too high in this case.
                             .interactions.withbounds.formula.part.components <- c();
                             for( .interactor.i in 1:length( .interactors ) ) {
                                 .interactor <- .interactors[ .interactor.i ];
                                 .interactions.withbounds.formula.part.components <-
                                     c( .interactions.withbounds.formula.part.components,
                                       paste( .interactor, .interactees.withbounds, collapse = "+", sep = ":" ) );
                             } # End foreach .interactor.i
                                                                      
                             .interactions.withbounds.formula.part <-
                                 paste( .interactions.withbounds.formula.part.components, collapse =" + " );
                             .formula.withbounds <-
                                 paste( .formula.withbounds, .interactions.withbounds.formula.part, sep = " + " );
                         } # End if there's any "interactees.withbounds"
                         ## Also add interactions with the estimate colname and everything included so far.
                         if( .estimate.colname != "none" ) {
                             .everything.included.so.far.in.formula.withbounds <- 
                                 strsplit( .formula.withbounds, split = "\\s*\\+\\s*" )[[1]];

                             # Add interactions.
                             .new.interactions.with.estimator.in.formula.withbounds <-
                                 sapply( .everything.included.so.far.in.formula.withbounds[ -1 ], function ( .existing.part ) {
                                     paste( .existing.part, .estimate.colname, sep = ":" )
                                 } );

                             .formula.withbounds <- paste( c( .everything.included.so.far.in.formula.withbounds, .new.interactions.with.estimator.in.formula.withbounds ), collapse = " + " );
                         } # End if .estimate.colname != "none"
                     } # End if there are any interactors.
                     
                     # To this point these are still strings:
                     #.formula <- as.formula( .formula );
                     #.formula.withbounds <- as.formula( .formula.withbounds );
                     
                    .df <- regression.df.without.ptid.i[ , .retained.covars, drop = FALSE ];
                    .glm.result <- lm( .formula, data = regression.df.without.ptid.i );
                      
                    for( .row.i in the.rows.for.ptid ) {
                      .pred.value.glm <- predict( .glm.result, regression.df[ .row.i, , drop = FALSE ] );
                      glm.validation.results.per.person[ .row.i, .col.i ] <- 
                          .pred.value.glm;
                    }
                      
                    # glm.withbounds:
                    .glm.withbounds.result <- lm( .formula.withbounds, data = regression.df.without.ptid.i );
                    for( .row.i in the.rows.for.ptid ) {
                      .pred.value.glm.withbounds <- predict( .glm.withbounds.result, regression.df[ .row.i, , drop = FALSE ] );
                      glm.withbounds.validation.results.per.person[ .row.i, .col.i ] <- 
                          .pred.value.glm.withbounds;
                    }
                    
                    ## # glm.nointercept:
                    ## .glm.nointercept <- lm( .formula.nointercept, data = regression.df.without.ptid.i );
                    ## glm.nointercept.validation.results.per.person.adjR2[ .ptid.i, .col.i ] <- 
                    ##     ( summary( .glm.nointercept ) )$adj.r.squared;
                    ## glm.nointercept.validation.results.per.person.pcoef[ .ptid.i, .col.i ] <- 
                    ##     coef( summary( .glm.nointercept ) )[ 1, "Pr(>|t|)" ];
                    ## for( .row.i in the.rows.for.ptid ) {
                    ##   .pred.value.glm.nointercept <- predict( .glm.nointercept, regression.df[ .row.i, , drop = FALSE ] );
                    ##   glm.nointercept.validation.results.per.person[ .row.i, .col.i ] <- 
                    ##       .pred.value.glm.nointercept;
                    ## }
                    ## 
                    ## # glm.withbounds.nointercept:
                    ## .glm.withbounds.nointercept <- lm( .formula.withbounds.nointercept, data = regression.df.without.ptid.i );
                    ## glm.withbounds.nointercept.validation.results.per.person.adjR2[ .ptid.i, .col.i ] <- 
                    ##     ( summary( .glm.withbounds.nointercept ) )$adj.r.squared;
                    ## glm.withbounds.nointercept.validation.results.per.person.pcoef[ .ptid.i, .col.i ] <- 
                    ##     coef( summary( .glm.withbounds.nointercept ) )[ 1, "Pr(>|t|)" ];
                    ## for( .row.i in the.rows.for.ptid ) {
                    ##   .pred.value.glm.withbounds.nointercept <- predict( .glm.withbounds.nointercept, regression.df[ .row.i, , drop = FALSE ] );
                    ##   glm.withbounds.nointercept.validation.results.per.person[ .row.i, .col.i ] <- 
                    ##       .pred.value.glm.withbounds.nointercept;
                    ## }
                  } # If the estimate can be used
    
                } # End if use.glm.validate
    
                if( use.lasso.validate ) {
                  # covariates for lasso
                  .covariates.lasso <- 
                      intersect( colnames( regression.df.without.ptid.i ), unique( c( helpful.additional.cols, all.additional.cols ) ) );
                  
                  # covariates for lasso.withbounds
                  .covariates.lasso.withbounds <-
                      intersect( colnames( regression.df.without.ptid.i ), unique( c( helpful.additional.cols, .lower.bound.colname, .upper.bound.colname, all.additional.cols ) ) );

                  # lasso:
                  if( .estimate.colname == "none" ) {
                    .lasso.mat <-
                        as.matrix( regression.df.without.ptid.i[ , .covariates.lasso, drop = FALSE ] );
                  } else {
                    .lasso.mat <-
                        as.matrix( regression.df.without.ptid.i[ , c( .covariates.lasso, .estimate.colname ), drop = FALSE ] );
                  }
                  mode( .lasso.mat ) <- "numeric";
                  
                  # exclude any rows with any NAs.
                  .retained.rows <-
                      which( apply( .lasso.mat, 1, function( .row ) { !any( is.na( .row ) ) } ) );
                  .lasso.mat <-
                      .lasso.mat[ .retained.rows, , drop = FALSE ];
                  .out <- regression.df.without.ptid.i[[ "days.since.infection" ]][ .retained.rows ];
                  .covars.to.exclude <- apply( .lasso.mat, 2, function ( .col ) {
                        return( ( sum( !is.na( .col ) ) <= 1 ) || ( var( as.numeric( .col ), na.rm = TRUE ) == 0 ) );
                   } );
                    .covars.to.exclude <- names( which( .covars.to.exclude ) );
                  
                  # At least one DF is needed for the variance estimate (aka the error term), and one for leave-one-out xvalidation.
                  MINIMUM.DF <- 2; # how much more should nrow( .lasso.mat ) be than ncol( .lasso.mat ) at minimum?
                  if( include.intercept ) {
                      MINIMUM.DF <- MINIMUM.DF + 1;
                  }

                  MINIMUM.CORRELATION.WITH.OUTCOME <- 0.1;
                  cors.with.the.outcome <-
                    sapply( setdiff( colnames( .lasso.mat ), c( .covars.to.exclude, .estimate.colname ) ), function( .covar.colname ) {
                          .cor <- 
                              cor( .out, .lasso.mat[ , .covar.colname ], use = "pairwise" );
                            return( .cor );
                        } );
                  sorted.cors.with.the.outcome <-
                    cors.with.the.outcome[ order( abs( cors.with.the.outcome ) ) ];
                  # Sort the columns of .lasso.mat by their correlation with the outcome. This is to ensure that more-relevant columns get selected when removing columns due to pairwise correlation among them.  We want to keep the best one.
                  if( .estimate.colname == "none" ) {
                    .lasso.mat <- .lasso.mat[ , rev( names( sorted.cors.with.the.outcome ) ), drop = FALSE ];
                  } else {
                    .lasso.mat <- .lasso.mat[ , c( .estimate.colname, rev( names( sorted.cors.with.the.outcome ) ) ), drop = FALSE ];
                  }
                  # Exclude covars that are not sufficiently correlated with the outcome.
                      .covars.to.exclude <- c( .covars.to.exclude,
                          names( which( sapply( setdiff( colnames( .lasso.mat ), c( .covars.to.exclude, .estimate.colname ) ), function( .covar.colname ) {
                            #print( .covar.colname );
                          .cor <- 
                              cor( .out, .lasso.mat[ , .covar.colname ], use = "pairwise" );
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
                      .covars.to.exclude <- c( .covars.to.exclude,
                          names( which( sapply( setdiff( colnames( .lasso.mat ), c( .covars.to.exclude, .estimate.colname ) ), function( .covar.colname ) {
                            #print( .covar.colname );
                          .cor <- 
                              cor( .lasso.mat[ , .estimate.colname ], .lasso.mat[ , .covar.colname ], use = "pairwise" );
                            #print( .cor );
                          if( is.na( .cor ) || ( abs( .cor ) >= COR.THRESHOLD ) ) {
                            TRUE
                          } else {
                            FALSE
                          }
                        } ) ) ) );
                    }
                  
                  # Exclude covars that are too highly correlated with each other.
                  .covars.to.consider <-
                    setdiff( colnames( .lasso.mat ), c( .covars.to.exclude, .estimate.colname, .lower.bound.colname, .upper.bound.colname ) );
                  ## Process these in reverse order to ensure that we prioritize keeping those towards the top.
                  .covars.to.consider <- rev( .covars.to.consider );
                  .new.covars.to.exclude <- rep( FALSE, length( .covars.to.consider ) );
                  names( .new.covars.to.exclude ) <- .covars.to.consider;
                  for( .c.i in 1:length( .covars.to.consider ) ) {
                    .covar.colname <- .covars.to.consider[ .c.i ];
                    #print( .covar.colname );
                    # Only consider those not already excluded.
                    .cor <- 
                      cor( .lasso.mat[ , .covar.colname ], .lasso.mat[ , names( which( !.new.covars.to.exclude ) ) ], use = "pairwise" );
                    #print( .cor );
                        if( ( length( .cor ) > 0 ) && any( .cor[ !is.na( .cor ) ] < 1 & .cor[ !is.na( .cor ) ] >= COR.THRESHOLD ) ) {
                          .new.covars.to.exclude[ .c.i ] <- TRUE;
                        } else {
                          .new.covars.to.exclude[ .c.i ] <- FALSE;
                        }
                  } # End foreach of the .covars.to.consider

                  .covars.to.exclude <- c( .covars.to.exclude, names( which( .new.covars.to.exclude ) ) );
                  .retained.covars <- setdiff( colnames( .lasso.mat ), .covars.to.exclude );

                  ## MARK
                  .covariates.lasso <- intersect( .retained.covars, .covariates.lasso );

                     if( .estimate.colname == "none" ) {
                        .cv.lasso <- intersect( .retained.covars, .covariates.lasso );
                        if( length( .cv.lasso ) == 0 ) {
                          .cv.lasso <- 1;
                      }
                      ## SEE NOTE BELOW ABOUT USE OF AN INTERCEPT. WHEN we're using "none", then we always do use an intercept (only for the non-bounded result).
                        if( include.intercept ) {
                            .formula <- paste( "days.since.infection ~", paste( .cv.lasso, collapse = "+" ) );
                        } else {
                            .cv.lasso.nointercept <- intersect( .retained.covars, .covariates.lasso );
                            if( length( .cv.lasso.nointercept ) == 0 ) {
                                .formula.nointercept <- "days.since.infection ~ 1";  # _only_ intercept!
                            } else {
                                .formula.nointercept <- paste( "days.since.infection ~ 0 + ", paste( .cv.lasso.nointercept, collapse = "+" ) );
                            }
                            .formula <- .formula.nointercept;
                        }
                    } else {
                        ## NOTE WE MUST BE CAREFUL ABOUT USING AN INTERCEPT, BECAUSE THEN WE JUST CHEAT BY PREDICTING EVERYTHING IS CENTERED AT the tight true bounds on the days.since.infection in our training data (when looking at one time at a time, and now with 6m.not.1m this should be true for both times as well.).
                      if( include.intercept ) {
                          .formula <- paste( "days.since.infection ~", paste( intersect( .retained.covars, c( .covariates.lasso, .estimate.colname ) ), collapse = "+" ) );
                      } else {
                          .covariates.lasso.nointercept <- .covariates.lasso;
                          .formula.nointercept <- paste( "days.since.infection ~ 0 + ", paste( intersect( .retained.covars, c( .covariates.lasso.nointercept, .estimate.colname ) ), collapse = "+" ) );
                          .formula <- .formula.nointercept;
                      }
                    }
                  
                     .interactors <-
                         intersect( .retained.covars, helpful.additional.cols.with.interactions );
                     # Add interactions among the interactors.
                     if( length( .interactors ) > 1 ) {
                         ## NOTE that for now it is up to the caller to ensure that the df is not too high in this case.
                         .interactions.among.interactors <- c();
                         for( .interactor.i in 1:( length( .interactors ) - 1 ) ) {
                             .interactor <- .interactors[ .interactor.i ];
                             .remaining.interactors <-
                                 .interactors[ ( .interactor.i + 1 ):length( .interactors ) ];
                             .interactions.among.interactors <-
                                 c( .interactions.among.interactors,
                                   paste( .interactor, .remaining.interactors, collapse = "+", sep = ":" ) );
                         } # End foreach .interactor.i
                         .interactors <- c( .interactors, .interactions.among.interactors );
                     } # End if there's more than one "interactor"
                     if( length( .interactors ) > 0 ) {
                         .interactees <-
                             setdiff( .retained.covars, .interactors );
                         if( length( .interactees ) > 0 ) {
                             ## NOTE that for now it is up to the caller to ensure that the df is not too high in this case.
                             .interactions.formula.part.components <- c();
                             for( .interactor.i in 1:length( .interactors ) ) {
                                 .interactor <- .interactors[ .interactor.i ];
                                 .interactions.formula.part.components <-
                                     c( .interactions.formula.part.components,
                                       paste( .interactor, .interactees, collapse = "+", sep = ":" ) );
                             } # End foreach .interactor.i
                                                                      
                             .interactions.formula.part <-
                                 paste( .interactions.formula.part.components, collapse =" + " );
                             .formula <-
                                 paste( .formula, .interactions.formula.part, sep = " + " );
                         } # End if there's any "interactees"
                         ## Also add interactions with the estimate colname and everything included so far.
                         if( .estimate.colname != "none" ) {
                             .everything.included.so.far.in.formula <- 
                                 strsplit( .formula, split = "\\s*\\+\\s*" )[[1]];

                             # Add interactions.
                             .new.interactions.with.estimator.in.formula <-
                                 sapply( .everything.included.so.far.in.formula[ -1 ], function ( .existing.part ) {
                                     paste( .existing.part, .estimate.colname, sep = ":" )
                                 } );

                             .formula <- paste( c( .everything.included.so.far.in.formula, .new.interactions.with.estimator.in.formula ), collapse = " + " );
                         } # End if .estimate.colname != "none"


                         .mf <- stats::model.frame( as.formula( .formula ), data = regression.df.without.ptid.i );
                         .lasso.mat <- model.matrix(as.formula( .formula ), .mf);
                         
                         .covars.to.exclude <- apply( .lasso.mat, 2, function ( .col ) {
                             return( ( sum( !is.na( .col ) ) <= 1 ) || ( var( .col, na.rm = TRUE  ) == 0 ) );
                         } );
                         .covars.to.exclude <- names( which( .covars.to.exclude ) );
                         
                         cors.with.the.outcome <-
                             sapply( setdiff( colnames( .lasso.mat ), c( .covars.to.exclude, .estimate.colname ) ), function( .covar.colname ) {
                          .cor <- 
                              cor( .out, .lasso.mat[ , .covar.colname ], use = "pairwise" );
                            return( .cor );
                         } );
                         sorted.cors.with.the.outcome <-
                             cors.with.the.outcome[ order( abs( cors.with.the.outcome ) ) ];
                         # Sort the columns of .lasso.mat by their correlation with the outcome. This is to ensure that more-relevant columns get selected when removing columns due to pairwise correlation among them.  We want to keep the best one.
                         if( .estimate.colname == "none" ) {
                             .lasso.mat <- .lasso.mat[ , rev( names( sorted.cors.with.the.outcome ) ), drop = FALSE ];
                         } else {
                             .lasso.mat <- .lasso.mat[ , c( .estimate.colname, rev( names( sorted.cors.with.the.outcome ) ) ), drop = FALSE ];
                         }
                         .retained.covars <- colnames( .lasso.mat );
                     } # End if length( .interactors ) > 0
                  .needed.df <-
                    ( length( .retained.covars ) - ( nrow( .lasso.mat ) - MINIMUM.DF ) );
                  if( .needed.df > 0 ) {
                    # Then remove some covars until there is at least MINIMUM.DF degrees of freedom.
                    # They are in order, so just chop them off the end.
                    .retained.covars <-
                      .retained.covars[ 1:( length( .retained.covars ) - .needed.df ) ];
                  }
                  if( ( .estimate.colname == "none" ) || ( .estimate.colname %in% .retained.covars ) ) {
                    
                      .lasso.mat <- .lasso.mat[ , .retained.covars, drop = FALSE ];

                    # penalty.factor = 0 to force the .estimate.colname variable.
    
                    .penalty.factor <-
                        as.numeric( colnames( .lasso.mat ) != .estimate.colname )
                                
                    tryCatch(
                        {
                            .cv.glmnet.fit <- cv.glmnet( .lasso.mat, .out, intercept = include.intercept,
                                                        penalty.factor = .penalty.factor, grouped = FALSE, nfold = length( .out ) ); # grouped = FALSE to avoid the warning.;
                            .lasso.validation.results.per.person.coefs.cell <-
                                coef( .cv.glmnet.fit, s = "lambda.min" );
                            if( return.lasso.coefs ) {
                                .lasso.validation.results.per.person.coefs.row[[ .col.i ]] <-
                                    .lasso.validation.results.per.person.coefs.cell;
                            }
                            .mf.ptid <- stats::model.frame( as.formula( .formula ), data = regression.df );
                            
                            .lasso.mat.ptid <-
                                model.matrix(as.formula( .formula ), .mf.ptid[ the.rows.for.ptid, , drop = FALSE ] );
                            for( .row.i.j in 1:length(the.rows.for.ptid) ) {
                                .row.i <- the.rows.for.ptid[ .row.i.j ];
                                .newx <-
                                    c( 1, .lasso.mat.ptid[ .row.i.j, rownames( as.matrix( .lasso.validation.results.per.person.coefs.cell ) )[-1], drop = FALSE ] );
                                names( .newx ) <- c( "(Intercept)", rownames( as.matrix( .lasso.validation.results.per.person.coefs.cell ) )[-1] );
                                        # Note that the intercept is always the first one, even when include.intercept == FALSE
                                .pred.value.lasso <-
                                    sum( as.matrix( .lasso.validation.results.per.person.coefs.cell ) * .newx );
                                lasso.validation.results.per.person[ .row.i, .col.i ] <- 
                                    .pred.value.lasso;
                            }
                        },
                    error = function( e )
                    {
                        if( .col.i == 1 ) {
                            warning( paste( "ptid", .ptid.i, "col", .col.i, "lasso failed with error", e, "\nReverting to simple regression with only an intrecept." ) );
                            .formula <- as.formula( "days.since.infection ~ 1" );
                        } else {
                            warning( paste( "ptid", .ptid.i, "col", .col.i, "lasso failed with error", e, "\nReverting to simple regression vs", .estimate.colname ) );
                            if( include.intercept ) {
                                .formula <- as.formula( paste( "days.since.infection ~ 1 + ", .estimate.colname ) );
                            } else {
                                .formula <- as.formula( paste( "days.since.infection ~ 0 + ", .estimate.colname ) );
                            }
                        }
                        .lm <- lm( .formula, data = regression.df.without.ptid.i );
                        if( return.lasso.coefs ) {
                          ## It's nothing, no coefs selected, so leave it as NA.
                        }
                        for( .row.i in the.rows.for.ptid ) {
                          .pred.value.lasso <-
                              predict( .lm, regression.df[ .row.i, , drop = FALSE ] );
                          lasso.validation.results.per.person[ .row.i, .col.i ] <<- 
                             .pred.value.lasso;
                        }
                   },
                   finally = {}
                   );
                  } # End if the lasso estimate variable is usable
                  
                  # lasso.withbounds:
                  if( .estimate.colname == "none" ) {
                      .lasso.withbounds.mat <-
                        as.matrix( regression.df.without.ptid.i[ , .covariates.lasso.withbounds ] );
                  } else {
                      .lasso.withbounds.mat <-
                        as.matrix( regression.df.without.ptid.i[ , c( .covariates.lasso.withbounds, .estimate.colname ) ] );
                  }
                  mode( .lasso.withbounds.mat ) <- "numeric";
                  
                  # exclude any rows with any NAs.
                  .retained.rows <-
                      which( apply( .lasso.withbounds.mat, 1, function( .row ) { !any( is.na( .row ) ) } ) );
                  .lasso.withbounds.mat <-
                      .lasso.withbounds.mat[ .retained.rows, , drop = FALSE ];
                  .out <- regression.df.without.ptid.i[[ "days.since.infection" ]][ .retained.rows ];
                  .covars.to.exclude <- apply( .lasso.withbounds.mat, 2, function ( .col ) {
                      return( ( sum( !is.na( .col ) ) <= 1 ) || ( var( .col, na.rm = TRUE  ) == 0 ) );
                  } );
                  .covars.to.exclude <- names( which( .covars.to.exclude ) );

                  cors.with.the.outcome <-
                    sapply( setdiff( colnames( .lasso.withbounds.mat ), c( .covars.to.exclude, .estimate.colname ) ), function( .covar.colname ) {
                          .cor <- 
                              cor( .out, .lasso.withbounds.mat[ , .covar.colname ], use = "pairwise" );
                            return( .cor );
                        } );
                  sorted.cors.with.the.outcome <-
                    cors.with.the.outcome[ order( abs( cors.with.the.outcome ) ) ];
                  # Sort the columns of .lasso.withbounds.mat by their correlation with the outcome. This is to ensure that more-relevant columns get selected when removing columns due to pairwise correlation among them.  We want to keep the best one.
                  if( .estimate.colname == "none" ) {
                    .lasso.withbounds.mat <- .lasso.withbounds.mat[ , rev( names( sorted.cors.with.the.outcome ) ), drop = FALSE ];
                  } else {
                    .lasso.withbounds.mat <- .lasso.withbounds.mat[ , c( .estimate.colname, rev( names( sorted.cors.with.the.outcome ) ) ), drop = FALSE ];
                  }
                  
                  # Exclude covars that are not sufficiently correlated with the outcome.
                      .covars.to.exclude <- c( .covars.to.exclude,
                          names( which( sapply( setdiff( colnames( .lasso.withbounds.mat ), c( .covars.to.exclude, .lower.bound.colname, .upper.bound.colname, .estimate.colname ) ), function( .covar.colname ) {
                            #print( .covar.colname );
                          .cor <- 
                              cor( .out, .lasso.withbounds.mat[ , .covar.colname ], use = "pairwise" );
                            #print( .cor );
                          if( is.na( .cor ) || ( abs( .cor ) <= MINIMUM.CORRELATION.WITH.OUTCOME ) ) {
                            TRUE
                          } else {
                            FALSE
                          }
                        } ) ) ) );
                  
                   # Exclude covars that are too highly correlated with the estimate.  Here we can exclude the upper or lower bound too if it is too highly correlated with the estimate.
                    if( .estimate.colname != "none" ) {
                      .covars.to.exclude <- c( .covars.to.exclude,
                          names( which( sapply( setdiff( colnames( .lasso.withbounds.mat ), c( .covars.to.exclude, .estimate.colname ) ), function( .covar.colname ) {
                            #print( .covar.colname );
                          .cor <- 
                              cor( .lasso.withbounds.mat[ , .estimate.colname ], .lasso.withbounds.mat[ , .covar.colname ], use = "pairwise" );
                            #print( .cor );
                          if( is.na( .cor ) || ( abs( .cor ) >= COR.THRESHOLD ) ) {
                            TRUE
                          } else {
                            FALSE
                          }
                        } ) ) ) );
                    }

                  # Exclude covars that are too highly correlated with the upper bound.
                  if( !( .upper.bound.colname %in% .covars.to.exclude ) ) {
                      .covars.to.exclude <- c( .covars.to.exclude,
                          names( which( sapply( setdiff( colnames( .lasso.withbounds.mat ), c( .covars.to.exclude, .estimate.colname, .upper.bound.colname ) ), function( .covar.colname ) {
                            #print( .covar.colname );
                          .cor <- 
                              cor( .lasso.withbounds.mat[ , .upper.bound.colname ], .lasso.withbounds.mat[ , .covar.colname ], use = "pairwise" );
                            #print( .cor );
                          if( is.na( .cor ) || ( abs( .cor ) >= COR.THRESHOLD ) ) {
                            TRUE
                          } else {
                            FALSE
                          }
                        } ) ) ) );
                    }
                  
                  # Exclude covars that are too highly correlated with each other (not including bounds)
                  .covars.to.consider <-
                    setdiff( colnames( .lasso.withbounds.mat ), c( .covars.to.exclude, .estimate.colname, .lower.bound.colname, .upper.bound.colname ) );
                  ## Process these in reverse order to ensure that we prioritize keeping those towards the top.
                  .covars.to.consider <- rev( .covars.to.consider );
                  .new.covars.to.exclude <- rep( FALSE, length( .covars.to.consider ) );
                  names( .new.covars.to.exclude ) <- .covars.to.consider;
                  for( .c.i in 1:length( .covars.to.consider ) ) {
                    .covar.colname <- .covars.to.consider[ .c.i ];
                    #print( .covar.colname );
                    # Only consider those not already excluded.
                    .cor <- 
                      cor( .lasso.withbounds.mat[ , .covar.colname ], .lasso.withbounds.mat[ , names( which( !.new.covars.to.exclude ) ) ], use = "pairwise" );
                    #print( .cor );
                        if( ( length( .cor ) > 0 ) && any( .cor[ !is.na( .cor ) ] < 1 & .cor[ !is.na( .cor ) ] >= COR.THRESHOLD ) ) {
                          .new.covars.to.exclude[ .c.i ] <- TRUE;
                        } else {
                          .new.covars.to.exclude[ .c.i ] <- FALSE;
                      }
                   } # End foreach of the .covars.to.consider
                  .covars.to.exclude <- c( .covars.to.exclude, names( which( .new.covars.to.exclude ) ) );
                  
                  .retained.covars <-
                      setdiff( colnames( .lasso.withbounds.mat ), .covars.to.exclude );

                  .covariates.lasso.withbounds <- intersect( .retained.covars, .covariates.lasso.withbounds );
                  if( .estimate.colname == "none" ) {
                      ## SEE NOTE BELOW ABOUT USE OF AN INTERCEPT. WHEN we're using "none", then we always do use an intercept (only for the non-bounded result).
                        if( include.intercept ) {
                            .formula.withbounds <- paste( "days.since.infection ~", paste( intersect( .retained.covars, .covariates.lasso.withbounds ), collapse = "+" ) );
                        } else {
                            .covariates.lasso.withbounds.nointercept <- .covariates.lasso.withbounds;
                            .formula.withbounds.nointercept <- paste( "days.since.infection ~ 0 + ", paste( intersect( .retained.covars, .covariates.lasso.withbounds.nointercept ), collapse = "+" ) );
                            .formula.withbounds <- .formula.withbounds.nointercept;
                        }
                    } else {
                        ## NOTE WE MUST BE CAREFUL ABOUT USING AN INTERCEPT, BECAUSE THEN WE JUST CHEAT BY PREDICTING EVERYTHING IS CENTERED AT the tight true bounds on the days.since.infection in our training data (when looking at one time at a time, and now with 6m.not.1m this should be true for both times as well.).
                      if( include.intercept ) {
                          .formula.withbounds <- paste( "days.since.infection ~", paste( intersect( .retained.covars, c( .covariates.lasso.withbounds, .estimate.colname ) ), collapse = "+" ) );
                      } else {
                          .covariates.lasso.withbounds.nointercept <- .covariates.lasso.withbounds;
                          .formula.withbounds.nointercept <-
                              paste( "days.since.infection ~ 0 + ", paste( intersect( .retained.covars, c( .covariates.lasso.withbounds.nointercept, .estimate.colname ) ), collapse = "+" ) );
                          .formula.withbounds <- .formula.withbounds.nointercept;
                      }
                    }
                
                .interactors <-
                         intersect( .retained.covars, helpful.additional.cols.with.interactions );
                     # Add interactions among the interactors.
                     if( length( .interactors ) > 1 ) {
                         ## NOTE that for now it is up to the caller to ensure that the df is not too high in this case.
                         .interactions.among.interactors <- c();
                         for( .interactor.i in 1:( length( .interactors ) - 1 ) ) {
                             .interactor <- .interactors[ .interactor.i ];
                             .remaining.interactors <-
                                 .interactors[ ( .interactor.i + 1 ):length( .interactors ) ];
                             .interactions.among.interactors <-
                                 c( .interactions.among.interactors,
                                   paste( .interactor, .remaining.interactors, collapse = "+", sep = ":" ) );
                         } # End foreach .interactor.i
                         .interactors <- c( .interactors, .interactions.among.interactors );
                     } # End if there's more than one "interactor"
                     if( length( .interactors ) > 0 ) {
                         .interactees <-
                             setdiff( .retained.covars, .interactors );
                         if( length( .interactees ) > 0 ) {
                             ## NOTE that for now it is up to the caller to ensure that the df is not too high in this case.
                             .interactions.formula.part.components <- c();
                             for( .interactor.i in 1:length( .interactors ) ) {
                                 .interactor <- .interactors[ .interactor.i ];
                                 .interactions.formula.part.components <-
                                     c( .interactions.formula.part.components,
                                       paste( .interactor, .interactees, collapse = "+", sep = ":" ) );
                             } # End foreach .interactor.i
                                                                      
                             .interactions.formula.part <-
                                 paste( .interactions.formula.part.components, collapse =" + " );
                             .formula.withbounds <-
                                 paste( .formula.withbounds, .interactions.formula.part, sep = " + " );
                         } # End if there's any "interactees"
                         ## Also add interactions with the estimate colname and everything included so far.
                         if( .estimate.colname != "none" ) {
                             .everything.included.so.far.in.formula.withbounds <- 
                                 strsplit( .formula.withbounds, split = "\\s*\\+\\s*" )[[1]];

                             # Add interactions.
                             .new.interactions.with.estimator.in.formula.withbounds <-
                                 sapply( .everything.included.so.far.in.formula.withbounds[ -1 ], function ( .existing.part ) {
                                     paste( .existing.part, .estimate.colname, sep = ":" )
                                 } );

                             .formula.withbounds <- paste( c( .everything.included.so.far.in.formula.withbounds, .new.interactions.with.estimator.in.formula.withbounds ), collapse = " + " );
                         } # End if .estimate.colname != "none"

                         .mf.withbounds <- stats::model.frame( as.formula( .formula.withbounds ), data = regression.df.without.ptid.i );
                         .lasso.withbounds.mat <- model.matrix(as.formula( .formula.withbounds ), .mf.withbounds);
                         
                         .covars.to.exclude <- apply( .lasso.withbounds.mat, 2, function ( .col ) {
                             return( ( sum( !is.na( .col ) ) <= 1 ) || ( var( .col, na.rm = TRUE  ) == 0 ) );
                         } );
                         .covars.to.exclude <- names( which( .covars.to.exclude ) );
                         
                         cors.with.the.outcome <-
                             sapply( setdiff( colnames( .lasso.withbounds.mat ), c( .covars.to.exclude, .estimate.colname ) ), function( .covar.colname ) {
                          .cor <- 
                              cor( .out, .lasso.withbounds.mat[ , .covar.colname ], use = "pairwise" );
                            return( .cor );
                         } );
                         sorted.cors.with.the.outcome <-
                             cors.with.the.outcome[ order( abs( cors.with.the.outcome ) ) ];
                         # Sort the columns of .lasso.withbounds.mat by their correlation with the outcome. This is to ensure that more-relevant columns get selected when removing columns due to pairwise correlation among them.  We want to keep the best one.
                         if( .estimate.colname == "none" ) {
                             .lasso.withbounds.mat <- .lasso.withbounds.mat[ , rev( names( sorted.cors.with.the.outcome ) ), drop = FALSE ];
                         } else {
                             .lasso.withbounds.mat <- .lasso.withbounds.mat[ , c( .estimate.colname, rev( names( sorted.cors.with.the.outcome ) ) ), drop = FALSE ];
                         }
                         .retained.covars <- colnames( .lasso.withbounds.mat );
                     } # End if length( .interactors ) > 0
                  
                  .needed.df <-
                    ( length( .retained.covars ) - ( nrow( .lasso.withbounds.mat ) - MINIMUM.DF ) );
                  if( .needed.df > 0 ) {
                    # Then remove some covars until there is at least MINIMUM.DF degrees of freedom.
                    # They are in order, so just chop them off the end.
                    .retained.covars <-
                      .retained.covars[ 1:( length( .retained.covars ) - .needed.df ) ];
                  }
                
                  if( ( .estimate.colname == "none" ) || ( .estimate.colname %in% .retained.covars ) ) {
                    .lasso.withbounds.mat <-
                      .lasso.withbounds.mat[ , .retained.covars, drop = FALSE ];
                    # penalty.factor = 0 to force the .estimate.colname variable.
    
                    tryCatch(
                    {
                      .cv.glmnet.fit.withbounds <-
                        cv.glmnet( .lasso.withbounds.mat, .out, intercept = include.intercept,
                                  penalty.factor = as.numeric( colnames( .lasso.withbounds.mat ) != .estimate.colname ), grouped = FALSE, nfold = length( .out ) ); # grouped = FALSE to avoid the warning.
                      .lasso.withbounds.validation.results.per.person.coefs.cell <-
                          coef( .cv.glmnet.fit.withbounds, s = "lambda.min" );
                      if( return.lasso.coefs ) {
                        .lasso.withbounds.validation.results.per.person.coefs.row[[ .col.i ]] <-
                            .lasso.withbounds.validation.results.per.person.coefs.cell;
                      }

                      .mf.ptid <- stats::model.frame( as.formula( .formula.withbounds ), data = regression.df );
                            
                      .lasso.withbounds.mat.ptid <- model.matrix(as.formula( .formula.withbounds ), .mf.ptid[ the.rows.for.ptid, , drop = FALSE ] );
                      for( .row.i.j in 1:length(the.rows.for.ptid) ) {
                          .row.i <- the.rows.for.ptid[ .row.i.j ];
                        .newx <-
                            c( 1, .lasso.withbounds.mat.ptid[ .row.i.j, rownames( .lasso.withbounds.validation.results.per.person.coefs.cell )[ -1 ], drop = FALSE ] );
                          names( .newx ) <- c( "(Intercept)", rownames( as.matrix( .lasso.withbounds.validation.results.per.person.coefs.cell ) )[-1] );

                        .pred.value.lasso.withbounds <-
                          sum( as.matrix( .lasso.withbounds.validation.results.per.person.coefs.cell ) * .newx );
                          #predict( .cv.glmnet.fit.withbounds, newx = .newx, s = "lambda.min" );
                        lasso.withbounds.validation.results.per.person[ .row.i, .col.i ] <- 
                          .pred.value.lasso.withbounds;
                      }
                    },
                    error = function( e )
                    {
                        warning( paste( "ptid", .ptid.i, "col", .col.i, "lasso withbounds failed with error", e, "\nReverting to simple regression vs", .estimate.colname ) );
                        if( include.intercept ) {
                            .formula <- as.formula( paste( "days.since.infection ~ 1 + ", .estimate.colname ) );
                        } else {
                            .formula <- as.formula( paste( "days.since.infection ~ 0 + ", .estimate.colname ) );
                        }
                        for( .row.i in the.rows.for.ptid ) {
                          .pred.value.lasso.withbounds <-
                              predict( lm( .formula, data = regression.df.without.ptid.i ), regression.df[ .row.i, , drop = FALSE ] );
                          lasso.withbounds.validation.results.per.person[ .row.i, .col.i ] <- 
                              .pred.value.lasso.withbounds;
                        }
                   },
                   finally = {}
                   );
                  } # End if the lasso.withbounds estimate variable is usable

                } # End if use.lasso.validate
    
            } # End foreach .col.i
            if( use.lasso.validate && return.lasso.coefs ) {
                lasso.validation.results.per.person.coefs[[ .ptid.i ]] <- 
                    .lasso.validation.results.per.person.coefs.row;
                lasso.withbounds.validation.results.per.person.coefs[[ .ptid.i ]] <- 
                    .lasso.withbounds.validation.results.per.person.coefs.row;
            } # End if use.lasso.validate
        } # End foreach .ptid.i
        if( use.glm.validate ) {
          ## glm:
            colnames( glm.validation.results.per.person ) <-
                paste( "glm.validation.results", estimate.cols, sep = "." );
            rownames( glm.validation.results.per.person ) <-
                rownames( regression.df );
            results.per.person <-
                cbind( results.per.person,
                      glm.validation.results.per.person );
          ## glm.withbounds:
            colnames( glm.withbounds.validation.results.per.person ) <-
                paste( "glm.withbounds.validation.results", estimate.cols, sep = "." );
            rownames( glm.withbounds.validation.results.per.person ) <-
                rownames( regression.df );
            results.per.person <-
                cbind( results.per.person,
                      glm.withbounds.validation.results.per.person );
        }
        if( use.lasso.validate ) {
          ## lasso:
          colnames( lasso.validation.results.per.person ) <-
            paste( "lasso.validation.results", estimate.cols, sep = "." );
          rownames( lasso.validation.results.per.person ) <-
            rownames( regression.df );
          results.per.person <-
            cbind( results.per.person,
                  lasso.validation.results.per.person );
          ## lasso.withbounds:
          colnames( lasso.withbounds.validation.results.per.person ) <-
            paste( "lasso.withbounds.validation.results", estimate.cols, sep = "." );
          rownames( lasso.withbounds.validation.results.per.person ) <-
            rownames( regression.df );
          results.per.person <-
            cbind( results.per.person,
                  lasso.withbounds.validation.results.per.person );
        }
      } # End if use.glm.validate || use.lasso.validate
      
      ## For fairness in evaluating when some methods
      ## completely fail to give a result, (so it's NA
      ## presently), we change all of these estimates
      ## from NA to 0.  When we put bounds on the
      ## results, below, they will be changed from 0 to
      ## a boundary endpoint if 0 is outside of the
      ## bounds.
      results.per.person.zeroNAs <- apply( results.per.person, 1:2, function( .value ) {
          if( is.na( .value ) ) {
              0
          } else {
              .value
          }
      } );
      
       ## unbounded results:
      diffs.by.stat.zeroNAs <- compute.diffs.by.stat( results.per.person.zeroNAs, days.since.infection );
      diffs.by.stat <- compute.diffs.by.stat( results.per.person, days.since.infection );

      get.results.list.for.bounds.type <- function ( diffs.by.stat, diffs.by.stat.zeroNAs, bounds.type ) {
          .results <- 
              list( bias = lapply( diffs.by.stat, mean, na.rm = T ), se = lapply( diffs.by.stat, sd, na.rm = T ), rmse = lapply( diffs.by.stat, rmse, na.rm = T ), n = lapply( diffs.by.stat, function( .vec ) { sum( !is.na( .vec ) ) } ), bias.zeroNAs = lapply( diffs.by.stat.zeroNAs, mean, na.rm = T ), se.zeroNAs = lapply( diffs.by.stat.zeroNAs, sd, na.rm = T ), rmse.zeroNAs = lapply( diffs.by.stat.zeroNAs, rmse, na.rm = T ), n.zeroNAs = lapply( diffs.by.stat.zeroNAs, function( .vec ) { sum( !is.na( .vec ) ) } ) );
       if( use.lasso.validate && return.lasso.coefs ) {
         #.results <- c( .results, list( lasso.coefs = list( lasso = lasso.validation.results.per.person.coefs, lasso.withbounds = lasso.withbounds.validation.results.per.person.coefs, lasso.nointercept = lasso.nointercept.validation.results.per.person.coefs, lasso.nointercept.withbounds = lasso.nointercept.withbounds.validation.results.per.person.coefs ) ) );
         .results <- c( .results, list( lasso.coefs = list( lasso = lasso.validation.results.per.person.coefs, lasso.withbounds = lasso.withbounds.validation.results.per.person.coefs ) ) );
       }
       
          results.list <- list();
          results.list[[ bounds.type ]] <- .results;

       ## If there are multi-timepoint and multi-region predictors, also include per-region, per-timepoint diffs.by.stat results.
       if( !is.null( ppt.suffices ) ) {
           unique.ppt.suffices <- unique( ppt.suffices );
           .results.by.suffix <- lapply( unique.ppt.suffices, function ( .ppt.suffix ) {
               ## TODO: REMOVE
               #print( paste( "suffix:", .ppt.suffix ) );
               .diffs.by.stat <-
                   lapply( diffs.by.stat, function( .diffs.for.stat ) { .diffs.for.stat[ ppt.suffices == .ppt.suffix ]; } );
               .diffs.by.stat.zeroNAs <-
                   lapply( diffs.by.stat.zeroNAs, function( .diffs.for.stat ) { .diffs.for.stat[ ppt.suffices == .ppt.suffix ]; } );
               list( bias = lapply( .diffs.by.stat, mean, na.rm = T ), se = lapply( .diffs.by.stat, sd, na.rm = T ), rmse = lapply( .diffs.by.stat, rmse, na.rm = T ), n = lapply( .diffs.by.stat, function( .vec ) { sum( !is.na( .vec ) ) } ), bias.zeroNAs = lapply( .diffs.by.stat.zeroNAs, mean, na.rm = T ), se.zeroNAs = lapply( .diffs.by.stat.zeroNAs, sd, na.rm = T ), rmse.zeroNAs = lapply( .diffs.by.stat.zeroNAs, rmse, na.rm = T ), n.zeroNAs = lapply( .diffs.by.stat.zeroNAs, function( .vec ) { sum( !is.na( .vec ) ) } ) );
           } );
           names( .results.by.suffix ) <- paste( bounds.type, unique.ppt.suffices, sep = "" );
           results.list <- c( results.list, .results.by.suffix );
       } # End if there are suffices, also include results by suffix.

        return( results.list );
      } # get.results.list.for.bounds.type ( diffs.by.stat, diffs.by.stat.zeroNAs );
      results.list <- get.results.list.for.bounds.type( diffs.by.stat, diffs.by.stat.zeroNAs, "unbounded" );
       
       #if( use.glm.validate ) {
       #  results.list <- c( results.list, list( glm.fit.statistics = glm.fit.statistics ) );
       #}
       
      ## bounded results:
      if( use.bounds ) {
    
        ## Ok, well, now we can also evaluate versions
        ## of variants of each method, but bounding the
        ## results.  Ie if the estimated time is within
        ## the bounds, that time is used, otherwise,
        ## the boundary.  Note we don't do this with
        ## the deterministic bounds, and we only do it
        ## for the time corresponding to the sample (
        ## mtn003 for "1m" and hvtn502 for "6m" ) [each
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
           #grep( ifelse( the.time == "6m", ifelse( the.time == "1m.6m", "1mmtn003_6msixmonths", "sixmonths" ), "mtn003" ), grep( "deterministic", names( the.artificial.bounds ), invert = TRUE, value = TRUE ), value = TRUE );
            grep( ifelse( the.time == "6m", ifelse( the.time == "1m.6m", "1mmtn003_6mhvtn502", "hvtn502" ), "mtn003" ), grep( "deterministic", names( the.artificial.bounds ), invert = TRUE, value = TRUE ), value = TRUE );
        ## Also note that NAs are bounded too, replaced by the lower bound.
        results.per.person.bounded <-
          lapply( .artificial.bounds.to.use, function ( .artificial.bounds.name ) {
            .mat <-
            apply( results.per.person, 2, function ( .results.column ) {
            sapply( names( .results.column ), function ( .ppt ) {
              .value <- .results.column[ .ppt ];
              if( !( .ppt %in% rownames( the.artificial.bounds[[ .artificial.bounds.name ]] ) ) ) {
                  stop( paste( .ppt, "has no bounds!" ) );
              }
              if( is.na( .value ) ) {
                return( the.artificial.bounds[[ .artificial.bounds.name ]][ .ppt, "lower" ] );
              } else {
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
            rownames( .mat ) <- rownames( results.per.person );
            return( .mat );
        } );
        names( results.per.person.bounded ) <- .artificial.bounds.to.use;
    
        bounded.results.by.bounds.type <- lapply( .artificial.bounds.to.use, function ( .bounds.type ) {
            .results.per.person <- results.per.person.bounded[[ .bounds.type ]];
            .diffs.by.stat <- compute.diffs.by.stat( .results.per.person, days.since.infection );
            .results.per.person.zeroNAs <-
                apply( .results.per.person, 1:2, function( .value ) {
                    if( is.na( .value ) ) {
                        0
                    } else {
                        .value
                    }
                } );

            .diffs.by.stat.zeroNAs <- compute.diffs.by.stat( .results.per.person.zeroNAs, days.since.infection );
            .results.list <-
                get.results.list.for.bounds.type( .diffs.by.stat, .diffs.by.stat.zeroNAs, .bounds.type );

          return( .results.list );
        } );
        return( c( results.list, do.call( c, bounded.results.by.bounds.type ) ) );
      } else {
        return( results.list );
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

        # Note that we use the pvl at the earliest time ie for "1m6m" we use timepoint 2.
        if( the.region == "nflg" || ( length( grep( "rv217", the.region ) ) > 0 ) ) {
            pvl.at.the.time <- sapply( rownames( results.in ), function( .ptid ) { as.numeric( as.character( rv217.pvl.in[ ( rv217.pvl.in[ , "ptid" ] == .ptid ) & ( rv217.pvl.in[ , "timepoint" ] == ifelse( the.time == "6m", 3, 2 ) ), "viralload" ] ) ) } );
        } else {
            stopifnot( the.region == "v3" );
            pvl.at.the.time <- sapply( rownames( results.in ), function( .ptid ) { as.numeric( as.character( caprisa002.pvl.in[ ( caprisa002.pvl.in[ , "ptid" ] == .ptid ) & ( caprisa002.pvl.in[ , "timepoint" ] == ifelse( the.time == "6m", 3, 2 ) ), "viralload" ] ) ) } );
        }
        ## Add log plasma viral load (lPVL).
        results.with.lPVL <- cbind( results.in, log( pvl.at.the.time ) );
        colnames( results.with.lPVL )[ ncol( results.with.lPVL ) ] <- "lPVL";

        results.covars.per.person.with.extra.cols <-
          summarizeCovariatesOnePerParticipant( results.with.lPVL );

        # #### TODO: REMOVE. TESTING VALUE OF KNOWING IS.MULTIPLE
        # if( the.region == "nflg" || ( length( grep( "rv217", the.region ) ) > 0 ) ) {
        #     gold.is.multiple <- rv217.gold.is.multiple[ rownames( results.covars.per.person.with.extra.cols ) ];
        # } else {
        #     stopifnot( the.region == "v3" );
        #     gold.is.multiple <- caprisa002.gold.is.multiple[ rownames( results.covars.per.person.with.extra.cols ) ];
        # }
        # ## Add gold.is.multiple
        # results.covars.per.person.with.extra.cols <- cbind( results.covars.per.person.with.extra.cols, gold.is.multiple );
        # colnames( results.covars.per.person.with.extra.cols )[ ncol( results.covars.per.person.with.extra.cols ) ] <-
        #   "gold.is.multiple";
        
        days.colnames <- c( grep( "time", colnames( results.in ), value = T ), grep( "days", colnames( results.in ), value = T ) );
        
        days.est.colnames <- grep( "est", days.colnames, value = TRUE );
        days.est <- results.in[ , days.est.colnames, drop = FALSE ];
        lambda.est.colnames <-
            gsub( "PFitter\\.lambda\\.est", "PFitter.lambda", gsub( "(?:days|time|fits)", "lambda", days.est.colnames, perl = TRUE ) );
        stopifnot( all( lambda.est.colnames %in% colnames( results.in ) ) );
        days.est.colnames.nb <- gsub( "[^\\.]+\\.Star[Pp]hy", "PFitter", gsub( "(?:days|time|fits).*$", "nbases", days.est.colnames, perl = TRUE ) );
        days.est.nb <- results.in[ , days.est.colnames.nb, drop = FALSE ];
        days.est.colnames.nseq <- gsub( "[^\\.]+\\.Star[Pp]hy", "PFitter", gsub( "(?:days|time|fits).*$", "nseq", days.est.colnames, perl = TRUE ) );
        days.est.nseq <- results.in[ , days.est.colnames.nseq, drop = FALSE ];
        
        results <- results.with.lPVL[ , days.est.colnames, drop = FALSE ];
        
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
        } # End if use.anchre and the.time is 1m6m, add anchre results too.
       
        if( is.na( partition.size ) ) {
          ## Now the issue is that there are multiple input files per ppt, eg for the NFLGs ther are often "LH" and "RH" files.  What to do?  The number of sequences varies.  Do a weighted average.
          .weights <- days.est.nseq*days.est.nb;
          results.per.person <-
              compute.results.per.person( results, .weights );
        
          if( use.bounds ) {
              the.artificial.bounds <- getArtificialBounds( the.region, the.time, results.dirname );
              # Only keep the "sampledwidth" bounds.
              the.artificial.bounds <-
                the.artificial.bounds[ grep( "sampledwidth", names( the.artificial.bounds ), value = TRUE ) ];
              
              center.of.bounds.table <- sapply( names( the.artificial.bounds ), function ( .artificial.bounds.name ) {
                round( apply( the.artificial.bounds[[ .artificial.bounds.name ]], 1, mean ) )
              } );
              colnames( center.of.bounds.table ) <-
                paste( "COB", gsub( "_", ".", colnames( center.of.bounds.table ) ), "time.est", sep = "." );
              
              results.per.person <-
                cbind( results.per.person, center.of.bounds.table[ rownames( results.per.person ), , drop = FALSE ] );

              ## TODO: REMOVE. TEMPORARY HACK TO REMOVE OLD RESULTS. Remove "onemonth" and "sixmonths" results.
              results.per.person <-
                results.per.person[ , grep( "uniform\\.(one|six)month", colnames( results.per.person ), invert = TRUE ), drop = FALSE ];
              
              return( list( results.per.person = results.per.person, days.since.infection = days.since.infection, results.covars.per.person.with.extra.cols = results.covars.per.person.with.extra.cols, bounds = the.artificial.bounds, evaluated.results = bound.and.evaluate.results.per.ppt( results.per.person, days.since.infection, results.covars.per.person.with.extra.cols, the.time, the.artificial.bounds ) ) );
            } else {
              return( list( results.per.person = results.per.person, days.since.infection = days.since.infection, results.covars.per.person.with.extra.cols = results.covars.per.person.with.extra.cols, evaluated.results = bound.and.evaluate.results.per.ppt( results.per.person, days.since.infection, results.covars.per.person.with.extra.cols, the.time ) ) );
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
                .thissample.results.per.person <-
                    sapply( colnames( results ), function( est.name ) {
                        sapply( names( .thissample.the.partition.ids ), function( .ppt ) {
                            unname( results.per.ppt.by.id[[ est.name ]][[ .ppt ]][ .thissample.the.partition.ids[ .ppt ] ] )
                        } ) } );
                return( list( results.per.person = .thissample.results.per.person, evaluated.results = bound.and.evaluate.results.per.ppt( .thissample.results.per.person, days.since.infection, results.covars.per.person.with.extra.cols, the.time, the.artificial.bounds ) ) );
          } # do.one.sample (..)
            
          set.seed( partition.bootstrap.seed );
          bootstrap.results <- 
              mclapply( 1:partition.bootstrap.samples, do.one.sample, mc.cores = partition.bootstrap.num.cores );
        
            matrix.of.unbounded.results.rmses <- sapply( bootstrap.results, function( .results.for.bootstrap ) { .results.for.bootstrap[[ "evaluated.results" ]][[ "unbounded" ]]$rmse } );
            mode( matrix.of.unbounded.results.rmses ) <- "numeric";
            ## This uses the second position, which is the first of the unbounded ones, and for now the only one.  It's "mtn003" unless the.time is "1m" in which case it is "hvtn502".
            matrix.of.bounded.results.rmses <- sapply( bootstrap.results, function( .results.for.bootstrap ) { .results.for.bootstrap[[ "evaluated.results" ]][[ 2 ]]$rmse } );
            mode( matrix.of.bounded.results.rmses ) <- "numeric";
            
            #hist( apply( matrix.of.unbounded.results.rmses, 1, diff ) )
            
          return( list( summary = the.summary, mean.bias.per.ppt.by.est = mean.bias.per.ppt, sd.bias.per.ppt.by.est = sd.bias.per.ppt, num.partitions.per.ppt = num.partitions.per.ppt, days.since.infection = days.since.infection, results.covars.per.person.with.extra.cols = results.covars.per.person.with.extra.cols, bounds = the.artificial.bounds, bootstrap.results = bootstrap.results, bootstrap.unbounded.rmse = matrix.of.unbounded.results.rmses, bootstrap.bounded.rmse = matrix.of.bounded.results.rmses ) );
        } # End if is.na( partition.size ) .. else ..
        
    } # get.timings.results.for.region.and.time (..)
    
    getTimingsResultsByRegionAndTime <- function ( partition.size = NA ) {
      getResultsByRegionAndTime( gold.standard.varname = "days.since.infection", get.results.for.region.and.time.fn = get.timings.results.for.region.and.time, evaluate.results.per.person.fn = bound.and.evaluate.results.per.ppt, partition.size = partition.size, regions = regions, times = times )
    } # getTimingsResultsByRegionAndTime (..)
    
    if( force.recomputation || !file.exists( results.by.region.and.time.Rda.filename ) ) {
        results.by.region.and.time <- getTimingsResultsByRegionAndTime();
        save( results.by.region.and.time, file = results.by.region.and.time.Rda.filename );
    } else {
        # loads results.by.region.and.time
        load( file = results.by.region.and.time.Rda.filename );
    }

    if( include.intercept ) {
        writeResultsTables( results.by.region.and.time, "_evaluateTimings_include_intercept.tab", regions = regions, results.are.bounded = TRUE );
    } else {
        writeResultsTables( results.by.region.and.time, "_evaluateTimings.tab", regions = regions, results.are.bounded = TRUE );
    }
    
    ## TODO
    ##  *) select a set of best predictions from isMultiple to use here (instead of gold.standard.varname) -- if there's any benefit to using gold.standard.varname.
    ##  *) check out results of the partitions of evaluateTimings.

    if( FALSE ) {
        get.rmses <- function ( evaluate.regions = train.regions, evaluate.times = train.times, train.regions = c( "nflg", "v3" ), train.times = c( "1m", "6m" ) ) {
            if( length( train.times ) == 2 ) {
                train.time <- "1m.6m";
                the.bound <- "sampledwidth_uniform_1mmtn003_6mhvtn502";
            } else {
                train.time <- train.times;
                if( train.time == "6m" ) {
                    the.bound <- "sampledwidth_uniform_hvtn502";
                    stopifnot( length( evaluate.times ) == 1 );
                    stopifnot( evaluate.times == "6m" );
                } else {
                    the.bound <- "sampledwidth_uniform_mtn003";
                    stopifnot( length( evaluate.times ) == 1 );
                    stopifnot( evaluate.times == "1m" );
                }
            }
            if( length( evaluate.times ) == 2 ) {
                # Leave the.bound alone.
            } else if( length( train.times ) == 2 ) {
                the.bound <- paste( the.bound, evaluate.times, sep = "." );
            }
            if( length( train.regions ) == 2 ) {
                .results.for.region <- results.by.region.and.time[[3]][[1]][[1]];
              if( length( evaluate.regions ) == 2 ) {
                  # Leave the.bound alone, then.
              } else {
                the.bound <- paste( the.bound, evaluate.regions, sep = "." );
              }
            } else {
              if( length( evaluate.regions ) == 2 ) {
                  stop( "can't evaluate more regions than trained" );
              } else if( evaluate.regions != train.regions ) {
                  stop( "can't evaluate a different region than trained" );
              }
              .results.for.region <- results.by.region.and.time[[ train.regions ]];
            }
            .lst <- sort( unlist( .results.for.region[[ train.time ]][[ "evaluated.results" ]][[ the.bound ]][[ "rmse.zeroNAs" ]] ), decreasing = T );
            ## TODO: REMOVE. Temporary.
            .lst <- .lst[ grep( "(one|six)month", names( .lst ), value = TRUE, invert = TRUE ) ];
            return( .lst );
        } # get.rmses (..)

        # This returns biases in order of rmses.
        get.biases <- function ( evaluate.regions = train.regions, evaluate.times = train.times, train.regions = c( "nflg", "v3" ), train.times = c( "1m", "6m" ) ) {
            if( length( train.times ) == 2 ) {
                train.time <- "1m.6m";
                the.bound <- "sampledwidth_uniform_1mmtn003_6mhvtn502";
            } else {
                train.time <- train.times;
                if( train.time == "6m" ) {
                    the.bound <- "sampledwidth_uniform_hvtn502";
                    stopifnot( length( evaluate.times ) == 1 );
                    stopifnot( evaluate.times == "6m" );
                } else {
                    the.bound <- "sampledwidth_uniform_mtn003";
                    stopifnot( length( evaluate.times ) == 1 );
                    stopifnot( evaluate.times == "1m" );
                }
            }
            if( length( evaluate.times ) == 2 ) {
                # Leave the.bound alone.
            } else if( length( train.times ) == 2 ) {
                the.bound <- paste( the.bound, evaluate.times, sep = "." );
            }
            if( length( train.regions ) == 2 ) {
                .results.for.region <- results.by.region.and.time[[3]][[1]][[1]];
              if( length( evaluate.regions ) == 2 ) {
                  # Leave the.bound alone, then.
              } else {
                the.bound <- paste( the.bound, evaluate.regions, sep = "." );
              }
            } else {
              if( length( evaluate.regions ) == 2 ) {
                  stop( "can't evaluate more regions than trained" );
              } else if( evaluate.regions != train.regions ) {
                  stop( "can't evaluate a different region than trained" );
              }
              .results.for.region <- results.by.region.and.time[[ train.regions ]];
            }
            .lst <- unlist( .results.for.region[[ train.time ]][[ "evaluated.results" ]][[ the.bound ]][[ "bias.zeroNAs" ]] )[ order( unlist( .results.for.region[[ train.time ]][[ "evaluated.results" ]][[ the.bound ]][[ "rmse.zeroNAs" ]] ), decreasing = T ) ];
            ## TODO: REMOVE. Temporary.
            .lst <- .lst[ grep( "(one|six)month", names( .lst ), value = TRUE, invert = TRUE ) ];
            return( .lst );
        } # get.biases (..)

        get.bias.and.rmse <- function ( ... ) { cbind( bias = get.biases( ... ), rmse = get.rmses( ... ) ) }
        
        get.uses <- function ( .varname = "none", withbounds = TRUE, regions = c( "nflg", "v3" ), times = c( "1m", "6m" ) ) {
            if( withbounds ) {
                .withbounds.string <- "lasso.withbounds";
            } else {
                .withbounds.string <- "lasso";
            }
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
                .results.for.region[[ the.time ]][[ "evaluated.results" ]][["unbounded"]][[ "lasso.coefs" ]][[ .withbounds.string ]];
            .uses <- sapply( 1:length( .results.by.removed.ptid ), function( .i ) { .dgCMatrix <- .results.by.removed.ptid[[.i]][[.varname]]; .rv <- as.logical( .dgCMatrix ); names( .rv ) <- rownames( .dgCMatrix ); return( .rv ); } );
            if( is.null( dim( .uses ) ) ) {
              table( names( which( unlist( .uses ) ) ) );
            } else {
              apply( .uses, 1, sum );
            }
        } # get.uses (..)

        ## ERE I AM. NOTE that include.intercept requires loading a separate set of data.
        evaluate.specific.model <-
          function ( model.vars, .include.intercept = include.intercept, step = FALSE, reload.data = ( .include.intercept != include.intercept ) ) {
            if( reload.data ) {
              if( include.intercept ) {
                  results.by.region.and.time.Rda.filename <-
                      paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/Timings.results.by.region.and.time.include.intercept.Rda", sep = "" );
              } else {
                  results.by.region.and.time.Rda.filename <-
                      paste( "/fh/fast/edlefsen_p/bakeoff_analysis_results/", results.dirname, "/Timings.results.by.region.and.time.Rda", sep = "" );
              }
              load( results.by.region.and.time.Rda.filename );
             } # End if( reload.data )
        results.covars.per.person.with.extra.cols <-
            results.by.region.and.time[[3]][[1]][[1]][[1]][["results.covars.per.person.with.extra.cols"]];
       results.covars.per.person.df <-
           data.frame( results.covars.per.person.with.extra.cols );
        ## DO NOT Undo conversion of the colnames (X is added before "6m.not.1m").  We want it to be called "X6m.not.1m" so it can work in the regression formulas.
        #colnames( results.covars.per.person.df ) <- colnames( results.covars.per.person.with.extra.cols );

        regression.df <- cbind( data.frame( days.since.infection = results.by.region.and.time[[3]][[1]][[1]][[1]][["days.since.infection" ]][ rownames( results.covars.per.person.df ) ] ), results.by.region.and.time[[3]][[1]][[1]][[1]][["results.per.person"]][ rownames( results.covars.per.person.df ), , drop = FALSE ], results.covars.per.person.df, results.by.region.and.time[[3]][[1]][[1]][[1]][["bounds" ]] );

        if( include.intercept ) {
            .formula <- as.formula( paste( "days.since.infection ~ ", paste( model.vars, collapse = "+" ) ) );
        } else {
            .formula <- as.formula( paste( "days.since.infection ~ 0 + ", paste( model.vars, collapse = "+" ) ) );
        }
        .lm <-
            suppressWarnings( lm( .formula, data = regression.df ) );
        if( step ) {
            .step.rv <- step( .lm ); # Stepwise regression, both forward and backward.
            return( summary( .step.rv ) );            
        }
       return( summary( .lm ) );
     } # evaluate.specific.model
                                        #        evaluate.specific.model( model.vars = c( "COB.sampledwidth.uniform.1mmtn003.6mhvtn502.time.est", "sampledwidth_uniform_1mmtn003_6mhvtn502.lower", "sampledwidth_uniform_1mmtn003_6mhvtn502.upper", "lPVL" ), step = TRUE );
     evaluate.specific.model( model.vars = c( "X6m.not.1m", "lPVL", "v3_not_nflg:lPVL","X6m.not.1m:v3_not_nflg:lPVL"  ), step = TRUE );
  # Start:  AIC=555.72
  # days.since.infection ~ 0 + X6m.not.1m + lPVL + v3_not_nflg:lPVL + 
  #     X6m.not.1m:v3_not_nflg:lPVL
  # 
  #                               Df Sum of Sq   RSS    AIC
  # <none>                                     17216 555.72
  # - X6m.not.1m:lPVL:v3_not_nflg  1    4306.3 21523 577.83
  # 
  # Call:
  # lm(formula = days.since.infection ~ 0 + X6m.not.1m + lPVL + v3_not_nflg:lPVL + 
  #     X6m.not.1m:v3_not_nflg:lPVL, data = regression.df)
  # 
  # Residuals:
  #     Min      1Q  Median      3Q     Max 
  # -33.702  -4.799   2.636  10.244  31.015 
  # 
  # Coefficients:
  #                             Estimate Std. Error t value Pr(>|t|)    
  # X6m.not.1m                  145.9451     2.8690  50.870  < 2e-16 ***
  # lPVL                          4.0772     0.1986  20.528  < 2e-16 ***
  # lPVL:v3_not_nflg              1.3970     0.3222   4.336 3.36e-05 ***
  # X6m.not.1m:lPVL:v3_not_nflg  -2.3826     0.4671  -5.100 1.53e-06 ***
  # ---
  # Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
  # 
  # Residual standard error: 12.87 on 104 degrees of freedom
  # Multiple R-squared:  0.9909,	Adjusted R-squared:  0.9906 
  # F-statistic:  2834 on 4 and 104 DF,  p-value: < 2.2e-16

        ## ERE I AM.STILL WORKING ON THAT . CHECK THIS. lPVL is up early down late, not unexpected!
# > plot( regression.df$days.since.infection[ regression.df$days.since.infection < 100 ], regression.df$lPVL[ regression.df$days.since.infection < 100 ] )
# > cor( regression.df$days.since.infection[ regression.df$days.since.infection < 100 ], regression.df$lPVL[ regression.df$days.since.infection < 100 ] )
# [1] 0.2237065
# > plot( regression.df$days.since.infection[ regression.df$days.since.infection > 100 ], regression.df$lPVL[ regression.df$days.since.infection > 100 ] )
# > cor( regression.df$days.since.infection[ regression.df$days.since.infection > 100 ], regression.df$lPVL[ regression.df$days.since.infection > 100 ] )
# [1] -0.2758797
        
        
    } # End if FALSE
    
    if( FALSE ) {
      ## For partition size == 10
      ## NOTE that this part is presently broken. TODO: FIX IT.
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
    } # End if FALSE
    
    return( invisible( NULL ) );
} # evaluateTimings (..)

## Here is where the action is.
evaluateTimings();
