
    missing.column.safe.rbind <- function ( matA, matB, matA.name, matB.name ) {
        all.colnames <- union( colnames( matA ), colnames( matB ) );
        .rv <- matrix( NA, nrow = nrow( matA ) + nrow( matB ), ncol = length( all.colnames ) );
        colnames( .rv ) <- all.colnames;
        rownames( .rv ) <-
            c( paste( rownames( matA ), matA.name, sep = "." ),
              paste( rownames( matB ), matB.name, sep = "." ) );
        .rv[ 1:nrow( matA ), colnames( matA ) ] <-
            as.matrix( matA );
        .rv[ nrow( matA ) + 1:nrow( matB ), colnames( matB ) ] <-
            as.matrix( matB );
        
        return( .rv );
    } # missing.column.safe.rbind (..)


getResultsByRegionAndTime <- function ( gold.standard.varname, get.results.for.region.and.time.fn, evaluate.results.per.person.fn, partition.size = NA, regions = c( "nflg", "v3", "rv217_v3" ), times = c( "1m", "6m", "1m6m" ) ) {
        if( !is.na( partition.size ) ) {
            regions <- "v3"; # Only v3 has partition results at this time.
        }
        results.by.region.and.time <-
            lapply( regions, function( the.region ) {
                ## TODO: REMOVE
                cat( the.region, fill = T );
           results.by.time <-
               lapply( times, function( the.time ) {
                   ## TODO: REMOVE
                   cat( the.time, fill = T );
                   get.results.for.region.and.time.fn( the.region, the.time, partition.size );               
               } );
           names( results.by.time ) <- times;

           .vars <- setdiff( names( results.by.time[[1]] ), "evaluated.results" );
           results.1m.6m <- lapply( .vars, function ( .varname ) {
               #print( .varname );
               if( .varname == "bounds" ) {
                   .rv <- 
                   lapply( names( results.by.time[[ "1m" ]][[ .varname ]] ), function( .bounds.type ) {
                       #print( .bounds.type );
                       missing.column.safe.rbind(
                           results.by.time[[ "1m" ]][[ .varname ]][[ .bounds.type ]],
                           results.by.time[[ "6m" ]][[ .varname ]][[ .bounds.type ]],
                           "1m",
                           "6m"
                       )
                   } );
                   names( .rv ) <-
                       names( results.by.time[[ "1m" ]][[ .varname ]] );
                   ## Add a new bounds.type called "uniform_1m5weeks_6m30weeks"
                   new.bounds.table <-
                       missing.column.safe.rbind(
                           results.by.time[[ "1m" ]][[ .varname ]][[ "uniform_5weeks" ]],
                           results.by.time[[ "6m" ]][[ .varname ]][[ "uniform_30weeks" ]],
                           "1m",
                           "6m"
                       );
                   .rv <- c( list( "uniform_1m5weeks_6m30weeks" = new.bounds.table ), .rv );
                   ## Add a new bounds.type called "exponentialwidth_uniform_1m5weeks_6m30weeks"
                   another.new.bounds.table <-
                       missing.column.safe.rbind(
                           results.by.time[[ "1m" ]][[ .varname ]][[ "exponentialwidth_uniform_5weeks" ]],
                           results.by.time[[ "6m" ]][[ .varname ]][[ "exponentialwidth_uniform_30weeks" ]],
                           "1m",
                           "6m"
                       );
                   .rv <- c( list( "exponentialwidth_uniform_1m5weeks_6m30weeks" = another.new.bounds.table ), .rv );
                   return( .rv );
               } else if( .varname == gold.standard.varname ) {
                   # one dimensional
                   .rv <- c( 
                             results.by.time[[ "1m" ]][[ .varname ]],
                             results.by.time[[ "6m" ]][[ .varname ]]
                       );
                     names( .rv ) <-
                         c( paste( names( results.by.time[[ "1m" ]][[ .varname ]] ), "1m", sep = "." ),
                           paste( names( results.by.time[[ "6m" ]][[ .varname ]] ), "6m", sep = "." ) );
                   return( .rv );
               } else {
                     .rv <- 
                       missing.column.safe.rbind(
                           results.by.time[[ "1m" ]][[ .varname ]],
                           results.by.time[[ "6m" ]][[ .varname ]],
                           "1m",
                           "6m"
                       );
                     return( .rv );
                 }
           } );
           names( results.1m.6m ) <-
             .vars;

           time.dependent.estimate.colname.roots <- c( "COB", "Infer" );
           for( .colname.root in time.dependent.estimate.colname.roots ) {
             if( length( grep( .colname.root, colnames( results.by.time[[ "1m" ]][[ "results.per.person" ]] ) ) ) > 0 ) {
               ## Add a new center-of-bounds result called "COB.uniform.1m5weeks.6m30weeks.time.est"
               new.estimates.table <-
                 rbind(
                   results.by.time[[ "1m" ]][[ "results.per.person" ]][ , paste( .colname.root, "uniform.5weeks.time.est", sep = "." ), drop = FALSE ],
                   results.by.time[[ "6m" ]][[ "results.per.person" ]][ , paste( .colname.root, "uniform.30weeks.time.est", sep = "." ), drop = FALSE ]
                 );
               rownames( new.estimates.table ) <-
                 c(
                   paste( rownames( results.by.time[[ "1m" ]][[ "results.per.person" ]] ), "1m", sep = "." ),
                   paste( rownames( results.by.time[[ "6m" ]][[ "results.per.person" ]] ), "6m", sep = "." )
                 );
               colnames( new.estimates.table ) <- paste( .colname.root, "uniform.1m5weeks.6m30weeks.time.est", sep = "." );
               results.1m.6m[[ "results.per.person" ]] <-
                   cbind( new.estimates.table, results.1m.6m[[ "results.per.person" ]] );
               
               ## Add a new center-of-bounds result called "COB.exponentialwidth.uniform.1m5weeks.6m30weeks.time.est"
               another.new.estimates.table <-
                 rbind(
                   results.by.time[[ "1m" ]][[ "results.per.person" ]][ , paste( .colname.root, "exponentialwidth.uniform.5weeks.time.est", sep = "." ), drop = FALSE ],
                   results.by.time[[ "6m" ]][[ "results.per.person" ]][ , paste( .colname.root, "exponentialwidth.uniform.30weeks.time.est", sep = "." ), drop = FALSE ]
                 );
               rownames( another.new.estimates.table ) <-
                 c(
                   paste( rownames( results.by.time[[ "1m" ]][[ "results.per.person" ]] ), "1m", sep = "." ),
                   paste( rownames( results.by.time[[ "6m" ]][[ "results.per.person" ]] ), "6m", sep = "." )
                 );
               colnames( another.new.estimates.table ) <- paste( .colname.root, "exponentialwidth.uniform.1m5weeks.6m30weeks.time.est", sep = "." );
               results.1m.6m[[ "results.per.person" ]] <-
                 cbind( another.new.estimates.table, results.1m.6m[[ "results.per.person" ]] );
               
             }
           } # End foreach .colname.root
           results.1m.6m <-
               c( results.1m.6m,
                 list( evaluated.results = evaluate.results.per.person.fn( results.1m.6m[[ "results.per.person" ]], results.1m.6m[[ gold.standard.varname ]], results.1m.6m[[ "results.covars.per.person.with.extra.cols" ]], the.time = "1m.6m", results.1m.6m[[ "bounds" ]] ) ) );
           return( c( list( "1m.6m" = results.1m.6m ), results.by.time ) );
       } ); # End foreach the.region
        names( results.by.region.and.time ) <- regions;

      # We now evaluate pooled results for every pair of regions, except nflg&rv217_v3.
      results.across.regions.by.time <-
        lapply( 1:( length( regions ) - 1), function( from.region.i ) {
          from.region <- regions[ from.region.i ];
        .rv.from.region.i <- 
        lapply( ( from.region.i + 1 ):length( regions ), function( to.region.j ) {
            to.region <- regions[ to.region.j ];
            if( ( from.region == "nflg" ) && ( to.region == "rv217_v3" ) ) {
                return( NA );
            }
          .times <- names( results.by.region.and.time[[ from.region ]] );
          .rv.from.region.i.to.region.j <- 
          lapply( .times, function ( the.time ) {
            ## TODO: REMOVE
            print( paste( "Pooling regions", from.region, "and", to.region, "at time", the.time ) );
              .vars <-
                  setdiff( names( results.by.region.and.time[[ from.region ]][[ the.time ]] ), "evaluated.results" );
              .rv.for.time <- 
              lapply( .vars, function ( .varname ) {
               #print( .varname );
               if( .varname == "bounds" ) {
                   .rv <- 
                   lapply( names( results.by.region.and.time[[ from.region ]][[ the.time ]][[ .varname ]] ), function( .bounds.type ) {
                       missing.column.safe.rbind(
                           results.by.region.and.time[[ from.region ]][[ the.time ]][[ .varname ]][[ .bounds.type ]],
                           results.by.region.and.time[[ to.region ]][[ the.time ]][[ .varname ]][[ .bounds.type ]],
                           from.region,
                           to.region
                       )
                   } );
                   names( .rv ) <-
                       names( results.by.region.and.time[[ from.region ]][[ the.time ]][[ .varname ]] );
                   return( .rv );
               } else if( .varname == gold.standard.varname ) {
                   # one dimensional
                   .rv <- c( 
                             results.by.region.and.time[[ from.region ]][[ the.time ]][[ .varname ]],
                             results.by.region.and.time[[ to.region ]][[ the.time ]][[ .varname ]]
                       );
                     names( .rv ) <-
                         c( paste( names( results.by.region.and.time[[ from.region ]][[ the.time ]][[ .varname ]] ), from.region, sep = "." ),
                           paste( names( results.by.region.and.time[[ to.region ]][[ the.time ]][[ .varname ]] ), to.region, sep = "." ) );
                   return( .rv );
               } else {
                     .rv <-
                         missing.column.safe.rbind(
                             results.by.region.and.time[[ from.region ]][[ the.time ]][[ .varname ]],
                             results.by.region.and.time[[ to.region ]][[ the.time ]][[ .varname ]],
                           from.region,
                           to.region
                       );
                     return( .rv );
                 }
           } );
           names( .rv.for.time ) <- .vars;





              
           # Add the evaluated.results:
           .evaluated.results <-
               evaluate.results.per.person.fn( .rv.for.time[[ "results.per.person" ]], .rv.for.time[[ gold.standard.varname ]], .rv.for.time[[ "results.covars.per.person.with.extra.cols" ]], the.time = the.time, .rv.for.time[[ "bounds" ]] );
           .rv.for.time <- c( .rv.for.time, 
                             list( evaluated.results = .evaluated.results ) )
           return( .rv.for.time );
          } );
        names( .rv.from.region.i.to.region.j ) <- .times;
          return( .rv.from.region.i.to.region.j );
        } );
        names( .rv.from.region.i ) <- regions[ ( from.region.i + 1 ):length( regions ) ];
          return( .rv.from.region.i );
      } );
      names( results.across.regions.by.time ) <- regions[ 1:( length( regions ) - 1 ) ];
      return( c( results.by.region.and.time, list( results.across.regions.by.time = results.across.regions.by.time ) ) );
} # getResultsByRegionAndTime ( .. )

