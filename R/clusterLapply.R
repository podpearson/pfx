# clusterLapply.R
# 
# Package: malariagen
# Author: Richard D Pearson, based on work by Richard and Gavin Band
# Email: richard.pearson@well.ox.ac.uk
#
###############################################################################


# Version of lapply that parallelises on the cluster
# This is done via the cluster and involves writing scripts to the script dirs.
# The getScriptDir argument should be a function which takes the jobSpec and returns
# an absolute directory for the scripts.
clusterLapply <- function(
		listOfObjects,
		functionName,
		runname 		                = "unnamed",
		getScriptDir                = file.path(getOption("malariagen.homeDir"), "clusterLapplyScripts"),
		SGEqueue                    = NULL,
		SGE_P                       = NULL,
#		SGEqueue                    = getOption("malariagen.SGEqueue"),
#		SGE_P                       = getOption("malariagen.SGE_P"),
#		SGE_P                       = NULL,
		secondsToWaitBetweenQsubs   = 1,
		masterScriptFilename        = NA,
		really.apply                = TRUE,
		wait.for.results            = TRUE,
    maxRAMinGB                  = 1,
		coresPerJob                 = getOption("malariagen.coresPerJob"),
		shouldJoinResultsOnly       = getOption("malariagen.shouldJoinResultsOnly"),
		Rexecutable                 = getOption("malariagen.Rexecutable"),
		initialRcommands            = "library(malariagen)\n",
    shouldUseLmn                = FALSE,
    shouldUseSmn                = FALSE,
    additionalSGEcommands       = NULL,
		...
) {
  if(!is.null(maxRAMinGB)) {
    if(is.null(SGEqueue)) {
      additionalSGEcommands <- paste(additionalSGEcommands, RAMlimitSGEcommands(maxRAMinGB))
    } else {
      additionalSGEcommands <- paste(additionalSGEcommands, RAMlimitSGEcommands(maxRAMinGB, SGEqueue=SGEqueue))
    }
  }
	if(!file.exists(getScriptDir)) {
		dir.create(getScriptDir, recursive=TRUE)
	}
	
	# If we have an empty list, just return now.
	if( length( listOfObjects ) == 0 ) {
		return( "clusterApply: nothing to do.  (List of objects was empty.)" )
	}
	
	# If we are passed something other than a function for getScriptDir, turn it into 
	# a function returning the script dir.
	if( typeof( getScriptDir ) != 'closure' ) {
		constantScriptDir = getScriptDir
		getScriptDir <- function( jobSpec ) {
			return( constantScriptDir )
		}
	}
	
	# If the list of objects is not named, number them "job1"...."jobN"
	if( length( listOfObjects ) > 0 && is.null( names( listOfObjects ))) {
		cat( 'clusterApply: list of objects is not named, so I will number them.\n' )
		names( listOfObjects ) <- sapply(
				seq( from = 1, to = length( listOfObjects )),
				function( i ) { paste( 'job', i, sep='' ) }
		)
	}
	
	# Construct functions which return the filenames
	getBashScriptFilename <- function( jobSpec, jobSpecName ) {
		file.path(
				getScriptDir( jobSpec ),
				paste( functionName, ".runall.sh", sep='' )
		)
	}
	
	getJobSpecFilename <- function( jobSpec, jobSpecName ) {
		file.path(
				getScriptDir( jobSpec ),
				paste(
						functionName,
						jobSpecName,
						'jobSpec',
						sep='.'
				)
		)
	}
	
	getRscriptFilename <- function( jobSpec, jobSpecName ) {
		file.path(
				getScriptDir( jobSpec ),
				paste(
						functionName,
						jobSpecName,
						'R',
						sep='.'
				)
		)
	}
	
	getQsubScriptFilename <- function( jobSpec, jobSpecName ) {
		file.path(
				getScriptDir( jobSpec ),
				paste(
						functionName,
						jobSpecName,
						'qsub',
						sep='.'
				)
		)
	}
	
	getResultsFilename <- function( jobSpec, jobSpecName ) {
		file.path(
				getScriptDir( jobSpec ),
				paste(
						functionName,
						jobSpecName,
						'results.rda',
						sep='.'
				)
		)
	}
	
# First sanity check script dir and remove any existing scripts...
	if(!shouldJoinResultsOnly) {
		for( i in seq(1, length( listOfObjects ))) {
			jobSpec = listOfObjects[[ i ]]
			jobSpecName = names( listOfObjects )[ i ]
			cat( jobSpecName, '\n', sep='' )
			if( substring( getScriptDir( jobSpec ), 1, 1 ) != "/" ) {
				stop( "clusterApply: I only accept absolute paths in getScriptDir( jobSpec ) as it could be a bit dodgy otherwise." )
			}
			
			if( !file.exists( getScriptDir( jobSpec ))) {
				cat( "clusterApply: the requested script dir \"", getScriptDir( jobSpec ), "\" does not exist.\n", sep='' )
				stop( "clusterApply: the requested script dir does not exist")
			}
			
			for( filename in c(
					getBashScriptFilename( jobSpec, jobSpecName ),
					getJobSpecFilename( jobSpec, jobSpecName ),
					getRscriptFilename( jobSpec, jobSpecName ),
					getQsubScriptFilename( jobSpec, jobSpecName ),
					getResultsFilename( jobSpec, jobSpecName )
			)) {
				if( file.exists( filename )) {
					cat( 'clusterApply: Removing existing script file \"', filename, '\"...\n', sep='' )
					file.remove( filename )
				}
			}
			
			logPath = file.path( getScriptDir( jobSpec ), "log" ) 
			if( !file.exists( logPath )) {
				dir.create( logPath )
			}
		}
		
		if( !is.na( masterScriptFilename ) && file.exists( masterScriptFilename )) {
			cat( 'clusterApply: removing master script \"', masterScriptFilename, '\"...\n', sep='' )
			file.remove( masterScriptFilename )
		}
		
		
		# Our job name will be a randomly chosen 8-character number.
		# However, we make sure it doesn't accidentally match any job that's currently running.
		jobName = paste( 'cA', sample( 1:(1000000 - 100000), size=1 ) + 99999, sep='' )
		while( length( getClusterJobsWithName( jobName )) > 0 ) {
			jobName = paste( 'cA', sample( 1:(1000000 - 100000), size=1 ) + 99999, sep='' )
		}
		
		# Now construct the new scripts.
		bashScriptFilenames <- c()
		
		for( i in seq(1, length( listOfObjects ))) {
			jobSpec = listOfObjects[[ i ]]
			jobSpecName = names( listOfObjects )[ i ]
			
			bashScriptFilename <- getBashScriptFilename( jobSpec, jobSpecName )
			jobSpecFilename <- getJobSpecFilename( jobSpec, jobSpecName )
			RscriptFilename <- getRscriptFilename( jobSpec, jobSpecName )
			qsubScriptFilename <- getQsubScriptFilename( jobSpec, jobSpecName )
			resultsFilename <- getResultsFilename( jobSpec, jobSpecName )
			
			# Write job spec...
			cat( '--> Writing scripts for \"', jobSpecName, '\"...\n', sep='' ) 
			save( jobSpec, file = jobSpecFilename )
			
			# Write R script...
			cat( "# ", functionName, "( ", jobSpecName, " )\n", sep="", file = RscriptFilename )
			cat( "# This file created automatically by function clusterLapply from package t2d\n", file=RscriptFilename, append=TRUE )
			cat( "# File created ", date(), "\n", file=RscriptFilename, append=TRUE )
			cat( "# File last modified ", date(), "\n\n", file=RscriptFilename, append=TRUE )
			cat(initialRcommands, file=RscriptFilename, append=TRUE )
			if(length(as.character(match.call(expand.dots=FALSE))) %in% grep("^list", as.character(match.call(expand.dots=FALSE)))) { # i.e. ... arguments have been supplied
				extraArguments <- as.character(match.call(expand.dots=FALSE))[length(as.character(match.call(expand.dots=FALSE)))]
				extraArguments <- sub("^list\\(", "", extraArguments)
				extraArguments <- sub("\\)$", "", extraArguments)
				extraArguments <- strsplit(extraArguments, ", ")[[1]]
				extraArguments <- sub(" = \\.\\..*$", "", extraArguments)
				for(j in 1:(length(extraArguments))) {
					if(class(eval(paste("..", j, sep=""))) == "character") {
						if(length(eval(parse(text=paste("..", j, sep="")))) == 1) {
							extraArguments[j] <- paste(extraArguments[j], "=\"", eval(parse(text=paste("..", j, sep=""))), "\"", sep="")
						} else {
							tempArguments <- paste(eval(parse(text=paste("..", j, sep=""))), collapse="\",\"")
							extraArguments[j] <- paste(extraArguments[j], "=c(\"", tempArguments, "\")", sep="")
						}
					} else {
						if(length(eval(parse(text=paste("..", j, sep="")))) == 1) {
							extraArguments[j] <- paste(extraArguments[j], "=", eval(parse(text=paste("..", j, sep=""))), sep="")
						} else {
							tempArguments <- paste(eval(parse(text=paste("..", j, sep=""))), collapse=",")
							extraArguments[j] <- paste(extraArguments[j], "=c(", tempArguments, ")", sep="")
						}
					}
				}
				extraArguments <- paste(extraArguments, collapse=", ")

				cat( "results <- clusterApplyFunctionProxy( ", functionName, ", jobSpecFilename = \"", jobSpecFilename, '\", ', extraArguments, ')\n', sep="", file=RscriptFilename, append=TRUE )
			} else {
				cat( "results <- clusterApplyFunctionProxy( ", functionName, ", jobSpecFilename = \"", jobSpecFilename, '\" )\n', sep="", file=RscriptFilename, append=TRUE )
			}
			cat( "save(results, file=\"", resultsFilename, "\")\n", sep="", file=RscriptFilename, append=TRUE )
			
			# Write bash script...
			cat( "#!/bin/bash\n", file = qsubScriptFilename )
			logPath = file.path( getScriptDir( jobSpec ), "log ")
			cat( "#$ -N ",  jobName, "\n",          sep="",   file = qsubScriptFilename, append = TRUE )
      if(!is.null(SGE_P) && !is.na(SGE_P)) {
  			cat( "#$ -P",  SGE_P, "\n",        sep=" ",  file = qsubScriptFilename, append = TRUE )
      }
      if(!is.null(SGEqueue) && !is.na(SGEqueue)) {
  			cat( "#$ -pe",  SGEqueue, coresPerJob, "\n",        sep=" ",  file = qsubScriptFilename, append = TRUE )
      }
      if(!is.null(additionalSGEcommands)) {
  			cat( "#$", additionalSGEcommands, "\n",        sep=" ",  file = qsubScriptFilename, append = TRUE )
      }
      if(shouldUseLmn) {
  			cat( "#$ -l lmn\n",                     file = qsubScriptFilename, append = TRUE )
      }
      if(shouldUseSmn) {
  			cat( "#$ -l smn\n",                     file = qsubScriptFilename, append = TRUE )
      }
			cat( "#$ -e",   logPath, "\n",            file = qsubScriptFilename, append = TRUE )
			cat( "#$ -o",   logPath, "\n",            file = qsubScriptFilename, append = TRUE )
			cat( "#$ -V\n",            file = qsubScriptFilename, append = TRUE )
			cat( "echo \"Running job ", jobName, "...\"\n",   sep="",   file = qsubScriptFilename, append = TRUE )
			cat(
					Rexecutable,
					" CMD BATCH --vanilla \"",
					RscriptFilename,
					"\" \"",
					paste( RscriptFilename, "out", sep="" ),
					"\"\n",
					sep = "",
					file = qsubScriptFilename,
					append = TRUE
			)
			cat( "echo \"Finished job ", jobName, ".\"\n",  sep="",   file = qsubScriptFilename, append = TRUE )
			
			# Invoke this from the bash script for this jobSpec

		        cat( "echo \"Submitting job ", qsubScriptFilename, "...\"\n", sep="", file = bashScriptFilename, append = TRUE )
			# cat( "qsub /users/jeffchen/tools/scripts/checkIB.qsub \n", sep="", file = bashScriptFilename, append = TRUE )
			cat( "qsub \"", qsubScriptFilename, "\"\n",           sep="", file = bashScriptFilename, append = TRUE )
			if(secondsToWaitBetweenQsubs > 0) {
				cat("sleep ", secondsToWaitBetweenQsubs, "\n",        sep="", file = bashScriptFilename, append = TRUE)
			}
			
			bashScriptFilenames <- c( bashScriptFilenames, bashScriptFilename )
			
			cat(
					'--> Wrote \"',
					qsubScriptFilename,
					'\".\n',
					sep=''
			)
		}
		
		sort( bashScriptFilenames )
		bashScriptFilenames = unique( bashScriptFilenames )
		
		cat( "clusterApply: I wrote the following bash scripts:\n" )
		for( filename in bashScriptFilenames ) {
			cat( "--> \"", filename, "\".\n", sep='' )
		}
		
		# Construct the master script if there is one...
		if( !is.na( masterScriptFilename )) {
			cat( "clusterApply: writing master script...\n" )
			for( bashScriptFilename in bashScriptFilenames ) {
				cat( "/bin/bash \"", bashScriptFilename, "\"\n", sep='', file = masterScriptFilename, append = TRUE )
			}
		}
		
		masterScriptFilename <- bashScriptFilename
		
		# Run the master script if requested...
		if( really.apply && !is.na( masterScriptFilename )) {
			cat( 'clusterApply: invoking \"', masterScriptFilename, '\"...\n', sep='' )
			output = system(
					paste(".", masterScriptFilename),
					intern = TRUE
			)
			
			# Wait for all jobs to finish if requested...
			if( wait.for.results ) {
				cat( 'clusterApply: waiting for jobs to complete...\n')
				Sys.sleep(20)
				L = length( getClusterJobsWithName( jobName ))
				while( L > 0 ) {
					Sys.sleep(20)
					L = length( getClusterJobsWithName( jobName ))
					cat( '...', L, ' jobs left', sep='' )
				}
				cat( '\nclusterApply: all jobs finished.\n' )
			
  			output <- list()
  			for( i in seq(1, length( listOfObjects ))) {
  				jobSpec = listOfObjects[[ i ]]
  				jobSpecName = names( listOfObjects )[ i ]
  				
  				resultsFilename <- getResultsFilename( jobSpec, jobSpecName )
  				
  				load(resultsFilename)
  				output <- c(output, list(results))
  			}
  			if(!is.null( names( listOfObjects ))) {
  				names(output) <- names(listOfObjects)
  			}
			} else {
        output <- "Jobs left running"
      }
			
		}
		else if( !is.na( masterScriptFilename )) {
			cat( 'To apply these jobs, invoke \"', masterScriptFilename, '\"...\n', sep='' )
			output = masterScriptFilename
		}
		else {
			cat( 'All job scripts written (no master script).\n')
			output = 'Success'
		}
	} else { # shouldJoinResultsOnly
		output <- list()
		for( i in seq(1, length( listOfObjects ))) {
			jobSpec = listOfObjects[[ i ]]
			jobSpecName = names( listOfObjects )[ i ]
			
			resultsFilename <- getResultsFilename( jobSpec, jobSpecName )
			
			load(resultsFilename)
			output <- c(output, list(results))
		}
		if(!is.null( names( listOfObjects ))) {
			names(output) <- names(listOfObjects)
		}
	}
	
	cat( 'clusterLapply() complete.\n' )
	return( output )
}

