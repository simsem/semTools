### Sunthud Pornprasertmanit & Terrence D. Jorgensen
### Last updated: 11 April 2017
### Copy or save each aspect of the lavaan object into a clipboard or a file


#' Copy or save the result of \code{lavaan} or \code{FitDiff} objects into a
#' clipboard or a file
#'
#' Copy or save the result of \code{lavaan} or \code{\linkS4class{FitDiff}}
#' object into a clipboard or a file. From the clipboard, users may paste the
#' result into the Microsoft Excel or spreadsheet application to create a table
#' of the output.
#'
#'
#' @aliases clipboard saveFile
#' @param object The \code{lavaan} or \code{\linkS4class{FitDiff}} object
#' @param what The attributes of the \code{lavaan} object to be copied in the
#' clipboard. \code{"summary"} is to copy the screen provided from the
#' \code{summary} function. \code{"mifit"} is to copy the result from the
#' \code{\link{miPowerFit}} function. Other attributes listed in the
#' \code{inspect} method in the \link[lavaan]{lavaan-class} could also be used,
#' such as \code{"coef"}, \code{"se"}, \code{"fit"}, \code{"samp"}, and so on.
#' For the The \code{\linkS4class{FitDiff}} object, this argument is not active
#' yet.
#' @param file A file name used for saving the result
#' @param tableFormat If \code{TRUE}, save the result in the table format using
#'   tabs for seperation. Otherwise, save the result as the output screen
#'   printed in the R console.
#' @param fit.measures \code{character} vector specifying names of fit measures
#'   returned by \code{\link[lavaan]{fitMeasures}} to be copied/saved.  Only
#'   relevant if \code{object} is class \code{\linkS4class{FitDiff}}.
#' @param writeArgs \code{list} of additional arguments to be passed to
#'   \code{\link[utils]{write.table}}
#' @param \dots Additional argument listed in the \code{\link{miPowerFit}}
#'   function (for \code{lavaan} object only).
#' @return The resulting output will be saved into a clipboard or a file. If
#'   using the \code{clipboard} function, users may paste it in the other
#'   applications.
#' @author
#'  Sunthud Pornprasertmanit (\email{psunthud@@gmail.com})
#'
#'  Terrence D. Jorgensen (University of Amsterdam; \email{TJorgensen314@gmail.com})
#' @examples
#'
#' \dontrun{
#' library(lavaan)
#' HW.model <- ' visual  =~ x1 + c1*x2 + x3
#'               textual =~ x4 + c1*x5 + x6
#'                speed   =~ x7 + x8 + x9 '
#'
#' fit <- cfa(HW.model, data=HolzingerSwineford1939, group="school", meanstructure=TRUE)
#'
#' # Copy the summary of the lavaan object
#' clipboard(fit)
#'
#' # Copy the modification indices and the model fit from the miPowerFit function
#' clipboard(fit, "mifit")
#'
#' # Copy the parameter estimates
#' clipboard(fit, "coef")
#'
#' # Copy the standard errors
#' clipboard(fit, "se")
#'
#' # Copy the sample statistics
#' clipboard(fit, "samp")
#'
#' # Copy the fit measures
#' clipboard(fit, "fit")
#'
#' # Save the summary of the lavaan object
#' saveFile(fit, "out.txt")
#'
#' # Save the modification indices and the model fit from the miPowerFit function
#' saveFile(fit, "out.txt", "mifit")
#'
#' # Save the parameter estimates
#' saveFile(fit, "out.txt", "coef")
#'
#' # Save the standard errors
#' saveFile(fit, "out.txt", "se")
#'
#' # Save the sample statistics
#' saveFile(fit, "out.txt", "samp")
#'
#' # Save the fit measures
#' saveFile(fit, "out.txt", "fit")
#' }
#'
#' @export
clipboard <- function(object, what = "summary", ...) {
	if (.Platform$OS.type == "windows") {
		saveFile(object, file="clipboard-128", what=what, tableFormat=TRUE, ...)
		cat("File saved in the clipboard; please paste it in any program you wish.\n")
	} else {
		if (system("pbcopy", ignore.stderr = TRUE) == 0) {
			saveFile(object, file=pipe("pbcopy", "w"), what=what, tableFormat=TRUE, ...)
			cat("File saved in the clipboard; please paste it in any program you wish. If you cannot paste it, it is okay because this function works for some computers, which I still have no explanation currently. Please consider using the 'saveFile' function instead.\n")
		} else if (system("xclip -version", ignore.stderr = TRUE) == 0) {
			saveFile(object, file=pipe("xclip -i", "w") , what=what, tableFormat=TRUE, ...)
			cat("File saved in the xclip; please paste it in any program you wish. If you cannot paste it, it is okay because this function works for some computers, which I still have no explanation currently. Please consider using the 'saveFile' function instead.\n")
		} else {
			stop("For Mac users, the 'pbcopy' command in the shell file does not work. For linux users, this function depends on the 'xclip' application. Please install and run the xclip application before using this function in R (it does not guarantee to work though). Alternatively, use the 'saveFile' function to write the output into a file.")
		}
	}
}

#' @rdname clipboard
#' @export
saveFile <- function(object, file, what = "summary", tableFormat = FALSE,
                     fit.measures = "default", writeArgs = list(), ...) {
	# Check whether the object is in the lavaan class
	if (is(object, "lavaan")) {
		saveFileLavaan(object, file, what = what, tableFormat = tableFormat,
		               writeArgs = writeArgs, ...)
	} else if (is(object, "FitDiff")) {
		saveFileFitDiff(object, file, what = what, tableFormat = tableFormat,
		                fit.measures = fit.measures, writeArgs = writeArgs)
	} else {
		stop("The object must be in the `lavaan' output or the",
		     " output from the compareFit function.")
	}
}



## ----------------
## Hidden functions
## ----------------

#' @importFrom lavaan lavInspect
saveFileLavaan <- function(object, file, what = "summary", tableFormat = FALSE,
                           writeArgs = list(), ...) {
	if (length(what) > 1) message("only the first `what' option is used")
  # be case insensitive
  what <- tolower(what[1])

  writeArgs$file <- file
  if (is.null(writeArgs$sep)) writeArgs$sep <- "\t"
  if (is.null(writeArgs$quote)) writeArgs$quote <- FALSE

	if (what == "summary") {
		if (tableFormat) {
			copySummary(object, file = file, writeArgs = writeArgs)
		} else {
			write(paste(utils::capture.output(summary(object, rsquare = TRUE, fit = TRUE,
			                                          standardize = TRUE)),
			            collapse = "\n"), file = file)
		}
	} else if (what == "mifit") {
		if (tableFormat) {
		  writeArgs$x <- miPowerFit(object, ...)
		  if (is.null(writeArgs$row.names)) writeArgs$row.names <- FALSE
		  if (is.null(writeArgs$col.names)) writeArgs$col.names <- TRUE
		} else {
			write(paste(utils::capture.output(miPowerFit(object, ...)),
			            collapse = "\n"), file = file)
		}
	} else {
		target <- lavInspect(object, what=what)
		if (tableFormat) {
			if (is(target, "lavaan.data.frame") || is(target, "data.frame")) {
			  writeArgs$x <- target
			  if (is.null(writeArgs$row.names)) writeArgs$row.names <- FALSE
			  if (is.null(writeArgs$col.names)) writeArgs$col.names <- TRUE
			} else if (is(target, "list")) {
				if (is(target[[1]], "list")) {
					target <- lapply(target, listToDataFrame)
					target <- mapply(function(x, y) rbind(rep("", ncol(y)), c(x, rep("", ncol(y) - 1)), y),
					                 names(target), target, SIMPLIFY = FALSE)
					writeArgs$x <- do.call(rbind, target)
					if (is.null(writeArgs$row.names)) writeArgs$row.names <- FALSE
					if (is.null(writeArgs$col.names)) writeArgs$col.names <- FALSE
				} else {
					writeArgs$x <- listToDataFrame(target)
					if (is.null(writeArgs$row.names)) writeArgs$row.names <- FALSE
					if (is.null(writeArgs$col.names)) writeArgs$col.names <- FALSE
				}
			} else {
			  writeArgs$x <- target
			  if (is.null(writeArgs$row.names)) writeArgs$row.names <- TRUE
			  if (is.null(writeArgs$col.names)) writeArgs$col.names <- TRUE
			}
		} else {
			write(paste(utils::capture.output(target), collapse = "\n"), file = file)
		}
	}
  do.call("write.table", writeArgs)
}


# copySummary: copy the summary of the lavaan object into the clipboard and
# potentially be useful if users paste it into the Excel application
# object = lavaan object input
copySummary <- function(object, file, writeArgs = list()) {
	# Capture the output of the lavaan class
	outputText <- utils::capture.output(lavaan::summary(object, rsquare=TRUE, standardize=TRUE, fit.measure=TRUE))

	# Split the text by two spaces
	outputText <- strsplit(outputText, "  ")

	# Trim and delete the "" elements
	outputText <- lapply(outputText, function(x) x[x != ""])
	outputText <- lapply(outputText, trim)
	outputText <- lapply(outputText, function(x) x[x != ""])

	# Group the output into three sections: fit, parameter estimates, and r-squared
	cut1 <- grep("Estimate", outputText)[1]
	cut2 <- grep("R-Square", outputText)[1]
	set1 <- outputText[1:(cut1 - 1)]
	set2 <- outputText[cut1:(cut2 - 1)]
	set3 <- outputText[cut2:length(outputText)]

	# Assign the number of columns in the resulting data frame and check whether the output contains any labels
	numcol <- 7
	test <- set2[-grep("Estimate", set2)]
	test <- test[sapply(test, length) >= 2]
	if (any(sapply(test, function(x) is.na(suppressWarnings(as.numeric(x[2])))))) numcol <- numcol + 1

	# A function to parse the fit-measures output
	set1Parse <- function(x, numcol) {
		if (length(x) == 0) {
			return(rep("", numcol))
		} else if (length(x) == 1) {
			return(c(x, rep("", numcol - 1)))
		} else if ((length(x) >= 2) & (length(x) <= numcol)) {
			return(c(x[1], rep("", numcol - length(x)), x[2:length(x)]))
		} else {
			stop("Cannot parse text")
		}
	}
	set1 <- t(sapply(set1, set1Parse, numcol))

	# A function to parse the parameter-estimates output
	set2Parse <- function(x, numcol) {
		if (length(x) == 0) return(rep("", numcol))
		if (any(grepl("Estimate", x))) return(c(rep("", numcol-length(x)), x))
		if (length(x) == 1) {
			return(c(x, rep("", numcol-1)))
		} else {
			group1 <- x[1]
			group2 <- x[2:length(x)]
			if (is.na(suppressWarnings(as.numeric(x[2])))) {
				group1 <- x[1:2]
				group2 <- x[3:length(x)]
			} else if (numcol == 8) {
				group1 <- c(group1, "")
			}
			if (length(group2) == 1) {
				group2 <- c(group2, rep("", 6 - length(group2)))
			} else if (length(group2) == 4) {
				group2 <- c(group2, rep("", 6 - length(group2)))
			} else {
				group2 <- c(group2[1], rep("", 6 - length(group2)), group2[2:length(group2)])
			}
			return(c(group1, group2))
		}
	}
	set2 <- t(sapply(set2, set2Parse, numcol))

	# A function to parse the r-squared output
	set3Parse <- function(x, numcol) {
		if (length(x) == 0) {
			return(rep("", numcol))
		} else {
			return(c(x, rep("", numcol - length(x))))
		}
	}
	set3 <- t(sapply(set3, set3Parse, numcol))

	# Copy the output into the clipboard
	writeArgs$x <- rbind(set1, set2, set3)
	writeArgs$file <- file
	if (is.null(writeArgs$quote)) writeArgs$quote <- FALSE
	if (is.null(writeArgs$sep)) writeArgs$sep <- "\t"
	if (is.null(writeArgs$row.names)) writeArgs$row.names <- FALSE
	if (is.null(writeArgs$col.names)) writeArgs$col.names <- FALSE
	do.call("write.table", writeArgs)
}

# trim function from the R.oo package
trim <- function(object) {
	s <- sub("^[\t\n\f\r ]*", "", as.character(object));
	s <- sub("[\t\n\f\r ]*$", "", s);
	s;
}

# listToDataFrame: Change a list with multiple elements into a single data.frame
listToDataFrame <- function(object) {
	name <- names(object)

	# Count the maximum number of column (+1 is for the column for row name)
	numcol <- max(sapply(object, function(x) ifelse(is(x, "lavaan.matrix") || is(x, "lavaan.matrix.symmetric") || is(x, "matrix") || is(x, "data.frame"), return(ncol(x)), return(1)))) + 1

	# Change all objects in the list into a data.frame with the specified column
	target <- lapply(object, niceDataFrame, numcol)

	# Paste the name of each object into each data.frame
	target <- mapply(function(x, y) rbind(rep("", ncol(y)), c(x, rep("", ncol(y) - 1)), y), name, target, SIMPLIFY=FALSE)

	# Combine into a single data.frame
	target <- do.call(rbind, target)
	target[-1,]
}

# niceDataFrame: Change an object into a data.frame with a specified number of columns and the row and column names are included in the data.frame
niceDataFrame <- function(object, numcol) {
	temp <- NULL
	if (is(object, "lavaan.matrix.symmetric")) {
		# save only the lower diagonal of the symmetric matrix
		temp <- matrix("", nrow(object), ncol(object))
		for (i in 1:nrow(object)) {
			temp[i, 1:i] <- object[i, 1:i]
		}
	} else if (is(object, "data.frame") || is(object, "matrix") || is(object, "lavaan.matrix")) {
		# copy the matrix
		temp <- object
	} else if (is(object, "vector") || is(object, "lavaan.vector")) {
		# transform a vector into a matrix
		object <- as.matrix(object)
		temp <- object
	} else {
		stop("The 'niceDataFrame' function has a bug. Please contact the developer.")
	}

	# Transform into the result with a specified number of columns, excluding the row name
	result <- matrix("", nrow(temp), numcol - 1)

	# Parse the column names
	result[,1:ncol(temp)] <- temp
	firstRow <- colnames(object)
	ifelse(is.null(firstRow), firstRow <- rep("", ncol(result)), firstRow <- c(firstRow, rep("", numcol - length(firstRow) - 1)))

	# Parse the row names
	result <- rbind(firstRow, result)
	firstCol <- rownames(object)
	ifelse(is.null(firstCol), firstCol <- rep("", nrow(result)), firstCol <- c("", firstCol))
	result <- cbind(firstCol, result)
	dimnames(result) <- NULL
	result
}
