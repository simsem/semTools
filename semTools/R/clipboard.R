### Title: Copy or save each aspect of the lavaan object into a clipboard or a file
### Author: Sunthud Pornprasertmanit
### Last updated: 14 October 2016
### Description: Copy or print each aspect of the lavaan object into a clipboard or a file

# Clipboard: copy each aspect of the lavaan object into a clipboard; this function will be compatible with lavaan::lavInspect
clipboard <- function(object, what="summary", ...) {
	if(.Platform$OS.type == "windows") {
		saveFile(object, file="clipboard-128", what=what, tableFormat=TRUE, ...)
		cat("File saved in the clipboard; please paste it in any program you wish.\n")
	} else {
		if(system("pbcopy", ignore.stderr = TRUE) == 0) {
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

# saveFile: save each aspect of the lavaan object into a file; this function will be compatible with lavaan::lavInspect
saveFile <- function(object, file, what="summary", tableFormat=FALSE, ...) {
	# Check whether the object is in the lavaan class
	if(is(object, "lavaan")) {
		saveFileLavaan(object, file, what=what, tableFormat=tableFormat, ...)
	} else if(is(object, "FitDiff")) {
		saveFileFitDiff(object, file, what=what, tableFormat=tableFormat)
	} else {
		stop("The object must be in the `lavaan' output or the output from the compareFit function.")
	}
}

saveFileLavaan <- function(object, file, what="summary", tableFormat=FALSE, ...) {
	if(length(what) > 1) {
        stop("`what' arguments contains multiple arguments; only one is allowed")
    } 
    # be case insensitive
    what <- tolower(what)
	
	if(what == "summary") {
		if(tableFormat) {
			copySummary(object, file=file)
		} else {
			write(paste(capture.output(summary(object, rsquare=TRUE, standardize=TRUE, fit.measure=TRUE)), collapse="\n"), file=file)
		}
	} else if (what == "mifit") {
		if(tableFormat) {
			write.table(miPowerFit(object, ...), file=file, sep="\t", row.names=FALSE, col.names=TRUE)
		} else {
			write(paste(capture.output(miPowerFit(object, ...)), collapse="\n"), file=file)
		}
	} else {
		target <- lavaan::lavInspect(object, what=what)
		if(tableFormat) {
			if(is(target, "lavaan.data.frame") || is(target, "data.frame")) {
				utils::write.table(target, file=file, sep="\t", row.names=FALSE, col.names=TRUE)
			} else if (is(target, "list")) {
				if(is(target[[1]], "list")) {
					target <- lapply(target, listToDataFrame)
					target <- mapply(function(x, y) rbind(rep("", ncol(y)), c(x, rep("", ncol(y) - 1)), y), names(target), target, SIMPLIFY=FALSE)
					target <- do.call(rbind, target)
					utils::write.table(target[-1,], file=file, sep="\t", row.names=FALSE, col.names=FALSE)
				} else {
					target <- listToDataFrame(target)
					utils::write.table(target, file=file, sep="\t", row.names=FALSE, col.names=FALSE)
				}
			} else {
				utils::write.table(target, file=file, sep="\t", row.names=TRUE, col.names=TRUE)
			}
		} else {
			write(paste(utils::capture.output(target), collapse="\n"), file=file)
		}
	}
}


# copySummary: copy the summary of the lavaan object into the clipboard and potentially be useful if users paste it into the excel application
# object = lavaan object input
copySummary <- function(object, file) {
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
	test <- test[sapply(test, length) >=2]
	if(any(sapply(test, function(x) is.na(suppressWarnings(as.numeric(x[2])))))) numcol <- numcol + 1

	# A function to parse the fit-measures output
	set1Parse <- function(x, numcol) {
		if(length(x) == 0) { 
			return(rep("", numcol))
		} else if(length(x) == 1) {
			return(c(x, rep("", numcol - 1)))
		} else if((length(x) >= 2) & (length(x) <= numcol)) {
			return(c(x[1], rep("", numcol - length(x)), x[2:length(x)]))
		} else {
			stop("Cannot parse text")
		}
	}
	set1 <- t(sapply(set1, set1Parse, numcol))

	# A function to parse the parameter-estimates output
	set2Parse <- function(x, numcol) {
		if(length(x) == 0) return(rep("", numcol))
		if(any(grepl("Estimate", x))) return(c(rep("", numcol-length(x)), x))
		if(length(x) == 1) {
			return(c(x, rep("", numcol-1)))
		} else {
			group1 <- x[1]
			group2 <- x[2:length(x)]
			if(is.na(suppressWarnings(as.numeric(x[2])))) {
				group1 <- x[1:2]
				group2 <- x[3:length(x)]
			} else if (numcol == 8) {
				group1 <- c(group1, "")
			}
			if(length(group2) == 1) {
				group2 <- c(group2, rep("", 6 - length(group2)))
			} else if(length(group2) == 4) {
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
		if(length(x) == 0) { 
			return(rep("", numcol))
		} else {
			return(c(x, rep("", numcol - length(x))))
		} 
	}
	set3 <- t(sapply(set3, set3Parse, numcol))

	# Copy the output into the clipboard
	utils::write.table(rbind(set1, set2, set3), file=file, sep="\t", row.names=FALSE, col.names=FALSE)
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
	if(is(object, "lavaan.matrix.symmetric")) {
		# save only the lower diagonal of the symmetric matrix
		temp <- matrix("", nrow(object), ncol(object))
		for(i in 1:nrow(object)) {
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
