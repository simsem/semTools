### Corbin Quick
### Last updated: 4 April 2017


#' Randomly Split a Data Set into Halves
#' 
#' This function randomly splits a data set into two halves, and saves the
#' resulting data sets to the same folder as the original.
#' 
#' This function randomly orders the rows of a data set, divides the data set
#' into two halves, and saves the halves to the same folder as the original
#' data set, preserving the original formatting. Data set type (*.csv or *.dat)
#' and formatting (headers) are automatically detected, and output data sets
#' will preserve input type and formatting unless specified otherwise. Input
#' can be in the form of a file path (*.dat or *.csv), or an R object (matrix or
#' dataframe). If input is an R object and path is default, output data sets
#' will be returned as a list object.
#' 
#' 
#' @importFrom stats runif
#' 
#' @param dataset The original data set to be divided. Can be a file path to a
#' *.csv or *.dat file (headers will automatically be detected) or an R object
#' (matrix or dataframe). (Windows users: file path must be specified using
#' FORWARD SLASHES (\code{/}) ONLY.)
#' @param path File path to folder for output data sets. NOT REQUIRED if
#' dataset is a filename. Specify ONLY if dataset is an R object, or desired
#' output folder is not that of original data set. If path is specified as
#' "object", output data sets will be returned as a list, and not saved to hard
#' drive.
#' @param div Number of output data sets. NOT REQUIRED if default, 2 halves.
#' @param type Output file format ("dat" or "csv"). NOT REQUIRED unless desired
#' output formatting differs from that of input, or dataset is an R object and
#' csv formatting is desired.
#' @param name Output file name. NOT REQUIRED unless desired output name
#' differs from that of input, or input dataset is an R object. (If input is an
#' R object and name is not specified, name will be "splitSample".)
#' @return If \code{path = "object"}, \code{list} of output data sets. 
#'   Otherwise, output will saved to hard drive in the same format as input.
#' @author Corbin Quick (University of Michigan; \email{corbinq@@umich.edu})
#' @examples
#' 
#' #### Input is .dat file
#' #splitSample("C:/Users/Default/Desktop/MYDATA.dat")
#' #### Output saved to "C:/Users/Default/Desktop/" in .dat format
#' #### Names are "MYDATA_s1.dat" and "MYDATA_s2.dat"
#' 
#' #### Input is R object
#' ## Split C02 dataset from the datasets package
#' library(datasets)
#' splitMyData <- splitSample(CO2, path = "object")
#' summary(splitMyData[[1]])
#' summary(splitMyData[[2]])
#' #### Output object splitMyData becomes list of output data sets
#' 
#' #### Input is .dat file in "C:/" folder
#' #splitSample("C:/testdata.dat", path = "C:/Users/Default/Desktop/", type = "csv")
#' #### Output saved to "C:/Users/Default/Desktop/" in *.csv format
#' #### Names are "testdata_s1.csv" and "testdata_s2.csv"
#' 
#' #### Input is R object
#' #splitSample(myData, path = "C:/Users/Default/Desktop/", name = "splitdata")
#' #### Output saved to "C:/Users/Default/Desktop/" in *.dat format
#' #### Names are "splitdata_s1.dat" and "splitdata_s2.dat"
#' 
#' @export
splitSample <- function(dataset, path = "default", div = 2,
                        type = "default", name = "splitSample") {
  
  type1 <- type
  hea = FALSE
  file <- dataset
  
  if (is.character(file)) {
    temp <- strsplit(file, "/", fixed = TRUE)
    if (path == "default") {
      path <- paste(temp[[1]][1:(length(temp[[1]]) - 1)], "/",
                    sep = "", collapse = "")
    }
    fileN <- temp[[1]][length(temp[[1]])]
    temp <- strsplit(fileN, ".", fixed = TRUE)
    type <- temp[[1]][2]
    name <- temp[[1]][1]
    if (type == "dat") {
      if (is.numeric(as.matrix(utils::read.table(file, nrows = 1))) == FALSE) {
        data <- as.matrix(utils::read.table(file, header = TRUE))
        hea = TRUE
      } else {
        data <- as.matrix(utils::read.table(file))
      }
    }
    if (type == "csv") {
      if (is.numeric(as.matrix(utils::read.table(file, nrows = 1))) == FALSE) {
        data <- as.matrix(utils::read.csv(file, header = TRUE))
        hea = TRUE
      } else {
        data <- as.matrix(utils::read.csv(file))
      }
    }
  } else {
    if (is.matrix(file) | is.data.frame(file)) {
      data <- as.matrix(file)
    } else {
      stop("Provide data in *.dat or *.csv format")
    }
  }
  
  if (type1 != "default") {
    type <- type1
  }
  
  if (is.character(colnames(data))) {
    hea = TRUE
  }
  
  random <- runif(nrow(data), 1, nrow(data))
  data <- cbind(random, data)
  data <- data[order(random), ]
  data <- data[, 2:ncol(data)]
  
  size <- split((1:nrow(data)), cut((1:nrow(data)), div, labels = FALSE))
  size <- as.matrix(as.data.frame(lapply(size, length)))
  dataL <- list()
  dataL[[1]] <- data[1:size[1, 1], ]
  for (i in 2:div) {
    size[1, i] <- size[1, (i - 1)] + size[1, i]
    dataL[[i]] <- data[(size[1, (i - 1)] + 1):size[1, i], ]
  }
  
  if (path == "default") {
    return(dataL)
  } else {
    if (path == "object") {
      return(dataL)
    } else {
      for (i in 1:div) {
        if (type == "dat") {
          utils::write.table(dataL[[i]],
                             paste(path, name, "_s", i, ".dat", sep = ""),
                             sep = "  ", row.names = FALSE, col.names = hea)
        }
        if (type == "csv") {
          utils::write.table(dataL[[i]],
                             paste(path, name, "_s", i, ".csv", sep = ""),
                             sep = ",", row.names = FALSE, col.names = hea)
        }
        if (type == "default") {
          utils::write.table(dataL[[i]],
                             paste(path, name, "_s", i, ".dat", sep = ""),
                             sep = "  ", row.names = FALSE, col.names = hea)
        }
      }
    }
  }
}


