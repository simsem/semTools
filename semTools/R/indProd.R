## Title: Orthogonalize data for 2-way and 3-way interaction in SEM
## Author: Sunthud Pornprasertmanit and Alexander M. Schoemann
## Description: Orthogonalize data for 2-way and 3-way interaction in SEM
##----------------------------------------------------------------------------##

# indProd: Make a product of indicators using mean centering, double-mean centering, or residual centering

indProd <- function(data, var1, var2, var3=NULL, match = TRUE, meanC = TRUE, residualC = FALSE, doubleMC = TRUE, namesProd = NULL) {
	# Get all variable names
    if (all(is.numeric(var1))) 
        var1 <- colnames(data)[var1]
    if (all(is.numeric(var2))) 
        var2 <- colnames(data)[var2]
	if (!is.null(var3) && all(is.numeric(var3))) var3 <- colnames(data)[var3]
    dat1 <- data[, var1]
    dat2 <- data[, var2]
    dat3 <- NULL
	if (!is.null(var3)) dat3 <- data[, var3]
	
	# Mean centering on the original indicators
    if (meanC) {
        dat1 <- scale(dat1, scale = FALSE)
        dat2 <- scale(dat2, scale = FALSE)
		if (!is.null(dat3)) dat3 <- scale(dat3, scale = FALSE)
    }
    if (match) {
		# Check whether the number of variables are equal across variable sets
        if (length(var1) != length(var2)) 
            stop("If the match-paired approach is used, the number of variables in all sets must be equal.")
		if (!is.null(var3) && (length(var1) != length(var3)))
            stop("If the match-paired approach is used, the number of variables in all three sets must be equal.")
		datProd <- NULL
		if (is.null(var3)) {
			# Two-way interaction
			datProd <- dat1 * dat2
			if (residualC) {
				notmissing <- which(!apply(datProd, 1, function(x) any(is.na(x))))
				colnames(datProd) <- paste("interactionProduct", 1:ncol(datProd), sep = "")
				# Write the expression for linear model and residualize the products
				temp <- data.frame(datProd, dat1, dat2)
				express <- paste("cbind(", paste(colnames(datProd), collapse = ", "), ") ~ ", paste(c(colnames(dat1), colnames(dat2)), collapse = " + "), 
					sep = "")
				datProd[notmissing,] <- lm(express, data = temp)$residuals
			}
		} else {
			# Three-way interaction
			datProd2way <- cbind(dat1 * dat2, dat1 * dat3, dat2 * dat3)
			datProd3way <- dat1 * dat2 * dat3
			if (residualC) {
				notmissing2way <- which(!apply(datProd2way, 1, function(x) any(is.na(x))))
				colnames(datProd2way) <- paste("interaction2Product", 1:ncol(datProd2way), sep = "")
				
				# Write the expression for linear model and residualize the two-way products
				temp2 <- data.frame(datProd2way, dat1, dat2, dat3)
				express2 <- paste("cbind(", paste(colnames(datProd2way), collapse = ", "), ") ~ ", paste(c(colnames(dat1), colnames(dat2), colnames(dat3)), collapse = " + "), sep = "")
				datProd2way[notmissing2way,] <- lm(express2, data = temp2)$residuals
				
				# Making all possible products to residualize the 3-way interaction
				datProd2wayFull <- matrix(0, nrow(data), 1)
				for (i in 1:length(var1)) datProd2wayFull <- data.frame(datProd2wayFull, matrix(rep(dat1[, i], length(var2)), ncol = length(var2)) * dat2)
				for (i in 1:length(var1)) datProd2wayFull <- data.frame(datProd2wayFull, matrix(rep(dat1[, i], length(var3)), ncol = length(var3)) * dat3)
				for (i in 1:length(var2)) datProd2wayFull <- data.frame(datProd2wayFull, matrix(rep(dat2[, i], length(var3)), ncol = length(var3)) * dat3)
				datProd2wayFull <- datProd2wayFull[, -1]
				colnames(datProd2wayFull) <- paste("interaction2Product", 1:ncol(datProd2wayFull), sep = "")
				
				notmissing3way <- which(!apply(datProd3way, 1, function(x) any(is.na(x))))
				colnames(datProd3way) <- paste("interaction3Product", 1:ncol(datProd3way), sep = "")
				# Write the expression for linear model and residualize the three-way products
				temp3 <- data.frame(datProd3way, dat1, dat2, dat3, datProd2wayFull)
				express3 <- paste("cbind(", paste(colnames(datProd3way), collapse = ", "), ") ~ ", paste(c(colnames(dat1), colnames(dat2), colnames(dat3), colnames(datProd2wayFull)), collapse = " + "), sep = "")
				datProd3way[notmissing3way,] <- lm(express3, data = temp3)$residuals
			}			
			datProd <- cbind(datProd2way, datProd3way)
		}
		# Mean-centering the final product
        if (doubleMC) 
            datProd <- scale(datProd, scale = FALSE)
			
		# Rename the obtained product terms
        if (is.null(namesProd)) {
			if (is.null(var3)) {
				colnames(datProd) <- paste(var1, var2, sep = ".")
			} else {
				colnames(datProd) <- c(paste(var1, var2, sep = "."), paste(var1, var3, sep = "."), paste(var2, var3, sep = "."), paste(var1, var2, var3, sep = "."))
			}
        } else {
            colnames(datProd) <- namesProd
        }
    } else {
		datProd <- NULL
		if (is.null(var3)) {
			# Create all possible combinations of the products of indicators
			datProd <- matrix(0, nrow(data), 1)
			for (i in 1:length(var1)) datProd <- data.frame(datProd, matrix(rep(dat1[, i], length(var2)), ncol = length(var2)) * dat2)
			datProd <- datProd[, -1]
			if (residualC) {
				notmissing <- which(!apply(datProd, 1, function(x) any(is.na(x))))
				colnames(datProd) <- paste("interactionProduct", 1:ncol(datProd), sep = "")
				# Write the expression for linear model and residualize the two-way products
				temp <- data.frame(datProd, dat1, dat2)
				express <- paste("cbind(", paste(colnames(datProd), collapse = ", "), ") ~ ", paste(c(colnames(dat1), colnames(dat2)), collapse = " + "), 
					sep = "")
				datProd[notmissing,] <- lm(express, data = temp)$residuals
			}
		} else {
			# Create all possible combinations of the products of indicators
			datProd2way <- matrix(0, nrow(data), 1)
			for (i in 1:length(var1)) datProd2way <- data.frame(datProd2way, matrix(rep(dat1[, i], length(var2)), ncol = length(var2)) * dat2)
			for (i in 1:length(var1)) datProd2way <- data.frame(datProd2way, matrix(rep(dat1[, i], length(var3)), ncol = length(var3)) * dat3)
			for (i in 1:length(var2)) datProd2way <- data.frame(datProd2way, matrix(rep(dat2[, i], length(var3)), ncol = length(var3)) * dat3)
			datProd3way <- matrix(0, nrow(data), 1)
			for (i in 1:length(var1)) {
				for(j in 1:length(var2)) {
					datProd3way <- data.frame(datProd3way, matrix(rep(dat1[, i], length(var3)), ncol = length(var3)) * matrix(rep(dat2[, j], length(var3)), ncol = length(var3)) * dat3)
				}
			}
			datProd2way <- datProd2way[, -1]
			datProd3way <- datProd3way[, -1]
			if (residualC) {
				notmissing2way <- which(!apply(datProd2way, 1, function(x) any(is.na(x))))
				colnames(datProd2way) <- paste("interaction2Product", 1:ncol(datProd2way), sep = "")
				# Write the expression for linear model and residualize the two-way products
				temp2 <- data.frame(datProd2way, dat1, dat2, dat3)
				express2 <- paste("cbind(", paste(colnames(datProd2way), collapse = ", "), ") ~ ", paste(c(colnames(dat1), colnames(dat2), colnames(dat3)), collapse = " + "), sep = "")
				datProd2way[notmissing2way,] <- lm(express2, data = temp2)$residuals
				notmissing3way <- which(!apply(datProd3way, 1, function(x) any(is.na(x))))
				colnames(datProd3way) <- paste("interaction3Product", 1:ncol(datProd3way), sep = "")
				# Write the expression for linear model and residualize the three-way products
				temp3 <- data.frame(datProd3way, dat1, dat2, dat3, datProd2way)
				express3 <- paste("cbind(", paste(colnames(datProd3way), collapse = ", "), ") ~ ", paste(c(colnames(dat1), colnames(dat2), colnames(dat3), colnames(datProd2way)), collapse = " + "), sep = "")
				datProd3way[notmissing3way,] <- lm(express3, data = temp3)$residuals
			}
			datProd <- cbind(datProd2way, datProd3way)
		}
		# Double-mean centering
        if (doubleMC) 
            datProd <- scale(datProd, scale = FALSE)
			
		# Name the resulting product terms
        if (is.null(namesProd)) {
			temp <- NULL
			if (is.null(var3)) {
				for (i in 1:length(var1)) temp <- c(temp, paste(var1[i], var2, sep = "."))
			} else {
				for (i in 1:length(var1)) temp <- c(temp, paste(var1[i], var2, sep = "."))
				for (i in 1:length(var1)) temp <- c(temp, paste(var1[i], var3, sep = "."))
				for (i in 1:length(var2)) temp <- c(temp, paste(var2[i], var3, sep = "."))
				for (i in 1:length(var1)) {
					for(j in 1:length(var2)) {
						temp <- c(temp, paste(var1[i], var2[j], var3, sep = "."))
					}
				}
			}
			colnames(datProd) <- temp
        } else {
            colnames(datProd) <- namesProd
        }
    }
	# Bind the products back to the original data
    data <- data.frame(data, datProd)
    return(data)
} 

# orthogonalize: the shortcut for residual centering
orthogonalize <- function(data, var1, var2, var3=NULL, match=TRUE, namesProd=NULL) {
	indProd(data=data, var1=var1, var2=var2, var3=var3, match=match, meanC=FALSE, residualC=TRUE, doubleMC=FALSE, namesProd=namesProd)
}
