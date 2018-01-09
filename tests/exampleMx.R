# library(tools)
# dirMan <- "C:/Users/student/Dropbox/semTools/semTools/man/runMI.Rd"
# showNonASCIIfile(dirMan)

sourceDir <- function(path, trace = TRUE, ...) {
     for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
		if(nm != "AllClass.R" & nm != "AllGenerics.R") {
        if(trace) cat(nm,":") 
        source(file.path(path, nm), ...)
        if(trace) cat("\n")
		}
     }
}

sourceDirData <- function(path, trace = TRUE) {
     for (nm in list.files(path, pattern = "\\.[Rr]da$")) {
        if(trace) cat(nm,":") 
        load(paste0(path, nm), envir = .GlobalEnv)
        if(trace) cat("\n")
	}
}

library(lavaan)
library(OpenMx)

dir <- "C:/Users/Sunthud/Dropbox/semTools/semTools/R/"
sourceDir(dir)

# dir2 <- "C:/Users/sunthud/Desktop/multcomp/R"
# sourceDir(dir2)

dirData <- "C:/Users/Sunthud/Dropbox/semTools/semTools/data/"
sourceDirData(dirData)

################# nullMX

data(demoOneFactor)
nullModel <- nullMx(demoOneFactor)

################# saturatedMx

data(demoOneFactor)
satModel <- saturateMx(demoOneFactor)

manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("One Factor", 
    type="RAM",
    manifestVars=manifests, 
    latentVars=latents,
    mxPath(from=latents, to=manifests),
	mxPath(from='one', to=manifests),
    mxPath(from=manifests, arrows=2),
    mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
    mxData(observed=demoOneFactor, type="raw")
)
factorFit <- mxRun(factorModel)
summary(factorFit, SaturatedLikelihood=satModel, IndependenceLikelihood=nullModel)

################# StandardizeMx

data(myFADataRaw)
myFADataRaw <- myFADataRaw[,c("x1","x2","x3","x4","x5","x6")]
oneFactorModel <- mxModel("Common Factor Model Path Specification", 
	type="RAM",
	mxData(
		observed=myFADataRaw, 
		type="raw"
	),
	manifestVars=c("x1","x2","x3","x4","x5","x6"),
	latentVars="F1",
	mxPath(from=c("x1","x2","x3","x4","x5","x6"),
		arrows=2,
		free=TRUE,
		values=c(1,1,1,1,1,1),
		labels=c("e1","e2","e3","e4","e5","e6")
	), 
	# residual variances
	# -------------------------------------
	mxPath(from="F1",
		arrows=2,
		free=TRUE,
		values=1,
		labels ="varF1"
	), 
	# latent variance
	# -------------------------------------
	mxPath(from="F1",
		to=c("x1","x2","x3","x4","x5","x6"),
		arrows=1,
		free=c(FALSE,TRUE,TRUE,TRUE,TRUE,TRUE),
		values=c(1,1,1,1,1,1),
		labels =c("l1","l2","l3","l4","l5","l6")
	), 
	# factor loadings
	# -------------------------------------
	mxPath(from="one",
		to=c("x1","x2","x3","x4","x5","x6","F1"),
		arrows=1,
		free=c(TRUE,TRUE,TRUE,TRUE,TRUE,TRUE,FALSE),
		values=c(1,1,1,1,1,1,0),
		labels =c("meanx1","meanx2","meanx3","meanx4","meanx5","meanx6",NA)
	) 
	# means
	# -------------------------------------
) # close model
# Create an MxModel object
# -----------------------------------------------------------------------------
oneFactorFit <- mxRun(oneFactorModel)      
standardizeMx(oneFactorFit)

################## Fit Measures

data(demoOneFactor)
manifests <- names(demoOneFactor)
latents <- c("G")
factorModel <- mxModel("One Factor", 
    type="RAM",
    manifestVars=manifests, 
    latentVars=latents,
    mxPath(from=latents, to=manifests),
    mxPath(from=manifests, arrows=2),
    mxPath(from=latents, arrows=2, free=FALSE, values=1.0),
    mxData(observed=cov(demoOneFactor), type="cov", numObs=500)
)
factorFit <- mxRun(factorModel)
round(fitMeasuresMx(factorFit), 3)
