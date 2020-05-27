.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file = system.file("DESCRIPTION", package = pkgname), fields = "Version")
    packageStartupMessage(" ")
    packageStartupMessage("###############################################################################")
    packageStartupMessage("This is ", paste(pkgname, version))
    packageStartupMessage("All users of R (or SEM) are invited to submit functions or ideas for functions.")
    packageStartupMessage("###############################################################################")
}

.onLoad <- function(libname, pkgname) {
    if (requireNamespace("emmeans", quietly = TRUE)){
        emmeans::.emm_register("lavaan", pkgname)
        message('semTools now supports using the emmeans package on lavaan ',
                'objects. For multigroup models, remember to set ',
                '"nesting = NULL" in order to treat the grouping variable as ',
                'another design factor in the ref_grid.')
    }
}

