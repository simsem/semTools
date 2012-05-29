.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file = system.file("DESCRIPTION", package = pkgname), fields = "Version")
    packageStartupMessage(" ")
    packageStartupMessage("###############################################################################")
    packageStartupMessage("This is ", paste(pkgname, version))
    packageStartupMessage("All users of R (or SEM) are invited to submit functions or ideas for functions.")
    packageStartupMessage("###############################################################################")
}
 
