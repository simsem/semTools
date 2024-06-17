.onAttach <- function(libname, pkgname) {
    version <- read.dcf(file = system.file("DESCRIPTION", package = pkgname), fields = "Version")
    packageStartupMessage(" ")
    packageStartupMessage("###############################################################################")
    packageStartupMessage("This is ", paste(pkgname, version))
    packageStartupMessage("All users of R (or SEM) are invited to submit functions or ideas for functions.")
    packageStartupMessage("###############################################################################")

    ## if lavaan.mi is already loaded, warn users to use :: or switch order
    if ("lavaan.mi" %in% names(utils::sessionInfo()$otherPkgs)) {
      packageStartupMessage(" ")
      packageStartupMessage("The lavaan.mi package was already attached when semTools was attached.")
      packageStartupMessage("To access lavaan.mi functions masked by semTools, use the double-colon:")
      packageStartupMessage(" ")
      packageStartupMessage("\t lavaan.mi::cfa.mi() or lavaan.mi::lavTestLRT.mi()\n")
      packageStartupMessage("It is preferable to first attach semTools, then attach lavaan.mi, but")
      packageStartupMessage("you can change the path order by detaching and reloading lavaan.mi,")
      packageStartupMessage("using this syntax:\n")
      packageStartupMessage('\t detach("package:lavaan.mi", unload = TRUE)')
      packageStartupMessage('\t library(lavaan.mi) \n')
      packageStartupMessage("Then the deprecated semTools::lavaan.mi functionality will be masked.")
    }
}

.onLoad <- function(libname, pkgname) {
  ## "register" emmeans functionality
  emInstalled <- try(loadNamespace("emmeans"), silent = TRUE)
  if (!inherits(emInstalled, "try-error")) {
      emmeans::.emm_register("lavaan", pkgname)
  }
}

