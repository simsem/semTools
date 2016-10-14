########### Mauricio Garnier Villarreal (mgv@ku.edu)
### Last updated: 14 October 2016
######This function estimates the Fraction of Missing Information for the variance and mean of each variable in a list of multiple imputed data sets
#### dat.imp is a list of the imputed data sets
#### method is the model used for the estimation
#### varnames is used to select a subset of variables
#### digits is the number of decimals
#### group is the grouping variable, in case you want to get the fmi for each group
#### exclude are the variables that you wnat to exclude from the analysis

fmi <- function(dat.imp, method="saturated", varnames=NULL, group=NULL, exclude=NULL, digits=3){ 
  
  if(is.character(varnames)){
    vars <- varnames
  } else {
    vars <- colnames(dat.imp[[1]])
  }
  
  if(!is.null(group)){
    vars <- vars[vars!=group]
  }
  
  if(!is.null(exclude)){
    vars <- vars[vars!=exclude]
  }
  
  if(method == "saturated" | method == "sat"){
    par.tab <- satParFMI(dat.imp, var.names=vars, groups=group) 
  }
  if(method == "null"){
    par.tab <- nullParFMI(dat.imp, var.names=vars, groups=group)
  }
    
  comb.results1 <- cfa.mi(par.tab, dat.imp, chi="none", meanstructure = TRUE, group = group)
  
  comb.results <- inspect(comb.results1, "impute")[[2]]  ## FIXME: can't just be lavInspect because it is a lavaanStar
  
  comb.results <- data.frame(comb.results[,c("lhs","op","rhs","group")], 
                             round(lavaan::parameterEstimates(comb.results1)[,"est"], digits), 
                             round(comb.results[,c("fmi1","fmi2")], digits))
  
  colnames(comb.results) <- c('lhs', 'op', 'rhs', 'group', 'coef', 'fmi.1', 'fmi.2')
  
  variances <- comb.results[comb.results$lhs==comb.results$rhs,]
  
  variances <- data.frame(variances[,"lhs"], variances[,"group"], variances[,"coef"],
                          variances[,"fmi.1"], variances[,"fmi.2"])
  
  colnames(variances) <- c('var', 'group', 'coef', 'fmi.1', 'fmi.2')
  
  var.means <- comb.results[comb.results$op=="~1",]
  
  var.means <- data.frame(var.means[,"lhs"], var.means[,"group"], var.means[,"coef"],
                          var.means[,"fmi.1"], var.means[,"fmi.2"])
  
  colnames(var.means) <- c('var', 'group', 'coef', 'fmi.1', 'fmi.2')
  
  if(method == "null"){
    mes <- "These estimates used the null model, they may not be as precise as the saturated model estimates"
    results<-list(Variances=variances, Means=var.means, Message=mes)
  } else {
    results<-list(Variances=variances, Means=var.means)
  }
  
  return(results)
}

#### function to produce a parameter table for the saturated model
satParFMI <- function(dat.imp, var.names=NULL, groups=NULL){
  
  if(!is.null(groups)){
    ngroups <- length(table(dat.imp[[1]][,groups]))
  } else {
    ngroups <- 1
  }
  
  # gets the parameter table from the null model
  par.null <- nullParFMI(dat.imp, var.names, groups=groups)
  lhs.diag <- par.null$lhs
  op.diag <- par.null$op
  rhs.diag <- par.null$rhs
  gnull <- par.null$group
  #combine the variable names to set al the covariances
  combs <- t(combn(var.names, 2))
  lhs.up <- rep(combs[, 1],times=ngroups)
  op.up <- rep("~~", length(lhs.up))
  rhs.up <- rep(combs[, 2],times=ngroups)
  galt <- sort(rep(1:ngroups,times=length(lhs.up)/ngroups))
  #put together the null table and the covariances
  lhs.all <- c(lhs.up, lhs.diag)
  id <- seq(1:length(lhs.all))
  op.all <- c(op.up, op.diag)
  rhs.all <- c(rhs.up, rhs.diag)
  user <- rep(1,length(lhs.all))
  group <- as.integer(c(galt,gnull))
  free <- as.integer(id)
  ustart <- rep(NA, length(lhs.all))
  exo <- rep(0, length(lhs.all))
  label <- rep("", length(lhs.all))
  plabel <- rep("", length(lhs.all))
  par.sat <- list(id, lhs.all, op.all, rhs.all, user, group,
                  free, ustart, exo, label, plabel)
  names(par.sat) <- c("id", "lhs", "op", "rhs", "user", "group", "free", "ustart", "exo", "label", "plabel")
  return(par.sat)
}

#### function to produce a parameter table for the null model
nullParFMI <- function(dat.imp, var.names=NULL, groups=NULL){
  
  if(!is.null(groups)){
    ngroups <- length(table(dat.imp[[1]][,groups]))
  } else {
    ngroups <- 1
  }
  
  lhs.diag1 <- rep(c(var.names),times=ngroups)
  op.diag1 <- rep("~~",ngroups*(length(var.names)))
  rhs.diag1 <- rep(var.names,times=ngroups)
  group1 <- sort(rep(1:ngroups,times=length(lhs.diag1)/ngroups))
  
  lhs.diag2 <- rep(c(var.names),times=ngroups)
  op.diag2 <- rep("~1",ngroups*(length(var.names)))
  rhs.diag2 <- rep("",ngroups*length(var.names))
  group2 <- sort(rep(1:ngroups,times=length(lhs.diag2)/ngroups))
  
  lhs.diag <- c(lhs.diag1, lhs.diag2)
  op.diag <- c(op.diag1, op.diag2)
  rhs.diag <- c(rhs.diag1, rhs.diag2)
  group <- c(group1, group2)
  first <- data.frame(lhs.diag,op.diag,rhs.diag,group)
  first <- first[order(first$group),]
  id <- seq(1:length(lhs.diag))
  user <- rep(1,length(lhs.diag))
  free <- as.integer(id)
  ustart <- rep(NA, length(lhs.diag))
  exo <- rep(0, length(lhs.diag))
  label <- rep("", length(lhs.diag))
  plabel <- rep("", length(lhs.diag))
  null.sat.fmi <- list(id, as.character(first$lhs.diag), as.character(first$op.diag), 
                       as.character(first$rhs.diag), user, first$group,
                       free, ustart, exo, label, plabel)
  names(null.sat.fmi) <- c("id","lhs","op","rhs","user","group","free","ustart","exo","label","plabel")
  return(null.sat.fmi)
}
