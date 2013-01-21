#############################################
## tukeySEM -- a function to compute       ##
## tukey's post-hoc test with unequal      ##
## sample sizes and variances              ##
##                                         ##
##  Alexander M. Schoemann                 ##
##  Last edited on 01/16/2013              ##
#############################################

##inputs: mean of group 1, mean of group2, variance of group 1, variance of group 2
## sample size or group1, sample size of group2, number of groups in the ANOVA
##Output: vector containing the q statistic, degrees of freedom, and associated p value

tukeySEM <- function(m1, m2, var1, var2, n1, n2, ng){
qNum <- abs(m1 - m2)
qDenom <- sqrt(((var1/n1) + (var2/n2))/2)
Tukeyq <- qNum/qDenom
Tukeydf <- ((var1/n1) + (var2/n2))^2 / (((var1/n1)^2/(n1-1)) + ((var2/n2)^2/(n2-2)))
p <- 1- ptukey(Tukeyq, ng, Tukeydf)
cols <- c("q", "df", "p")
res <- c(Tukeyq, Tukeydf, p)
names(res) <- cols
res
}

##Example from Schoemann (2013)
##Bio vs. policial science on evo misconceptions
#tukeySEM(3.91, 3.96,.46, .62, 246, 425,3)
