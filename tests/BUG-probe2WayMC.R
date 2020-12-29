
#double mean centering

dbc<-read.csv('dogintol2016.dmc.csv',header=TRUE,sep=",")
names(dbc)

library(semTools)

model<-'
intol=~1*intol1+intol2+intol3+intol4+intol5
dogma=~1*dog1c+dog2c+dog3c
agree=~1*agree1c+agree2c+agree3c
interact=~1*d1a1dmc+d1a2dmc+d1a3dmc+d2a1dmc+d2a2dmc+d2a3dmc+d3a1dmc+d3a2dmc+d3a3dmc
intol~dogma+agree+interact
d1a1dmc~~d1a2dmc+d1a3dmc+d2a1dmc+d3a1dmc
d1a2dmc~~d1a3dmc+d2a2dmc+d3a2dmc
d1a3dmc~~d2a3dmc+d3a3dmc
d2a1dmc~~d2a2dmc+d2a3dmc+d3a1dmc
d2a2dmc~~d2a3dmc+d3a2dmc
d2a3dmc~~d3a3dmc
d3a1dmc~~d3a2dmc+d3a3dmc
d3a2dmc~~d3a3dmc
dogma~~agree
agree~~interact
dogma~~interact
'

fit<-sem(model,dbc,meanstructure=TRUE)
summary(fit,fit.measures=TRUE,stand=TRUE)

## manually calculate simple slopes as user-defined parameters
model1<-'
intol=~1*intol1+intol2+intol3+intol4+intol5
dogma=~1*dog1c+dog2c+dog3c
agree=~1*agree1c+agree2c+agree3c
interact=~1*d1a1dmc+d1a2dmc+d1a3dmc+d2a1dmc+d2a2dmc+d2a3dmc+d3a1dmc+d3a2dmc+d3a3dmc
# labels added to first order effect of dogmatism and to the interaction term to compute a+b(-/+ sd for agree)
intol~d*dogma+a*agree+b*interact
d1a1dmc~~d1a2dmc+d1a3dmc+d2a1dmc+d3a1dmc
d1a2dmc~~d1a3dmc+d2a2dmc+d3a2dmc
d1a3dmc~~d2a3dmc+d3a3dmc
d2a1dmc~~d2a2dmc+d2a3dmc+d3a1dmc
d2a2dmc~~d2a3dmc+d3a2dmc
d2a3dmc~~d3a3dmc
d3a1dmc~~d3a2dmc+d3a3dmc
d3a2dmc~~d3a3dmc
dogma~~agree
agree~~interact
dogma~~interact

# simple slopes of agreeableness
ssa_low := a+b*(-.356)
ssa_med := a
ssa_high := a+b*(.356)

#simple slopes of dogma
ssd_low := d+b*(-.356)
ssd_med := d
ssd_high := d+b*(.356)
'

fit1<-sem(model1,dbc,meanstructure=TRUE)
summary(fit1,fit.measures=TRUE,stand=TRUE)

## compare to semTools
valProbe <- c(-.356,0,.356)

## effect of agree -- list moderator first
probe2WayMC(fit, nameX=c("dogma","agree","interact"), nameY="intol", modVar="dogma", valProbe = valProbe)
## effect of Dogma -- list moderator first
probe2WayMC(fit, nameX=c("agree","dogma","interact"), nameY="intol", modVar="agree", valProbe = valProbe)

## slopes and SE match, but they are assigned to the wrong focal variable!
## now list the focal variable first, then the moderator:

## effect of agree -- list moderator second (SEs are WRONG)
probe2WayMC(fit, nameX=c("agree","dogma","interact"), nameY="intol", modVar="dogma", valProbe = valProbe)
## effect of Dogma -- list moderator second (SEs are WRONG)
probe2WayMC(fit, nameX=c("dogma","agree","interact"), nameY="intol", modVar="agree", valProbe = valProbe)

## the estimates are now appropriate for the focal predictor,
## but the SEs are switched!


                 