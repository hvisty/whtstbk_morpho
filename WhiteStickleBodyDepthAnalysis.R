##Initial White Stickleback analysis and graphs##


data<-read.csv("EastCoastMorphometrics2015.csv")
data<-as.data.frame(data) #it was displaying as a function instead of a data frame?
class(data)
View(data)

species<-data$species
population<-data$population
depth.1<-data$depth.1
depth.2<-data$depth.2
SL<-data$std.length

library(ggplot2)


##GRAPHS##
#Body Depth by Population
BDP<-ggplot(data, aes(population, depth.1, colour=factor(species))) + xlab("Population") + ylab("Body Depth (cm)")
BDP + geom_boxplot()
                   
BDP + geom_boxplot(aes(population, depth.2)) #Both depth measurements show the same pattern.
  #First measurement is likely more reliable, as ventral spines (or whatever they're called) sometimes obscured the body measurement from the second spine

#Standard Length by Population
SLP<-ggplot(data, aes(population, SL, colour=factor(species))) + xlab("Population") + ylab("Standard Length (cm)")
SLP + geom_boxplot()

#Body Depth by White or Common individuals
BDW<-ggplot(data, aes(species, depth.1, colour=factor(species))) + xlab("Species") + ylab("Body Depth (cm)")
BDW + geom_boxplot() #lookin' good

#Standard Length by White or Common individuals
SLW<-ggplot(data, aes(species, SL, colour=factor(species))) + xlab("Species") + ylab("Standard Length (cm)")
SLW + geom_boxplot() #no difference...


#Bimodal (violin distribution) visual with Body Depth for each population
BDVP<-ggplot(data, aes(population, depth.1, color=factor(species))) + xlab("Population") + ylab("Body Depth (cm)")
BDVP + geom_violin()

#Bimodal (violin distribution) visual with Standard Length for each population
SLVP<-ggplot(data, aes(population, SL, color=factor(species))) + xlab("Population") + ylab("Standard Length (cm)")
SLVP + geom_violin()

#Violin plot for Body Depth by species (extra, out of interest)
BDVW<-ggplot(data, aes(species, depth.1, color=factor(species))) + xlab("Species") + ylab("Body Depth (cm)")
BDVW + geom_violin() #cool because BOTH shows bimodal distribution, others not so much



##HYPOTHESES TESTING## 

#constructing a data frame of averages...
popsummary<-aggregate(SL, by=list(population), FUN='mean')
names(popsummary)[names(popsummary)=="Group.1"]<-"Population"
names(popsummary)[names(popsummary)=="x"]<-"Avg.SL"
depth.means<-aggregate(depth.1, by=list(population), FUN="mean")
popsummary[,"Avg.depth"]<-depth.means$x
species.agg<-aggregate(species, by=list(population), FUN="unique")
popsummary[,"Species"]<-species.agg$x
popsummary<-popsummary[c(1,4,2,3)]
View(popsummary) #woot woot 

#one way ANOVA: mean by population Standard Length for White/Common/Both
SLANOVA.AVGS<-aov(popsummary$Avg.SL~popsummary$Species)
summary(SLANOVA.AVGS) #F=0.156, p>0.856, DF=2

#one way ANOVA: mean by pop. Body Depth for White/Common/Both
BDANOVA.AVGS<-aov(popsummary$Avg.depth~popsummary$Species)
summary(BDANOVA.AVGS) #F=5.556, p>0.012*, DF=2
TukeyHSD(BDANOVA.AVGS) #Common and Both are not significantly different, but White species is significantly different from both others


#test for bimodaltiy by population: normal with p value and bimodal test if not normal
#done with body depth since this value was the one with significant differences 

#Use Shapiro Wilk Test for normality (adding values to data sheet)
normaltest<-aggregate(depth.1, by=list(population), FUN="shapiro.test")
normaltest<-normaltest[,2]
normaltest<-as.data.frame(normaltest)
normaltest$statistic<-as.numeric(normaltest$statistic)
normaltest$p.value<-as.numeric(normaltest$p.value)

popsummary[,"W"]<-normaltest$statistic
popsummary[,"W.p.value"]<-normaltest$p.value
popsummary[,"Normal"]<-NULL

for(i in 1:nrow(popsummary)) {
  x1<-popsummary[i,"W.p.value"]
  if((x1<0.05)==TRUE){
    popsummary[i,"Normal"]<-"Y"
  } else {
    popsummary[i,"Normal"]<-"N"
  }}

View(popsummary)

#Test bimodality for values failing the Normality test: Diptest for Unimodality 
library(diptest) 

diptest<-aggregate(depth.1, by=list(population), FUN="dip.test")
diptest<-diptest[,2]
diptest<-as.data.frame(diptest)
diptest$statistic<-as.numeric(diptest$statistic)
diptest$p.value<-as.numeric(diptest$statistic) #changing class from list to numeric changes p value??

popsummary[,"D"]<-diptest$statistic
popsummary[,"D.p.value"]<-diptest$p.value
View(popsummary)

popsummary[,"Unimodal"]<-NULL

for(i in 1:nrow(popsummary)) {
  x2<-popsummary[i,"D.p.value"]
  x3<-popsummary[i,"Normal"]
  if(x3=="N"){
  if((x2<0.05)==TRUE){
    popsummary[i,"Unimodal"]<-"N"
  }}}
for(i in 1:nrow(popsummary)) {
  x2<-popsummary[i,"D.p.value"]
  x3<-popsummary[i,"Normal"]
  if(x3=="N"){
    if((x2>0.05)==TRUE){
      popsummary[i,"Unimodal"]<-"Y"
    }}}

write.csv(popsummary,file="WhiteStickelBackAnalysisSummary.csv")

remove(list=ls()) 
