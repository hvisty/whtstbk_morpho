##landmark analysis in Geomorph##

#### LIBRARIES ####

# required libraries
install.packages("geomorph")
install.packages("candisc")
install.packages("MASS")
install.packages("candisc")
install.packages("devtools")


# google sheets API
library(devtools)
devtools::install_github("jennybc/googlesheets")

# load libraries
library("geomorph")
library("ggplot2")
library("candisc")
library("MASS")
library("dplyr")
library("googlesheets")

#### LIBRARIES ####



#### LOAD RAW DATA ####

# google sheets method

## list your sheets and authenticate with server (will open browswer)
gs_ls()

## load in sheet from google sheets
morphodat.gs <- gs_title("East Coast Morphometrics - 2015")
morphodat <- data.frame(get_via_csv(morphodat.gs))


## write csv to file (for working offline)
write.csv(morphodat, file="east_coast_morphometrics_2015.csv")
# morphodat <- read.csv("east_coast_morphometrics_2015.csv")

# build 'corrected' landmark file
# pop, id, sex, species, x.1,y.1,x.2,y.2, etc.

## select landmark columns, remove nas
landmarks.full <- morphodat %>%
  select(population, individual, sex, species, correction, contains("lndmrk")) %>%
  filter(!is.na(lndmrk.x1))

## extract ids
landmarks.ids <- landmarks.full %>%
  select(population, individual, sex, species)

## fix names
names(landmarks.full) <- gsub("lndmrk.","",names(landmarks.full)) %>% gsub ("^x", "x.", .) %>% gsub ("^y", "y.", .)

## strip whitespace (...) and apply correction
landmarks.matrix <- landmarks.full %>%
  select(matches("^[x,y]{1}.")) %>%
  as.list %>%
  lapply(function(x)gsub(" ","",x)) %>%
  lapply(as.numeric) %>%
  as.data.frame %>%
  as.matrix

## apply correction
landmarks.matrix <- data.frame(sweep(landmarks.matrix, MARGIN=1, landmarks.full$correction, `/`))

## rearrange so x.1,y.1, not x.1,x.2
landmarks.order <- landmarks.matrix  %>%
  names %>%
  gsub ("^x.", "", .) %>%
  gsub ("^y.", "", .) %>%
  as.numeric %>%
  order

landmarks.raw <- cbind(landmarks.ids, landmarks.matrix[,landmarks.order])

## manual method
## landmarks <- read.csv("Corrected Landmark Data.csv")

#### LOAD RAW DATA ####

##Removing outliers
#run initial analysis, remove outliers from this, then re-do data sheet for final landmarks

landmarks.outliers <- landmarks.raw %>%
  filter(!is.na(x.1))
landmarks.outliers$ID <- paste(landmarks.outliers$population, landmarks.outliers$individual, sep="_")


landmarks.data <- landmarks.outliers %>%
  select(matches("^[x,y]{1}."))

landmark.array <- arrayspecs(landmarks.data, 19, 2) #3D array

GPA.landmarks <- gpagen(landmark.array) #Procrustes analysis

GPA.2D <- two.d.array(GPA.landmarks$coords) #2D Data frame of procrustes coordinates

#ouliers found from principle component stuff. 

PCA<- plotTangentSpace(GPA.landmarks$coords, groups= factor(landmarks.outliers$population), verbose=T)
PC.scores <- as.data.frame(PCA$pc.scores) #gives the PCA scores
PC.scores$ID <- paste(landmarks.outliers$population, landmarks.outliers$individual, sep="_")


#determined outliers from PC plot....
PC.scores <- filter(PC.scores, PC1> -0.10 & PC1< 0.10 & PC2> -0.06 & PC2< 0.06)
ggplot(PC.scores, aes(x=PC1, y=PC2, label=ID)) +geom_point() +geom_text() #filter works

#create updated landmarks data frame, outliers removed 
landmarks <- merge(landmarks.outliers, PC.scores, by=c("ID"))
landmarks <- landmarks[!duplicated(landmarks[,"ID"]),] #duplicated IDs upon merging..
landmarks <- subset(landmarks, select=-c(PC1:PC38))


#NOW start real analysis

#### GEOMORPH ANALYSIS ####
#Make a 3D array; only x/y cordinates

#reclassify SF and SR as 'W' populations

for (i in 1:nrow(landmarks)){
  if(landmarks$population[i] =="SF2014" | landmarks$population[i] == "SR2014") {
    landmarks$species[i] <- "W"
  }
}

landmarks.data <- landmarks %>%
  select(matches("^[x,y]{1}."))

landmark.array <- arrayspecs(landmarks.data, 19, 2) #3D array

GPA.landmarks <- gpagen(landmark.array) #Procrustes analysis

GPA.2D <- two.d.array(GPA.landmarks$coords) #2D Data frame of procrustes coordinates


#principle component stuff. 
PCA<- plotTangentSpace(GPA.landmarks$coords, groups= factor(landmarks$population), verbose=T)
PCA$pc.scores #gives the PCA scores

#unbending; just show me PC2 vs PC3
plotTangentSpace(GPA.landmarks$coords, groups=factor(landmarks$species), axis1=2, axis2=3)

#make a prettier plot! 
PC.scores<-as.data.frame(PCA$pc.scores)
ggplot(PC.scores, aes(x=PC3, y=PC4, color=landmarks$species)) + geom_point(size=3) 

ggplot(PC.scores, aes(x=PC5, y=PC6, color=landmarks$species)) + geom_point(size=3) 
#...nothin to show for shape differences. Dang. 

#MANOVA 
PC1 <- as.vector(PCA$pc.scores[,1]) #get a vector of PC1 to use as a covariate

procD.lm(GPA.2D~landmarks$population*landmarks$species, iter=99) #without PC1 bending accounted for
procD.lm(GPA.2D~landmarks$population+landmarks$species*PC1, iter=99) #if PC1 is a covariate...

#Linear discriminant function analysis 

lda.species <- lda(GPA.2D, landmarks$species) #what does this mean? 

# the default plot, looks meh
plot(lda.species)

# project the original values in lda space
plda <- predict(object = lda.species,
                newdata = GPA.2D)

# the percent of variance explained by the LD funcitons
prop.lda <- lda.species$svd^2/sum(lda.species$svd^2) 
prop.lda <- round(prop.lda*100)

# data frame of projected data with species names
lda.project <- data.frame(species = landmarks$species, 
                          population = landmarks$population, 
                          individual = landmarks$individual,
                          ID = landmarks$ID,
                          ld1 = plda$x[,1], 
                          ld2 = plda$x[,2])

# plot LDA without "Both" locations
lda.project %>%
  filter(species != "B") %>%
  ggplot(aes(color = species, x = ld1,y = ld2))+
  geom_point(size = 3) +
  labs(x = paste0("LD1 (", prop.lda[1], "%)"),
       y = paste0("LD2 (", prop.lda[2], "%)"))


# the loadings of the LDA axes
lda.species$scaling

# make a data frame with their original names (x.1, y.1 ,etc.)
landmark.names <- names(landmarks)[-c(1:4)][-39]

lm.loadings <- data.frame(lm.name = landmark.names, 
                          ld1.loading = lda.species$scaling[,1],
                          ld2.loading = lda.species$scaling[,2])

head(lm.loadings)


####REGRESSION WORK####

##make a data frame to look at lda values versus procrustes x/y

GPA <- as.data.frame(GPA.2D)
GPA.X <- GPA[,seq(from=1, to=38, by=2)] #data was in x,y not x(1:19)...
GPA.X$ID <- lda.project$ID
GPA.Y <- GPA[,seq( from=2, to=38, by=2)]
GPA.Y$ID <- lda.project$ID
GPA <- merge(GPA.X, GPA.Y, by="ID")
names(GPA)[c(2:39)] <- c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10",
                         "X11","X12","X13","X14","X15","X16","X17","X18","X19",
                         "Y1","Y2","Y3","Y4","Y5","Y6","Y7","Y8","Y9","Y10",
                         "Y11","Y12","Y13","Y14","Y15","Y16","Y17","Y18","Y19")

lda.xy <- merge(lda.project, GPA, by="ID")


#create a linear regression coefficient vector for each 1-19 x:ld1 and y:ld2
lda.lm<-data.frame(matrix(NA, nrow = nrow(lda.xy) , ncol = length(PC.scores)))
colnames(lda.lm)<-c("LMx.1","LMx.2","LMx.3","LMx.4","LMx.5","LMx.6","LMx.7","LMx.8",
                    "LMx.9","LMx.10","LMx.11","LMx.12","LMx.13","LMx.14","LMx.15",
                    "LMx.16","LMx.17","LMx.18","LMx.19",
                    "LMy.1","LMy.2","LMy.3","LMy.4","LMy.5","LMy.6","LMy.7","LMy.8",
                    "LMy.9","LMy.10","LMy.11","LMy.12","LMy.13","LMy.14","LMy.15",
                    "LMy.16","LMy.17","LMy.18","LMy.19", "species")
lda.lm$species <- lda.xy$species
ld1 <- lda.xy$ld1
ld2 <- lda.xy$ld2

# loop the regressions. Separate for X and Y cause my looping skillz are weak

for (l in 1:19){
  pie <- as.matrix(lda.xy[paste("X", l, sep="")])
    lm.1<-lm(ld1 ~ pie)
    lda.lm[,l] <- lm.1$fitted.values
}

for (l in 1:19){
  cake <- as.matrix(lda.xy[paste("Y", l, sep="")])
    lm.2<-lm(ld2 ~ cake)
    lda.lm[,(l+19)] <- lm.2$fitted.values
}

plot(ld1 ~ lda.lm$LMx.4)

#now, get the LM averages and combine with the avg procrustes values.... 

library(reshape2)



#get averages of lm values for each landmark, by species
lm.avg <- aggregate(lda.lm[,c(1:38)], by=list(lda.lm$species), FUN='mean', na.rm=T)
lm.avg$species <- c("B","C","W")
lm.avg <- lm.avg[,-1]
lm.avg <- melt(lm.avg, id.vars="species")
colnames(lm.avg) <- c("species", "landmark", "Avg LM fit")

#reorganize data frame for easier plotting
lm.avg$landmark <- sub('LMx.', 'X_', lm.avg$landmark, fixed=T)
lm.avg$landmark <- sub('LMy.', 'Y_', lm.avg$landmark, fixed=T)

lm.avg$XorY <- sapply(strsplit(lm.avg$landmark,'_'), "[", 1)
lm.avg$landmark.num <- sapply(strsplit(lm.avg$landmark,'_'), "[", 2)
lm.avg <- as.data.frame(split(lm.avg, lm.avg$XorY))
lm.avg <-lm.avg[order(lm.avg[,1]), ]
lm.split <- as.data.frame(split(lm.avg, lm.avg$X.species))

lm.avg <- data.frame(number= c(1:19),
                     LM.X.B=lm.split$B.X.Avg.LM.fit,
                     LM.X.C=lm.split$C.X.Avg.LM.fit,
                     LM.X.W=lm.split$W.X.Avg.LM.fit,
                     LM.Y.B=lm.split$B.Y.Avg.LM.fit,
                     LM.Y.C=lm.split$C.Y.Avg.LM.fit,
                     LM.Y.W=lm.split$W.Y.Avg.LM.fit)


#get avg. procrustes coordinates by landmark
GPA.avgs <- as.data.frame(colMeans(GPA[,2:39],))
GPA.avgs$number <- c(1:19)
GPA.avgs$XorY <- c(rep("X", 19), rep("Y", 19))
names(GPA.avgs)[1] <- "Avg.GPA"
GPA.avgs<- as.data.frame(split(GPA.avgs, GPA.avgs$XorY))
GPA.avgs<- data.frame(number=c(1:19),
                      GPA.x = GPA.avgs$X.Avg.GPA,
                      GPA.y = GPA.avgs$Y.Avg.GPA)

all.avg <- merge(GPA.avgs, lm.avg, by="number")

#differene between displacement values by procrustes avgs
all.avg$LM.X.B <- all.avg$GPA.x - lm.avg$LM.X.B
all.avg$LM.X.C <- all.avg$GPA.x - lm.avg$LM.X.C 
all.avg$LM.X.W <- all.avg$GPA.x - lm.avg$LM.X.W 

all.avg$LM.Y.B <- all.avg$GPA.y - lm.avg$LM.Y.B 
all.avg$LM.Y.C <- all.avg$GPA.y - lm.avg$LM.Y.C 
all.avg$LM.Y.W <- all.avg$GPA.y - lm.avg$LM.Y.W 

##testing for errors
lm.avg.test <- (lm.avg$LM.X.C * lm.avg$LM.X.W)/2
all.avg$LM.X.test<- all.avg$GPA.x - lm.avg.test
lm.avg.testy <- (lm.avg$LM.Y.C * lm.avg$LM.Y.W)/2
all.avg$LM.Y.test<- all.avg$GPA.y - lm.avg.testy #avg of avgs... does overall avg work too??

#overall averages (really small values)
avg.test <- lda.lm[-39]
avg.test <- as.data.frame(colMeans(avg.test))
names(avg.test)[1] <- "avgs"
avg.test2 <- data.frame(number=rep(c(1:19), 2),
                      XorY= c(rep("X", 19), rep("Y", 19)),
                      avgs= avg.test$avgs)
avg.test <- as.data.frame(split(avg.test2, avg.test2$XorY))
avg.test$X.avgs <- all.avg$GPA.x - avg.test$X.avgs
avg.test$Y.avgs <- all.avg$GPA.y - avg.test$Y.avgs #displacements are tiny.


#plot avg. procrustes and their displacements 
ggplot(all.avg, aes(x=LM.X.C, y=LM.Y.C))+
  geom_text(aes(label=number)) +
  geom_text(data=all.avg, aes(x=LM.X.W, y=LM.Y.W, label=number), color="blue") 

#plot of C*W averages 
ggplot(all.avg, aes(x=LM.X.test, y=LM.Y.test))+
   geom_text(aes(label=number)) +
  geom_text(data=all.avg, aes(x=GPA.x, y=GPA.y, label=number), color="blue") 
 

#plot overall averages as a check. 
ggplot(all.avg, aes(x=GPA.x, y=GPA.y))+
  geom_text(aes(label=number), color="blue") +
  geom_text(data=avg.test, aes(x=X.avgs, y=Y.avgs, label=X.number), color="red")   


######


##### CROSS VALIDATION of LDA #####

#take half of the values and create the LDA, then plot the rest to check model 


#LDA with CV function
GPA.CW <- GPA
GPA.CW$species <- landmarks$species
GPA.CW$centroid = GPA.landmarks$Csize

GPA.CW<- GPA.CW %>%
  filter(species != "B") #remove "Both" category

landmarks.CW<- landmarks %>%
  filter(species != "B")

GPA.CW <- GPA.CW[-c(39,40)]

lda.species.CV <- lda(GPA.CW, landmarks.CW$species, CV=T)

# percent correct (accuracy of prediction)
ct <- table(landmarks.CW$species, lda.species.CV$class)
diag(prop.table(ct, 1))
# total percent correct
sum(diag(prop.table(ct)))

#visually plot this 
ct.plot <- as.data.frame(ct)
ggplot(ct.plot, aes(color=Var1, x=Var2, y=Freq, shape=Var2)) + 
  geom_point(size=6)
  #ehhh not that useful

#manual: sample a random half of the samples, 
  #with equal representation of W and C 

GPA.50<- data.frame(species = landmarks$species, 
      population = landmarks$population, 
      individual = landmarks$individual,
      ID = landmarks$ID,
      centroid = GPA.landmarks$Csize)

GPA.50<-merge(GPA.50, GPA, by.x="ID")

GPA.parse<- GPA.50 %>% 
  filter(species != "B") %>%
  group_by(species) 

GPA.parse <- sample_n(GPA.parse, 403)
GPA.parse1<- GPA.parse[c(1:201, 403:604),]
GPA.parse2 <- GPA.parse[c(202:402,605:806),]

GPA.1 <- GPA.parse1[-c(1:4)]
GPA.2 <- GPA.parse2[-c(1:4)]


#run LDA for each 1/2 sample
lda.CV1 <- lda(GPA.1, as.character(GPA.parse1$species))
lda.CV2 <- lda(GPA.2, as.character(GPA.parse2$species))

#subtract scaling between the two halves and look at differences 

lda.CV <- as.vector(lda.CV1$scaling - lda.CV2$scaling)
plot(lda.CV) #weird LDA units 

#lets look at the plda and overplotting for each half
#project in lda space
plda.CV1 <- predict(object = lda.CV1,
                newdata = GPA.1)

plda.CV2 <- predict(object = lda.CV2,
                    newdata = GPA.2)

# data frames of projected data with species names
lda.project.CV1 <- data.frame(species = GPA.parse1$species, 
                          population = GPA.parse1$population, 
                          individual = GPA.parse1$individual,
                          ID = GPA.parse1$ID,
                          ld1 = plda.CV1$x)

lda.project.CV2 <- data.frame(species = GPA.parse2$species, 
                              population = GPA.parse2$population, 
                              individual = GPA.parse2$individual,
                              ID = GPA.parse2$ID,
                              ld1 = plda.CV2$x)

##plot both halves on one plot... looks like they match 
ggplot(lda.project.CV1,aes(fill = species, x = LD1))+
  geom_histogram() +
  geom_histogram(aes(x = lda.project.CV2$LD1, fill = species, size=2))
 
  


#other option: project the remaining half in lda space. Use one half LDA and 
#  predict the other half GPA, then plot it
plda.CV <- predict(object = lda.CV1,
                    newdata = GPA.2)

lda.project.CV <- data.frame(species = GPA.parse1$species, 
                              population = GPA.parse1$population, 
                              individual = GPA.parse1$individual,
                              ID = GPA.parse1$ID,
                              ld1 = plda.CV$x)

#plot predicted lda with matching half and cross-fit half
ggplot(lda.project.CV1, aes(x = LD1, colour = species))+
  geom_histogram() +
  geom_histogram( aes(fill = species, x = lda.project.CV$LD1)) 
  #they're identical (nearly)


######## GENETICS VS. MORPHOLOGY ########

# extract centroids
landmark.centroids <- data.frame(population = landmarks$population, 
                                 id = landmarks$ID, 
                                 species = landmarks$species, 
                                 centroid = GPA.landmarks$Csize)
# plot centroid sizes
landmark.centroids %>%
  ggplot(aes(x = population, y = centroid, color=species)) +
  geom_point()+
  facet_grid(~ species, scale = "free", space = "free")

#read in genetic cluster data

cluster.dat <- read.table(file=list.files(pattern="structure",full.names = TRUE), stringsAsFactors = FALSE)
cluster.dat.reduced <- cluster.dat %>%
  filter(k.value.run == 2) %>%
  select(id,X1,X2,membership)

landmark.centroids.rename <- landmark.centroids
landmark.centroids.rename$id <- landmark.centroids$id %>% gsub("2014","",.) %>% gsub("_","",.)

cluster.morph <- left_join(cluster.dat.reduced, landmark.centroids.rename)

cluster.morph  %>%
  filter(!is.na(population)) %>%
  ggplot(aes(x = population, y = centroid, color=factor(membership))) +
  geom_point(size=3)+
  facet_grid(~ species, scale = "free", space = "free")


