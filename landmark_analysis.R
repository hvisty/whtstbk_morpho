##landmark analysis in Geomorph##

install.packages("geomorph")
install.packages("candisc")
install.packages("MASS")
install.packages("candisc")

library("geomorph")
library("ggplot2")
library("ggthemes")
library("candisc") # will use this later
library("MASS")
library("dplyr")
library("cats")

landmarks <- read.csv("Corrected Landmark Data.csv")

#Used to find outliers, plot any x/y points with indiv. ID#
#####
landmarks$ID <- paste(landmarks$population, landmarks$individual, sep="_")
ggplot(landmarks, aes(x=x.4, y=y.4, label=ID, color=population)) +geom_point() +geom_text()

#####

#Make a 3D array; only x/y cordinates, outliers removed...

landmarks <- landmarks[-c(10,125,271,253,327,360,448,544,630,641,665,877,882,887,898:900,947,959:961,997, 1061, 1172),]
landmarks <-landmarks[complete.cases(landmarks[,"x.1"]),]
landmarks.data <- landmarks[, c(5:42)]

landmark.array <- arrayspecs(landmarks.data, 19, 2) #3D array

GPA.landmarks <- gpagen(landmark.array) #Procrustes analysis

GPA.2D <- two.d.array(GPA.landmarks$coords) #2D Data frame of procrustes coordinates
  
#principle component stuff. 
PCA<- plotTangentSpace(GPA.landmarks$coords, groups= landmarks$population, verbose=T)
PCA$pc.scores #gives the PCA scores

#unbending; just show me PC2 vs PC3
col.gpa <- rainbow(length(levels(landmarks$species)))
names(col.gpa) <- levels(landmarks$species)
col.gpa <- col.gpa[match(landmarks$species, names(col.gpa))]

plotTangentSpace(GPA.landmarks$coords, groups=landmarks$ID, axis1=1, axis2=2)
plotTangentSpace(GPA.landmarks$coords, groups=col.gpa, axis1=2, axis2=3)
plotTangentSpace(GPA.landmarks$coords, groups=col.gpa, axis1=3, axis2=4)


#make a prettier plot! 
PC.scores<-as.data.frame(PCA$pc.scores)

ggplot(PC.scores, aes(x=PC2, y=PC3, color=landmarks$species, label=landmarks$ID))+ 
  geom_point() +
  geom_text()

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
                          ld2 = plda$x[,2], 
                          X = landmarks$X)
########
##use this data frame to look at lda values versus original x/y

lda.xy <- merge(lda.project, landmarks, by=c("population", "individual", "species", "ID", "X"))
lda.xy <- subset(lda.xy, select=-c(X)) 

#create a linear regression coefficient for each x:ld1 and y:ld2
lda.lm<-data.frame(matrix(NA, nrow = nrow(lda.xy) , ncol = length(PC.scores)))
colnames(lda.lm)<-c("LMx.1","LMx.2","LMx.3","LMx.4","LMx.5","LMx.6","LMx.7","LMx.8",
                    "LMx.9","LMx.10","LMx.11","LMx.12","LMx.13","LMx.14","LMx.15",
                    "LMx.16","LMx.17","LMx.18","LMx.19",
                    "LMy.1","LMy.2","LMy.3","LMy.4","LMy.5","LMy.6","LMy.7","LMy.8",
                    "LMy.9","LMy.10","LMy.11","LMy.12","LMy.13","LMy.14","LMy.15",
                    "LMy.16","LMy.17","LMy.18","LMy.19")
ld1 <- lda.xy$ld1
ld2 <- lda.xy$ld2

#separate loops for X and Y cause my looping skillz are weak

for (l in 1:19){
  pie <- as.matrix(lda.xy[paste("x", l, sep=".")])
    for(i in 1:length(pie)){ 
      lm.1<-lm(ld1 ~ pie)
      lda.lm[,l] <- lm.1$fitted.values
  }
}

for (l in 1:19){
  cake <- as.matrix(lda.xy[paste("y", l, sep=".")])
  for(i in 1:length(cake)){ 
    lm.2<-lm(ld2 ~ cake)
    lda.lm[,(l+19)] <- lm.2$fitted.values
  }
}


#######

# plot with ggplot
lda.project %>%
  ggplot(aes(color = species, x = ld1,y = ld2))+
    geom_point(size = 3) +
    labs(x = paste0("LD1 (", prop.lda[1], "%)"),
         y = paste0("LD2 (", prop.lda[2], "%)"))

# no "Both" locations
# ooo pretty
lda.project %>%
  filter(species != "B") %>%
  ggplot(aes(color = species, x = ld1,y = ld2))+
    geom_point(size = 3) +
    labs(x = paste0("LD1 (", prop.lda[1], "%)"),
         y = paste0("LD2 (", prop.lda[2], "%)")) +
  theme_classic() +
  guides(color=guide_legend(title="Species")) +
  scale_colour_manual(values=c("firebrick1", "cornflower blue"))+
  theme(axis.text=element_text(size=20), axis.title=element_text(size=22), 
        legend.text=element_text(size=20), legend.title=element_text(size=20))







####working on getting lda values
# so what do the lda axes mean?


# the loadings of the axes
lda.species$scaling

# make a data frame with their original names (x.1, y.1 ,etc.)
landmark.names <- names(landmarks)[-c(1:4)][-39]

lm.loadings <- data.frame(lm.name = landmark.names, 
           ld1.loading = lda.species$scaling[,1],
           ld2.loading = lda.species$scaling[,2])

head(lm.loadings)

######## GENETICS VS. MORPHOLOGY ########

# extract centroids
landmark.centroids <- data.frame(population = landmarks$population, 
                                 id = landmarks$ID, 
                                 species = landmarks$species, 
                                 centroid = GPA.landmarks$Csize)
# plot centroid sizes
landmark.centroids %>%
  ggplot(aes(x = population, y = centroid, color=species)) +
  geom_boxplot()+ xlab("Population") +ylab("Centroid") +
  scale_colour_manual(values=c("darkorchid4", "firebrick1", "cornflower blue")) +
  theme_classic() +
  guides(color=guide_legend(title="Species")) +
  theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=90, vjust=1),
        axis.title=element_text(size=22), strip.text=element_text(size=20), 
        legend.text=element_text(size=20), legend.title=element_text(size=20))+
  facet_grid(~ species, scale = "free", space = "free")

#read in genetic cluster data

cluster.dat <- read.table(file=list.files(pattern="structure",full.names = TRUE), stringsAsFactors = FALSE)
cluster.dat.reduced <- cluster.dat %>%
  filter(k.value.run == 2) %>%
  select(id,X1,X2,membership)

landmark.centroids.rename <- landmark.centroids
landmark.centroids.rename$id <- landmark.centroids$id %>% gsub("2014","",.) %>% gsub("_","",.)

# read in sex data

sex.dat <- read.csv(file=list.files(pattern="sex",full.names = TRUE), stringsAsFactors = FALSE, na.strings = "")
names(sex.dat) <- c("population","id","sex")
sex.dat$id <- sex.dat$id %>% gsub("2014","",.) 
sex.dat$population <- sex.dat$population %>% gsub("2014","",.) 
sex.dat$id <- paste0(sex.dat$population,sex.dat$id)

cluster.morph <- left_join(cluster.dat.reduced, landmark.centroids.rename)
cluster.morph <- left_join(cluster.morph, sex.dat)

cluster.morph  %>%
  filter(!is.na(population)) %>%
  ggplot(aes(x = population, y = centroid, color=factor(membership))) +
  geom_point(size=3)+
  theme_classic() +
  theme(axis.text=element_text(size=20), axis.text.x = element_text(angle=90, vjust=1),
        axis.title=element_text(size=22), strip.text=element_text(size=20), 
        legend.text=element_text(size=20), legend.title=element_text(size=20))+
  guides(color=guide_legend(title="Species")) +
  scale_colour_manual(values=c("firebrick1", "cornflower blue")) +
  facet_grid(~species, scale = "free", space = "free")


