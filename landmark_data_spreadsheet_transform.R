##Putting Landmark Data into master spreadsheet##

#get a list of files
landmark.data <- list.files("F2 landmark data", full.names=T) 


#looping basics

landmark.data.frame <- data.frame(row.names=landmark.data) 

for (i in 1:length(landmark.data)){
  
  file.x <- read.table(file= landmark.data[i], colClasses = c("NULL", NA, "NULL", "NULL", "NULL", "NULL"), header = TRUE)
  file.y <- read.table(file= landmark.data[i], colClasses = c("NULL", "NULL", NA, "NULL", "NULL", "NULL"), header = TRUE)
  landmark.row <- c(t(file.x), t(file.y))
    landmark.data.frame <- rbind(landmark.data.frame, landmark.row)
  
}


#name the columns, split the ID into same format as master spreadsheet. 

landmark.data.frame$ID <- list.files("F2 landmark data")
names(landmark.data.frame) <- c("lndmrk.x1", "lndmrk.x2", "lndmrk.x3", "lndmrk.x4", "lndmrk.x5", "lndmrk.x6", "lndmrk.x7", "lndmrk.x8", "lndmrk.x9", "lndmrk.x10", "lndmrk.x11", "lndmrk.x12", "lndmrk.x13", "lndmrk.x14", "lndmrk.x15", "lndmrk.x16", "lndmrk.x17", "lndmrk.x18", "lndmrk.x19",
                                "lndmrk.y1", "lndmrk.y2", "lndmrk.y3", "lndmrk.y4", "lndmrk.y5", "lndmrk.y6", "lndmrk.y7", "lndmrk.y8", "lndmrk.y9", "lndmrk.y10", "lndmrk.y11", "lndmrk.y12", "lndmrk.y13", "lndmrk.y14", "lndmrk.y15", "lndmrk.y16", "lndmrk.y17", "lndmrk.y18", "lndmrk.y19",
                                "ID")       

landmark.data.frame <- landmark.data.frame[,c(39, 1:38)]

landmark.data.frame$ID <- gsub("_.txt", "", landmark.data.frame$ID)

landmark.data.frame$family <- as.character(lapply(strsplit(landmark.data.frame$ID, split="-"), "[", 2))
landmark.data.frame$individual <- as.character(lapply(strsplit(landmark.data.frame$ID, split="-"), "[", 3))

landmark.data.frame <- landmark.data.frame[,c(1, 40, 41, 2:39)]


#Merge with master data sheet, re-order columns 

master.data <- read.csv("F2 white stickleback morphometrics.csv")

data <- merge(landmark.data.frame, master.data, by="ID", all.y=T)



#make a landmark data sheet that is corrected for standard pixels and ready for Geomorph

data$x.1 <- data$lndmrk.x1 / data$pixel.correction
data$y.1 <- data$lndmrk.y1 / data$pixel.correction

data$x.2 <- data$lndmrk.x2 / data$pixel.correction
data$y.2 <- data$lndmrk.y2 / data$pixel.correction

data$x.3 <- data$lndmrk.x3 / data$pixel.correction
data$y.3 <- data$lndmrk.y3 / data$pixel.correction

data$x.4 <- data$lndmrk.x4 / data$pixel.correction
data$y.4 <- data$lndmrk.y4 / data$pixel.correction

data$x.5 <- data$lndmrk.x5 / data$pixel.correction
data$y.5 <- data$lndmrk.y5 / data$pixel.correction

data$x.6 <- data$lndmrk.x6 / data$pixel.correction
data$y.6 <- data$lndmrk.y6 / data$pixel.correction

data$x.7 <- data$lndmrk.x7 / data$pixel.correction
data$y.7 <- data$lndmrk.y7 / data$pixel.correction

data$x.8 <- data$lndmrk.x8 / data$pixel.correction
data$y.8 <- data$lndmrk.y8 / data$pixel.correction

data$x.9 <- data$lndmrk.x9 / data$pixel.correction
data$y.9 <- data$lndmrk.y9 / data$pixel.correction

data$x.10 <- data$lndmrk.x10 / data$pixel.correction
data$y.10 <- data$lndmrk.y10 / data$pixel.correction

data$x.11 <- data$lndmrk.x11 / data$pixel.correction
data$y.11 <- data$lndmrk.y11 / data$pixel.correction

data$x.12 <- data$lndmrk.x12 / data$pixel.correction
data$y.12 <- data$lndmrk.y12 / data$pixel.correction

data$x.13 <- data$lndmrk.x13 / data$pixel.correction
data$y.13 <- data$lndmrk.y13 / data$pixel.correction

data$x.14 <- data$lndmrk.x14 / data$pixel.correction
data$y.14 <- data$lndmrk.y14 / data$pixel.correction

data$x.15 <- data$lndmrk.x15 / data$pixel.correction
data$y.15 <- data$lndmrk.y15 / data$pixel.correction

data$x.16 <- data$lndmrk.x16 / data$pixel.correction
data$y.16 <- data$lndmrk.y16 / data$pixel.correction

data$x.17 <- data$lndmrk.x17 / data$pixel.correction
data$y.17 <- data$lndmrk.y17 / data$pixel.correction

data$x.18 <- data$lndmrk.x18 / data$pixel.correction
data$y.18 <- data$lndmrk.y18 / data$pixel.correction

data$x.19 <- data$lndmrk.x19 / data$pixel.correction
data$y.19 <- data$lndmrk.y19 / data$pixel.correction


data <- data[,-c(4:41 )]
data <- data[,-c(4:9 )]

write.csv(data, "F2 white stickleback morphometrics.csv")



