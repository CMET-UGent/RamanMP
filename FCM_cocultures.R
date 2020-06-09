####################################################
########### FCM calculations "Cocultures" ##########
########## Cristina García-Timermans ###############
####################################################

## Load packages
require('flowFDA')
require("vegan")
require("MESS")
require("flowFDA")
require("gridExtra")
require("boot")
library(Phenoflow)

#Upload fcs files
baseFolder <- "/media/projects2/CristinaGT/RamanMP/fcs_cocultures/"
labels = list.files(baseFolder, pattern=".fcs")

Strain_plus<- unlist(lapply(strsplit(labels,split="_"),function(x) (x[1])))
Treat<- unlist(lapply(strsplit(Strain_plus,split=" "),function(x) (x[2])))
Dil <- unlist(lapply(strsplit(labels,split="_"),function(x) (x[2])))
Strain<- unlist(lapply(strsplit(Strain_plus,split="_"),function(x) (x[1])))


cell.name=rep("",length(labels))
for (i in 1:length(cell.name)) {
  #cell.name[i] <-paste(Donor[i], Treatment[i], Antibiotic[i], Timepoint[i],  sep= " ", collapse=NULL)
  cell.name[i] <-paste(Strain[i], sep= " ", collapse=NULL)
}


## Fingerprint
set.seed(777)

path = "/media/projects2/CristinaGT/Xiaona_FCMfingerprint_HOBhet/fcs/"
flowData <- read.flowSet(path = "/media/projects2/CristinaGT/Xiaona_FCMfingerprint_HOBhet/fcs/", transformation = FALSE, pattern=".fcs")

flowData_transformed <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))

param=c("FL1-H", "FL3-H","SSC-H","FSC-H")
flowData_transformed = flowData_transformed[,param]
#remove(flowData)

#### Date and time of the measurements ####
Date <- c()
Volume <- c()
for(i in 1:length(flowData)){
  Date[i] <- flowData[[i]]@description$`$DATE`
  Volume[i] <- as.numeric(flowData[[i]]@description$`$VOL`)/1000 #(in µL)
}


### Create a PolygonGate for denoising the dataset
### Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(8.25,8.25,16,16,3.5,10.5,16.7,10),ncol=2, nrow=4)

colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")
#polyGate1 <- polygonGate

####  Gating quality check ####


library(flowViz)
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[2],filter=polyGate1,
       scales=list(y=list(limits=c(0.2,16)),
                   x=list(limits=c(0.3,16))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)


#### Isolate only the cellular information based on the polyGate1####
flowData_transformed <- Subset(flowData_transformed, polyGate1)

#Multiply the dilution

Dilution <- mapvalues(Dil, from= c("x1E3", "x1E4"), to = as.numeric(c("1000", "10000")))
Dilution_2 <- as.numeric(Viacounts$true/as.numeric(Volume))
Dilution_2 <- Dilution_2 * as.numeric(Dilution)

results_counts <- data.frame(Samples=Treat, 
                             Total.cells = Dilution_2 *1000) #to mL