####################################################
########## FCM calculations "Carbon sources" #######
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


baseFolder <- "/media/projects2/CristinaGT/RamanMP/fcs_carbonsource/"

a = list.files(baseFolder, pattern=".fcs")
Media_code<- unlist(lapply(strsplit(a,split="_"),function(x) (x[1])))
Media<- unlist(lapply(strsplit(Media_code,split=" "),function(x) (x[2])))
Time<- unlist(lapply(strsplit(a,split="_"),function(x) (x[2])))

cell.name=rep("",length(a))
for (i in 1:length(cell.name)) {
  cell.name[i] <-paste(Media_code[i],Time[i],Rep[i], sep= " ", collapse=NULL)
}

media_rep=rep("",length(a))
for (i in 1:length(media_rep)) {
  media_rep[i] <-paste(Media_code[i],Rep[i], sep= " ", collapse=NULL)
}

media_time=rep("",length(a))
for (i in 1:length(media_rep)) {
  media_time[i] <-paste(Media_code[i],Time[i], sep= " ", collapse=NULL)
}


## Fingerprint
set.seed(777)
flowData <- read.flowSet(path = "/media/projects2/CristinaGT/CollabMyrsini_April2019/FCS/", transformation = FALSE, pattern=".fcs")


flowData_transformed <- transform(flowData,`FL1-H`=asinh(`FL1-H`), 
                                  `SSC-H`=asinh(`SSC-H`), 
                                  `FL3-H`=asinh(`FL3-H`), 
                                  `FSC-H`=asinh(`FSC-H`))
param=c("FL1-H", "FL3-H","SSC-H","FSC-H")

#### Date and time of the measurements ####
Date <- c()
Volume <- c()
for(i in 1:length(flowData)){
  Date[i] <- flowData[[i]]@description$`$DATE`
  Volume[i] <- as.numeric(flowData[[i]]@description$`$VOL`)/1000 #volume in µL
}



### Create a PolygonGate for denoising the dataset
### Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
sqrcut1 <- matrix(c(7.75,7.75,14,14,3,8.15,16,3),ncol=2, nrow=4)
colnames(sqrcut1) <- c("FL1-H","FL3-H")
polyGate1 <- polygonGate(.gate=sqrcut1, filterId = "Total Cells")

###  Gating quality check
xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[83], filter=polyGate1,
       scales=list(y=list(limits=c(0,14)),
                   x=list(limits=c(6,16))),
       axis = axis.default, nbin=125, 
       par.strip.text=list(col="white", font=2, cex=2), smooth=FALSE)


##Isolate gate
Viacounts <- flowCore::filter(flowData_transformed, polyGate1)
Viacounts <- toTable(summary(Viacounts))

results_counts <- data.frame(Samples=flowCore::sampleNames(flowData_transformed), 
                             Total.cells = (Viacounts$true/as.numeric(Volume)))  #uL and no dilution

results_counts_FCM <- data.frame(Samples=flowCore::sampleNames(flowData_transformed), 
                                 Total.cells = (Viacounts$true/as.numeric(Volume)) *100 *10000)#mL and dilution
results_counts_FCM$media <- Media_code
results_counts_FCM$time <- Time
results_counts_FCM$rep <- Rep
results_counts_FCM$media_time <-media_time

