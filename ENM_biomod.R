library("sp")
library("raster")
library("dplyr")
library("rgdal")
library("fuzzySim")
library("rgeos")
library("magrittr")
library("tools")
library("readr")
library("dismo")
library("biomod2")

setwd("C:/")

####DataFormating ####
shapePath <- 'C:/WWF_ecoregions'
shapeLayer <- "wwf_terr_ecos"
regionalizacion <- rgdal::readOGR(shapePath, shapeLayer)

#Current climate
bioclima<-stack(list.files(path="C:/", pattern = "*.tif$", full.names=TRUE)) 

#Future climate
#--------------------------------- BCC_SSP245
covarDataFolder_T1<-stack(list.files(path = "C:/T1", pattern="*.tif$", full.names=TRUE))
covarDataFolder_T2<-stack(list.files(path = "C:/T2", pattern="*.tif$", full.names=TRUE))
covarDataFolder_T3<-stack(list.files(path = "C:/T3", pattern="*.tif$", full.names=TRUE))

#--------------------------------- Can_SSP245
covarDataFolder_T1<-stack(list.files(path = "C:/T1", pattern="*.tif$", full.names=TRUE))
covarDataFolder_T2<-stack(list.files(path = "C:/T2", pattern="*.tif$", full.names=TRUE))
covarDataFolder_T3<-stack(list.files(path = "C:/T3", pattern="*.tif$", full.names=TRUE))

#--------------------------------- BCC_SSP585
covarDataFolder_T1<-stack(list.files(path = "C:/T1", pattern="*.tif$", full.names=TRUE))
covarDataFolder_T2<-stack(list.files(path = "C:/T2", pattern="*.tif$", full.names=TRUE))
covarDataFolder_T3<-stack(list.files(path = "C:/T2", pattern="*.tif$", full.names=TRUE))

#--------------------------------- Can_SSP585
covarDataFolder_T1<-stack(list.files(path = "C:/T1", pattern="*.tif$", full.names=TRUE))
covarDataFolder_T2<-stack(list.files(path = "C:/T2", pattern="*.tif$", full.names=TRUE))
covarDataFolder_T3<-stack(list.files(path = "C:/T3", pattern="*.tif$", full.names=TRUE))

#
args <- list.files("C:/", pattern = "*.csv$",full.names = TRUE) #A folder with all species csv (We run this code in a cluster with slurm)


inputDataFile <- args[3]
outputFolder1 <- inputDataFile %>%
  basename %>%
  file_path_sans_ext;outputFolder1

outputFolder<- gsub(" ", ".", outputFolder1);outputFolder

if (!dir.exists(outputFolder)) {
  dir.create(outputFolder, recursive = TRUE)
}

# Extract envorimental varibales with species occurrences
crs.wgs84 <- sp::CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
occsData <- readr::read_csv(inputDataFile)#working directory must be the same of your csv

#Filter species with more than 250 (this could be change if your region is bigger)
if (dim(occsData)[1] > 250) {occsData <- occsData[sample(1:dim(occsData)[1], 250, replace=FALSE),] }

sp::coordinates(occsData) <- c("x", "y") #Change if the columns are name differently
sp::proj4string(occsData) <- crs.wgs84

covarData <- raster::extract(bioclima, occsData)
covarData <- cbind(occsData, covarData)

completeDataCases <- covarData@data %>%
  dplyr::select(.dots=names(bioclima)) %>%
  complete.cases
covarData <- covarData[completeDataCases, ]

####Variables selection####
speciesCol <- match(outputFolder, names(occsData))
varCols <- ncol(occsData) + 1

correlacion <- corSelect(
  data = covarData@data,
  sp.cols = speciesCol,
  var.cols = varCols:ncol(covarData),
  cor.thresh = 0.8,
  use = "pairwise.complete.obs"
)

select_var <- correlacion$selected.vars
write(select_var, file = file.path(outputFolder, "selected_variables.txt"))

#Raster covariables selected for model calibration
enviromentalVariables <- bioclima[[select_var]]

# Raster covariables selected for model calibration
selectedVariables <- enviromentalVariables[[select_var]]

# Raster covariables selected for model transfer
#seleccion variables CC

#--------------------------------- 1_45: BCC_SSP245
env_1_45_T1t<-covarDataFolder_T1[[select_var]]
env_1_45_T2t<-covarDataFolder_T2[[select_var]]
env_1_45_T3t<-covarDataFolder_T3[[select_var]]
#--------------------------------- 2_45: Can_SSP245
env_2_45_T1t<-covarDataFolder_T1[[select_var]]
env_2_45_T2t<-covarDataFolder_T2[[select_var]]
env_2_45_T3t<-covarDataFolder_T3[[select_var]]
#--------------------------------- 1_85: BCC_SSP585
env_1_85_T1t<-covarDataFolder_T1[[select_var]]
env_1_85_T2t<-covarDataFolder_T2[[select_var]]
env_1_85_T3t<-covarDataFolder_T3[[select_var]]
#--------------------------------- 2_85: Can_SSP585
env_2_85_T1t<-covarDataFolder_T1[[select_var]]
env_2_85_T2t<-covarDataFolder_T2[[select_var]]
env_2_85_T3t<-covarDataFolder_T3[[select_var]]

# M####
# Intersects the occurrence data with polygons
ecoregionsOfInterest <- sp::over(occsData, regionalizacion) %>%
  filter(!is.na(ECO_ID))

idsEcoRegions <- unique(ecoregionsOfInterest$ECO_ID)
polygonsOfInterest <- regionalizacion[regionalizacion$ECO_ID %in% idsEcoRegions, ]
pts_b <- gBuffer(occsData, width=4)#This can be change 
pts_b <- as(pts_b, 'SpatialPolygonsDataFrame')

#Poligono area de calibracion
polygonsOfInterest<-gIntersection(pts_b, polygonsOfInterest, drop_lower_td = T)
polygonsOfInterest<-gBuffer(polygonsOfInterest, width=0.5)
plot(polygonsOfInterest)

# Variables ambientales
# Cortar bioclimas con mascara M
selectedVariablesCrop <- raster::crop(selectedVariables, polygonsOfInterest)
myExpl <- raster::mask(selectedVariablesCrop,polygonsOfInterest) #Species variables delimited by M
myExpl<-stack(myExpl)

# Mask future raster with ecoregions of interest
#BCC---------------------------------------------
#1_45_T1
env_1_45_T1a <- raster::crop(env_1_45_T1t, polyTransferencia)
env_1_45_T1 <- raster::mask(env_1_45_T1a ,  polyTransferencia)
env_1_45_T1<-stack(env_1_45_T1)
rm(env_1_45_T1t, env_1_45_T1a)
#1_45_T2
env_1_45_T2a <- raster::crop(env_1_45_T2t, polyTransferencia)
env_1_45_T2 <- raster::mask(env_1_45_T2a ,  polyTransferencia)
env_1_45_T2<-stack(env_1_45_T2)
rm(env_1_45_T2t, env_1_45_T2a)
#1_45_T3
env_1_45_T3a <- raster::crop(env_1_45_T3t, polyTransferencia)
env_1_45_T3 <- raster::mask(env_1_45_T3a ,  polyTransferencia)
env_1_45_T3<-stack(env_1_45_T3)
rm(env_1_45_T3t, env_1_45_T3a)

#1_85_T1
env_1_85_T1a <- raster::crop(env_1_85_T1t, polyTransferencia)
env_1_85_T1 <- raster::mask(env_1_85_T1a ,  polyTransferencia)
env_1_85_T1<-stack(env_1_85_T1)
rm(env_1_85_T1t, env_1_85_T1a)
#1_85_T2
env_1_85_T2a <- raster::crop(env_1_85_T2t, polyTransferencia)
env_1_85_T2 <- raster::mask(env_1_85_T2a ,  polyTransferencia)
env_1_85_T2<-stack(env_1_85_T2)
rm(env_1_85_T2t, env_1_85_T2a)
#1_85_T3
env_1_85_T3a <- raster::crop(env_1_85_T3t, polyTransferencia)
env_1_85_T3 <- raster::mask(env_1_85_T3a ,  polyTransferencia)
env_1_85_T3<-stack(env_1_85_T3)
rm(env_1_85_T3t, env_1_85_T3a)

#Can---------------------------------------------
#2_45_T1
env_2_45_T1a <- raster::crop(env_2_45_T1t, polyTransferencia)
env_2_45_T1 <- raster::mask(env_2_45_T1a ,  polyTransferencia)
env_2_45_T1<-stack(env_2_45_T1)
rm(env_2_45_T1t, env_2_45_T1a)
#2_45_T2
env_2_45_T2a <- raster::crop(env_2_45_T2t, polyTransferencia)
env_2_45_T2 <- raster::mask(env_2_45_T2a ,  polyTransferencia)
env_2_45_T2<-stack(env_2_45_T2)
rm(env_2_45_T2t, env_2_45_T2a)
#2_45_T3
env_2_45_T3a <- raster::crop(env_2_45_T3t, polyTransferencia)
env_2_45_T3 <- raster::mask(env_2_45_T3a ,  polyTransferencia)
env_2_45_T3<-stack(env_2_45_T3)
rm(env_2_45_T3t, env_2_45_T3a)

#2_85_T1
env_2_85_T1a <- raster::crop(env_2_85_T1t, polyTransferencia)
env_2_85_T1 <- raster::mask(env_2_85_T1a ,  polyTransferencia)
env_2_85_T1<-stack(env_2_85_T1)
rm(env_2_85_T1t, env_2_85_T1a)
#2_85_T2
env_2_85_T2a <- raster::crop(env_2_85_T2t, polyTransferencia)
env_2_85_T2 <- raster::mask(env_2_85_T2a ,  polyTransferencia)
env_2_85_T2<-stack(env_2_85_T2)
rm(env_2_85_T2t, env_2_85_T2a)
#2_85_T3
env_2_85_T3a <- raster::crop(env_2_85_T3t, polyTransferencia)
env_2_85_T3 <- raster::mask(env_2_85_T3a ,  polyTransferencia)
env_2_85_T3<-stack(env_2_85_T3)
rm(env_2_85_T3t, env_2_85_T3a)

#DataFormating####
presencias<-data.frame(occsData)
names(presencias)[names(presencias)=="outputFolder"] <- outputFolder
presencias<-dplyr::select(presencias, c("x", "y",outputFolder))
names(presencias)

## BIOMOD####
DataSpecies<-presencias

myRespName <- outputFolder
myResp <- as.numeric(DataSpecies[,myRespName])
myRespCoord = DataSpecies[c("x", "y")]
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespCoord,
                                     resp.name = myRespName,
                                     PA.nb.rep = 5,
                                     PA.nb.absences = 10000,
                                     PA.strategy = 'random')

myBiomodData
#plot(myBiomodData)

### Parametrizacion
myBiomodOption <-BIOMOD_ModelingOptions()

### Modelling
myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('GLM', 'GBM', 'GAM', 'CTA', 'ANN','RF'),
  models.options = myBiomodOption,
  NbRunEval=10,
  DataSplit=70,
  models.eval.meth = c('KAPPA','TSS','ROC'),
  SaveObj = TRUE,
  rescal.all.models = FALSE,
  do.full.models = FALSE,
  modeling.id = paste(outputFolder))

#myBiomodModelOut
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
write.csv(myBiomodModelEval, file = file.path(outputFolder, "myBiomodModelEval.csv"),
          row.names = FALSE)

### Hacer predicciones sobre el raster
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl,
  proj.name = outputFolder,
  selected.models = 'all',
  binary.meth = "TSS",
  compress = 'xz',
  build.clamping.mask = TRUE,
  output.format = '.grd')

plot(myBiomodProj)

### ensemble_modeling#####
myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  em.by='all',
  eval.metric = c('ROC'),
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  eval.metric.quality.threshold = c(0.7),
  prob.mean.weight = T, 
  VarImport = 10)

myVarImportEM<-data.frame(get_variables_importance(myBiomodEM))
myVarImportEM<-myVarImportEM[5]
write.csv(myVarImportEM, file = file.path(outputFolder, "myVarImportEM.csv"),
          row.names = T)

# print summary
#myBiomodEM

# get evaluation scores
myBiomodEMEval<-get_evaluations(myBiomodEM)
write.csv(myBiomodEMEval, file = file.path(outputFolder, "myBiomodEMEval.csv"),
          row.names = FALSE)

####Current####
# Creating ensembles projections
myBiomodEM_proj <-BIOMOD_EnsembleForecasting(EM.output  = myBiomodEM,
                                             projection.output = myBiomodProj,
                                             selected.models = 'all',
                                             proj.name = outputFolder,
                                             binary.meth = "TSS")

outpath<-file.path(outputFolder,paste0("proj_",outputFolder))
currentPred <- stack(file.path(outputFolder,paste("proj_",outputFolder,"/proj_",outputFolder,"_",outputFolder,"_ensemble_TSSbin.grd",sep="")))
writeRaster(currentPred,
            file.path(outpath,"TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

currentPred_c <- stack(file.path(outputFolder,paste("proj_",outputFolder,"/proj_",outputFolder,"_",outputFolder,"_ensemble.grd",sep="")))
writeRaster(currentPred_c,
            file.path(outpath,"c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)


## Prohection to future scenarios ####
###_1_45_T1: BCC####
### Projection to future rasters
myBiomodEM_env_1_45_T1<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                        new.env = env_1_45_T1,
                                                        selected.models = 'all',
                                                        proj.name = "env_1_45_T1",
                                                        binary.meth = "TSS")

myBiomodEM_env_1_45_T1
plot(myBiomodEM_GreenControl_T1)

FutureProj <- stack(file.path(outputFolder, paste("proj_env_1_45_T1/proj_env_1_45_T1_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_env_1_45_T1/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_env_1_45_T1/proj_env_1_45_T1_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_env_1_45_T1/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_env_1_45_T1)

###_1_45_T2: BCC####
myBiomodEM_env_1_45_T2<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_1_45_T2,
                                                    selected.models = 'all',
                                                    proj.name = "env_1_45_T2",
                                                    binary.meth = "TSS")

myBiomodEM_env_1_45_T2
plot(myBiomodEM_GreenControl_T2)

FutureProj <- stack(file.path(outputFolder, paste("proj_env_1_45_T2/proj_env_1_45_T2_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_env_1_45_T2/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_env_1_45_T2/proj_env_1_45_T2_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_env_1_45_T2/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_env_1_45_T2)

###_1_45_T3: BCC####
myBiomodEM_env_1_45_T3<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_1_45_T3,
                                                    selected.models = 'all',
                                                    proj.name = "env_1_45_T3",
                                                    binary.meth = "TSS")

myBiomodEM_env_1_45_T3
plot(myBiomodEM_GreenControl_T3)

FutureProj <- stack(file.path(outputFolder, paste("proj_env_1_45_T3/proj_env_1_45_T3_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_env_1_45_T3/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_env_1_45_T3/proj_env_1_45_T3_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_env_1_45_T3/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_env_1_45_T3)


###_1_85_T1: BCC####
myBiomodEM_env_1_85_T1<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_1_85_T1,
                                                    selected.models = 'all',
                                                    proj.name = "env_1_85_T1",
                                                    binary.meth = "TSS")

myBiomodEM_env_1_85_T1
plot(myBiomodEM_GreenControl_T1)

FutureProj <- stack(file.path(outputFolder, paste("proj_env_1_85_T1/proj_env_1_85_T1_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_env_1_85_T1/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_env_1_85_T1/proj_env_1_85_T1_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_env_1_85_T1/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_env_1_85_T1)

###_1_85_T2: BCC####
myBiomodEM_env_1_85_T2<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_1_85_T2,
                                                    selected.models = 'all',
                                                    proj.name = "env_1_85_T2",
                                                    binary.meth = "TSS")

myBiomodEM_env_1_85_T2
plot(myBiomodEM_GreenControl_T2)

FutureProj <- stack(file.path(outputFolder, paste("proj_env_1_85_T2/proj_env_1_85_T2_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_env_1_85_T2/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_env_1_85_T2/proj_env_1_85_T2_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_env_1_85_T2/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_env_1_85_T2)

###_1_85_T3: BCC####
myBiomodEM_env_1_85_T3<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_1_85_T3,
                                                    selected.models = 'all',
                                                    proj.name = "env_1_85_T3",
                                                    binary.meth = "TSS")

myBiomodEM_env_1_85_T3
plot(myBiomodEM_GreenControl_T3)

FutureProj <- stack(file.path(outputFolder, paste("proj_env_1_85_T3/proj_env_1_85_T3_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_env_1_85_T3/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_env_1_85_T3/proj_env_1_85_T3_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_env_1_85_T3/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_env_1_85_T3)

###_2_45_T1: CAN####
myBiomodEM_env_2_45_T1<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_2_45_T1,
                                                    selected.models = 'all',
                                                    proj.name = "env_2_45_T1",
                                                    binary.meth = "TSS")

myBiomodEM_env_2_45_T1
plot(myBiomodEM_GreenControl_T1)

FutureProj <- stack(file.path(outputFolder, paste("proj_env_2_45_T1/proj_env_2_45_T1_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_env_2_45_T1/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_env_2_45_T1/proj_env_2_45_T1_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_env_2_45_T1/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_env_2_45_T1)

###_2_45_T2: CAN####
myBiomodEM_env_2_45_T2<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_2_45_T2,
                                                    selected.models = 'all',
                                                    proj.name = "env_2_45_T2",
                                                    binary.meth = "TSS")

myBiomodEM_env_2_45_T2
plot(myBiomodEM_GreenControl_T2)

FutureProj <- stack(file.path(outputFolder, paste("proj_env_2_45_T2/proj_env_2_45_T2_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_env_2_45_T2/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_env_2_45_T2/proj_env_2_45_T2_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_env_2_45_T2/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_env_2_45_T2)

###_2_45_T3: CAN####
myBiomodEM_env_2_45_T3<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_2_45_T3,
                                                    selected.models = 'all',
                                                    proj.name = "env_2_45_T3",
                                                    binary.meth = "TSS")

myBiomodEM_env_2_45_T3
plot(myBiomodEM_GreenControl_T3)

FutureProj <- stack(file.path(outputFolder, paste("proj_env_2_45_T3/proj_env_2_45_T3_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_env_2_45_T3/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_env_2_45_T3/proj_env_2_45_T3_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_env_2_45_T3/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_env_2_45_T3)


###_2_85_T1: CAN####
myBiomodEM_env_2_85_T1<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_2_85_T1,
                                                    selected.models = 'all',
                                                    proj.name = "env_2_85_T1",
                                                    binary.meth = "TSS")

myBiomodEM_env_2_85_T1
plot(myBiomodEM_GreenControl_T1)

FutureProj <- stack(file.path(outputFolder, paste("proj_env_2_85_T1/proj_env_2_85_T1_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_env_2_85_T1/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_env_2_85_T1/proj_env_2_85_T1_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_env_2_85_T1/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_env_2_85_T1)

###_2_85_T2: CAN####
myBiomodEM_env_2_85_T2<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_2_85_T2,
                                                    selected.models = 'all',
                                                    proj.name = "env_2_85_T2",
                                                    binary.meth = "TSS")

myBiomodEM_env_2_85_T2
plot(myBiomodEM_GreenControl_T2)

FutureProj <- stack(file.path(outputFolder, paste("proj_env_2_85_T2/proj_env_2_85_T2_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_env_2_85_T2/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_env_2_85_T2/proj_env_2_85_T2_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_env_2_85_T2/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_env_2_85_T2)

###_2_85_T3: CAN####
myBiomodEM_env_2_85_T3<-BIOMOD_EnsembleForecasting( EM.output  = myBiomodEM,
                                                    new.env = env_2_85_T3,
                                                    selected.models = 'all',
                                                    proj.name = "env_2_85_T3",
                                                    binary.meth = "TSS")

myBiomodEM_env_2_85_T3
plot(myBiomodEM_GreenControl_T3)

FutureProj <- stack(file.path(outputFolder, paste("proj_env_2_85_T3/proj_env_2_85_T3_",outputFolder,"_ensemble_TSSbin.grd", sep="")))
writeRaster(FutureProj,
            file.path(outputFolder,"proj_env_2_85_T3/TSSbin.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

FutureProj_c <- stack(file.path(outputFolder, paste("proj_env_2_85_T3/proj_env_2_85_T3_",outputFolder,"_ensemble.grd", sep="")))
writeRaster(FutureProj_c,
            file.path(outputFolder,"proj_env_2_85_T3/c.tif"),
            suffix='names',
            bylayer=TRUE,
            overwrite= TRUE)

rm(FutureProj, FutureProj_c, myBiomodEM_env_2_85_T3)
