######################
## Predicting Soil Property Mapping Workflow
## Random Forest script that includes:
## Extraction of covariates to points
## Point depth training set creation
## Prediction interval creation
## Relative Precition interval analysis
## Cross Validation
## Most steps parallelized
######################
## Authors: Travis Nauman...

# Workspace setup
# Install packages if not already installed
required.packages <- c("raster", "sp", "rgdal", "randomForest", "snow", "snowfall", "quantregForest","dplyr", "ggplot2","hexbin","parallel")# might need snowfall
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase actuve memory useable by raster package: Windows only
memory.limit(500000)
## Raster settings: adjust based on system
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08, memfrac = 0.8)

## Key Folder Locations
predfolder <- "/home/tnaum/data/DSM_DepthCompare"
covfolder <- "/home/tnaum/OneDrive/USGS/BLM_projects/BLM_CO_ESGs/ESGs_UCRB/ESG_30mCovsInt"

######## Load shapefile (if needed) ##############
setwd(predfolder)## FOlder with points
shp.pts <-readOGR(".", "ec_12pre_ncss_LIMS_UPCO")
point.proj <- projection(shp.pts)
# If no prj file and you know proj, can specify by name
temp.proj <- CRS("+proj=longlat +datum=WGS84")
projection(shp.pts) <- temp.proj


######## Get points for extraction if in table form ###########
setwd(predfolder)
pts <- read.delim("NCSS17_PSDA_rkFrags_ttab.txt") # If in delimited file other than csv
# pts <- readRDS("nasispts_gSSURGO18hor_ucrb_final.rds") # If in rds
### Weed out points with imprecise coordinates (if needed) ###
pts$latnchar <- nchar(abs(pts$ywgs84)) # update with correct field name
pts$longnchar <- nchar(abs(pts$xwgs84))
ptsc <- subset(pts, pts$latnchar > 5 & pts$longnchar > 6)
### Turn into spatial file
shp.pts <- ptsc
coordinates(shp.pts) <- ~ xwgs84 + ywgs84 # modify with field names in table
temp.proj <- CRS("+proj=longlat +datum=WGS84") ## specify projection if needed
projection(shp.pts) <- temp.proj


######## Load map clip boundary (if needed) ###########
setwd("/home/tnaum/Dropbox/USGS/BLM_projects/Utah_BLM_Salinity/Huc6_boundary")
polybound <- readOGR(".", "CO_River_watershed_Meade_alb")
polybound <- spTransform(polybound, temp.proj)
## Now clip points and check with visualization
shp.pts = shp.pts[polybound,]#clip by outer extent of all polybound features
plot(polybound)
plot(shp.pts, add=TRUE)

######### Grid Prep #################
## Make list of grids
setwd(covfolder)
cov.grids <- list.files(pattern=".tif$")
## If points need to be matched up to grids ###
projgrid <- raster(cov.grids[1])
cov.proj <- projection(projgrid)
shp.pts <- spTransform(shp.pts, CRS(cov.proj)) # project to match rasters

## Plot to ensure alignment bw points and rasters
plot(projgrid)
plot(shp.pts, add=TRUE)

## Parallelized extract: (larger datasets)
cpus <- detectCores(all.tests = FALSE, logical = TRUE)-1
sfInit(parallel=TRUE, cpus=cpus)
sfExport("shp.pts", "cov.grids")
sfLibrary(raster)
sfLibrary(rgdal)
ov.lst <- sfLapply(cov.grids, function(i){try( raster::extract(raster(i), shp.pts) )})
snowfall::sfStop()
ov.lst <- as.data.frame(ov.lst)
names(ov.lst) = tools::file_path_sans_ext(basename(cov.grids))
ov.lst$DID <- seq.int(nrow(ov.lst))
shp.pts$DID <- seq.int(nrow(shp.pts))
pts.ext <- merge(as.data.frame(shp.pts),ov.lst, by="DID")

## Save points with extracted covariates (if desired)
setwd(predfolder)
saveRDS(pts.ext, "pts_ext_covs.rds")

## Prep training data for Random Forest
pts.ext$prop <- pts.ext$sandtotal_r ## UPDATE EVERY TIME with proper response field
prop <- "sandtotal_r" ## Dependent variable

## Create location IDs for extra duplicate removal step
pts.ext$LocID <- paste(pts.ext$xwgs84, pts.ext$ywgs84, sep = "")

##### Loop to train and predict properties for all depths
depths <- c(0,5,15,30,60,100,150)
for(d in depths){
pts.extc <- subset(pts.ext, as.numeric(pts.ext$hzdept_r) <= d & as.numeric(pts.ext$hzdepb_r) > d) # subset to chosen depth
pedonLocations <- unique(pts.extc$LocID) # if length differs from # of rows in pts, there are duplicates
pts.extc <- subset(pts.extc, !duplicated(pts.extc[c("LocID")])) #removes duplicates
ptspred.list <- gsub(".tif","", cov.grids)# Take .tif off of the grid list to just get the variable names
ptspred.list <- c(ptspred.list,"prop") #Add dependent variable
pts.extc <- pts.extc[c(ptspred.list)]## Or create a specific list of dependent variable and covariate names to use 
pts.extcc <- na.omit(pts.extc)# Remove any record with NA's (in any column - be careful)
xtrain <- as.matrix(pts.extcc[c(gsub(".tif","", cov.grids))])
ytrain <- c(as.matrix(pts.extcc[c("prop")]))
# logytrain <- log(ytrain+0.1)
# sqrtytrain <- sqrt(ytrain)
varrange <- as.numeric(quantile(pts.extcc$prop, probs=c(0.975), na.rm=T)-quantile(pts.extcc$prop, probs=c(0.025),na.rm=T)) # For RPI

############### Build quantile Random Forest
Qsoiclass <- quantregForest(x=xtrain, y=ytrain, importance=TRUE, ntree=120, keep.forest=TRUE) # If dataset is ~>300: , nthreads = cpus)
soiclass <- Qsoiclass
class(soiclass) <- "randomForest"
## Linear Adjustment for bias
pts.extcc$trainpreds <- predict(soiclass, newdata=xtrain)
pts.extcc$prop_t <- pts.extcc$prop ## TRANSFORM IF NEEDED
attach(pts.extcc)
rf_lm_adj <- lm(prop_t ~ trainpreds)
detach(pts.extcc)
pts.extcc$trainpredsadj <- predict(rf_lm_adj, newdata=pts.extcc)
## plot linear adjustment
# x1 <-c(-100,0,100,10000,100000000)
# y1 <-c(-100,0,100,10000,100000000)
# lines(x1,y1, col = 'red')#1:1 line
# plot(ytrain~predict(soiclass, newdata=xtrain)) #Fit plot
# lines(x1,y1, col = 'red')#1:1 line
# plot(ytrain~pts.extcc$trainpredsadj) #Fit plot
# lines(x1,y1, col = 'red')#1:1 line
varImpPlot(soiclass) # look at variable imporance
## Save models
setwd(predfolder)
saveRDS(Qsoiclass, paste("Qsoiclass_RFmodel", prop, d, "cm_nasisSSURGO_ART_SG100.rds",sep="_"))
saveRDS(rf_lm_adj, paste("rflmadj_RFmodel",prop, d, "cm_nasisSSURGO_ART_SG100.rds",sep="_"))


## Reference covar rasters to use in prediction
setwd(covfolder)
rasterOptions(maxmemory = 1e+09,chunksize = 2e+08)# maxmemory = 1e+09,chunksize = 1e+08 for soilmonster
rasters <- stack(cov.grids)

## Predict onto covariate grid
setwd(predfolder)
## Parallelized predict
beginCluster(cpus,type='SOCK')
predl <- clusterR(rasters, predict, args=list(model=Qsoiclass,what=c(0.025)),progress="text")
predh <- clusterR(rasters, predict, args=list(model=Qsoiclass,what=c(0.975)),progress="text")
Sys.time()
pred <- clusterR(rasters, predict, args=list(model=soiclass),progress="text")
Sys.time()
names(pred) <- "trainpreds"
## Linear Adjustment
predlm <- clusterR(pred, predict, args=list(model=rf_lm_adj),progress="text")
## PI widths
s <- stack(predh,predl)
PIwidth.fn <- function(a,b) {
  ind <- a-b
  return(ind)
}
PIwidth <- clusterR(s, overlay, args=list(fun=PIwidth.fn),progress = "text")
# Determine 95% interquantile range of original training data for horizons that include the depth being predicted
PIrelwidth.fn <- function(a,b) {
  ind <- (a-b)/varrange
  return(ind)
}
PIrelwidth <- clusterR(s, overlay, args=list(fun=PIrelwidth.fn),progress = "text",export='varrange')
## Back transformation workflow 
# bt.fn <- function(x) {
#   ind <- (exp(x))-0.1 #If a backtransform is needed 10^(x) or exp(x) or ^2
#   return(ind)
# }
# predh_bt <- clusterR(predh, calc, args=list(fun=bt.fn),progress='text')
# predl_bt <- clusterR(predl, calc, args=list(fun=bt.fn),progress='text')
# pred_bt <- clusterR(pred, calc, args=list(fun=bt.fn),progress='text')
# s_bt <- stack(predh_bt,predl_bt)
# PIwidth_bt.fn <- function(a,b) {
#   ind <- a-b
#   return(ind)
# }
# PIwidth_bt <- clusterR(s_bt, overlay, args=list(fun=PIwidth_bt.fn),progress = "text")
# ## If transformed, use the following code for PI width prep steps
# PIrelwidth_bt.fn <- function(a,b) {
#   ind <- (a-b)/varrange
#   return(ind)
# }
# PIrelwidth_bt <- clusterR(s_bt, overlay, args=list(fun=PIrelwidth_bt.fn),progress = "text", export='varrange')

endCluster()
## Write new geotiff files
setwd(predfolder)
## Untranformed code block
writeRaster(predlm, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
writeRaster(predl, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_l.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
writeRaster(predh, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_h.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
writeRaster(PIrelwidth, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_relwidth.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
# writeRaster(PIwidth, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_width.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
# ## Transformed code block
# writeRaster(pred_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_bt_ART_SG100covs.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
# writeRaster(predl_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_l_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
# writeRaster(predh_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_h_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
# writeRaster(PIwidth_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_width_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
# writeRaster(PIrelwidth_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_relwidth_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")



################### Manual Cross validation ################################
pts.extcvm <- pts.extcc
nfolds <- 10
pts.extcvm$folds <- sample.int(nfolds,size =length(pts.extcvm[,1]),replace=T)
pts.extcvm$prop_t <- pts.extcvm$prop ## UPDATE: tranform if needed else just create new version of prop
formulaStringCVm <- as.formula(paste('prop_t ~', paste(gsub(".tif","", cov.grids), collapse="+")))
#for (g in seq(nfolds)){
CV_factorRF <- function(g,pts.extcvm, formulaStringCVm){
  traindf <- subset(pts.extcvm, pts.extcvm$folds != g)
  testdf <- subset(pts.extcvm, pts.extcvm$folds == g)
  xtrain.t <- as.matrix(traindf[c(gsub(".tif","", cov.grids))])
  ytrain.t <- c(as.matrix(traindf$prop_t))
  rf.pcv <- quantregForest(x=xtrain.t, y=ytrain.t, importance=TRUE, ntree=120, keep.forest=TRUE,nthreads = max(cpus/10,1))
  rf.pcvc <- rf.pcv
  class(rf.pcvc) <- "randomForest"
  traindf$pcvpredpre <- predict(rf.pcvc, newdata=traindf)
  testdf$pcvpredpre <- predict(rf.pcvc, newdata=testdf)
  #traindf$pcvpredpre <- predict(rf.pcv, newdata=traindf, what=c(0.5)) ## If median is desired
  #testdf$pcvpredpre <- predict(rf.pcv, newdata=testdf,, what=c(0.5)) ## If median is desired
  testdf$pcvpredpre.025 <- predict(rf.pcv, newdata=testdf, what=c(0.025))
  testdf$pcvpredpre.975 <- predict(rf.pcv, newdata=testdf, what=c(0.975))
  attach(traindf)
  lm.pcv <- lm(prop_t~pcvpredpre)
  detach(traindf)
  testdf$pcvpred <- predict(lm.pcv, newdata=testdf)
  return(testdf)
}
snowfall::sfInit(parallel=TRUE, cpus=nfolds)
snowfall::sfExport("pts.extcvm","formulaStringCVm","CV_factorRF","cov.grids")
snowfall::sfLibrary(randomForest)
snowfall::sfLibrary(quantregForest)
pts.extpcv <- snowfall::sfLapply(1:nfolds, function(g){CV_factorRF(g, pts.extcvm=pts.extcvm,formulaStringCVm=formulaStringCVm)})
snowfall::sfStop()
pts.extpcv <- plyr::rbind.fill(pts.extpcv)
pts.extpcv$pcvpred <- as.numeric(pts.extpcv$pcvpred)
## PCV statistics
cvp.RMSE <- sqrt(mean((pts.extpcv$prop_t - pts.extpcv$pcvpred)^2, na.rm=TRUE))
cvp.Rsquared <- 1-var(pts.extpcv$prop_t - pts.extpcv$pcvpred, na.rm=TRUE)/var(pts.extpcv$prop_t, na.rm=TRUE)
## Back transformed: create pcvpred_bt even if not tranformed for cv.depth function
pts.extpcv$pcvpred_bt <- pts.extpcv$pcvpred ## ALWAYS Update with backtransformation (if needed)
cvp.RMSE_bt <- sqrt(mean((pts.extpcv$prop - pts.extpcv$pcvpred_bt)^2, na.rm=TRUE))
cvp.Rsquared_bt <- 1-var(pts.extpcv$prop - pts.extpcv$pcvpred_bt, na.rm=TRUE)/var(pts.extpcv$prop, na.rm=TRUE)
## RPI
pts.extpcv$prop_bt <- pts.extpcv$prop_t # UPDATE: backtransform if necessary. Used for PICP and to characterize backtransformation bias
pts.extpcv$pcvpredpre.025_bt <- pts.extpcv$pcvpredpre.025 # UPDATE: backtransform if necessary
pts.extpcv$pcvpredpre.975_bt <- pts.extpcv$pcvpredpre.975 # UPDATE: backtransform if necessary
pts.extpcv$abs.resid <- abs(pts.extpcv$prop - pts.extpcv$pcvpred_bt)
pts.extpcv$RPI <- (pts.extpcv$pcvpredpre.975_bt - pts.extpcv$pcvpredpre.025_bt)/varrange
# plot(pts.extpcv$abs.resid~pts.extpcv$RPI) # Quick look at relationship
## Summarize RPI and residuals
pts.extpcv$rel.abs.resid <- pts.extpcv$abs.resid/varrange
RPI.cvave <- mean(pts.extpcv$RPI)
RPI.cvmed <- median(pts.extpcv$RPI)
rel.abs.res.ave <- mean(pts.extpcv$rel.abs.resid)
rel.abs.res.med <- median(pts.extpcv$rel.abs.resid)
pts.extpcv$BTbias <- pts.extpcv$prop_bt - pts.extpcv$prop
BTbias.abs.max <- max(abs(pts.extpcv$BTbias))
BTbias.ave <- mean(pts.extpcv$BTbias)
PICP <- sum(ifelse(pts.extpcv$prop_bt <= pts.extpcv$pcvpredpre.975_bt & pts.extpcv$prop_bt >= pts.extpcv$pcvpredpre.025_bt,1,0))/length(pts.extpcv[,1])
## Summarize RPI in full raster prediction using sample for speed (tested against full average)
rpi_samp <- sampleRegular(relpi,size=200000,useGDAL=T)
rpi_samp <- na.omit(rpi_samp)
gRPI.ave <- mean(rpi_samp)
gRPI.med <- median(rpi_samp)
## Create PCV table
CVdf <- data.frame(cvp.RMSE, cvp.Rsquared, cvp.RMSE_bt, cvp.Rsquared_bt, n_scd,RPI.cvave,RPI.cvmed,gRPI.ave,gRPI.med,PICP,rel.abs.res.ave,rel.abs.res.med,BTbias.abs.max,BTbias.ave)
names(CVdf) <- c("cvp.RMSE","cvp.Rsquared","cvp.RMSE_bt", "cvp.Rsquared_bt","RPI.CVave","RPI.CVmed","gRPI.ave","gRPI.med","PICP","rel.abs.res.ave","rel.abs.res.med","BTbias.abs.max","BTbias.ave")
setwd(predfolder)
write.table(CVdf, paste("PCVstats", prop, d, "cm_nasisSSURGO_ART_SG100.txt",sep="_"), sep = "\t", row.names = FALSE)
## CV 1:1 plot
viri <- c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF", "#FDE725FF") # color ramp
gplt.dcm.2D.CV <- ggplot(data=pts.extpcv, aes(prop_t, pcvpred)) +
  stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1)  + xlim(0,100) + ylim(0,100) +
  theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
  xlab("Measured") + ylab("CV Prediction") + scale_fill_gradientn(name = "log(Count)", trans = "log", colours = rev(viri)) +
  ggtitle(paste("Cross val", prop, d, "cm",sep=" "))
gplt.dcm.2D.CV
## Save Cross validation graph and data for future plotting
setwd(predfolder)
ggsave(paste(prop,'_CV_plot.tif',sep=""), plot = gplt.dcm.2D.CV, device = "tiff", dpi = 600, limitsize = TRUE, width = 4, height = 4, units = 'in')
saveRDS(pts.extpcv, paste(prop, "cvlm_preds_2D", d, "cm_nasisSSURGO_ART_SG100.rds", sep="_"))
} # End of depth loop

############# Masking water pixels out ############
## Code for producting water mask from NLCD (if necessary)
# nlcd <- raster("/home/tnaum/data/UCRB_Covariates/NLCDcl.tif")
# beginCluster(30,type='SOCK')
# # Make a mask raster
# mask_fn <- function(nlcd){ind <- ifelse(nlcd!=11,1,NA)
#   return(ind)
# }
# mask <- clusterR(nlcd, calc, args=list(fun=mask_fn),progress='text')
# endCluster()
# plot(mask)
# writeRaster(mask, overwrite=TRUE,filename="/home/tnaum/data/BLMsoils/nlcd_watermask.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text",datatype='INT1U')
# rm(mask)
## Now set up a list of rasters and function to mask out water
rasterOptions(maxmemory = 1e+09,chunksize = 1e+08)
setwd("/home/tnaum/data/BLMsoils/Sand_2D_NASIS_SSURGO_SCD")
grids <- list.files(pattern=".tif$")
mskfn <- function(rast,mask){
  ind <- rast*mask
  ind[ind<0]<-0 # to bring the slighly negative predictions back to zero (if appropriate for property)
  return(ind)
}
## par list apply fn
watermask_fn <- function(g){
  setwd(predfolder)
  rast <- raster(g)
  names(rast) <- "rast"
  mask <- raster("/home/tnaum/data/BLMsoils/nlcd_watermask.tif") # path to water mask
  h2ostk <- stack(rast,mask)
  setwd(paste(predfolder,"/masked",sep=""))
  overlay(h2ostk,fun=mskfn,progress='text',filename=g, options=c("COMPRESS=DEFLATE", "TFW=YES"))
  gc()
}
snowfall::sfInit(parallel=TRUE, cpus=cpus)
snowfall::sfExport("watermask_fn","mskfn","grids")
snowfall::sfLibrary(raster)
Sys.time()
snowfall::sfLapply(grids, function(g){watermask_fn(g)})
Sys.time()
snowfall::sfStop()



