######################
## Predicting Soil Property Mapping Workflow
## Random Forest script that includes:
## Extraction of covariates to points
## Point depth training set creation
## Prediction interval creation
## Relative Predition interval analysis
## Cross Validation
## Most steps parallelized
######################
## Authors: Travis Nauman, Suzann Kienast-Brown

######## Workspace Setup ########
## install packages if not already installed
required.packages <- c("raster", "sp", "rgdal", "randomForest", "snow", "snowfall", "quantregForest","dplyr", "ggplot2","hexbin","parallel", "aqp", "tidyr")

new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## increase active memory useable by raster package: Windows only
memory.limit(500000)
## raster settings: adjust based on system
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08, memfrac = 0.9)

## key folder locations
# travis
#predfolder <- "/home/tnaum/data/BLMsoils/Sand_2D_NASIS_SSURGO_SCD"
#covfolder <- "/home/tnaum/data/UCRB_Covariates"
# suz
predfolder <- "L:/methodology/predict"
covfolder <- "L:/methodology/covariates"
#Jess
predfolder <- "E:/DSMFocus/Properties/Methodology/pt_depth_predict"
covfolder <- "E:/DSMFocus/Properties/Methodology/covariates"



######## Prepare Property Data ########
## load Rdata file with lab data - soil profile collection
load("L:/methodology/NCSS_Lab_Data_Mart_20180914.Rdata")
#load("E:/DSMFocus/Properties/Methodology/R/NCSS_Lab_Data_Mart_20180914.Rdata")

## check it out
aqp::rebuildSPC(spc_access)
spc_access
str(spc_access)

## subset SPC into needed site and horizon data
sites <- site(spc_access)[ , c('pedon_key', 'pedlabsampnum', 'latitude_decimal_degrees', 'longitude_decimal_degrees')]
str(sites)
summary(sites)

horizons <- horizons(spc_access)[ , c('pedon_key', 'labsampnum', 'hzn_top', 'hzn_bot', 'clay_tot_psa', 'silt_tot_psa', 'sand_tot_psa', 'ph_h2o', 'oc')]
str(horizons)
summary(horizons)

## remove duplicates
horizons$test <- paste0(spc_access$pedon_key, spc_access$hzn_desgn, spc_access$hzn_top, spc_access$hzn_bot, sep = "-")
horizons <-horizons[!duplicated(horizons$test), ]

## merge site and horizon data frames by pedon key
site_hzn <- merge(sites, horizons, by = "pedon_key", all=TRUE)
str(site_hzn)
summary(site_hzn)

## remove NA from location 
pts <- site_hzn %>% drop_na(latitude_decimal_degrees, hzn_top, hzn_bot)
str(pts)
summary(pts)

### make property point data spatial 
## weed out points with imprecise coordinates (if needed) 
pts$latnchar <- nchar(abs(pts$latitude_decimal_degrees)) # update with correct field name
pts$longnchar <- nchar(abs(pts$longitude_decimal_degrees))
ptsc <- subset(pts, pts$latnchar > 5 & pts$longnchar > 6)
summary(ptsc)

## turn into spatial file
shp.pts <- ptsc
coordinates(shp.pts) <- ~ longitude_decimal_degrees + latitude_decimal_degrees # modify with coordinate field names in pts 
temp.proj <- CRS("+proj=longlat +datum=WGS84") # specify projection if needed
projection(shp.pts) <- temp.proj



######## Load Map Clip Boundary (if point dataset extends beyond project area) ########
setwd("L:/methodology/boundary")
polybound <- readOGR(".", "CO_River_watershed_Meade_alb")
polybound <- spTransform(polybound, temp.proj)

## clip points and check with visualization
shp.pts = shp.pts[polybound,]# clip by outer extent of all polybound features
summary(shp.pts)
plot(polybound)
plot(shp.pts, add=TRUE)



######## Grid Prep and Point Extract ########
## make list of grids
setwd(covfolder)
cov.grids <- list.files(pattern=".tif$")
cov.grids

### if points need to be matched up to grids ###
projgrid <- raster(cov.grids[1])
cov.proj <- projection(projgrid)
shp.pts <- spTransform(shp.pts, CRS(cov.proj)) # project to match rasters

## plot to ensure alignment bw points and rasters
plot(projgrid)
plot(shp.pts, add=TRUE)

### parallelized extract: (larger datasets)
cpus <- detectCores(all.tests = FALSE, logical = TRUE)-1
sfInit(parallel=TRUE, cpus=cpus)
sfExport("shp.pts", "cov.grids")
sfLibrary(raster)
sfLibrary(rgdal)
ov.lst <- sfLapply(cov.grids, function(i){try( raster::extract(raster(i), shp.pts) )})
snowfall::sfStop()
#endCluster()
ov.lst <- as.data.frame(ov.lst)
names(ov.lst) = tools::file_path_sans_ext(basename(cov.grids))
ov.lst$DID <- seq.int(nrow(ov.lst))
shp.pts$DID <- seq.int(nrow(shp.pts))
pts.ext <- merge(as.data.frame(shp.pts),ov.lst, by="DID")
summary(pts.ext)

## save points with extracted covariates (if desired)
setwd(predfolder)
saveRDS(pts.ext, "pts_ext_covs.rds")



######## Prep Training Data for Random Forest ########
pts.ext$prop <- pts.ext$oc # UPDATE EVERY TIME with proper response field
prop <- "oc" # dependent variable

## create location IDs for extra duplicate removal step
#####pts.ext$LocID <- paste(pts.ext$longitude_decimal_degrees, pts.ext$latitude_decimal_degrees, sep = "")



######## Loop to Train and Predict Properties for All Depths ########
#depths <- c(2.5,10,22.5,45,80,125) 
d <- 2.5 # if running or checking for one depth

#for(d in depths){
pts.extc <- subset(pts.ext, as.numeric(pts.ext$hzn_top) <= d & as.numeric(pts.ext$hzn_bot) > d) # subset to chosen depth
#####pedonLocations <- unique(pts.extc$LocID) # if length differs from # of rows in pts, there are duplicates
#####pts.extc <- subset(pts.extc, !duplicated(pts.extc[c("LocID")])) #removes duplicates
ptspred.list <- gsub(".tif","", cov.grids) # take .tif off of the grid list to just get the variable names
ptspred.list <- c(ptspred.list,"prop") # add dependent variable
pts.extc <- pts.extc[c(ptspred.list)] ## or create a specific list of dependent variable and covariate names to use 
pts.extcc <- na.omit(pts.extc) # remove any record with NA's (in any column - be careful)
xtrain <- as.matrix(pts.extcc[c(gsub(".tif","", cov.grids))])
ytrain <- c(as.matrix(pts.extcc[c("prop")]))

### transform property data if needed (oc, clay) ###
#hist(ytrain) # examine distribution to determine transformation
ytrain <- log(ytrain+0.1) # oc
#ytrain <- sqrt(ytrain) # clay
#hist(ytrain) # check transformed distribution

### for RPI
varrange <- as.numeric(quantile(pts.extcc$prop, probs=c(0.975), na.rm=T)-quantile(pts.extcc$prop, probs=c(0.025),na.rm=T)) 

##### build quantile Random Forest
set.seed(915)
#set.seed(NULL)
Qsoiclass <- quantregForest(x=xtrain, y=ytrain, importance=TRUE, ntree=120, keep.forest=TRUE) # if dataset is ~>300, nthreads = cpus
soiclass <- Qsoiclass
class(soiclass) <- "randomForest"

### linear adjustment for bias ###
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

varImpPlot(soiclass) # look at variable importance

### Save models
setwd(predfolder)
saveRDS(Qsoiclass, paste("Qsoiclass_RFmodel", prop, d, "cm_nasisSSURGO_ART_SG100.rds",sep="_"))
saveRDS(rf_lm_adj, paste("rflmadj_RFmodel", prop, d, "cm_nasisSSURGO_ART_SG100.rds",sep="_"))


### reference covariate rasters to use in prediction
setwd(covfolder)
rasterOptions(maxmemory = 1e+09,chunksize = 2e+08) # maxmemory = 1e+09,chunksize = 1e+08 for soilmonster
rasters <- stack(cov.grids)

##### predict onto covariate grid
setwd(predfolder)
closeAllConnections()
### parallelized predict
beginCluster(cpus,type='SOCK')
####predl <- clusterR(rasters, predict, args=list(model=Qsoiclass,what=c(0.025)),progress="text")
####predh <- clusterR(rasters, predict, args=list(model=Qsoiclass,what=c(0.975)),progress="text")
Sys.time()
pred <- clusterR(rasters, predict, args=list(model=soiclass),progress="text")
Sys.time()
names(pred) <- "trainpreds"

## linear adjustment
predlm <- clusterR(pred, predict, args=list(model=rf_lm_adj),progress="text")
### prediction interval widths
s <- stack(predh,predl)

####PIwidth.fn <- function(a,b) {
  ####ind <- a-b
  ####return(ind)
####}
####PIwidth <- clusterR(s, overlay, args=list(fun=PIwidth.fn),progress = "text")

## Determine 95% interquantile range of original training data for horizons that include the depth being predicted
####PIrelwidth.fn <- function(a,b) {
  ####ind <- (a-b)/varrange
  ####return(ind)
####}
####PIrelwidth <- clusterR(s, overlay, args=list(fun=PIrelwidth.fn),progress = "text",export='varrange')

### back transformation workflow 
bt.fn <- function(x) {
  ind <- (exp(x))-0.1 #If a backtransform is needed 10^(x) or exp(x) or ^2
  return(ind)
}
####predh_bt <- clusterR(predh, calc, args=list(fun=bt.fn),progress='text')
####predl_bt <- clusterR(predl, calc, args=list(fun=bt.fn),progress='text')
pred_bt <- clusterR(pred, calc, args=list(fun=bt.fn),progress='text')
####s_bt <- stack(predh_bt,predl_bt)
####PIwidth_bt.fn <- function(a,b) {
####ind <- a-b
####return(ind)
####}
#PIwidth_bt <- clusterR(s_bt, overlay, args=list(fun=PIwidth_bt.fn),progress = "text")
# ## If transformed, use the following code for PI width prep steps
####PIrelwidth_bt.fn <- function(a,b) {
   ####ind <- (a-b)/varrange
   ####return(ind)
####}
####PIrelwidth_bt <- clusterR(s_bt, overlay, args=list(fun=PIrelwidth_bt.fn),progress = "text", export='varrange')

endCluster()

##### write prediction results to new geotiff files
### untransformed code block
####writeRaster(predlm, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
####writeRaster(predl, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_l.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
####writeRaster(predh, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_h.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
####writeRaster(PIrelwidth, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_relwidth.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
# writeRaster(PIwidth, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_width.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")

### back-transformed code block
writeRaster(pred_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_bt_ART_SG100covs.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
####writeRaster(predl_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_l_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
####writeRaster(predh_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_h_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
####writeRaster(PIwidth_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_width_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
####writeRaster(PIrelwidth_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_relwidth_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")



######## Manual Cross validation ########
pts.extcvm <- pts.extcc
nfolds <- 10
pts.extcvm$folds <- sample.int(nfolds,size =length(pts.extcvm[,1]),replace=T)
pts.extcvm$prop_t <- log(pts.extcvm$prop+0.1) # oc # transform if needed else just create new version of prop
formulaStringCVm <- as.formula(paste('prop_t ~', paste(gsub(".tif","", cov.grids), collapse="+")))
#for (g in seq(nfolds)){
CV_factorRF <- function(g,pts.extcvm, formulaStringCVm){
  traindf <- subset(pts.extcvm, pts.extcvm$folds != g)
  testdf <- subset(pts.extcvm, pts.extcvm$folds == g)
  xtrain.t <- as.matrix(traindf[c(gsub(".tif","", cov.grids))])
  ytrain.t <- c(as.matrix(traindf$prop_t))
  rf.pcv <- quantregForest(x=xtrain.t, y=ytrain.t, importance=TRUE, ntree=120, keep.forest=TRUE,nthreads = max(cpus/nfolds,1))
  rf.pcvc <- rf.pcv
  class(rf.pcvc) <- "randomForest"
  traindf$pcvpredpre <- predict(rf.pcvc, newdata=traindf)
  testdf$pcvpredpre <- predict(rf.pcvc, newdata=testdf)
  testdf$pcvpredpre.025 <- predict(rf.pcv, newdata=testdf, what=c(0.025))
  testdf$pcvpredpre.975 <- predict(rf.pcv, newdata=testdf, what=c(0.975))
  attach(traindf)
  lm.pcv <- lm(prop_t~pcvpredpre)
  detach(traindf)
  testdf$pcvpred <- predict(lm.pcv, newdata=testdf)
  return(testdf)
}
snowfall::sfInit(parallel=TRUE, cpus=nfolds)
snowfall::sfExport("pts.extcvm","formulaStringCVm","CV_factorRF","cov.grids","cpus","nfolds")
snowfall::sfLibrary(randomForest)
snowfall::sfLibrary(quantregForest)
pts.extpcv <- snowfall::sfLapply(1:nfolds, function(g){CV_factorRF(g, pts.extcvm=pts.extcvm,formulaStringCVm=formulaStringCVm)})
snowfall::sfStop()
pts.extpcv <- plyr::rbind.fill(pts.extpcv)
pts.extpcv$pcvpred <- as.numeric(pts.extpcv$pcvpred)

### PCV statistics - untransformed
cvp.RMSE <- sqrt(mean((pts.extpcv$prop_t - pts.extpcv$pcvpred)^2, na.rm=TRUE))
cvp.Rsquared <- 1-var(pts.extpcv$prop_t - pts.extpcv$pcvpred, na.rm=TRUE)/var(pts.extpcv$prop_t, na.rm=TRUE)

## Back transformed: create pcvpred_bt even if not tranformed for cv.depth function
pts.extpcv$pcvpred_bt <- (exp(pts.extpcv$pcvpred))-0.1 #oc ## ALWAYS update with backtransformation (if needed)
cvp.RMSE_bt <- sqrt(mean((pts.extpcv$prop - pts.extpcv$pcvpred_bt)^2, na.rm=TRUE))
cvp.Rsquared_bt <- 1-var(pts.extpcv$prop - pts.extpcv$pcvpred_bt, na.rm=TRUE)/var(pts.extpcv$prop, na.rm=TRUE)

### RPI
pts.extpcv$prop_bt <- (exp(pts.extpcv$prop_t))-0.1 # UPDATE: backtransform if necessary. Used for PICP and to characterize backtransformation bias
pts.extpcv$pcvpredpre.025_bt <- (exp(pts.extpcv$pcvpredpre.025))-0.1 # UPDATE: backtransform if necessary
pts.extpcv$pcvpredpre.975_bt <- (exp(pts.extpcv$pcvpredpre.975))-0.1 # UPDATE: backtransform if necessary
pts.extpcv$abs.resid <- abs(pts.extpcv$prop - pts.extpcv$pcvpred_bt)
pts.extpcv$RPI <- (pts.extpcv$pcvpredpre.975_bt - pts.extpcv$pcvpredpre.025_bt)/varrange
plot(pts.extpcv$abs.resid~pts.extpcv$RPI) # Quick look at relationship

## summarize RPI and residuals
pts.extpcv$rel.abs.resid <- pts.extpcv$abs.resid/varrange
RPI.cvave <- mean(pts.extpcv$RPI)
RPI.cvmed <- median(pts.extpcv$RPI)
rel.abs.res.ave <- mean(pts.extpcv$rel.abs.resid)
rel.abs.res.med <- median(pts.extpcv$rel.abs.resid)
pts.extpcv$BTbias <- pts.extpcv$prop_bt - pts.extpcv$prop
BTbias.abs.max <- max(abs(pts.extpcv$BTbias))
BTbias.ave <- mean(pts.extpcv$BTbias)
PICP <- sum(ifelse(pts.extpcv$prop_bt <= pts.extpcv$pcvpredpre.975_bt & pts.extpcv$prop_bt >= pts.extpcv$pcvpredpre.025_bt,1,0))/length(pts.extpcv[,1])

## summarize RPI in full raster prediction using sample for speed (tested against full average)
rpi_samp <- sampleRegular(PIrelwidth,size=200000,useGDAL=T)
rpi_samp <- na.omit(rpi_samp)
gRPI.ave <- mean(rpi_samp)
gRPI.med <- median(rpi_samp)

### Create PCV table
##if gRPI was created
#CVdf <- data.frame(cvp.RMSE, cvp.Rsquared, cvp.RMSE_bt, cvp.Rsquared_bt,RPI.cvave,RPI.cvmed,gRPI.ave,gRPI.med,PICP,rel.abs.res.ave,rel.abs.res.med,BTbias.abs.max,BTbias.ave, nrow(pts.extcvm))
##if gRPI was not created
CVdf <- data.frame(cvp.RMSE, cvp.Rsquared, cvp.RMSE_bt, cvp.Rsquared_bt,RPI.cvave,RPI.cvmed,PICP,rel.abs.res.ave,rel.abs.res.med,BTbias.abs.max,BTbias.ave, nrow(pts.extcvm))
##if gRPI was created
#names(CVdf) <- c("cvp.RMSE","cvp.Rsquared","cvp.RMSE_bt", "cvp.Rsquared_bt","RPI.CVave","RPI.CVmed","gRPI.ave","gRPI.med","PICP","rel.abs.res.ave","rel.abs.res.med","BTbias.abs.max","BTbias.ave","obs.pts.extcvm")
##if gRPI was not created
names(CVdf) <- c("cvp.RMSE","cvp.Rsquared","cvp.RMSE_bt", "cvp.Rsquared_bt","RPI.CVave","RPI.CVmed","PICP","rel.abs.res.ave","rel.abs.res.med","BTbias.abs.max","BTbias.ave","obs.pts.extcvm")
## write table
setwd(predfolder)
write.table(CVdf, paste("PCVstats", prop, d, "cm_ptdepth_DC_UCRB.txt",sep="_"), sep = "\t", row.names = FALSE)

### CV 1:1 plot
viri <- c("#440154FF", "#39568CFF", "#1F968BFF", "#73D055FF", "#FDE725FF") # color ramp
gplt.dcm.2D.CV <- ggplot(data=pts.extpcv, aes(prop_t, pcvpred)) +
  stat_binhex(bins = 30) + geom_abline(intercept = 0, slope = 1,lwd=1)  + xlim(0,4) + ylim(0,4) +
  theme(axis.text=element_text(size=8), legend.text=element_text(size=10), axis.title=element_text(size=10),plot.title = element_text(size=10,hjust=0.5)) +
  xlab("Measured") + ylab("CV Prediction") + scale_fill_gradientn(name = "log(Count)", trans = "log", colours = rev(viri)) +
  ggtitle(paste("Cross val", prop, d, "cm",sep=" "))
gplt.dcm.2D.CV
## Save Cross validation graph and data for future plotting
setwd(predfolder)
ggsave(paste(prop, d, 'cm_CV_plot.tif', sep="_"), plot = gplt.dcm.2D.CV, device = "tiff", dpi = 600, limitsize = TRUE, width = 4, height = 4, units = 'in')
saveRDS(pts.extpcv, paste(prop, "cvlm_preds_2D", d, "cm_ptdepth_DC_UCRB.rds", sep="_"))
} # End of depth loop 

######## masking water pixels out ########
## Code for producting water mask from NLCD (if necessary)
# nlcd <- raster("/home/tnaum/data/UCRB_Covariates/NLCDcl.tif")
# beginCluster(30,type='SOCK')
## make a mask raster
# mask_fn <- function(nlcd){ind <- ifelse(nlcd!=11,1,NA)
#   return(ind)
# }
# mask <- clusterR(nlcd, calc, args=list(fun=mask_fn),progress='text')
# endCluster()
# plot(mask)
# writeRaster(mask, overwrite=TRUE,filename="/home/tnaum/data/BLMsoils/nlcd_watermask.tif", options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text",datatype='INT1U')
# rm(mask)
## now set up a list of rasters and function to mask out water
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



