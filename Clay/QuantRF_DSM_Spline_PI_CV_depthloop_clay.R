######################
## Predicting Soil Property Mapping Workflow
## Random Forest script that includes:
## Extraction of covariates to points
## Spline depth interval normalization
## Prediction interval creation
## Relative Prediction Interval Analysis
## Cross Validation
## Most steps parallelized
######################
## Authors: Travis Nauman, Jessica Philippe

######## Workspace setup ########
# Install packages if not already installed

required.packages <- c("raster", "sp", "rgdal", "randomForest", "snow", "snowfall", "quantregForest","dplyr", "tidyr", "ggplot2","hexbin","plyr","aqp","GSIF","parallel","snow")
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(required.packages, require, character.only=T)
rm(required.packages, new.packages)
## Increase active memory useable by raster package: Windows only
memory.limit(500000)
## Raster settings: adjust based on system
rasterOptions(maxmemory = 1e+09, chunksize = 1e+08, memfrac = 0.9)

## Key Folder Locations
predfolder <- "E:/DSMFocus/Properties/Methodology/spline_predict"
covfolder <- "E:/DSMFocus/Properties/Methodology/covariates"
testfolder <- "E:/DSMFocus/Properties/Methodology/test"


######## Prepare Property Data ########

## Load RData file with lab data
load(file = "E:/DSMFocus/Properties/Methodology/R/NCSS_Lab_Data_Mart_20180914.RData")

## check it out
spc_access
str(spc_access)

## subset SPC into needed site and horizon data
site_sub <- site(spc_access)[ ,c('pedon_key','latitude_decimal_degrees','longitude_decimal_degrees')]
summary(site_sub)

horizon_sub <- horizons(spc_access)[ ,c('pedon_key','labsampnum','hzn_top','hzn_bot','clay_tot_psa','silt_tot_psa','sand_tot_psa','ph_h2o','oc')]
summary(horizon_sub)

## merge site and horizon data frames by pedon key
horz_site <- merge(horizon_sub, site_sub, by = "pedon_key", all=TRUE)
str(horz_site)
summary(horz_site) #check for NAs

## remove NAs from location
pts <- horz_site %>% drop_na(latitude_decimal_degrees, hzn_top, hzn_bot) #remove NA locations
summary(pts)


## make property point data spatial
## Weed out points with imprecise coordinates (if needed) 
pts$latnchar <- nchar(abs(pts$latitude_decimal_degrees)) # update with correct field name
pts$longnchar <- nchar(abs(pts$longitude_decimal_degrees))
ptsc <- subset(pts, pts$latnchar > 5 & pts$longnchar > 6)
summary(ptsc)

## Turn into spatial file
shp.pts <- ptsc
coordinates(shp.pts) <- ~ longitude_decimal_degrees + latitude_decimal_degrees # modify with field names in table
temp.proj <- CRS("+proj=longlat +datum=WGS84") ## specify projection if needed
projection(shp.pts) <- temp.proj

######## Load map clip boundary (if needed) ###########
setwd("E:/DSMFocus/Properties/Methodology/boundary")
polybound <- readOGR(".", "CO_River_watershed_Meade_alb")
polybound <- spTransform(polybound, temp.proj)
## clip points and check with visualization
shp.pts <- shp.pts[polybound,]#clip by outer extent of all polybound features
plot(polybound)
plot(shp.pts, add=TRUE)

######## Grid Prep and extract values to points ########
## Make list of grids
setwd(covfolder)
cov.grids <- list.files(pattern=".tif$")
cov.grids
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


######## Prep training data for Random Forest ########
pts.ext$prop <- pts.ext$clay_tot_psa ## UPDATE EVERY TIME with proper response field
prop <- "clay_tot_psa" ## Dependent variable

## Create pedon ids for spline
pts.ext$pid <- pts.ext$pedon_key
## Create location IDs for duplicate removal step
pts.ext$LocID <- paste(pts.ext$longitude_decimal_degrees, pts.ext$latitude_decimal_degrees, sep = "")

## Combine datasets and remove NAs for spline
ptspred.list <- gsub(".tif","", cov.grids)# Take .tif off of the grid list to just get the variable names
ptspred.list <- c(ptspred.list,"prop","pid","top","bottom") #Add dependent variable
pts.ext$top <- pts.ext$hzn_top
pts.ext$bottom <- pts.ext$hzn_bot
pts.extc <- pts.ext[c(ptspred.list)]## Or create a specific list of dependent variable and covariate names to use 
pts.ext.com <- na.omit(pts.extc)# Remove any record with NA's (in any column - be careful)
pts.ext.com$hzid <- paste(pts.ext.com$pid,pts.ext.com$top,sep="_")

## Run spline to create standard depth interval average
pts.soilpr <- pts.ext.com
# pts.soilpr <- ## NEED TO GET NAs out of prop field
depths(pts.soilpr)<- pid ~ top + bottom # turns object into soil profile collection
# coordinates(pts.soilpr) <- ~ coords.x1 +coords.x2 # Doesn't work without organizing by site first
## Run Spline
pts.soilpr.spl <- mpspline(pts.soilpr, var.name="prop",vlow = min(pts.ext.com$prop),vhigh = max(pts.ext.com$prop),d = t(c(0,5,15,30,60,100,150))) # set vlow and vhigh based on range of each property

#create dataframe with unique pedon ID for covariates
pts.pid <- pts.ext.com[,c("pid",gsub(".tif","", cov.grids))]
pts.pid <- unique(pts.pid)

######## Loop to train and predict properties for all depths ########
depths <- c("0-5 cm","5-15 cm","15-30 cm","30-60 cm","60-100 cm","100-150 cm")
d <- c("15-30 cm") #if running one depth

for(d in depths){
pts.extcc <- data.frame(pid=pts.soilpr.spl$idcol, prop=pts.soilpr.spl$var.std[,d], stringsAsFactors = F) #change object as needed

pts.extcc <- left_join(pts.extcc,pts.pid, by="pid") #left join to pts.extcc by pedon id


## Clean up from merge (if needed)
####pts.extcc$prop <- pts.extcc$prop.x ####commented out because it got rid of "prop" field
#pts.extcc$prop.1 <- NULL
#pts.extcc$prop.x <- NULL
#pts.extcc$pid.y <- NULL
pts.extcc <- na.omit(pts.extcc) #remove any record with NA's (in any column - be careful)
xtrain <- as.matrix(pts.extcc[c(gsub(".tif","", cov.grids))])
ytrain <- c(as.matrix(pts.extcc[c("prop")]))

## transform property data if needed ##
hist(ytrain) #examine distribution to determine transformation
#ytrain <- log(ytrain+0.1) #oc
ytrain <- sqrt(ytrain) #clay
hist(ytrain) # check transformed distribution

varrange <- as.numeric(quantile(pts.extcc$prop, probs=c(0.975), na.rm=T)-quantile(pts.extcc$prop, probs=c(0.025),na.rm=T)) # For RPI

## Build quantile Random Forest
set.seed(915)
Qsoiclass <- quantregForest(x=xtrain, y=ytrain, importance=TRUE, ntree=120, keep.forest=TRUE) # If dataset is ~>300: , nthreads = cpus)
soiclass = randomForest(prop ~ ., data = pts.extcc, importance=TRUE, proximity=FALSE, ntree=100, keep.forest=TRUE)
soiclass <- Qsoiclass
class(soiclass) <- "randomForest"
soiclass ## Get oob error

## Linear Adjustment for bias
pts.extcc$trainpreds <- predict(soiclass, newdata=xtrain)
pts.extcc$prop_t <- pts.extcc$prop ## TRANSFORM IF NEEDED
attach(pts.extcc)
rf_lm_adj <- lm(prop_t ~ trainpreds)
detach(pts.extcc)

pts.extcc$trainpredsadj <- predict(rf_lm_adj, newdata=pts.extcc)
##plot linear adjustment
 #x1 <-c(-100,0,100,10000,100000000)
 #y1 <-c(-100,0,100,10000,100000000)
 #lines(x1,y1, col = 'red')#1:1 line
 #plot(ytrain~predict(soiclass, newdata=xtrain)) #Fit plot
 #lines(x1,y1, col = 'red')#1:1 line
 #plot(ytrain~pts.extcc$trainpredsadj) #Fit plot
# lines(x1,y1, col = 'red')#1:1 line

varImpPlot(soiclass)# look at variable importance

##save models
setwd(predfolder)
saveRDS(Qsoiclass, paste("Qsoiclass_RFmodel", prop, d, "cm_nasisSSURGO_ART_SG100.rds",sep="_"))
saveRDS(rf_lm_adj, paste("rflmadj_RFmodel",prop, d, "cm_nasisSSURGO_ART_SG100.rds",sep="_"))

## Reference covariate rasters to use in prediction
setwd(covfolder)
rasterOptions(maxmemory = 1e+09,chunksize = 2e+08)# maxmemory = 1e+09,chunksize = 1e+08 for soilmonster
rasters<-stack(cov.grids)


## Predict onto covariate grid
setwd(testfolder)
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
## Prediction Interval widths
s <- stack(predh,predl)
PIwidth.fn <- function(a,b) {
  ind <- a-b
  return(ind)
}
PIwidth <- clusterR(s, overlay, args=list(fun=PIwidth.fn),progress = "text")
## Determine 95% interquantile range of original training data for horizons that include the depth being predicted
PIrelwidth.fn <- function(a,b) {
  ind <- (a-b)/varrange
  return(ind)
}
PIrelwidth <- clusterR(s, overlay, args=list(fun=PIrelwidth.fn),progress = "text",export='varrange')

## Back transformation workflow
bt.fn <- function(x) {
   ind <- (x)^2 #If a backtransform is needed 10^(x) or exp(x) or ^2
   return(ind)
 }
 predh_bt <- clusterR(predh, calc, args=list(fun=bt.fn),progress='text')
 predl_bt <- clusterR(predl, calc, args=list(fun=bt.fn),progress='text')
 pred_bt <- clusterR(pred, calc, args=list(fun=bt.fn),progress='text')
 s_bt <- stack(predh_bt,predl_bt)
 PIwidth_bt.fn <- function(a,b) {
   ind <- a-b
   return(ind)
 }
 PIwidth_bt <- clusterR(s_bt, overlay, args=list(fun=PIwidth_bt.fn),progress = "text")
 ## If transformed, use the following code for PI width prep steps
 PIrelwidth_bt.fn <- function(a,b) {
   ind <- (a-b)/varrange
   return(ind)
 }
 PIrelwidth_bt <- clusterR(s_bt, overlay, args=list(fun=PIrelwidth_bt.fn),progress = "text", export='varrange')

endCluster()

## Write predictions to new geotiff files
setwd(predfolder)

## Untransformed code block
#writeRaster(predlm, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
#writeRaster(predl, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_l.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
#writeRaster(predh, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_h.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
#writeRaster(PIrelwidth, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_relwidth.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
#writeRaster(PIwidth, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_width.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")

## Transformed code block
 writeRaster(pred_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_bt_ART_SG100covs.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
 writeRaster(predl_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_l_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
 writeRaster(predh_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_h_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
 writeRaster(PIwidth_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_width_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")
 writeRaster(PIrelwidth_bt, overwrite=TRUE,filename=paste(prop,d,"cm_2D_QRF_95PI_relwidth_bt.tif",sep="_"), options=c("COMPRESS=DEFLATE", "TFW=YES"), progress="text")



################### Manual Cross validation ################################
pts.extcvm <- pts.extcc #new copy of pts object
nfolds <- 10 # numer of folds for cross validation
pts.extcvm$folds <- sample.int(nfolds,size =length(pts.extcvm[,1]),replace=T) #random samples - random assignment of 1-10 for every point
pts.extcvm$prop_t <- pts.extcvm$prop ## UPDATE: transform if needed else just create new version of prop
formulaStringCVm <- as.formula(paste('prop_t ~', paste(gsub(".tif","", cov.grids), collapse="+"))) #cross-validation version of formula for random forest
#for (g in seq(nfolds)){
CV_factorRF <- function(g,pts.extcvm, formulaStringCVm){ #initiates cross-validation function - does the same thing for each 10 fold. Parallel list apply for each fold
    traindf <- subset(pts.extcvm, pts.extcvm$folds != g) #if g is fold 1, training will be other 9 sets, test will be fold 1
  testdf <- subset(pts.extcvm, pts.extcvm$folds == g)
  xtrain.t <- as.matrix(traindf[c(gsub(".tif","", cov.grids))]) #create covariate matrix
  ytrain.t <- c(as.matrix(traindf$prop_t)) #response matrix
  rf.pcv <- quantregForest(x=xtrain.t, y=ytrain.t, importance=TRUE, ntree=120, keep.forest=TRUE,nthreads = max(cpus/nfolds,1)) #quant RF model trained on training set
  rf.pcvc <- rf.pcv
  class(rf.pcvc) <- "randomForest"
  traindf$pcvpredpre <- predict(rf.pcvc, newdata=traindf) 
  testdf$pcvpredpre <- predict(rf.pcvc, newdata=testdf)
  testdf$pcvpredpre.025 <- predict(rf.pcv, newdata=testdf, what=c(0.025)) #creating prediction intervals
  testdf$pcvpredpre.975 <- predict(rf.pcv, newdata=testdf, what=c(0.975)) #creating prediction intervals
  attach(traindf) #train linear adjustment model on training set
  lm.pcv <- lm(prop_t~pcvpredpre)
  detach(traindf)
  testdf$pcvpred <- predict(lm.pcv, newdata=testdf) #apply linear adjustment model to test set
  return(testdf)
}
snowfall::sfInit(parallel=TRUE, cpus=nfolds) #run the function with parallel list apply
snowfall::sfExport("pts.extcvm","formulaStringCVm","CV_factorRF","cov.grids","nfolds","cpus")
snowfall::sfLibrary(randomForest)
snowfall::sfLibrary(quantregForest)
pts.extpcv <- snowfall::sfLapply(1:nfolds, function(g){CV_factorRF(g, pts.extcvm=pts.extcvm,formulaStringCVm=formulaStringCVm)})
snowfall::sfStop()
pts.extpcv <- plyr::rbind.fill(pts.extpcv) #pull rows of data out of list format into data frame format
pts.extpcv$pcvpred <- as.numeric(pts.extpcv$pcvpred) #turn prediction into numeric

## PCV statistics
cvp.RMSE <- sqrt(mean((pts.extpcv$prop_t - pts.extpcv$pcvpred)^2, na.rm=TRUE)) #create objects that can be used to display results
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
plot(pts.extpcv$abs.resid~pts.extpcv$RPI) # Quick look at relationship

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
#rpi values from cross validation vs summary of rendered map
rpi_samp <- sampleRegular(PIrelwidth,size=200000,useGDAL=T) #take sample of whole raster, summarize values, estimate global average
rpi_samp <- na.omit(rpi_samp)
gRPI.ave <- mean(rpi_samp) #global RPI - estimating map RPI - lots of computing
gRPI.med <- median(rpi_samp)

##Create PCV table 
##if gRPI was created
#CVdf <- data.frame(cvp.RMSE, cvp.Rsquared, cvp.RMSE_bt, cvp.Rsquared_bt, n_scd,RPI.cvave,RPI.cvmed,gRPI.ave,gRPI.med,PICP,rel.abs.res.ave,rel.abs.res.med,BTbias.abs.max,BTbias.ave)
#names(CVdf) <- c("cvp.RMSE","cvp.Rsquared","cvp.RMSE_bt", "cvp.Rsquared_bt","RPI.CVave","RPI.CVmed","gRPI.ave","gRPI.med","PICP","rel.abs.res.ave","rel.abs.res.med","BTbias.abs.max","BTbias.ave")
##if gRPI was not created
CVdf <- data.frame(cvp.RMSE, cvp.Rsquared, cvp.RMSE_bt, cvp.Rsquared_bt,RPI.cvave,RPI.cvmed,PICP,rel.abs.res.ave,rel.abs.res.med,BTbias.abs.max,BTbias.ave)
names(CVdf) <- c("cvp.RMSE","cvp.Rsquared","cvp.RMSE_bt", "cvp.Rsquared_bt","RPI.CVave","RPI.CVmed","PICP","rel.abs.res.ave","rel.abs.res.med","BTbias.abs.max","BTbias.ave")
##write table
setwd(predfolder)
write.table(CVdf, paste("PCVstats", prop, d, "cm_spline_DC_UCRB.txt",sep="_"), sep = "\t", row.names = FALSE)

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
ggsave(paste(prop, d,'_CV_plot.tif',sep="_"), plot = gplt.dcm.2D.CV, device = "tiff", dpi = 600, limitsize = TRUE, width = 4, height = 4, units = 'in')
saveRDS(pts.extpcv, paste(prop, "cvlm_preds_2D", d, "cm_spline_DC_UCRB.rds", sep="_"))
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




