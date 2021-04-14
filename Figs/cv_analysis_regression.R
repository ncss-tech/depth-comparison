##########
## Simple regression analysis
## of Cross Valid. Results
##########

library(dplyr)

# data
cvdf <- read.csv("cv_analysis2.csv",stringsAsFactors = F)
cvdf$sym <- as.numeric(as.factor(cvdf$Prop))

## Plot based on all data
dev.off() # Clear plot
tiff("CV_Rsq_vs_Ndiffall.tif",width = 5, height = 5, units = 'in', res = 200 )
par(mar = c(5, 4, 4, 2),mgp=c(2.5,1,0))
lm_all <- lm(Rsq_diff~N_diff,data=cvdf)
plot(Rsq_diff~N_diff,data=cvdf,pch = as.character(cvdf$sym),xlab="Fractional Difference in Sample size",ylab=expression(paste("Change in   ", R^2,sep="")),ylim=c(-0.12,0.12))
#grid(NULL,NULL, lty = 2,lwd=0.3, col = "black")
x <- data.frame(N_diff=c(-0.25,0.25))
pred <- unname(predict(lm_all,x))
y <- c(0,0)
lines(x$N_diff,pred,lwd=1,lty=2)
lines(x$N_diff,y,lwd=2)
legend(-0.005, 0.125,c("1 - Clay", "2 - OC", "3 - pH", "4 - Sand","5 - Silt"),title="Properties",cex=0.8)#T - Texture, F - Fragments, D - Dynamic, C- Chemistry, H - Hydrologic, R - Restriction Depth
dev.off()

## Plot based on groupings
# By property
dev.off() # Clear plot
tiff("CV_Rsq_vs_Ndiffall_byProp.tif",width = 6.5, height = 4.5, units = 'in', res = 200 )
par(mar = c(4, 4, 1.5, 1),mgp=c(2.5,1,0),mfrow=c(1,5),mfrow=c(2,3))
plot(Rsq_diff~N_diff,data=cvdf[cvdf$sym==1,],pch = as.character(cvdf$Layer),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Change in   ", R^2,sep="")),ylim=c(-0.12,0.12),main=cvdf[cvdf$sym==1,]$Prop[1])
lines(x$N_diff,y,lwd=2)
plot(Rsq_diff~N_diff,data=cvdf[cvdf$sym==2,],pch = as.character(cvdf$Layer),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Change in   ", R^2,sep="")),ylim=c(-0.12,0.12),main=cvdf[cvdf$sym==2,]$Prop[1])
lines(x$N_diff,y,lwd=2)
plot(Rsq_diff~N_diff,data=cvdf[cvdf$sym==3,],pch = as.character(cvdf$Layer),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Change in   ", R^2,sep="")),ylim=c(-0.12,0.12),main=cvdf[cvdf$sym==3,]$Prop[1])
lines(x$N_diff,y,lwd=2)
plot(Rsq_diff~N_diff,data=cvdf[cvdf$sym==4,],pch = as.character(cvdf$Layer),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Change in   ", R^2,sep="")),ylim=c(-0.12,0.12),main=cvdf[cvdf$sym==4,]$Prop[1])
lines(x$N_diff,y,lwd=2)
plot(Rsq_diff~N_diff,data=cvdf[cvdf$sym==5,],pch = as.character(cvdf$Layer),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Change in   ", R^2,sep="")),ylim=c(-0.12,0.12),main=cvdf[cvdf$sym==5,]$Prop[1])
lines(x$N_diff,y,lwd=2)
plot(Rsq_diff~N_diff,data=cvdf, type = "n", axes=FALSE,xlab="",ylab="")
legend(x="center",c("1 - 0-5", "2 - 5-15", "3 - 15-30", "4 - 30-60","5 - 60-100","6 - 100-150"),title="Depths (cm)",cex=1.3)
dev.off()

## Plot based on groupings: relative change in Rsq
# By property
dev.off() # Clear plot
tiff("CV_RelRsq_vs_Ndiffall_byProp.tif",width = 6.5, height = 4.5, units = 'in', res = 200 )
par(mar = c(4, 4, 1.5, 1),mgp=c(2.5,1,0),mfrow=c(1,5),mfrow=c(2,3))
plot(Rsq_reldiff~N_diff,data=cvdf[cvdf$sym==1,],pch = as.character(cvdf$Layer),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Relative change in   ", R^2,sep="")),ylim=c(-0.5,0.5),main=cvdf[cvdf$sym==1,]$Prop[1])
lines(x$N_diff,y,lwd=2)
plot(Rsq_reldiff~N_diff,data=cvdf[cvdf$sym==2,],pch = as.character(cvdf$Layer),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Relative change in   ", R^2,sep="")),ylim=c(-0.5,0.5),main=cvdf[cvdf$sym==2,]$Prop[1])
lines(x$N_diff,y,lwd=2)
plot(Rsq_reldiff~N_diff,data=cvdf[cvdf$sym==3,],pch = as.character(cvdf$Layer),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Relative change in   ", R^2,sep="")),ylim=c(-0.5,0.5),main=cvdf[cvdf$sym==3,]$Prop[1])
lines(x$N_diff,y,lwd=2)
plot(Rsq_reldiff~N_diff,data=cvdf[cvdf$sym==4,],pch = as.character(cvdf$Layer),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Relative change in   ", R^2,sep="")),ylim=c(-0.5,0.5),main=cvdf[cvdf$sym==4,]$Prop[1])
lines(x$N_diff,y,lwd=2)
plot(Rsq_reldiff~N_diff,data=cvdf[cvdf$sym==5,],pch = as.character(cvdf$Layer),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Relative change in   ", R^2,sep="")),ylim=c(-0.5,0.5),main=cvdf[cvdf$sym==5,]$Prop[1])
lines(x$N_diff,y,lwd=2)
plot(Rsq_diff~N_diff,data=cvdf, type = "n", axes=FALSE,xlab="",ylab="")
legend(x="center",c("1 - 0-5", "2 - 5-15", "3 - 15-30", "4 - 30-60","5 - 60-100","6 - 100-150"),title="Depths (cm)",cex=1.3)
dev.off()

# By depth
dev.off() # Clear plot
tiff("CV_Rsq_vs_Ndiffall_byDepth.tif",width = 6.5, height = 4.5, units = 'in', res = 200 )
par(mar = c(4, 4, 4, 2),mgp=c(2.5,1,0),mfrow=c(2,3))
plot(Rsq_diff~N_diff,data=cvdf[cvdf$Layer==1,],pch = as.character(cvdf[cvdf$Layer==1,]$sym),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Change in   ", R^2,sep="")),ylim=c(-0.12,0.12),main="0-5 cm")
lines(x$N_diff,y,lwd=2)
plot(Rsq_diff~N_diff,data=cvdf[cvdf$Layer==2,],pch = as.character(cvdf[cvdf$Layer==1,]$sym),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Change in   ", R^2,sep="")),ylim=c(-0.12,0.12),main="5-15 cm")
lines(x$N_diff,y,lwd=2)
plot(Rsq_diff~N_diff,data=cvdf[cvdf$Layer==3,],pch = as.character(cvdf[cvdf$Layer==1,]$sym),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Change in   ", R^2,sep="")),ylim=c(-0.12,0.12),main="15-30 cm")
lines(x$N_diff,y,lwd=2)
plot(Rsq_diff~N_diff,data=cvdf[cvdf$Layer==4,],pch = as.character(cvdf[cvdf$Layer==1,]$sym),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Change in   ", R^2,sep="")),ylim=c(-0.12,0.12),main="30-60 cm")
lines(x$N_diff,y,lwd=2)
plot(Rsq_diff~N_diff,data=cvdf[cvdf$Layer==5,],pch = as.character(cvdf[cvdf$Layer==1,]$sym),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Change in   ", R^2,sep="")),ylim=c(-0.12,0.12),main="60-100 cm")
lines(x$N_diff,y,lwd=2)
plot(Rsq_diff~N_diff,data=cvdf[cvdf$Layer==6,],pch = as.character(cvdf[cvdf$Layer==1,]$sym),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Change in   ", R^2,sep="")),ylim=c(-0.12,0.12),main="100-150 cm")
lines(x$N_diff,y,lwd=2)
legend(0.125, 0.13,c("1 - Clay", "2 - OC", "3 - pH", "4 - Sand","5 - Silt"),title="Properties",cex=0.6)
dev.off()

# By depth: Relative Rsq
dev.off() # Clear plot
tiff("CV_RelRsq_vs_Ndiffall_byDepth.tif",width = 6.5, height = 4.5, units = 'in', res = 200 )
par(mar = c(4, 4, 4, 2),mgp=c(2.5,1,0),mfrow=c(2,3))
plot(Rsq_reldiff~N_diff,data=cvdf[cvdf$Layer==1,],pch = as.character(cvdf[cvdf$Layer==1,]$sym),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Relative change in   ", R^2,sep="")),ylim=c(-0.5,0.5),main="0-5 cm")
lines(x$N_diff,y,lwd=2)
plot(Rsq_reldiff~N_diff,data=cvdf[cvdf$Layer==2,],pch = as.character(cvdf[cvdf$Layer==1,]$sym),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Relative change in   ", R^2,sep="")),ylim=c(-0.5,0.5),main="5-15 cm")
lines(x$N_diff,y,lwd=2)
plot(Rsq_reldiff~N_diff,data=cvdf[cvdf$Layer==3,],pch = as.character(cvdf[cvdf$Layer==1,]$sym),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Relative change in   ", R^2,sep="")),ylim=c(-0.5,0.5),main="15-30 cm")
lines(x$N_diff,y,lwd=2)
plot(Rsq_reldiff~N_diff,data=cvdf[cvdf$Layer==4,],pch = as.character(cvdf[cvdf$Layer==1,]$sym),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Relative change in   ", R^2,sep="")),ylim=c(-0.5,0.5),main="30-60 cm")
lines(x$N_diff,y,lwd=2)
plot(Rsq_reldiff~N_diff,data=cvdf[cvdf$Layer==5,],pch = as.character(cvdf[cvdf$Layer==1,]$sym),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Relative change in   ", R^2,sep="")),ylim=c(-0.5,0.5),main="60-100 cm")
lines(x$N_diff,y,lwd=2)
plot(Rsq_reldiff~N_diff,data=cvdf[cvdf$Layer==6,],pch = as.character(cvdf[cvdf$Layer==1,]$sym),xlab="Fractional Diff. in Sample size",ylab=expression(paste("Relative change in   ", R^2,sep="")),ylim=c(-0.5,0.5),main="100-150 cm")
lines(x$N_diff,y,lwd=2)
legend(0.125, 0.52,c("1 - Clay", "2 - OC", "3 - pH", "4 - Sand","5 - Silt"),title="Properties",cex=0.63)
dev.off()

## Try multiple linear regression
lm_multi <- lm(Rsq_diff~N_diff+Layer,data=cvdf)
summary(lm_multi)
boxplot(Rsq_diff~Prop,data=cvdf)
plot(Rsq_diff~Layer,data=cvdf)

