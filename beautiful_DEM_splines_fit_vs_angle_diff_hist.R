
#### Caroline Haslebacher
# 2025-12-03
# 
#%%
# algol:
# setwd('/d0/chaslebacher/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')
# DELL tower
setwd('E:/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')

source("Sunazimuth_IsotonicH.R")

# load libraries
library(raster)

source('AxialDataSphere.R')
source('GLM_3D.R')
source('Directional.R')




#%%
# We can use AxialAngularDensity(w, theta) and retrieve weights with fth[lineament_azimuth_math_crs]
# for this, let's load some data:
# ds <- read.csv("ULDR_children_for_vMF_spherical.csv",
# 			   header=TRUE,sep=',')
# children
# ds <- read.csv("dff_children_nofilter_for_vMF_spherical_phocube.csv",
# 			   header=TRUE,sep=',')
# PARENTS
ds <- read.csv("./data/manualsegs/dff_all_nofilter_for_vMF_spherical.csv",
 			   header=TRUE,sep=',')
str(ds)

X <- as.matrix(ds[,2:4])
V <- as.matrix(ds[,8:10])


sunazi <- pi/180*(ds$sun)

#
# NEW: use splines weights 1/h(t):
# first, we need to calculate the angle between the sun direction and the azimuth
# sun direction was transformed to mathematical crs in Generate_input_for_vMF.py,
# but alpha was not (therefore, you see pi/2-ds$alpha below)
# symmetry allows to make all sunazimuth values to be inside (0, pi)

# note that sunazi element (-3pi/2, -2pi) because of the transformation
sunazi[(sunazi >= -pi) & (sunazi < 0)] <- sunazi[(sunazi >= -pi) & (sunazi < 0)] + pi
sunazi[sunazi < -pi] <- 2*pi + sunazi[sunazi < -pi]
# make alphas in between 0 and pi (not /pi/2, pi/2)
alphs <- (pi/2 - ds$alpha) 
alphs[alphs<0] <- pi + alphs[alphs<0] # e.g. -60deg gets transformed to 180 + (-60deg) = 120deg

# angle difference:
angle_diff <- abs(sunazi - alphs)
angle_diff[angle_diff>pi/2] <- pi - angle_diff[angle_diff>pi/2]

# get weights by getting the nearest data value from dataframe in R
# with df[which.min(angle_diff[i],]
# retrieve the weights easily (pass path to dc):
# gsf <- GetSunFit('dff_children_nofilter_for_isotonic_spherical.csv', M=3)
# shift a little:
gsf <- GetSunFit('dff_children_nofilter_for_isotonic_spherical.csv', M=3, Jshifta = 0.1, Jshiftb = 0.4, Jshiftc=0.7)
# plot(gsf[,1], gsf[,2])
# to datafame:
dfgsf <- data.frame(gsf)
colnames(dfgsf) <- c('t','ht')
# retrieve weights:
customweights <- matrix(0, length(ds[,1]))
for (i in seq_along(ds[,1])){
    customweights[i] <- dfgsf[which.min(abs(angle_diff[i] - dfgsf$t)),]$ht
}

# look at it:
pibins <- seq(0, pi/2, length.out = 16)
adiffh <- hist(angle_diff, freq=FALSE, breaks=pibins, main=paste("Subsolar Sun Azimuth Bias Correction"), xlab="Difference of Sun Azimuth and Lineament Azimuth [rad]", ylab="Probability Density")
lines(dfgsf$t, 1/dfgsf$ht)
# lines(density(angle_diff, to=pi/2), col='red')

# show sunazimuth hist:
hist(sunazi)
hist(alphs)

#%% actually save it:
# visualisation of the calculated 5623 angle differences (histogram) and the fitting curve I got from the DEM study
# note how well they agree.
vis_savepath <- './results/'
png(filename=paste(vis_savepath, 'splines_DEM_WITHdffall_obs_3995PARENTS.png', sep=''), width=1500, height=2000, pointsize=25)
# 2. Adjust margins to fit ylabel
par(mar = c(5, 5, 1.5, 2) + 0.1) # mar= bottom, left, top, right (5, 4, 4, 2)+0.1 is the default
pibins <- seq(0, pi/2, length.out = 16)
# Aik is Difference of Sun Azimuth and Lineament Azimuth [rad]
adiffh <- hist(angle_diff, freq=FALSE, breaks=pibins, main=paste("Subsolar Sun Azimuth Bias Correction"), xlab="Aik", ylab="Probability Density",
 cex = 2,   
 cex.main = 2, # title   
     cex.lab  = 2,   # Increase axis label size
     cex.axis = 2)     # increase axis tick values
lines(dfgsf$t, 1/dfgsf$ht, lwd=2, col='red')
# Second x-axis (e.g., smaller interior/minor ticks without labels)
axis(side = 1, at =pi/180*c(0, 15, 30, 45, 60, 75, 90), labels =c("0°", "15°", "30°", "45°", "60°", "75°", "90°"), tcl = 0.4,
     col.ticks = "blue", # Ticks color
     col.axis = "blue",
	 cex.axis  = 1.8,   # Increase axis label size
	 mgp=c(3,-1.5,0))    # Labels color)
	 # par(mgp = c(title_line, label_line, axis_line)). The default is c(3, 1, 0).
# lines(density(angle_diff, to=pi/2), col='red')
dev.off()

vis_savepath <- './results/'
png(filename=paste(vis_savepath, 'splines_DEM_WITHdffall_obs_3995PARENTS_histonly.png', sep=''), width=1500, height=2000, pointsize=25)
# 2. Adjust margins to fit ylabel
par(mar = c(5, 5, 1.5, 2) + 0.1)
pibins <- seq(0, pi/2, length.out = 16)
# Aik is Difference of Sun Azimuth and Lineament Azimuth [rad]
adiffh <- hist(angle_diff, freq=FALSE, breaks=pibins, main=paste("Subsolar Sun Azimuth Bias Correction"), xlab="Aik", ylab="Probability Density",
 cex = 2,      
 cex.main = 2, # title  
     cex.lab  = 2,   # Increase axis label size
     cex.axis = 2)     # increase axis tick values
# Second x-axis (e.g., smaller interior/minor ticks without labels)
axis(side = 1, at =pi/180*c(0, 15, 30, 45, 60, 75, 90), labels =c("0°", "15°", "30°", "45°", "60°", "75°", "90°"), tcl = 0.4,
     col.ticks = "blue", # Ticks color
     col.axis = "blue",
	 cex.axis  = 1.8,   # Increase axis label size
	 mgp=c(3,-1.5,0))    # Labels color)
# lines(density(angle_diff, to=pi/2), col='red')
dev.off()



#%%
# and show extracted customweights:
plot(angle_diff, customweights, ylim=c(0,7))
lines(dfgsf$t, 1/dfgsf$ht)
lines(dfgsf$t, dfgsf$ht)


#%% lineament categories (parameter id_int)
# id_int = 1 are bands, but there are none in the DTM region
# id_int = 2 are double ridges
# id_int = 3 are ridge complexes
# id_int = 4 are undifferentiated lineae

# # for debugging, choose one, comment out the others:
# # all:
# savestr <- 'All lineaments'
# # double ridges
# savestr <- 'Double Ridges'
# # ridge complexes
# savestr <- 'Ridge Complexes'
# # undifferentiated lineae
# savestr <- 'Undifferentiated Lineae'

cat_names <- c('All lineaments', 'Double Ridges', 'Ridge Complexes', 'Undifferentiated Lineae')
# select rank of splines
Msel <- 3

for(savestr in cat_names){
	if(savestr == 'All lineaments'){
	csvp <- './isotonic_DTMdatasets/dff_children_nofilter_for_vMF_spherical.csv' # all
	ds <- read.csv("./data/manualsegs/dff_all_nofilter_for_vMF_spherical.csv",
				header=TRUE,sep=',')
	}else if (savestr == 'Double Ridges') {
		csvp <- './isotonic_DTMdatasets/DR_children_for_vMF_spherical.csv'
		ds <- read.csv("./data/manualsegs/DR_children_for_vMF_spherical.csv",
					header=TRUE,sep=',')
	}else if (savestr == 'Ridge Complexes') {
		csvp <- './isotonic_DTMdatasets/RC_children_for_vMF_spherical.csv'
		ds <- read.csv("./data/manualsegs/RC_children_for_vMF_spherical.csv",
					header=TRUE,sep=',')
	}else if (savestr == 'Undifferentiated Lineae') {
		csvp <- './isotonic_DTMdatasets/UL_children_for_vMF_spherical.csv'
		ds <- read.csv("./data/manualsegs/UL_children_for_vMF_spherical.csv",
					header=TRUE,sep=',')
	}

	str(ds)
	X <- as.matrix(ds[,2:4])
	V <- as.matrix(ds[,8:10])
	sunazi <- pi/180*(ds$sun)

	# note that sunazi element (-3pi/2, -2pi) because of the transformation
	sunazi[(sunazi >= -pi) & (sunazi < 0)] <- sunazi[(sunazi >= -pi) & (sunazi < 0)] + pi
	sunazi[sunazi < -pi] <- 2*pi + sunazi[sunazi < -pi]
	# make alphas in between 0 and pi (not /pi/2, pi/2)
	alphs <- (pi/2 - ds$alpha) 
	alphs[alphs<0] <- pi + alphs[alphs<0] # e.g. -60deg gets transformed to 180 + (-60deg) = 120deg

	# angle difference:
	angle_diff <- abs(sunazi - alphs)
	angle_diff[angle_diff>pi/2] <- pi - angle_diff[angle_diff>pi/2]

	gsf <- GetSunFit(csvp, M=Msel, Jshifta = 0.1, Jshiftb = 0.4, Jshiftc=0.7)
	# plot(gsf[,1], gsf[,2])
	# to datafame:
	dfgsf <- data.frame(gsf)
	colnames(dfgsf) <- c('t','ht')
	# retrieve weights:
	customweights <- matrix(0, length(ds[,1]))
	for (i in seq_along(ds[,1])){
		customweights[i] <- dfgsf[which.min(abs(angle_diff[i] - dfgsf$t)),]$ht
	}
	# for the first time, we initialize the data frame
	if (savestr == 'All lineaments') {
	   corr_funcs <- dfgsf
	}else{
		# we only append the ht function:
		corr_funcs <- cbind(corr_funcs, dfgsf$ht)
	}

	######### plot A,J original visibility study data as well
	dc <- read.csv(csvp)
	A <- c(dc$A00, dc$A45, dc$A90, dc$A135)
	J <- c(dc$J00, dc$J45, dc$J90, dc$J135)
	# modify so that Aik element [0, pi/2] instead of [0, pi]
	A[A>pi/2] <- pi - A[A>pi/2]

	# shift if desired:
	Jshifta <- 0.1 # former 0
	Jshiftb <- 0.4 # former 0.33
	Jshiftc <- 0.7 # former 0.67
	# no shfit for 1
	J[J==0] <- Jshifta
	J[J==0.33] <- Jshiftb
	J[J==0.67] <- Jshiftc
	####### ready to be plotted

	adiffh <- hist(angle_diff, freq=TRUE, breaks=pibins)
	# rescale to 0-1
	adiffh$counts <- adiffh$counts/max(adiffh$counts)

	# save it
	vis_savepath <- './results/'
	png(filename=paste(vis_savepath, 'splines_Msel', Msel ,'_DTM_', savestr, '.png', sep=''), width=1500, height=2000, pointsize=25)
	# 2. Adjust margins to fit ylabel
	par(mar = c(5, 5, 1.5, 2) + 0.1) # mar= bottom, left, top, right (5, 4, 4, 2)+0.1 is the default
	plot(adiffh, main=savestr, xlab="Aik", ylab="Probability Density (scaled to 0-1)", xlim=c(0, pi/2),
			cex = 2,      
		cex.main = 2, # title  
			cex.lab  = 2,   # Increase axis label size
			cex.axis = 2)     # increase axis tick values
	# Second x-axis (e.g., smaller interior/minor ticks without labels)
	axis(side = 1, at =pi/180*c(0, 15, 30, 45, 60, 75, 90), labels =c("0°", "15°", "30°", "45°", "60°", "75°", "90°"), tcl = 0.4,
     col.ticks = "blue", # Ticks color
     col.axis = "blue",
	 cex.axis  = 1.8,   # Increase axis label size
	 mgp=c(3,-1.5,0))    # Labels color)
	axis(2, at = c(0.2, 0.4, 0.6, 0.8), cex.axis = 2)
	points(A, J)
	pibins <- seq(0, pi/2, length.out = 16)
	lines(dfgsf$t, 1/dfgsf$ht, lwd=2, col='red')
	# lines(density(angle_diff, to=pi/2), col='red')
	dev.off()
}
colnames(corr_funcs) <- c('x', cat_names)
linestyles <- c('solid', 'dashed', "dotted", 'twodash')

png(filename=paste(vis_savepath, 'correction_functions.png', sep=''), width=1500, height=2000, pointsize=25)
# 2. Adjust margins to fit ylabel
par(mar = c(5, 5, 2, 2) + 0.1) # mar= bottom, left, top, right (5, 4, 4, 2)+0.1 is the default
# plot all correction functions:
plot(corr_funcs$x, corr_funcs[,2], type='l', lty=linestyles[1], ylim=c(0,max(corr_funcs)), lwd=3,
	main="Sun Bias Correction Functions", ylab="Weights (1/h(Aik))", xlab="Angular Difference (Aik)",
	cex = 2,      
	cex.main = 2, # title  
	cex.lab  = 2,   # Increase axis label size
	cex.axis = 2)     # increase axis tick values
lines(corr_funcs$x, corr_funcs[,3], lty=linestyles[2], lwd=3)
lines(corr_funcs$x, corr_funcs[,4], lty=linestyles[3], lwd=3)
lines(corr_funcs$x, corr_funcs[,5], lty=linestyles[4], lwd=3)
# Second x-axis (e.g., smaller interior/minor ticks without labels)
axis(side = 1, at =pi/180*c(0, 15, 30, 45, 60, 75, 90), labels =c("0°", "15°", "30°", "45°", "60°", "75°", "90°"), tcl = 0.4,
	col.ticks = "blue", # Ticks color
	col.axis = "blue",
	cex.axis  = 1.8,   # Increase axis label size
	mgp=c(3,-1.5,0))    # Labels color)
# Add the legend
legend("topright",               # Position keyword
       cat_names, # Vector of labels
       col = 'black',       # Vector of corresponding colors
       lty = linestyles,              # Vector of corresponding line types (1 = solid)
       lwd = c(3,3,3,3),              # Vector of corresponding line widths
       cex = 2,     # Increases text size
       pt.cex = 2.0)  # Forces legend symbols to stay normal size
dev.off()

# retrieve min/max of all weights:
# NOTE: we can simply exclude x
min(corr_funcs[,c(2,3,4,5)])
# 1
max(corr_funcs[,c(2,3,4,5)])
# 11.15
# > str(corr_funcs)
# 'data.frame':   181 obs. of  5 variables:
#  $ x                      : num  0 0.00873 0.01745 0.02618 0.03491 ...
#  $ All lineaments         : num  6.73 6.73 6.71 6.67 6.63 ...
#  $ Double Ridges          : num  5.37 5.37 5.36 5.34 5.31 ...
#  $ Ridge Complexes        : num  11.2 11.1 11.1 11.1 11 ...
#  $ Undifferentiated Lineae: num  7.34 7.33 7.3 7.26 7.21 ...


#%%
# plot bands without any DTM splines data (there aren't any bands in the DTM)
savestr <- 'Bands'
ds <- read.csv("./data/manualsegs/Band_children_for_vMF_spherical.csv",
				header=TRUE,sep=',')
str(ds)
X <- as.matrix(ds[,2:4])
V <- as.matrix(ds[,8:10])
sunazi <- pi/180*(ds$sun)

# note that sunazi element (-3pi/2, -2pi) because of the transformation
sunazi[(sunazi >= -pi) & (sunazi < 0)] <- sunazi[(sunazi >= -pi) & (sunazi < 0)] + pi
sunazi[sunazi < -pi] <- 2*pi + sunazi[sunazi < -pi]
# make alphas in between 0 and pi (not /pi/2, pi/2)
alphs <- (pi/2 - ds$alpha) 
alphs[alphs<0] <- pi + alphs[alphs<0] # e.g. -60deg gets transformed to 180 + (-60deg) = 120deg

# angle difference:
angle_diff <- abs(sunazi - alphs)
angle_diff[angle_diff>pi/2] <- pi - angle_diff[angle_diff>pi/2]

adiffh <- hist(angle_diff, freq=TRUE, breaks=pibins)
# rescale to 0-1
adiffh$counts <- adiffh$counts/max(adiffh$counts)

vis_savepath <- './results/'
png(filename=paste(vis_savepath, 'histogram_', savestr, '.png', sep=''), width=1500, height=2000, pointsize=25)
# 2. Adjust margins to fit ylabel
par(mar = c(5, 5, 2, 2) + 0.1) # mar= bottom, left, top, right (5, 4, 4, 2)+0.1 is the default
plot(adiffh, main=savestr, xlab="Aik", ylab="Probability Density (scaled to 0-1)", xlim=c(0, pi/2),
	cex = 2,      
	cex.main = 2, # title  
	cex.lab  = 2,   # Increase axis label size
	cex.axis = 2)     # increase axis tick values
axis(2, at = c(0.2, 0.4, 0.6, 0.8))
pibins <- seq(0, pi/2, length.out = 16)
axis(side = 1, at =pi/180*c(0, 15, 30, 45, 60, 75, 90), labels =c("0°", "15°", "30°", "45°", "60°", "75°", "90°"), tcl = 0.4,
	col.ticks = "blue", # Ticks color
	col.axis = "blue",
	cex.axis  = 1.8,   # Increase axis label size
	mgp=c(3,-1.5,0))    # Labels color)
dev.off()

#%%


