
#%%
# algol:
# setwd('/d0/chaslebacher/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')
# DELL tower
setwd('E:/Caroline/lineament_detection/galileo_manual_segmentation/azimuth_analysis/R_code/Europa_Bingham')
source("Sunazimuth_IsotonicH.R")


#%% CH check data
dc <- read.csv('dff_children_nofilter_for_isotonic_spherical.csv')
# > str(dc)
# 'data.frame':   242 obs. of  24 variables:
#  $ X               : int  0 1 2 3 4 5 6 7 8 9 ...
#  $ X1              : num  -0.565 -0.577 -0.573 -0.577 -0.566 ...
#  $ X2              : num  0.569 0.58 0.584 0.574 0.581 ...
#  $ X3              : num  0.598 0.574 0.575 0.581 0.584 ...
#  $ Y1              : num  -0.693 -0.813 -0.816 -0.695 0.597 ...
#  $ Y2              : num  0.0661 -0.4698 -0.3444 0.0283 0.7782 ...
#  $ Y3              : num  -0.718 -0.343 -0.464 -0.718 -0.196 ...
#  $ V1              : num  -0.594 -0.379 -0.444 -0.584 -0.815 ...
#  $ V2              : num  -0.783 -0.814 -0.811 -0.787 -0.292 ...
#  $ V3              : num  0.183 0.441 0.381 0.197 -0.5 ...
#  $ alpha           : num  1.34 1 1.09 1.33 2.23 ...
#  $ lon             : num  2.35 2.35 2.35 2.36 2.34 ...
#  $ lat             : num  0.641 0.612 0.613 0.62 0.624 ...
#  $ scalarproduct_XY: num  0 0 0 0 0 0 0 0 0 0 ...
#  $ sun             : num  -0.092 -0.092 -0.092 -0.092 -0.092 ...
#  $ cluster         : num  0 0 0 0 0 0 0 0 0 0 ...
#  $ J00             : num  1 0.67 1 0.67 0.33 1 0.33 0.33 1 0.33 ...
#  $ J45             : num  0.33 0 0.33 0.33 1 0.67 0 0.67 0.67 0 ...
#  $ J90             : num  0 0.33 0.67 0.67 0.67 0 0.67 1 0 0.67 ...
#  $ J135            : num  0.67 1 1 1 0.33 0.67 1 1 0.67 1 ...
#  $ A00             : num  1.34 1 1.09 1.33 2.23 ...
#  $ A45             : num  0.555 0.216 0.301 0.541 1.449 ...
#  $ A90             : num  0.231 0.569 0.484 0.245 0.664 ...
#  $ A135            : num  1.016 1.355 1.269 1.03 0.122 ...

# check
plot(dc$A00, dc$J00)
points(dc$A45, dc$J45)
points(dc$A90, dc$J90)
points(dc$A135, dc$J135)

A <- c(dc$A00, dc$A45, dc$A90, dc$A135)
J <- c(dc$J00, dc$J45, dc$J90, dc$J135)
# Check:
plot(A, J)



# modify so that Aik element [0, pi/2] instead of [0, pi]
A[A>pi/2] <- pi - A[A>pi/2]
# check:
plot(A, J)

# shift if desired:
Jshifta <- 0.1 # former 0
Jshiftb <- 0.4 # former 0.33
Jshiftc <- 0.7 # former 0.67
# no shfit for 1
J[J==0] <- Jshifta
J[J==0.33] <- Jshiftb
J[J==0.67] <- Jshiftc
# check:
plot(A, J, ylim=c(0,1))

# now, we split dataframe in preparation of whisker plots.
# Note that the A values have a mix of J values, because they were simply sorted by hillshade (00 sun azimuth, 45deg sun azimuth, etc)
J_vals <- c(0.1, 0.4, 0.7, 1)
tol <- 1e-8

A_list <- setNames(
  lapply(J_vals, function(val) A[abs(J - val) < tol]),
  paste0("A_J", J_vals)
)
# > str(A_list)
# List of 4
#  $ A_J0.1: num [1:171] 0.5018 0.0572 0.1309 0.0243 0.078 ...
#  $ A_J0.4: num [1:235] 0.907 0.664 0.383 0.622 0.615 ...
#  $ A_J0.7: num [1:289] 1.002 1.326 0.847 1.005 0.964 ...
#  $ A_J1  : num [1:273] 1.34 1.09 1.54 1.55 1.53 ...

N <- length(dc$J45)
m <- 4
AJ <- matrix(0,N,2*m) # 2*m, 1 for A, 1 for J
AJ[,seq(1,2*m-1,2)] <- A
AJ[,seq(2,2*m,2)] <- J

plot(AJ)

plot(A,J,xlab='t',ylab='h(t)', ylim=c(0,1))

# Isotonic fit:
TH <- IsotonicH(AJ)
lines(TH[,1],TH[,2],lwd=2,col='black')

# Parametric fit:
XH <- FittedH(AJ,M=6)
lines(XH[,1],XH[,2],lwd=2,col='blue')
XH <- FittedH(AJ,M=4)
lines(XH[,1],XH[,2],lwd=2,col='red')
XH <- FittedH(AJ,M=3)
lines(XH[,1],XH[,2],lwd=2,col='green')
XH <- FittedH(AJ,M=2)
lines(XH[,1],XH[,2],lwd=2,col='orange') # ends at 1, but starts above

xs <- seq(0, pi/2, 0.01)
lines(xs,sin(xs), lwd=2, col='purple')

# final fit:
plot(A,J,xlab='t',ylab='h(t)', ylim=c(0,1))
# Isotonic fit:
TH <- IsotonicH(AJ)
lines(TH[,1],TH[,2],lwd=2,col='black')
# Parametric fit:
XH <- FittedH(AJ,M=3)
lines(XH[,1],XH[,2],lwd=2,col='green')


## show final weighting function:
# 1/h(t)
XH <- FittedH(AJ,M=3)
plot(XH[,1],1/XH[,2],lwd=2,col='green')
points(A,J,xlab='t',ylab='h(t)')
# Isotonic fit:
TH <- IsotonicH(AJ)
lines(TH[,1],TH[,2],lwd=2,col='black')

## finally, I've implemented a function:
# retrieve the weights easily (pass path to dc):
source("Sunazimuth_IsotonicH.R")
gsf <- GetSunFit('dff_children_nofilter_for_isotonic_spherical.csv', M=3, Jshifta = 0.1, Jshiftb = 0.4, Jshiftc=0.7)
plot(gsf[,1], gsf[,2])
plot(gsf[,1], 1/gsf[,2], ylim=c(0,1)) #'original curve'

# retrieve min/max weight values
min(gsf[,2])
# 1
max(gsf[,2])
# 6.73

#%% figures for publication:

png(filename='./results/M3_Splines_fit_visibilitycurve.png', width = 1500, height = 1500, units = "px", pointsize = 32)
# 2. Adjust margins to fit ylabel
par(mar = c(5, 7, 4, 2) + 0.1)
# final fit:
# expression(A[i]^k)
# J <- c(dc$J00, dc$J45, dc$J90, dc$J135)
plot(A,J,xlab='Angular Difference Aik',ylab='Visibility h(Aik)', ylim=c(0,1), pch = 21, bg = rgb(0, 0.5, 1, 0.2), col = 'black', # Blue color with 50% transparency (alpha = 0.5)
 cex = 2,      
     cex.lab  = 2,   # Increase axis label size
     cex.axis = 2)     # increase axis tick values)
# add boxplots
boxplot(A_list, horizontal=TRUE, at=c(0.05, 0.35, 0.65, 0.95), boxwex = 0.05, col="bisque", 
        add=TRUE, # add to plot
        ann=FALSE, # do not annotate xlab and ylab
        xaxt = "n", yaxt = "n", # suppress xticks and yticks
        outline=FALSE # do not indicate outliers
        )
# Parametric fit:
XH <- FittedH(AJ,M=3)
lines(XH[,1],XH[,2],lwd=3,col='red')
legend("topleft", legend = c(expression(h(A[i]^k)), expression(J[i]^k), 'Boxplot'), col = c("red", 'black', 'bisque'), 
        lty = c(1, NA, 1), lwd = c(3, NA, 20), pch = c(NA, 21, NA),  pt.bg=rgb(0, 0.5, 1, 0.2),
       cex = 2,     # Increases text size
       pt.cex = 2.0)  # Forces legend symbols to stay normal size
# degree axis:
axis(side = 1, at =pi/180*c(0, 15, 30, 45, 60, 75, 90), labels =c("0°", "15°", "30°", "45°", "60°", "75°", "90°"), tcl = 0.8,
     col.ticks = "blue", # Ticks color
     col.axis = "blue",
	 cex.axis  = 2,   # Increase axis label size
	 mgp=c(3,-2,0))    # Labels color)
dev.off()


# 1/h(t)
png(filename='./results/M3_Splines_fit_sunbiasfunction.png', width = 1500, height = 1500, units = "px", pointsize = 32)
# 2. Adjust margins to fit ylabel
par(mar = c(5, 7, 4, 2) + 0.1)
# 1/h(t) final sun bias function
# Parametric fit:
plot(A,J,xlab='Angular Difference Aik',ylab='Weights 1/h(Aik)', ylim=c(0,7), pch = 21, bg = rgb(0, 0.5, 1, 0.2), col = 'black', # Blue color with 50% transparency (alpha = 0.5)
 cex = 2,      
     cex.lab  = 2,   # Increase axis label size
     cex.axis = 2)     # increase axis tick values)
XH <- FittedH(AJ,M=3)
points(XH[,1],1/XH[,2], type='l',lwd=3,col='#c49a31')
lines(XH[,1],XH[,2],lwd=3,col='red')
legend("topright", legend = c(expression(1/h(A[i]^k)), expression(h(A[i]^k)), expression(J[i]^k)), col = c('#c49a31', 'red', 'black'), 
        lty = c(1, 1, NA), lwd = c(3, 3, NA), pch = c(NA, NA, 21), pt.bg=rgb(0, 0.5, 1, 0.2),
        cex = 2,     # Increases text size
       pt.cex = 2.0)  # Forces legend symbols to stay normal size
# degree axis:
axis(side = 1, at =pi/180*c(0, 15, 30, 45, 60, 75, 90), labels =c("0°", "15°", "30°", "45°", "60°", "75°", "90°"), tcl = 0.8,
     col.ticks = "blue", # Ticks color
     col.axis = "blue",
	 cex.axis  = 2,   # Increase axis label size
	 mgp=c(3,-2,0))    # Labels color)
dev.off()

#%% Table as guideline

## finally, I've implemented a function:
# gsf are the weights, 1/gsf is the probability of visibility
source("Sunazimuth_IsotonicH.R")
gsf <- GetSunFit('dff_children_nofilter_for_isotonic_spherical.csv', M=3, Jshifta = 0.1, Jshiftb = 0.4, Jshiftc=0.7)
plot(gsf[,1], gsf[,2])
plot(gsf[,1], 1/gsf[,2], ylim=c(0,1)) #'original curve'

dfgsf <- data.frame(gsf)
colnames(dfgsf) <- c('t','ht')
# retrieve weights:
# for angles from 0 to 180deg in increments of 5deg
angles_list <- seq(0, 90, 5)
guide_weights <- matrix(0, length(angles_list))
guide_probs <- matrix(0, length(angles_list))
guide_missedprobs <- matrix(0, length(angles_list))
angles_rad <- matrix(0, length(angles_list))
for (i in seq_along(angles_list)){
        print(i)
        print(pi/180*angles_list[i])
        print(which.min(abs(pi/180*angles_list[i] - dfgsf$t)))
        angles_rad[i] <- pi/180*angles_list[i]
        # I am simply finding the correct index with abs(angle_diff[i] - dfgsf$t))
        guide_weights[i] <- dfgsf[which.min(abs(pi/180*angles_list[i] - dfgsf$t)),]$ht
        guide_probs[i] <- 1/guide_weights[i] 
        guide_missedprobs[i] <- 1 - guide_probs[i]
}

# to csv
df_guide <- data.frame(
  angles_list    = angles_list,
  angles_rad     = round(as.vector(angles_rad), 2),
  guide_probs    = 100*round(as.vector(guide_probs), 2), # to transform to percent
  guide_missedprobs = 100*round(as.vector(guide_missedprobs), 2), # to transform to percent
  guide_weights  = round(as.vector(guide_weights), 2)
)

write.csv(df_guide, "SunCorrection_Guideline.csv", row.names = FALSE)

#%% "integral"
trapezoid <- function(x, y) {
  n <- length(x)
  integral <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  return(integral)
}


angles_list <- seq(0, 90, 1)
guide_weights <- matrix(0, length(angles_list))
guide_probs <- matrix(0, length(angles_list))
angles_rad <- matrix(0, length(angles_list))
for (i in seq_along(angles_list)){
        angles_rad[i] <- pi/180*angles_list[i]
        # I am simply finding the correct index with abs(angle_diff[i] - dfgsf$t))
        guide_weights[i] <- dfgsf[which.min(abs(pi/180*angles_list[i] - dfgsf$t)),]$ht
        guide_probs[i] <- 1/guide_weights[i] 
}


full_trapz <- trapezoid(angles_rad, guide_probs)
total_area <- pi/2 # 1.57*1
full_perc <-  full_trapz / total_area # 0.61 = 61%
missed_perc <- 1 - full_perc # 0.387 = 39%

delta_angle <- angles_rad[2] - angles_rad[1] # same as pi/(2*18)
full_identified_sum <- sum(delta_angle * guide_probs) # this is basically a left-hand integration rule
full_percentage <- full_identified_sum/total_area

#%%
