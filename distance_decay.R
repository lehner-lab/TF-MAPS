### This is to fit the exponential decay based on the absolute binding residuals & to plot them. 

setwd("/PATH/to/base_dir") 
load("combined_muts_hnf1a.RData") #muts_hn
load("combined_muts_fg.RData") #muts_fg
load("combined_muts_fp.RData") #muts_fp 

###### fittting the decay curve

library(minpack.lm)

fit_exponential_decay <- function(dat, x, y) {
  # Check if the specified columns exist
  if (!(x %in% names(dat)) || !(y %in% names(dat))) {
    stop("Specified x or y column not found in the data frame.")
  }
  
  # Optional filtering if STOP column exists
  if ("STOP" %in% names(dat)) {
    dat <- dat[dat$STOP == "FALSE" & !is.na(dat[[y]]), c(x, y)]
  } else {
    dat <- dat[!is.na(dat[[y]]), c(x, y)]
  }
  
  # Extract variables
  x_vals <- dat[[x]]
  y_vals <- abs(dat[[y]])  # ensure positive
  
  # Initial parameter estimates
  start_vals <- list(
    a = max(y_vals, na.rm = TRUE),
    b = 0.5,
    c = min(y_vals, na.rm = TRUE)
  )
  
  # Fit exponential decay model
  fit <- nlsLM(y_vals ~ a * exp(-b * x_vals) + c, start = start_vals)
  
  return(fit)
}

###  fit to the data

decay_curves<- list(
  hn=fit_exponential_decay(muts_hn, "Min_Distance" ,"residual_scaled_BA"),
  fg=fit_exponential_decay(muts_fg, "Min_Distance" ,"residual_scaled_BA"),
  fp=fit_exponential_decay(muts_fp, "Min_Distance" ,"residual_scaled_BA")
)
 
muts<- list(
  hn=muts_hn, 
  fg=muts_fg, 
  fp=muts_fp
)

pos<- list(
  hn=pos_hn, 
  fg=pos_fg, 
  fp=pos_fp
)
 
names(pos[[1]])[names(pos[[1]]) == "abs_median_residual_scaled_BA"] <- "median_abs_residual"

for (i in 1:3) {
  muts[[i]]$scal_decay_residual <- NA
  valid_rows <- !is.na(muts[[i]]$Min_Distance)
  
  # make a data.frame for prediction with correct column name
  newdata <- data.frame(x_vals = muts[[i]]$Min_Distance[valid_rows])
  
  # predict using named argument
  preds <- predict(decay_curves[[i]], newdata = newdata)
  
  muts[[i]]$scal_decay_residual[valid_rows] <- preds
}

#
######## plot the decay curves 
# ######## plot the decay curves 

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
cbPalette <- c("#E69F00","#999999","#56B4E9",  "#0072B2","#009E73",    "#F0E442", "#D55E00", "#CC79A7", 
               'aquamarine' ,   'aquamarine1',  'aquamarine2',  'aquamarine3',  'aquamarine4') # first gray, second near yellow, third near blue


myplots<- list() 
myplots_pos <- list()  



for (i in 1:3) {
  if ("STOP" %in% names(muts[[i]])) {
      dat <- muts[[i]][muts[[i]]$STOP == "FALSE", ]
  } else {
    dat <- muts[[i]]
  }
    
  myplots[[i]] = ggplot() + 
    geom_point(
      data = dat, 
      aes(x = Min_Distance, y = abs(residual_scaled_BA), 
          col = surf_core ), 
      shape = 21,  # Circle with outline and fill
      alpha = 0.7, size = 2  # Adjust transparency and size
    ) + 
    theme_classic() + 
    labs(x = "Minimal Distance to DNA", y = "abs(binding-residual)", 
         title=names(muts)[i]) + 
    geom_line(
      data = dat, 
      aes(x = Min_Distance, y = scal_decay_residual), 
      col = "black", 
      size=1
    ) + 
    theme(
      legend.key.size = unit(0.3, "cm")  # Reduce the size of legend keys
    ) + 
    scale_colour_manual(values = cbPalette, name="")
  
  myplots_pos[[i]] = ggplot() + 
    geom_point(
      data = pos[[i]], 
      aes(x = Min_Distance, y = median_abs_residual, 
          col = surf_core ), 
      shape = 21,  # Circle with outline and fill
      alpha = 0.7, size = 2  # Adjust transparency and size
    ) + 
    theme_classic() + 
    labs(x = "Minimal Distance to DNA", y = "med.abs(binding-residual)", 
         title=names(muts)[i]) + 
    geom_line(
      data = dat, 
      aes(x = Min_Distance, y = scal_decay_residual), 
      col = "black", 
      size=1
    ) + 
    theme(
      legend.key.size = unit(0.3, "cm")  # Reduce the size of legend keys
    ) + 
    scale_colour_manual(values = cbPalette, name="") 
  
}







ggarrange(myplots[[1]], myplots[[2]],myplots[[3]],
       myplots_pos[[1]], myplots_pos[[2]],myplots_pos[[3]])

#





