# IC 50 fitting for abundance assay

setwd("/PATH/to/base_dir") 
load("dimsum_output_ss.RData") # ss_raw
load("combined_muts_hnf1a.RData") #muts_hn
load("combined_muts_fg.RData") #muts_fg
load("combined_muts_fp.RData") #muts_fp 

# biophysical classification:  

######## model ic50, and reorganise the data 
# try to fit the model: 
library(drc)
library(dplyr)
#### Fitting the model #### 

trial_frame<- for (i in 1:3){
  trial_frame[[i]]= unique(ss_raw[, c("aa_seq", "Nham_aa",  "STOP" )])
  trial_frame[[i]][, columns] <- NA 
} 

for (i in 1:3){ 
  for (j in 1:length(rownames(trial_frame[[i]]))) {
    dat= ss_raw[[i]][ss_raw[[i]]$aa_seq==trial_frame[[i]][j, "aa_seq"], c("SpecCon", "growthrate")]
    dat$SpecCon= log2(dat$SpecCon)
    tryCatch(
      
      {# Attempt to fit the model
        
        ## https://rpubs.com/TX-YXL/656451
        ## Fitting a four-parameter log-logistic model
        ## with user-defined parameter names
        # b = the Hill coefficient; c = min_value; d = max_value; e = EC50.
        mymod=drm(dat[, 2]~dat[,1], data = dat, fct = LL.4(names=c("Hill slope","Min","Max","IC50")),
                  start = c(b = 1, d = 0, e = max(dat[,2]), f = median(dat[,1]))) 
        
        st=summary(mymod)$coefficients
        trial_frame[i, columns] <- unname(c(st[, 1], st[, 2], st[, 3], st[, 4]))
      },
      
      error = function(e) {
        # Return NULL to indicate failure
        return(NA) }
    )
  }
}

library(drc)

# Example: Define the output columns you want to store
columns <- c("Hill_slope", "Min", "Max", "IC50", "SE_Hill", "SE_Min", "SE_Max", "SE_IC50")

# Initialize a list to store results
trial_frame <- vector("list", 3)

#### Fitting the model #### 
for (i in 1:3) {
  
  # Create the base table for this replicate
  trial_frame[[i]] <- unique(ss_raw[[i]][, c("aa_seq", "Nham_aa", "STOP")])
  trial_frame[[i]][, columns] <- NA
  
  # Iterate through each amino acid sequence
  for (j in seq_len(nrow(trial_frame[[i]]))) {
    
    dat <- ss_raw[[i]] %>%
      dplyr::filter(aa_seq == trial_frame[[i]][j, "aa_seq"]) %>%
      dplyr::select(SpecCon, growthrate)
    
  
    
    # Skip if insufficient data
    if (nrow(dat) < 4) next  
    
    dat$SpecCon <- log2(dat$SpecCon)
    
    # Fit model with error handling
    fit <- tryCatch({
      drm(
        growthrate ~ SpecCon,
        data = dat,
        fct = LL.4(names = c("Hill_slope", "Min", "Max", "IC50")),
        control = drmc(errorm = FALSE)
      )
    }, error = function(e) NULL)
    
    # Skip failed fits
    if (is.null(fit)) next
    
    # Extract summary statistics
    st <- summary(fit)$coefficients
    
    # Fill results: estimate + standard error
    vals <- as.numeric(t(st[, 1:2])) # flatten estimates and SEs
    names(vals) <- c("Hill_slope", "SE_Hill", "Min", "SE_Min", "Max", "SE_Max", "IC50", "SE_IC50")
    
    trial_frame[[i]][j, intersect(names(vals), columns)] <- vals[intersect(names(vals), columns)]
  }
}
#####
###### compare these values with mean of the three abundance enrichment 
tmp_list <- list(
  hn=muts_hn, 
  fg=muts_fg,
  fp=muts_fp
)
# Combine pairwise data frames in two lists by 'aa_seq' 
merged_list <- Map(function(df1, df2) {
  inner_join(df1, df2[, c("aa_seq","mean_scaled_ss_fitness"  ,    "mean_scaled_ss_se"  )], by = "aa_seq")
}, trial_frame, tmp_list)

names(merged_list) = names(tmp_list)

lapply(merged_list, function(v) print(summary(v))) 
plots<-list(lapply(seq_along(merged_list), function(i) {
  df <- merged_list[[i]]
  # Ensure no NA values (lm and smoothScatter cannot handle them)
  df <- df[complete.cases(df[, c("IC50", "mean_scaled_ss_fitness")]), ]
  
  smoothScatter( 
    x=df$IC50,
    y=df$mean_scaled_ss_fitness ,
      main = names(merged_list)[i],
      xlab = "IC50",
      ylab = "Mean Abundance"
    ) 
  
  
  # Add linear regression line
  model <- lm(mean_scaled_ss_fitness ~ IC50, data = df)
  abline(model, col = "blue", lwd = 2)
  
  # Add horizontal line at y = 0
  abline(h = 0, col = "grey40", lwd = 1.5, lty = 2)
  # Add vertical line at x = log2(2500) # limit of the assay
  abline(v = 11.28771, col = "red", lwd = 1.5,)
  
  
  # Optionally return the model if you want to store results
  return(model)
  
}))
##### any over 11.28771 in IC50 is not reliable as an IC50 estimate

