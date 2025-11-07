
setwd("/PATH/to/base_dir") 
load("combined_muts_hnf1a.RData") #muts_hn
load("combined_muts_fg.RData") #muts_fg
load("combined_muts_fp.RData") #muts_fp 

######
######## fitting the loess model 
# use only core and surface residues to fit the loess model 
library(dplyr)
library(rlang)
library(ggplot2)

cbPalette<- c( "#E69F00" ,    "#999999"    , "#56B4E9"   ,   "#D55E00" ,  "#009E73"    ,
 "#F0E442"   )

fit_loess_with_anchor <- function(df, binding_col, abundance_col, binding_se_col, abundance_se_col,
                                  span = 0.9, plot = TRUE) {
  
  # Convert to symbols
  binding_col_sym <- ensym(binding_col)
  abundance_col_sym <- ensym(abundance_col)
  binding_se_col_sym <- ensym(binding_se_col)
  abundance_se_col_sym <- ensym(abundance_se_col)
  
  # Filter valid data
  newdat <- df %>%
    filter(
      surf_core %in% c("core", "surf"),
      !is.na(!!binding_col_sym),
      !is.na(!!abundance_col_sym),
      !is.na(!!binding_se_col_sym),
      !is.na(!!abundance_se_col_sym)
    )
  
  # Create anchor points (to force curve through origin-like region)
  anchor_points <- tibble(
    !!as_string(abundance_col_sym) := c(0, -1),
    !!as_string(binding_col_sym)   := c(0, -1),
    !!as_string(binding_se_col_sym) := c(1e-5, 1e-5),
    !!as_string(abundance_se_col_sym) := c(1e-5, 1e-5),
    surf_core = "anchor"
  )
  
  # Combine real and anchor data
  newdat <- bind_rows(newdat, anchor_points)
  
  # Fit LOESS
  mylo <- loess(
    formula = as.formula(paste(as_string(binding_col_sym), "~", as_string(abundance_col_sym))),
    data = newdat,
    weights = 1 / sqrt((pull(newdat, !!binding_se_col_sym))^2),
    span = span,
    family = "symmetric"
  )
  
  # Predict and residuals
  df <- df %>%
    mutate(
      loess_model = predict(mylo, newdata = df),
      residual_scaled_BA = !!binding_col_sym - loess_model
    )
  
  # Optional diagnostic plots
  if (plot) {
    # Compute dynamic axis limits
    x_min <- min(df %>% pull(!!abundance_col_sym), na.rm = TRUE) - 0.1
    x_max <- max(df %>% pull(!!abundance_col_sym), na.rm = TRUE) + 0.1
    y_min <- min(df %>% pull(!!binding_col_sym), na.rm = TRUE) - 0.1
    y_max <- max(df %>% pull(!!binding_col_sym), na.rm = TRUE) + 0.1
    
    p1 <- ggplot(df, aes(x = !!abundance_col_sym, y = !!binding_col_sym)) +
      geom_point(alpha = 0.4, aes(col=surf_core)) +
      geom_line(aes(y = loess_model),  linewidth = 1,  col="blue") +
      labs(x = as_string(abundance_col_sym), y = as_string(binding_col_sym),
           title = "LOESS fit with anchor points") +
      theme_minimal() + 
      coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) + 
      geom_vline(xintercept = -0.8, linewidth=1, col="red", linetype = "dashed") + # highly destabilising variants < -0.8
      scale_colour_manual(values=cbPalette) 
    
    p2 <- ggplot(df, aes(x = !!abundance_col_sym, y = residual_scaled_BA)) +
      geom_point(alpha = 0.4, aes(col=surf_core)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      labs(x = as_string(abundance_col_sym), y = "Residuals (Binding - LOESS fit)",
           title = "Residual diagnostics") +
      theme_minimal() +
      scale_colour_manual(values=cbPalette) +
      geom_vline(xintercept = -0.8, linewidth=1, col="red", linetype = "dashed") # highly destabilising variants < -0.8
    
    print(p1)
    print(p2)
  }
  
  invisible(list(data = df, model = mylo))
  return(df)
}

#### example run 

hn_result <- fit_loess_with_anchor(
  muts_hn,
  binding_col = scaled_fitness_b1h,
  abundance_col = mean_scaled_ss_fitness,
  binding_se_col = scaled_se_b1h,
  abundance_se_col = mean_scaled_ss_se, 
  plot = TRUE
)


fg_result <- fit_loess_with_anchor(
  muts_fg,
  binding_col = scaled_b1h,
  abundance_col = mean_scaled_ss_fitness,
  binding_se_col = scaled_b1h_se,
  abundance_se_col = mean_scaled_ss_se, 
  plot = TRUE
)

fp_result <- fit_loess_with_anchor(
  muts_fp,
  binding_col = scaled_b1h,
  abundance_col = mean_scaled_ss_fitness,
  binding_se_col = scaled_b1h_se,
  abundance_se_col = mean_scaled_ss_se,
  plot = TRUE
)


# plot the loess 
# Optional diagnostic plots

 
}

invisible(list(data = df, model = mylo))

