
# Create, train, and test random forests to predict MESMA fractions from Landsat imagery

# Load library dependencies
library(terra)
library(tidyverse)
library(randomForest)
library(lubridate)
library(ggtext)

# Filepath inputs
oli_directory <- "G:/EMIT_MESMA/Landsat/OLI/"
etm_directory <- "G:/EMIT_MESMA/Landsat/ETM/"
emit_directory <- "G:/EMIT_MESMA/fractional_cover/"
emit_mask_directory <- "G:/EMIT_MESMA/envi/"
other_data_directory <- "G:/EMIT_MESMA/Landsat/Other/"
output_directory <- "D:/SERDP/scripts/asd_viper_tools_processing/"

# Get list of filepaths
# EMIT data
emit_images <- list.files(emit_directory, pattern="S$")
emit_datetime <- lapply(strsplit(emit_images, "_"), function(str){return(str[[5]])})
emit_dates <- paste(substr(emit_datetime,1,4), "_",
                    substr(emit_datetime,5,6), "_",
                    substr(emit_datetime,7,8), sep="")
names(emit_images) <- emit_dates
# EMIT cloud masks
emit_cloud_masks <- list.files(emit_mask_directory, pattern="_mask$")
emit_cloud_masks <- emit_cloud_masks[!grepl("band", emit_cloud_masks)]
emit_mask_datetime <- lapply(strsplit(emit_cloud_masks, "_"), function(str){return(str[[5]])})
emit_mask_dates <- paste(substr(emit_mask_datetime,1,4), "_",
                    substr(emit_mask_datetime,5,6), "_",
                    substr(emit_mask_datetime,7,8), sep="")
names(emit_cloud_masks) <- emit_mask_dates
# Other datasets
image_date_strings_oli <- list.files(oli_directory)
image_date_strings_oli <- substr(image_date_strings_oli, 5, nchar(image_date_strings_oli)-4)
image_date_strings_etm <- list.files(etm_directory)
image_date_strings_etm <- substr(image_date_strings_etm, 5, nchar(image_date_strings_etm)-4)



# Function to load and pre-process all imagery from a given date
readRasterData <- function(date_str, target_directory, satellite_prefix, return_image=FALSE){
  # Landsat spectral data
  landsat_data <- terra::rast(paste(target_directory, 
                                    satellite_prefix, "_", date_str, ".tif", sep=""))
  # EMIT data
  emit_data <- terra::rast(paste(emit_directory,
                                 emit_images[date_str],
                                 sep=""))
  fractions_fit <- emit_data[['GV_fraction']] + emit_data[['NPV_fraction']] + emit_data[['SOIL_fraction']] + emit_data[['shade_fraction']]
  fractions_fit <- (fractions_fit != 0)
  # Get EMIT cloud mask
  emit_mask <- terra::rast(paste(emit_mask_directory,
                                 emit_cloud_masks[date_str],
                                 sep=""))
  # Band 8 of the cloud mask is a binary image with 1 at sites with any cloud contamination
  # This includes some aggressive dilation around clouds to capture dilute cloud edges and shadows
  # So this should be conservative; in other applications maybe too much so, but probably fine here
  emit_data <- mask(emit_data, emit_mask[[8]], 
                    maskvalues=1, updatevalue=NA)
  # Other data - topography and land classes
  other_data <- terra::rast(paste(other_data_directory,
                                  "other_correlates_", date_str, ".tif", sep=""))
  # Get masks for natural riparian and upland areas
  natural_riparian <- (other_data[[4]]==5) * (other_data[[3]] <= 10)
  names(natural_riparian) <- "riparian_mask"
  natural_upland <- (other_data[[4]]==5) * (other_data[[3]] > 10)
  names(natural_upland) <- "upland_mask"
  # Assemble and return data
  output_raster <- c(landsat_data, other_data,
                     natural_riparian, natural_upland)
  # Aggregate all other datasets to the EMIT basis
  output_raster_reproj <- terra::project(output_raster, emit_data, method="average")
  
  # Mask to reject pixels that aren't natural riparian, or that have GV==0 (which should never happen when a model is fit) 
  masked_data <- mask(c(emit_data, output_raster_reproj[[1:11]]),
                      output_raster_reproj[['riparian_mask']]*fractions_fit,
                      maskvalues=c(0, NA), updatevalue=NA)
  
  all_raster_data <- as.data.frame(masked_data, xy=TRUE) %>%
    drop_na()
  
  if(!return_image)  
    return(all_raster_data)
  return(masked_data)
}

# ******************************************************
# ************* ETM Models (Landsat 5 & 7) *************

all_model_data_etm <- lapply(image_date_strings_etm,
                         function(date_str){
  print(paste("   Starting to work on imagery from date", date_str))
  newdata <- readRasterData(date_str, etm_directory, "ETM") %>%
    mutate(date_str = date_str,
           year = substr(date_str, 1,4),
           month = substr(date_str, 6,7),
           day = substr(date_str, 9,10))
})
all_model_data_etm <- bind_rows(all_model_data_etm)
all_model_data_etm$doy <- lubridate::yday(paste(all_model_data_etm$year, 
                                                all_model_data_etm$month,
                                                all_model_data_etm$day, 
                                                sep="-"))
all_model_data_etm$doy_sin <- sin(all_model_data_etm$doy/365*pi)
write_csv(all_model_data_etm, "G:/EMIT_MESMA/random_forest/etm_model_data.csv")
# all_model_data_etm <- read_csv("G:/EMIT_MESMA/random_forest/etm_model_data.csv")

# normalize GV fractions
all_model_data_etm <- all_model_data_etm %>%
  mutate(GV_fraction_raw = GV_fraction, 
         NPV_fraction_raw = NPV_fraction,
         SOIL_fraction_raw = SOIL_fraction) %>%
  mutate(GV_fraction = GV_fraction / (1-shade_fraction),
         NPV_fraction = NPV_fraction / (1-shade_fraction),
         SOIL_fraction = SOIL_fraction / (1-shade_fraction))

# Generate training and test datasets
set.seed(1)
etm_training_samples <- sample(1:nrow(all_model_data_etm), 1000)
etm_validation_samples <- sample(1:nrow(all_model_data_etm)[-etm_training_samples], 1000)
# ***** GV *****
# Get training data, stratified by GV level and month of year
etm_training_df_GV <- all_model_data_etm %>% 
  mutate(index = 1:n(),
         GV_bin = factor(floor(GV_fraction*10), levels=c(0:10)),
         month = factor(as.numeric(month), levels=1:12)) %>%
  group_by(month, GV_bin) %>%
  sample_n(size=12)
# Get validation data, holding out any values used in training
etm_validation_df_GV <- all_model_data_etm[-etm_training_df_GV$index,] %>%
  mutate(GV_bin = factor(floor(GV_fraction*10), levels=c(0:10)),
         month = factor(as.numeric(month), levels=1:12)) %>%
  group_by(month, GV_bin) %>%
  sample_n(size=12)

etm_model_GV <- randomForest(GV_fraction ~ SR_B1+SR_B2+SR_B3+SR_B4+SR_B5+SR_B7, data=etm_training_df_GV, ntree=1000, mtry=5)
etm_training_df_GV$GV_prediction <- predict(etm_model_GV, newdata=etm_training_df_GV)
etm_validation_df_GV$GV_prediction <- predict(etm_model_GV, newdata=etm_validation_df_GV)
summary(lm(data=etm_training_df_GV, GV_fraction~GV_prediction))
summary(lm(data=etm_validation_df_GV, GV_fraction~GV_prediction))
print(paste(" RMSE - Training: ",  
            sqrt(mean((etm_training_df_GV$GV_fraction-etm_training_df_GV$GV_prediction)^2)),
            "   Validation: ",
            sqrt(mean((etm_validation_df_GV$GV_fraction-etm_validation_df_GV$GV_prediction)^2)),
            sep=""))
etm_model_GV$importance


# ***** NPV *****
# Get training data, stratified by GV level and month of year
etm_training_df_NPV <- all_model_data_etm %>% 
  mutate(index = 1:n(),
         NPV_bin = factor(floor(NPV_fraction*10), levels=c(0:10)),
         month = factor(as.numeric(month), levels=1:12)) %>%
  group_by(month, NPV_bin) %>%
  slice_sample(n=12)
# Get validation data, holding out any values used in training
etm_validation_df_NPV <- all_model_data_etm[-etm_training_df_NPV$index,] %>%
  mutate(NPV_bin = factor(floor(NPV_fraction*10), levels=c(0:10)),
         month = factor(as.numeric(month), levels=1:12)) %>%
  group_by(month, NPV_bin) %>%
  slice_sample(n=12)

etm_model_NPV <- randomForest(NPV_fraction ~ SR_B1+SR_B2+SR_B3+SR_B4+SR_B5+SR_B7, data=etm_training_df_NPV, ntree=1000, mtry=5)
etm_training_df_NPV$NPV_prediction <- predict(etm_model_NPV, newdata=etm_training_df_NPV)
etm_validation_df_NPV$NPV_prediction <- predict(etm_model_NPV, newdata=etm_validation_df_NPV)
summary(lm(data=etm_training_df_NPV, NPV_fraction~NPV_prediction))
summary(lm(data=etm_validation_df_NPV, NPV_fraction~NPV_prediction))
print(paste(" RMSE - Training: ",  
            sqrt(mean((etm_training_df_NPV$NPV_fraction-etm_training_df_NPV$NPV_prediction)^2)),
            "   Validation: ",
            sqrt(mean((etm_validation_df_NPV$NPV_fraction-etm_validation_df_NPV$NPV_prediction)^2)),
            sep=""))
etm_model_NPV$importance


# ***** Soil *****
# Get training data, stratified by GV level and month of year
etm_training_df_soil <- all_model_data_etm %>% 
  mutate(index = 1:n(),
         soil_bin = factor(floor(SOIL_fraction*10), levels=c(0:10)),
         month = factor(as.numeric(month), levels=1:12)) %>%
  group_by(month, soil_bin) %>%
  slice_sample(n=12)
# Get validation data, holding out any values used in training
etm_validation_df_soil <- all_model_data_etm[-etm_training_df_soil$index,] %>%
  mutate(soil_bin = factor(floor(SOIL_fraction*10), levels=c(0:10)),
         month = factor(as.numeric(month), levels=1:12)) %>%
  group_by(month, soil_bin) %>%
  slice_sample(n=12)

etm_model_soil <- randomForest(SOIL_fraction ~ SR_B1+SR_B2+SR_B3+SR_B4+SR_B5+SR_B7, data=etm_training_df_soil, ntree=1000, mtry=5)
etm_training_df_soil$SOIL_prediction <- predict(etm_model_soil, newdata=etm_training_df_soil)
etm_validation_df_soil$SOIL_prediction <- predict(etm_model_soil, newdata=etm_validation_df_soil)
summary(lm(data=etm_training_df_soil, SOIL_fraction~SOIL_prediction))
summary(lm(data=etm_validation_df_soil, SOIL_fraction~SOIL_prediction))
print(paste(" RMSE - Training: ",  
            sqrt(mean((etm_training_df_soil$SOIL_fraction-etm_training_df_soil$SOIL_prediction)^2)),
            "   Validation: ",
            sqrt(mean((etm_validation_df_soil$SOIL_fraction-etm_validation_df_soil$SOIL_prediction)^2)),
            sep=""))
etm_model_soil$importance

# Save copies of training and test datasets
write_csv(etm_training_df_GV, paste(output_directory, "ETM_training_data_GV.csv"))
write_csv(etm_validation_df_GV, paste(output_directory, "ETM_validation_data_GV.csv"))
write_csv(etm_training_df_NPV, paste(output_directory, "ETM_training_data_NPV.csv"))
write_csv(etm_validation_df_NPV, paste(output_directory, "ETM_validation_data_NPV.csv"))
write_csv(etm_training_df_soil, paste(output_directory, "ETM_training_data_soil.csv"))
write_csv(etm_validation_df_soil, paste(output_directory, "ETM_validation_data_soil.csv"))

# Generate some example plots for error

# Generate some summary plots of error
generateSummaryPlots <- function(validation_df, variable, satellite)
{
  var_true <- sym(paste(variable, "_fraction", sep=""))
  var_pred <- sym(paste(variable, "_prediction", sep=""))
  linmod <- summary(lm(data=validation_df,
                       formula(paste(var_true, var_pred, sep="~"))))
  rmse <- sqrt(mean((as.vector(validation_df[,var_true] - validation_df[,var_pred])[[1]])^2, na.rm=TRUE))
  model_annotation <- paste("R<sup>2</sup> = ", round(linmod$adj.r.squared,3), "<br>",
                            "RMSE = ", round(rmse,3), sep="")
  error_plot <- ggplot() +
    geom_density_2d_filled(data=validation_df, aes(x=!!var_pred, y=!!var_true), n=9) + 
    geom_abline(intercept = linmod$coefficients[1,1], 
                slope = linmod$coefficients[2,1], col="red") + 
    geom_abline(intercept = 0, 
                slope = 1, col="black", linetype="dashed") +
    scale_x_continuous(limits=c(0,1), expand=c(0,0)) + 
    scale_y_continuous(limits=c(0,1), expand=c(0,0)) + 
    scale_fill_brewer(palette=1) + 
    theme_bw() + 
    xlab("GV - Landsat Prediction") + 
    ylab("GV - EMIT") + 
    geom_richtext(data=data.frame(text_out = model_annotation,
                                  x = 0.05, y=0.875),
                  aes(label=text_out, x=x, y=y),
                  hjust=0,
                  show.legend=FALSE) + 
    theme(legend.position="none")
  print(error_plot)
  ggsave(paste(output_directory, satellite, "_", variable, "_error.tif", sep=""),
         error_plot, width=4, height=4)
  return(error_plot)
}

#   GV Plot
generateSummaryPlots(etm_validation_df_GV, "GV", "ETM")
generateSummaryPlots(etm_validation_df_NPV, "NPV", "ETM")
generateSummaryPlots(etm_validation_df_soil, "SOIL", "ETM")



# Output an example predicted dataset
test_raster <- c(terra::rast("G:/EMIT_MESMA/Landsat/ETM/ETM_2023_06_30.tif"),
                 terra::rast("G:/EMIT_MESMA/Landsat/Other/other_correlates_2023_06_30.tif"))
gv_pred <- terra::predict(object=test_raster, model=etm_model_GV)
gv_pred_masked <- mask(gv_pred, 
                (test_raster[[10]] <= 10) * (test_raster[[11]]==5),
                maskvalues=c(0,NA), updatevalue=NA)
terra::writeRaster(gv_pred_masked, "G:/EMIT_MESMA/fractional_cover_predicted/GV_2023_06_30.tif", overwrite=TRUE)





# ******************************************************
# ************* OLI Models (Landsat 8 * 9) *************

all_model_data_oli <- lapply(image_date_strings_oli,
                             function(date_str){
                               print(paste("   Starting to work on imagery from date", date_str))
                               newdata <- readRasterData(date_str, oli_directory, "OLI") %>%
                                 mutate(date_str = date_str,
                                        year = substr(date_str, 1,4),
                                        month = substr(date_str, 6,7),
                                        day = substr(date_str, 9,10))
                             })
all_model_data_oli <- bind_rows(all_model_data_oli)
all_model_data_oli$doy <- lubridate::yday(paste(all_model_data_oli$year, 
                                                all_model_data_oli$month,
                                                all_model_data_oli$day, 
                                                sep="-"))
all_model_data_oli$doy_sin <- sin(all_model_data_oli$doy/365*pi)
write_csv(all_model_data_oli, "G:/EMIT_MESMA/random_forest/oli_model_data.csv")
# all_model_data_oli <- read_csv("G:/EMIT_MESMA/random_forest/oli_model_data.csv")

# normalize GV fractions
all_model_data_oli <- all_model_data_oli %>%
  mutate(GV_fraction_raw = GV_fraction, 
         NPV_fraction_raw = NPV_fraction,
         SOIL_fraction_raw = SOIL_fraction) %>%
  mutate(GV_fraction = GV_fraction / (1-shade_fraction),
         NPV_fraction = NPV_fraction / (1-shade_fraction),
         SOIL_fraction = SOIL_fraction / (1-shade_fraction))

# Generate training and test datasets
set.seed(1)
oli_training_samples <- sample(1:nrow(all_model_data_oli), 1000)
oli_validation_samples <- sample(1:nrow(all_model_data_oli)[-oli_training_samples], 1000)
# ***** GV *****
# Get training data, stratified by GV level and month of year
oli_training_df_GV <- all_model_data_oli %>% 
  mutate(index = 1:n(),
         GV_bin = factor(floor(GV_fraction*10), levels=c(0:10)),
         month = factor(as.numeric(month), levels=1:12)) %>%
  group_by(month, GV_bin) %>%
  sample_n(size=10)
# Get validation data, holding out any values used in training
oli_validation_df_GV <- all_model_data_oli[-oli_training_df_GV$index,] %>%
  mutate(GV_bin = factor(floor(GV_fraction*10), levels=c(0:10)),
         month = factor(as.numeric(month), levels=1:12)) %>%
  group_by(month, GV_bin) %>%
  sample_n(size=10)

oli_model_GV <- randomForest(GV_fraction ~ SR_B1+SR_B2+SR_B3+SR_B4+SR_B5+SR_B7, data=oli_training_df_GV, ntree=1000, mtry=5)
oli_training_df_GV$GV_prediction <- predict(oli_model_GV, newdata=oli_training_df_GV)
oli_validation_df_GV$GV_prediction <- predict(oli_model_GV, newdata=oli_validation_df_GV)
summary(lm(data=oli_training_df_GV, GV_fraction~GV_prediction))
summary(lm(data=oli_validation_df_GV, GV_fraction~GV_prediction))
print(paste(" RMSE - Training: ",  
            sqrt(mean((oli_training_df_GV$GV_fraction-oli_training_df_GV$GV_prediction)^2)),
            "   Validation: ",
            sqrt(mean((oli_validation_df_GV$GV_fraction-oli_validation_df_GV$GV_prediction)^2)),
            sep=""))
oli_model_GV$importance


# ***** NPV *****
# Get training data, stratified by GV level and month of year
oli_training_df_NPV <- all_model_data_oli %>% 
  mutate(index = 1:n(),
         NPV_bin = factor(floor(NPV_fraction*10), levels=c(0:10)),
         month = factor(as.numeric(month), levels=1:12)) %>%
  group_by(month, NPV_bin) %>%
  slice_sample(n=12)
# Get validation data, holding out any values used in training
oli_validation_df_NPV <- all_model_data_oli[-oli_training_df_NPV$index,] %>%
  mutate(NPV_bin = factor(floor(NPV_fraction*10), levels=c(0:10)),
         month = factor(as.numeric(month), levels=1:12)) %>%
  group_by(month, NPV_bin) %>%
  slice_sample(n=12)

oli_model_NPV <- randomForest(NPV_fraction ~ SR_B1+SR_B2+SR_B3+SR_B4+SR_B5+SR_B7, data=oli_training_df_NPV, ntree=1000, mtry=5)
oli_training_df_NPV$NPV_prediction <- predict(oli_model_NPV, newdata=oli_training_df_NPV)
oli_validation_df_NPV$NPV_prediction <- predict(oli_model_NPV, newdata=oli_validation_df_NPV)
summary(lm(data=oli_training_df_NPV, NPV_fraction~NPV_prediction))
summary(lm(data=oli_validation_df_NPV, NPV_fraction~NPV_prediction))
print(paste(" RMSE - Training: ",  
            sqrt(mean((oli_training_df_NPV$NPV_fraction-oli_training_df_NPV$NPV_prediction)^2)),
            "   Validation: ",
            sqrt(mean((oli_validation_df_NPV$NPV_fraction-oli_validation_df_NPV$NPV_prediction)^2)),
            sep=""))
oli_model_NPV$importance


# ***** Soil *****
# Get training data, stratified by GV level and month of year
oli_training_df_soil <- all_model_data_oli %>% 
  mutate(index = 1:n(),
         soil_bin = factor(floor(SOIL_fraction*10), levels=c(0:10)),
         month = factor(as.numeric(month), levels=1:12)) %>%
  group_by(month, soil_bin) %>%
  slice_sample(n=12)
# Get validation data, holding out any values used in training
oli_validation_df_soil <- all_model_data_oli[-oli_training_df_soil$index,] %>%
  mutate(soil_bin = factor(floor(SOIL_fraction*10), levels=c(0:10)),
         month = factor(as.numeric(month), levels=1:12)) %>%
  group_by(month, soil_bin) %>%
  slice_sample(n=12)

oli_model_soil <- randomForest(SOIL_fraction ~ SR_B1+SR_B2+SR_B3+SR_B4+SR_B5+SR_B7, data=oli_training_df_soil, ntree=1000, mtry=5)
oli_training_df_soil$SOIL_prediction <- predict(oli_model_soil, newdata=oli_training_df_soil)
oli_validation_df_soil$SOIL_prediction <- predict(oli_model_soil, newdata=oli_validation_df_soil)
summary(lm(data=oli_training_df_soil, SOIL_fraction~SOIL_prediction))
summary(lm(data=oli_validation_df_soil, SOIL_fraction~SOIL_prediction))
print(paste(" RMSE - Training: ",  
            sqrt(mean((oli_training_df_soil$SOIL_fraction-oli_training_df_soil$SOIL_prediction)^2)),
            "   Validation: ",
            sqrt(mean((oli_validation_df_soil$SOIL_fraction-oli_validation_df_soil$SOIL_prediction)^2)),
            sep=""))
oli_model_soil$importance

# Save copies of training and test datasets
write_csv(oli_training_df_GV, paste(output_directory, "OLI_training_data_GV.csv"))
write_csv(oli_validation_df_GV, paste(output_directory, "OLI_validation_data_GV.csv"))
write_csv(oli_training_df_NPV, paste(output_directory, "OLI_training_data_NPV.csv"))
write_csv(oli_validation_df_NPV, paste(output_directory, "OLI_validation_data_NPV.csv"))
write_csv(oli_training_df_soil, paste(output_directory, "OLI_training_data_soil.csv"))
write_csv(oli_validation_df_soil, paste(output_directory, "OLI_validation_data_soil.csv"))

#   GV Plot
generateSummaryPlots(oli_validation_df_GV, "GV", "OLI")
generateSummaryPlots(oli_validation_df_NPV, "NPV", "OLI")
generateSummaryPlots(oli_validation_df_soil, "SOIL", "OLI")



# Output an example predicted dataset
test_raster <- readRasterData("2023_06_30", oli_directory, "OLI", TRUE) #c(terra::rast("G:/EMIT_MESMA/Landsat/OLI/OLI_2023_06_30.tif"),
              #   terra::rast("G:/EMIT_MESMA/Landsat/Other/other_correlates_2023_06_30.tif"))
gv_pred <- terra::predict(object=test_raster, model=oli_model_GV)
gv_pred_masked <- mask(gv_pred, 
                       (test_raster[["relative_elevation"]] <= 10) * (test_raster[["cdl_class"]]==5),
                       maskvalues=c(0,NA), updatevalue=NA)
gv_pred_masked <- c(gv_pred_masked, test_raster)
rm(gv_pred, test_raster)
ndvi <- (gv_pred_masked[["SR_B4"]]-gv_pred_masked[["SR_B3"]])/(gv_pred_masked[["SR_B4"]]+gv_pred_masked[["SR_B3"]])
names(ndvi) <- "NDVI"
ndvi_gv_comparison <- as.data.frame(c(gv_pred_masked[["GV_fraction"]], ndvi), xy=TRUE)
ndvi_gv_model <- summary(lm(data = ndvi_gv_comparison, GV_fraction ~ NDVI))
gv_ndvi_pred <- ndvi_gv_model$coefficients[1,1] + ndvi_gv_model$coefficients[2,1]*gv_pred_masked[["GV_fraction"]]
names(gv_ndvi_pred) <- "GV_from_NDVI"

terra::writeRaster(gv_pred_masked, "G:/EMIT_MESMA/fractional_cover_predicted/GV_2023_06_30.tif", overwrite=TRUE)



