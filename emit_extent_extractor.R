
# Load all EMIT envi files in a folder, get their extent, and turn it into a shapefile
# Output is compressed to .zip format for upload to Google Earth Engine

library(terra)
library(sf)
library(utils)
library(zip)

target_directory <- "G:/EMIT_MESMA/envi/"

# Temporarily change the working directory to the target location
#    This is necessary because the Windows tools for zipping files are very tempermental
original_wd <- getwd()
setwd(target_directory)

emit_files <- list.files(target_directory, 
                         pattern = "*reflectance$")

extractExtent <- function(filename)
{
  # Filename format is:   EMIT_L2A_RFL_001_YYYYMMDDTHHMMSS_ddddddd_002_reflectance
  # We're interested in the 'YYYYMMDD' bit
  date_str <- strsplit(filename, "_")[[1]][[5]]
  date_str <- strsplit(date_str, "T")[[1]][[1]]
  year <- (substr(date_str, 1,4))
  month <- (substr(date_str, 5,6))
  day <- (substr(date_str, 7,8))
  # Generate output filename: extent_YYYY_MM_DD
  output_filename <- paste("extents/EMIT_extent_",
                           year, "_", month, "_", day, sep="")
  # Load the target image
  emit_data <- terra::rast(filename, lyrs=1)
  # Write the extent to disk
  st_write(obj=st_as_sfc(st_bbox(emit_data)),   # extent of the image
           dsn=paste(output_filename, ".shp", sep=""),
           delete_dsn=TRUE)
  # Compress .shp files into a .zip archive
  zip::zip(paste(output_filename,".zip",sep=""),
           c(paste(output_filename,".dbf",sep=""),
             paste(output_filename,".prj",sep=""),
             paste(output_filename,".shp",sep=""),
             paste(output_filename,".shx",sep="")))
  return(output_filename)
}

all_files <- lapply(emit_files, extractExtent)

setwd(original_wd)
