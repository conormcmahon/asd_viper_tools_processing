
# Read in ENVI format data from ASD, output by PRISM USGS software
# Convert to matrix format 
# Average across adjacent spectra and subset to only 'selected' data (after manual inspection with PRISM)
# Output averaged target spectra to Viper Tools .sli format
# Create associated metadata

# Load relevant packages
library(terra)
library(tidyverse)

# Input library in tif format (data file, not the .hdr file)
input_library_path <- "D:/SERDP/SRER/ASD_Data/ENVI/AZ_spectra_all_tif.tif"
# Spectral library header file
input_header_path <- "D:/SERDP/SRER/ASD_Data/ENVI/AZ_spectra_all.hdr"
# Metadata for input library
input_metadata_filepath <- "D:/SERDP/SRER/ASD_Data/ENVI/sli_metadata.csv"
# Output filepath for completed library
output_filepath <- "D:/SERDP/SRER/ASD_Data/averaged_data/az_2024_filtered_averaged"

# Read ENVI library
input_library <- terra::rast(input_library_path)
# Read header 
header <- readLines(input_header_path)
# Read metadata
input_metadata <- read_csv(input_metadata_filepath)

# Get list of spectral groups to consider
