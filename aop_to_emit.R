

# Convert the ASD spectra to use EMIT wavelengths
# ASD has 

# ************* Dependencies ************* 

# Load relevant packages
library(terra)
library(tidyverse)
library(janitor)

# ************* User Settings ************* 
# Filepaths 
aop_data_path <- "D:/SERDP/SRER/AOP_Spectra/spectra.tif"
aop_spectral_info_path <- "D:/SERDP/SRER/AOP_Spectra/aop_spectral_characteristics.csv"
aop_metadata_path <- "D:/SERDP/SRER/AOP_Spectra/metadata.csv"
emit_reflectance_path <- "D:/SERDP/SRER/EMIT/EMIT_L2A_RFLUNCERT_001_20231015T215344_2328814_020.nc"
output_spectrum_path <- "D:/SERDP/SRER/ASD_Data/EMIT/AZ_spectra_aop_emit"

# Variables for numeric integration in convolution
wavelength_samples <- 100

# ************* Load Data ************* 
aop_library <- terra::rast(aop_data_path)
emit_example <- ncdf4::nc_open(emit_reflectance_path)

# ************* Process Input ENVI Headers for ASD Data ************* 
aop_spectral_info <- read_csv(aop_spectral_info_path)
aop_fwhms <- aop_spectral_info$fwhm
aop_wavelengths <- aop_spectral_info$wavelength
aop_wavelengths <- aop_wavelengths / 1000
aop_fwhms <- aop_fwhms / 1000

# ************* Get Wavelengths and FWHMs from EMIT ************* 

# Load an EMIT scene from which to pull target wavelengths
emit_wavelengths <- ncvar_get(emit_example, emit_example$var[["sensor_band_parameters/wavelengths"]])
emit_fwhms <- ncvar_get(emit_example, emit_example$var[["sensor_band_parameters/fwhm"]])
# EMIT wavelengths are in nm - convert to um to match ASD
emit_wavelengths <- emit_wavelengths / 1000
emit_fwhms <- emit_fwhms / 1000


# ************* Get Weights Matrix ************* 
# Assumes all instrument responses are gaussian distributions
# Estimates EMIT response based on intersection of wavelength gaussians with equivalent
#   Gaussians from ASD distribution
# Estimate standard deviation for a Gaussian for each wavelength in each sensor, from FWHM
sd_emit <- emit_fwhms/2/sqrt(2*log(2))
sd_aop <- aop_fwhms/2/sqrt(2*log(2))
# Initialize empty weights matrix from ASD to EMIT
emit_weights <- matrix(0,
                       length(emit_fwhms),
                       length(aop_fwhms))
# Function to estimate weights from AOP to EMIT (for one cell of matrix)
getWeight <- function(emit_i, aop_j)
{
  integration_wavelengths <- seq(from = emit_wavelengths[emit_i] - 3*sd_emit[emit_i],
                                 to = emit_wavelengths[emit_i] + 3*sd_emit[emit_i],
                                 by = (6*sd_emit[emit_i])/wavelength_samples)
  emit_weights <- dnorm(integration_wavelengths, 
                        mean=emit_wavelengths[emit_i], 
                        sd=sd_emit[emit_i])
  aop_weights <- dnorm(integration_wavelengths, 
                       mean=aop_wavelengths[aop_j], 
                       sd=sd_aop[aop_j])
  return(sum(emit_weights*aop_weights))
}
# Apply function over all ASD and EMIT wavelength combinations to populate matrix
emit_weights <- unlist(lapply(1:length(aop_wavelengths), function(aop_ind){
  return(unlist(lapply(1:length(emit_wavelengths), getWeight, aop_j=aop_ind)))}))
# Convert from 1D vector to 2D matrix
emit_weights <- matrix(emit_weights, nrow=length(emit_fwhms), ncol=length(aop_fwhms))
# Scale each row (one EMIT wavelength) to have a summed value of 1
emit_weights <- emit_weights * (1/rowSums(emit_weights))




# ************* Apply weights to convolve ASD spectra to EMIT resolution ************* 

output_library <- terra::rast((as.matrix(aop_library, wide=TRUE)) %*% t(emit_weights))




# ************* Manually Filter Spectra ************* 

# After manually reviewing spectra, choosing a few problematic ones to drop (probable mixed pixels)

# Copy metadata file for ViperTools
aop_metadata_no_fid <- read_csv(aop_metadata_path) %>% dplyr::select(-fid)

# These spectra were manually selected
bad_spectra <- c("Celtis_1a",
                 "Fraxinus_velutina_2c",
                 "Fraxinus_velutina_3c",
                 "Pavementc",
                 "Populus_fremontii_4c",
                 "Populus_fremontii_11a",
                 "Populus_fremontii_19a",
                 "Populus_fremontii_19c")
bad_spectra_indices <- lapply(bad_spectra,
                              function(spectrum_name){
                                return(which(spectrum_name == aop_metadata_no_fid$Name))
                              })
bad_spectra_indices <- unlist(bad_spectra_indices)

# Also removing all the Fouquieria spectra, because these have a lot of mixing (odd crown geometry)
fouquieria_indices <- which(grepl("Fouquieria", aop_metadata_no_fid$Name))

# Last, removing all NPV and Mixture spectra EXCEPT these ones
npv_retained <- c("Populus_fremontii_16b", 
                  "Populus_fremontii_16c")
npv_spectral_names <- aop_metadata_no_fid[aop_metadata_no_fid$class %in% c("Mixture", "NPV"),]$Name
npv_spectral_names_removed <- npv_spectral_names[which(!(npv_spectral_names %in% npv_retained))]
npv_spectral_indices_removed <- lapply(npv_spectral_names_removed,
                                       function(spectrum_name){
                                         return(which(spectrum_name == aop_metadata_no_fid$Name))
                                       })
npv_spectral_indices_removed <- unlist(npv_spectral_indices_removed)

aop_metadata_no_fid <- aop_metadata_no_fid[-c(bad_spectra_indices, fouquieria_indices, npv_spectral_indices_removed),]
output_library <- terra::rast(as.matrix(output_library, wide=TRUE)[-c(bad_spectra_indices, fouquieria_indices, npv_spectral_indices_removed),])


# ************* Write Result to Disk ************* 

# Write new library data file to disk
#    Create ENVI file format and MESMA (Dar's) file format
terra::writeRaster(output_library/10000, filename=output_spectrum_path, filetype="ENVI", overwrite=TRUE)
# Create a copy of the original ENVI file in .sli format for Viper Tools
file.copy(output_spectrum_path, 
          paste(output_spectrum_path, "_sli.sli", sep=""), 
          overwrite=TRUE)
file.copy(paste(output_spectrum_path, ".hdr", sep=""), 
          paste(output_spectrum_path, "_sli.hdr", sep=""), 
          overwrite=TRUE)

# ************* Write Result to Disk ************* 

write_csv(aop_metadata_no_fid, paste(output_spectrum_path, "_sli.csv", sep=""))
# Read in spectrum names from metadata
spectra_names <- aop_metadata_no_fid$Name

# Add spectra names to header file (.hdr)
#    Read auto-generated header 
output_header <- readLines(paste(output_spectrum_path, "_sli.hdr", sep=""))
#    Generate well-formatted list of spectra names for output
spectra_names_lines <- paste("spectra names = {\n", paste(spectra_names, collapse=', \n'), "}", sep="")
spectra_wavelengths_lines <- paste("wavelength = {\n", paste(emit_wavelengths, collapse=', \n'), "}", sep="")
#    Correct a few header data lines to the Spectral Library format
output_header[8] <- "file type = ENVI Spectral Library"
output_header[12] <- "sensor_type = Unknown"
output_header[15] <- "wavelength units = Micrometers"
# NOTE for some reason BBL in ENVI are formatted with 0 for 'bad' and 1 for 'good' 
bad_bands <- rep(1, 285)
bad_bands[c(1:10, 131:141, 189:220, 280:285)] <- 0
bbl <- paste("bbl = {\n", paste(bad_bands, collapse=","), "}\n", sep="")
#    Add band names and wavelengths to header file
output_header <- c(output_header, 
                   spectra_names_lines,
                   spectra_wavelengths_lines,
                   bbl)
#    Add list of spectra names to end of text file and write to disk
writeLines(output_header, 
           paste(output_spectrum_path, "_sli.hdr", sep=""))



