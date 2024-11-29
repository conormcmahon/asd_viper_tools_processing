

# Convert the ASD spectra to use EMIT wavelengths
# ASD has 

# ************* Dependencies ************* 

# Load relevant packages
library(terra)
library(tidyverse)
library(janitor)
library(ncdf4)

# ************* User Settings ************* 
# Filepaths 
asd_header_path <- "D:/SERDP/SRER/ASD_Data/ENVI/AZ_spectra_all.hdr"
asd_data_path <- "D:/SERDP/SRER/ASD_Data/averaged_data/AZ_spectral_library"
emit_reflectance_path <- "D:/SERDP/SRER/EMIT/EMIT_L2A_RFLUNCERT_001_20231015T215344_2328814_020.nc"
output_spectrum_path <- "D:/SERDP/SRER/ASD_Data/EMIT/AZ_spectra_asd_emit"

# Variables for numeric integration in convolution
wavelength_samples <- 100

# ************* Load Data ************* 
asd_library <- terra::rast(asd_data_path)
emit_example <- ncdf4::nc_open(emit_reflectance_path)

# ************* Process Input ENVI Headers for ASD Data ************* 
asd_header <- readLines(asd_header_path)

# Get list of wavelengths of spectra from the header
#    Use regex to get start and end of list of spectra names
start_ind_asd_wavelengths <- which(grepl("wavelength = ", asd_header)) + 1
closing_bracket_indices <- which(grepl("}", asd_header))
last_ind_asd_wavelengths <- last_ind_asd_wavelengths <- closing_bracket_indices[(closing_bracket_indices >= start_ind_asd_wavelengths)][1]
#    Subset to those rows of the header file
asd_wavelengths_messy <- asd_header[start_ind_asd_wavelengths:last_ind_asd_wavelengths]
#    Reprocess header file text lines into short string names
asd_wavelengths_long <- substr(unlist(lapply(str_split(asd_wavelengths_messy, 
                                                       ".asd"), 
                                             function(str_list){
                                               return(str_list[[1]])})), 
                               2,1000)
#    Reprocess into a numeric format
asd_wavelengths <- substr(asd_wavelengths_long, 
                              1, 
                              ((nchar(asd_wavelengths_long))-1))
asd_wavelengths <- as.numeric(unlist(lapply(asd_wavelengths, 
                                            str_split, 
                                            pattern=", ")))
#    This text-filtering method adds some erroneous zeros on ends of scan lines - remove these
asd_wavelengths <- asd_wavelengths[asd_wavelengths != 0]
#    Clean up temporary variables
rm(asd_wavelengths_messy, asd_wavelengths_long, start_ind_asd_wavelengths, last_ind_asd_wavelengths)

# Get list of fwhm of spectra from the header
#    Use regex to get start and end of list of spectra names
start_ind_asd_fwhms <- which(grepl("fwhm = ", asd_header)) + 1
last_ind_asd_fwhms <- last_ind_asd_fwhms <- closing_bracket_indices[(closing_bracket_indices >= start_ind_asd_fwhms)][1]
#    Subset to those rows of the header file
asd_fwhms_messy <- asd_header[start_ind_asd_fwhms:last_ind_asd_fwhms]
#    Reprocess header file text lines into short string names
asd_fwhms_long <- substr(unlist(lapply(str_split(asd_fwhms_messy, 
                                                           ".asd"), 
                                                 function(str_list){
                                                   return(str_list[[1]])})), 
                                   2,1000)
#    Reprocess into a numeric format
asd_fwhms <- substr(asd_fwhms_long, 
                    1, 
                    ((nchar(asd_fwhms_long))-1))
asd_fwhms <- as.numeric(unlist(lapply(asd_fwhms, 
                                      str_split, 
                                      pattern=", ")))
#    Clean up temporary variables
rm(asd_fwhms_messy, asd_fwhms_long, start_ind_asd_fwhms, last_ind_asd_fwhms)

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
sd_asd <- asd_fwhms/2/sqrt(2*log(2))
# Initialize empty weights matrix from ASD to EMIT
emit_weights <- matrix(0,
                       length(emit_fwhms),
                       length(asd_fwhms))
# Function to estimate weights from ASD to EMIT (for one cell of matrix)
getWeight <- function(emit_i, asd_j)
{
  integration_wavelengths <- seq(from = emit_wavelengths[emit_i] - 3*sd_emit[emit_i],
                                 to = emit_wavelengths[emit_i] + 3*sd_emit[emit_i],
                                 by = (6*sd_emit[emit_i])/wavelength_samples)
  emit_weights <- dnorm(integration_wavelengths, 
                        mean=emit_wavelengths[emit_i], 
                        sd=sd_emit[emit_i])
  asd_weights <- dnorm(integration_wavelengths, 
                       mean=asd_wavelengths[asd_j], 
                       sd=sd_asd[asd_j])
  return(sum(emit_weights*asd_weights))
}
# Apply function over all ASD and EMIT wavelength combinations to populate matrix
emit_weights <- unlist(lapply(1:length(asd_wavelengths), function(asd_ind){
  return(unlist(lapply(1:length(emit_wavelengths), getWeight, asd_j=asd_ind)))}))
# Convert from 1D vector to 2D matrix
emit_weights <- matrix(emit_weights, nrow=length(emit_fwhms), ncol=length(asd_fwhms))
# Scale each row (one EMIT wavelength) to have a summed value of 1
emit_weights <- emit_weights * (1/rowSums(emit_weights))




# ************* Apply weights to convolve ASD spectra to EMIT resolution ************* 

output_library <- terra::rast((as.matrix(asd_library, wide=TRUE)) %*% t(emit_weights))



# ************* Write Result to Disk ************* 

# Write new library data file to disk
#    Create ENVI file format and MESMA (Dar's) file format
terra::writeRaster(output_library, filename=output_spectrum_path, filetype="ENVI", overwrite=TRUE)
# Create a copy of the original ENVI file in .sli format for Viper Tools
file.copy(output_spectrum_path, 
          paste(output_spectrum_path, "_sli.sli", sep=""), 
          overwrite=TRUE)
file.copy(paste(output_spectrum_path, ".hdr", sep=""), 
          paste(output_spectrum_path, "_sli.hdr", sep=""), 
          overwrite=TRUE)

# ************* Write Result to Disk ************* 

# Copy metadata file for ViperTools
file.copy(paste(asd_data_path, "_sli.csv", sep=""), paste(output_spectrum_path, "_sli.csv", sep=""))
# Read in spectrum names from metadata
spectra_names <- read_csv(paste(output_spectrum_path, "_sli.csv", sep=""))$Name

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


