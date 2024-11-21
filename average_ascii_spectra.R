
# Read in ENVI format data from ASD, output from PRISM USGS software
# Convert to matrix format 
# Average across adjacent spectra and subset to only 'selected' data (after manual inspection with PRISM)
# Output averaged target spectra to Viper Tools .sli format
# Create associated metadata


# ************* Dependencies ************* 

# Load relevant packages
library(terra)
library(tidyverse)
library(janitor)


# ************* Filepath Inputs ************* 

# Input library in tif format (data file, not the .hdr file)
input_library_path <- "D:/SERDP/SRER/ASD_Data/ENVI/AZ_spectra_all_tif.tif"
# Spectral library header file
input_header_path <- "D:/SERDP/SRER/ASD_Data/ENVI/AZ_spectra_all.hdr"
# Metadata for input library
input_metadata_filepath <- "D:/SERDP/SRER/ASD_Data/ENVI/sli_metadata.csv"
# Output filepath for completed library
output_filepath <- "D:/SERDP/SRER/ASD_Data/averaged_data"
# Output spectrum name prefix - add to start of each individual spectrum name in metadata
output_spectrum_name_prefix <- "AZ_"

# ************* Load Input Data ************* 

# Read ENVI library
input_library <- terra::rast(input_library_path)
# Read header 
header <- readLines(input_header_path)
# Read metadata
input_metadata <- read_csv(input_metadata_filepath) %>%
  janitor::clean_names()
input_metadata[input_metadata$day==25,]$series <- ""
input_metadata <- input_metadata %>%
  mutate(spectrum_name = paste(substr(year, 3,4),
                               sprintf("%02d", month),
                               day,
                               series, 
                               sep=""))


# ************* Spectrum Metadata ************* 

# Get list of names of spectra from the header
#    Use regex to get start and end of list of spectra names
start_ind_spectra_names <- which(grepl("spectra names", header)) + 1
closing_bracket_indices <- which(grepl("}", header))
last_ind_spectra_names <- last_ind_spectra_names <- closing_bracket_indices[(closing_bracket_indices >= start_ind_spectra_names)][1]
#    Subset to those rows of the header file
spectra_names_messy <- header[start_ind_spectra_names:last_ind_spectra_names]
#    Reprocess header file text lines into short string names
spectra_names_long <- substr(unlist(lapply(str_split(spectra_names_messy, 
                                                     ".asd"), 
                                           function(str_list){
                                             return(str_list[[1]])})), 
                             2,1000)
#    Get the datestamp and numeric index of each spectrum
spectra_names <- substr(spectra_names_long, 
                        1, 
                        ((nchar(spectra_names_long))-5))
spectra_numbers <- as.numeric(substr(spectra_names_long, 
                                     ((nchar(spectra_names_long))-4), 
                                     nchar(spectra_names_long)))

# Get list of wavelengths of spectra from the header
#    Use regex to get start and end of list of spectra names
start_ind_spectra_wavelengths <- which(grepl("wavelength = ", header)) + 1
last_ind_spectra_wavelengths <- last_ind_spectra_wavelengths <- closing_bracket_indices[(closing_bracket_indices >= start_ind_spectra_wavelengths)][1]
#    Subset to those rows of the header file
spectra_wavelengths_messy <- header[start_ind_spectra_wavelengths:last_ind_spectra_wavelengths]
#    Reprocess header file text lines into short string names
spectra_wavelengths_long <- substr(unlist(lapply(str_split(spectra_wavelengths_messy, 
                                                           ".asd"), 
                                                 function(str_list){
                                                   return(str_list[[1]])})), 
                                   2,1000)
#    Get the wavelengths for spectra
spectra_wavelengths <- substr(spectra_wavelengths_long, 
                              1, 
                              ((nchar(spectra_wavelengths_long))-1))

# ************* Filtering Bad Spectra ************* 

# Get list of input spectra which are marked 'keep' in metadata
input_metadata_retained <- input_metadata %>% 
  filter(consider == 1)
# Double-check that each retained data file has a unique name
print("")
print("Number of retained records by class: ")
print(input_metadata_retained %>% group_by(class) %>% tally() %>% arrange(-n))
print(paste("Number of unique spectrum names in retained list: ", length(unique(input_metadata_retained$name)), sep=""))
print(paste("Total number of spectra in retained list: ", length(input_metadata_retained$name), sep=""))


# ************* Extract Spectrum Data ************* 

# Function to retrieve a spectral group from metadata
getSpectra <- function(target_metadata, output_plot=FALSE)
{
  target_spectra_indices <- which(((spectra_names == target_metadata$spectrum_name) * 
                                     (spectra_numbers >= target_metadata$start) * 
                                     (spectra_numbers <= target_metadata$end)) == 1)
  target_spectra <- lapply(target_spectra_indices,
                           function(index){
                             as.vector(input_library[index,][[1]])
                           })
  target_spectra <- matrix(unlist(target_spectra), ncol=length(target_spectra_indices))
  
  if(output_plot)
  {
    plot(0,0, xlim=c(0,2151), ylim=c(0,1))
    for(ind in 1:length(target_spectra_indices))
      lines(target_spectra[,ind])
  }
  
  return(target_spectra)
}

# Run the function to extract averaged data for all spectra which are marked as 'keep'
#   Averages across each set of 5 adjacent samples which are collected together
selected_spectra_data <- lapply(1:nrow(input_metadata_retained),
                                function(ind){
                                  rowMeans(getSpectra(input_metadata_retained[ind,]))
                                })


# ************* Build Raster Dataset ************* 

# Generate new output library using the input as a template, subsetting based on new number of records
output_library <- input_library[1:length(selected_spectra_data),,drop=FALSE]
# Replace each record with correct data, spectrum by spectrum
for(ind in 1:length(selected_spectra_data)) {
  output_library[ind,] <- selected_spectra_data[[ind]] 
}
# Replace impossible reflectances (vary from 0 to 1)
output_library[output_library < 0] <- 0
output_library[output_library > 1] <- 1

# Write new library data file to disk
#    Create ENVI file format and MESMA (Dar's) file format
terra::writeRaster(output_library, paste(output_filepath, "/AZ_spectral_library", sep=""), filetype="ENVI", overwrite=TRUE)
# Create a copy of the original ENVI file in .sli format for Viper Tools
file.copy(paste(output_filepath, "/AZ_spectral_library", sep=""), paste(output_filepath, "/AZ_spectral_library_sli.sli", sep=""), overwrite=TRUE)
file.copy(paste(output_filepath, "/AZ_spectral_library.hdr", sep=""), paste(output_filepath, "/AZ_spectral_library_sli.hdr", sep=""), overwrite=TRUE)


# Create metadata file for ViperTools
output_metadata <- input_metadata_retained %>%
  mutate(Name = paste(output_spectrum_name_prefix, year, "_", month, "_", day, "_", name, sep="")) %>%
  dplyr::select(Name, class, subclass, life_history, binomial, subject, photograph)
write_csv(output_metadata, paste(output_filepath, "/AZ_spectral_library_sli.csv", sep=""))

# Add spectra names to header file (.hdr)
#    Read auto-generated header 
output_header <- readLines(paste(output_filepath, "/AZ_spectral_library_sli.hdr", sep=""))
#    Generate well-formatted list of spectra names for output
spectra_names_lines <- paste("spectra names = {\n", paste(output_metadata$Name, collapse=', \n'), "}", sep="")
spectra_wavelengths_lines <- paste("wavelength = {\n", paste(spectra_wavelengths, collapse=', \n'), "}", sep="")
#    Correct a few header data lines to the Spectral Library format
output_header[8] <- "file type = ENVI Spectral Library"
output_header[12] <- "sensor_type = Unknown"
output_header[15] <- "wavelength units = Micrometers"
#    Add band names and wavelengths to header file
output_header <- c(output_header, 
                   spectra_names_lines,
                   spectra_wavelengths_lines)
#    Add list of spectra names to end of text file and write to disk
writeLines(output_header, 
           paste(output_filepath, "/AZ_spectral_library_sli.hdr", sep=""))



