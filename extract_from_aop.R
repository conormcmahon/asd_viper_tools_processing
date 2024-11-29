
# Extract endmember spectra from AOP spectra at target polygons (for field-surveyed vegetation, etc.)

# Required libraries
library(rhdf5)
library(terra)
library(sf)
library(tidyverse)



# ******************** User Parameters ********************

# NEON Sitename
neon_sitename <- "SRER"
# AOP Spectra Directory
dir_aop <- "D:/SERDP/SRER/AOP_Spectra/NEON_refl-surf-dir-ortho-mosaic/NEON.D14.SRER.DP3.30006.001.2021-09.basic.20241121T232913Z.RELEASE-2024/"
# Field polygon filepath
filepath_polygons <- "D:/SERDP/SRER/veg_polygons/tree_polygons.gpkg"
# Internal Buffer - buffer polygons inward by this amount to avoid edge effects 
internal_buffer <- 1 # units are meters
# Output Directory
output_dir <- "D:/SERDP/SRER/AOP_Spectra/"


# ******************** Load source data ********************

# Get all .h5 data files in the target directory (tiled up)
aop_h5_files <- list.files(dir_aop, ".h5$") # the $ character specifies an end-of-line, excluding .xml files

# Get metadata from one example file...
# Specify a target .h5 file
aop_h5_file_example <- paste(dir_aop, aop_h5_files[[1]], sep="")
# View file structure
View(h5ls(aop_h5_file_example,all=T))

# Get reflectance data dimensions
reflectance_info <- h5readAttributes(aop_h5_file_example,
                   paste(neon_sitename, "/Reflectance/Reflectance_Data", sep=""))
n_rows <- reflectance_info$Dimensions[1]
n_cols <- reflectance_info$Dimensions[2]
n_bands <- reflectance_info$Dimensions[3]
# Get EPSG code for the AOP data
epsg_AOP <- as.numeric(h5read(aop_h5_file_example, 
                              "/SRER/Reflectance/Metadata/Coordinate_System")$'EPSG Code')
# Get extent of raster
x_min <- reflectance_info$Spatial_Extent_meters[1]
x_max <- reflectance_info$Spatial_Extent_meters[2]
y_min <- reflectance_info$Spatial_Extent_meters[3]
y_max <- reflectance_info$Spatial_Extent_meters[4]
tile_width <- x_max - x_min
tile_height <- y_max - y_min
raster_extent <- ext(x_min, x_max, y_min, y_max)

# Load polygons
polygons <- st_read(filepath_polygons)
# Reproject to same basis as AOP data
polygons <- st_transform(polygons, 
                         crs=epsg_AOP) %>%
  mutate(fid = 1:nrow(polygons))
# Buffer polygons inward to prevent edge effects  
buffered_polygons <- st_buffer(polygons, -internal_buffer)
buffered_polygons$area <- as.numeric(st_area(buffered_polygons))
buffered_polygons <- buffered_polygons %>%
  filter(area > 5)

# Load manually-generated metadata file
polygon_metadata <- read_csv("D:/SERDP/SRER/veg_polygons/tree_polygons_metadata.csv") %>% 
  filter(fid %in% buffered_polygons$fid)


# ******************** Get geospatial information from example AOP tile ********************

# Get minimum X and Y values of polygons
min_x_polygons <- (ext(polygons)$xmin %/% tile_width) * tile_width
max_x_polygons <- (ext(polygons)$xmax %/% tile_width + 1) * tile_width
min_y_polygons <- (ext(polygons)$ymin %/% tile_width) * tile_width
max_y_polygons <- (ext(polygons)$ymax %/% tile_width + 1) * tile_width
tile_x_range <- max_x_polygons - min_x_polygons
tile_y_range <- max_y_polygons - min_y_polygons
# Get large raster of tile objects
tile_layout <- rast(matrix(1:((tile_x_range)/tile_width * (tile_y_range)/tile_height), 
                           ncol=(tile_x_range)/tile_width,
                           nrow=(tile_y_range)/tile_height))
ext(tile_layout) <- ext(min_x_polygons, max_x_polygons, min_y_polygons, max_y_polygons)

# Get intersection of polygons with hyperspectral raster
tile_intersection_indices <- unique(terra::extract(tile_layout, polygons)$lyr.1) %>% sort()
tile_intersections <- tile_layout %in% tile_intersection_indices
tile_intersections_x <- min_x_polygons + (((tile_intersection_indices-1)*tile_height) %/% tile_y_range) * tile_height
tile_intersections_y <- max_y_polygons - (((tile_intersection_indices-1)*tile_height) %% tile_y_range) - tile_height



# ******************** Output example imagery from AOP ********************

# Function to build selected bands from a file into a Terra object
buildAndStack <- function(img_filepath, bands)
{
  img <- lapply(1:length(bands), function(ind)
  {
    return(terra::rast(t(matrix(h5read(img_filepath,
                                       "/SRER/Reflectance/Reflectance_Data",
                                       index=list(bands[ind], 1:tile_width, 1:tile_height)),
                                nrow=tile_width, ncol=tile_height)),
                       crs=crs(paste("+init=epsg:", epsg_AOP, sep=""))))
  })
  img <- terra::rast(img)
  ext(img) <- raster_extent
  return(img)
}

# Generate a SWIR-NIR-R falsecolor
falsecolor_example <- buildAndStack(paste(dir_aop, aop_h5_files[[2]], sep=""), c(61, 101, 241))
falsecolor_example[falsecolor_example > 5000] <- 0
plotRGB(falsecolor_example, r=3, g=2, b=1, scale=5000)

# Write an example tile with all wavelengths as a geotiff
full_spectrum_example <- buildAndStack(paste(dir_aop, aop_h5_files[[2]], sep=""), 1:n_bands)
full_spectrum_example[full_spectrum_example > 5000] <- 0
ext(full_spectrum_example) <- ext(x_min+1000, x_max+1000, y_min, y_max)
terra::writeRaster(full_spectrum_example, 
                   paste(dir_aop, str_split(aop_h5_files[[2]], "\\.")[[1]][[1]], "_out.tif"),
                   overwrite=TRUE)




# ******************** Sample spectra from within polygons ********************

sampleSpectra <- function(target_index)
{
  tile_x <- tile_intersections_x[[target_index]]
  tile_y <- tile_intersections_y[[target_index]]
  file_matches <- grepl(paste("*", tile_x, "_", tile_y, "_reflectance.h5", sep=""), 
                    aop_h5_files)
  if(sum(file_matches) == 0)
  {
    print(paste("The combination ", tile_x, "_", tile_y, " has no associated spectra.", sep=""))
    return(-1)
  }
  
  print(paste("Found a file for combination ", tile_x, "_", tile_y, sep=""))
  
  aop_file <- paste(dir_aop, aop_h5_files[which(file_matches)[[1]]], sep="")
  
  # Image where each value is the column index, with extent and size of the spectral imagery 
  col_img <- rast(t(matrix(rep(1:tile_width, tile_height),
                         ncol=tile_width, nrow=tile_height)))
  ext(col_img) <- ext(tile_x, tile_x+tile_height, tile_y, tile_y+tile_height)
  names(col_img) <- "col"
  # Image where each value is the row index, with extent and size of the spectral imagery 
  row_img <- rast(matrix(rep(1:tile_width, tile_height),
                         ncol=tile_width, nrow=tile_height))
  ext(row_img) <- ext(tile_x, tile_x+tile_height, tile_y, tile_y+tile_height)
  names(row_img) <- "row"

  # Extract indices of target cells which are overlapped by polygons
  indices <- (terra::extract(c(col_img, row_img), 
                             buffered_polygons[,'fid'])) %>%
    drop_na(row, col) %>%
    mutate(tile_index = ind,
           tile_x = tile_x,
           tile_y = tile_y)
  
  # Get values of spectra at each index
  output_spectra <- lapply(1:nrow(indices),
                           function(ind){
                             h5read(aop_file,
                                    "/SRER/Reflectance/Reflectance_Data",
                                    index=list(1:n_bands,indices[ind,]$col,indices[ind,]$row))[,1,1]
                           })
  return(list(output_spectra, indices))
}

# Get all spectra
all_spectral_samples <- lapply(1:length(aop_h5_files), sampleSpectra)
# Drop cases where there were no spectra (a few polygons are outside AOP extent)
all_spectral_samples <- all_spectral_samples[unlist(lapply(all_spectral_samples, length)) != 1]
# Combine everything into one dataframe of spectra
all_spectra <- bind_rows(lapply(all_spectral_samples,
                                function(newdata){ return(as.data.frame(do.call(rbind, newdata[[1]]))) }))
all_metadata <- bind_rows(lapply(all_spectral_samples,
                                 function(newdata){ return(as.data.frame(newdata[[2]])) }))




# ******************** Extract clusters within each polygon's set of spectra ********************

# List of the polygon indices which will be used in output (which intersect at least some spectral tiles)
output_IDs <- unique(all_metadata$ID) %>% sort()

# Number of clusters to extract center values for within each polygon
clusters_per_polygon <- 3

# Get list of bands which are contaminated with water vapor effects - exclude these from cluster analysis
bad_indices <- c(195:207, 287:312, 420:426)

# For each polygon, extract a list of cluster centers
getCenters <- function(ind)
{
  # List of all spectrum indices which correspond to the target polygon index
  spectral_indices <- which(all_metadata$ID == ind)
  # Extract those spectra and change to matrix
  target_spectra <- as.matrix(all_spectra[spectral_indices,])
  target_spectra[,bad_indices] <- 0
  # Apply clustering algorithm to those spectra
  clusters <- kmeans(target_spectra, clusters_per_polygon)
  return(clusters$centers)
}
# Function to create a simple plot of multiple spectra from one polygon
plotCenters <- function(centers)
{
  plot(1:ncol(centers), centers[1,], ylim=c(0, 5000))
  if(nrow(centers) > 1)
    for(ind in 1:(nrow(centers)-1))
      points(1:ncol(centers), centers[ind,])
}
# Generate all cluster spectra
spectral_clusters <- lapply(unique(all_metadata$ID) %>% sort(),
                            getCenters) %>%
  do.call(what=rbind)

# Function to triple metadata and add a new identifier for each cluster center
tripleMetadata <- function(metadata_index)
{
  if(clusters_per_polygon > 26)
    warning("NOTE this function will probably not label files correctly because the number of clusters is > 26 (more than letters in alphabet)")
  newdata <- polygon_metadata %>% 
    filter(fid == buffered_polygons[metadata_index,]$fid)
  newdata$Name <- paste(newdata$Name, "a", sep="")
  for(ind in 2:clusters_per_polygon)
  {
    newdata <- rbind(newdata, 
                     polygon_metadata %>% 
                       filter(fid == buffered_polygons[metadata_index,]$fid) %>%
                       mutate(Name = paste(Name, letters[ind], sep="")))
  }
  return(newdata)
}
spectral_metadata <- lapply(unique(all_metadata$ID) %>% 
                              sort(),
                            tripleMetadata) %>%
  bind_rows()

# Write output raster and metadata
terra::writeRaster(terra::rast(spectral_clusters),
                   paste(output_dir, "spectra.tif", sep=""), 
                   overwrite=TRUE)
write_csv(spectral_metadata,
          paste(output_dir, "metadata.csv", sep=""))




