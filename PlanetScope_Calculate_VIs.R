###########################################################################################
#
#  this script calculate core vegetation indices from plantscope surface reflectance data 

#  this is the second step of planetscope data processing

#    --- Last updated:  2024.10.01 By Daryl Yang <yangd@ornl.gov>
###########################################################################################

#******************** close all devices and delete all variables *************************#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#*****************************************************************************************#

#****************************** load required libraries **********************************#
### install and load required R packages
list.of.packages <- c("ggplot2","raster", "readr", "TDPanalysis")  
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
#*****************************************************************************************#

#*************************************** load data ***************************************#
# define the directory to extracted or mosaicked planetscope surface reflectance data, 
# processed using 'PlanetScope_Clip_and_Mosaic.R
planetDIR <- ""
# search all hdf files stored in the directory
planet.LIST <- list.dirs(planetDIR, full.names = TRUE, recursive = FALSE)
# check the number of hdf files
print(paste0('number of images to process: ', length(planet.LIST)))
# print out the first a couple of hdf files to check
str(planet.LIST)
#*****************************************************************************************#

#****************************** calculate vegetation index *******************************#
for (folder in planet.LIST)
{
  print(folder)
  # search for the planetscape reflectance file
  planet.refl.DIR <- list.files(folder, pattern = 'mosaicqc.tif$', full.names = TRUE)
  if (length(planet.refl.DIR) > 0)
  {
    ### extract day of year from planetscope file name
    date <- strsplit(basename(planet.refl.DIR), "_")[[1]][2] # this need to revise accordingly
    date <- as.Date(date, "%Y%m%d")
    date <- gsub("-", '/', date)
    # convert date to doy
    doy <- date.to.DOY(date, format = 'yyyy/mm/dd')
    
    ### read in planetscope raster file
    refl.RST <- brick(planet.refl.DIR)
    
    ###test
    planet.cal.stack <- refl.RST/10000.0
    names(planet.cal.stack) <- c("blue", "green", "red", 'nir')
    
    ### calculate ndvi
    ndvi <- (planet.cal.stack$nir-planet.cal.stack$red)/
      (planet.cal.stack$nir+planet.cal.stack$red)
    outname <- gsub("mosaicqc", "ndvi", planet.refl.DIR)
    writeRaster(ndvi, outname, format="GTiff", overwrite=TRUE)
    
    ### calculate NIRv
    NIRv <- ndvi*planet.cal.stack$nir
    outname <- gsub("mosaicqc", "NIRv", planet.refl.DIR)
    writeRaster(NIRv, outname, format="GTiff", overwrite=TRUE)
    
    ### calculate EVI
    evi <- 2.5*(planet.cal.stack$nir-planet.cal.stack$red)/
      (planet.cal.stack$nir + 6*planet.cal.stack$red - 7.5*planet.cal.stack$blue + 1)
    outname <- gsub("mosaicqc", "evi", planet.refl.DIR)
    writeRaster(evi, outname, format="GTiff", overwrite=TRUE)
    
    ### calculate snow index
    snow <- (planet.cal.stack$red+planet.cal.stack$green+planet.cal.stack$blue)/3
    outname <- gsub("mosaicqc", "snow", planet.refl.DIR)
    writeRaster(snow, outname, format="GTiff", overwrite=TRUE)
    
    ### calculate NDVI
    ndgi <- (0.69*planet.cal.stack$green + 0.21*planet.cal.stack$nir - planet.cal.stack$red)/
      (0.69*planet.cal.stack$green + 0.21*planet.cal.stack$nir + planet.cal.stack$red)
    outname <- gsub("mosaicqc", "ndgi", planet.refl.DIR)
    writeRaster(ndgi, outname, format="GTiff", overwrite=TRUE)
  }
}
#*****************************************************************************************#




