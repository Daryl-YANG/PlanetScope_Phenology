###########################################################################################
#
#  this script creates time series vegetation index from planetscope.

#  this is the third step of planetscope data processing

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

#************************************ user parameters ************************************#
# define an output directory to store outputs
out.dir <- file.path("")
# create output directory if not exist
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)

# define the type of vegetation index to use for generating time series
veg.index <- 'ndgi'

#*****************************************************************************************#

#************************************ user define funs ***********************************#
# function that extract YYYY/MM/DD from file name list and convert to DOY
extr_doy <- function(dir.list)
{
  doy.year.target <- c()
  for(dir in dir.list)
  {
    date <- strsplit(basename(dir), "_")[[1]][2] # this need to revise accordingly
    date <- as.Date(date, "%Y%m%d")
    date <- gsub("-", '/', date)
    # convert date to doy
    doy <- TDPanalysis::date.to.DOY(date, format = 'yyyy/mm/dd')
    
    doy.year.target <- c(doy.year.target, doy)
  }
  return(doy.year.target)
}
#*****************************************************************************************#

#*************************************** load data ***************************************#
# define the directory to daily vegetation index processed using 'PlanetScope_Calculate_VIs.R'
planetDIR <- "/Volumes/Jupyter/Yang_etal_ERE/Scripts/KFC_PlanetScope/day_mosaic"
vi.DIRs <- list.files(planetDIR, pattern = paste0(veg.index, '.tif$'), 
                      full.names = TRUE, recursive = TRUE)
# check the number of hdf files
print(paste0('number of images to process: ', length(vi.DIRs)))
# print out the first a couple of hdf files to check
str(vi.DIRs)
#*****************************************************************************************#

#******************** create vegetation index time series stack **************************#
###### calculate ndvi time series
print("generating vegetation index time series")

vi.RSTs <- list()
for(i in 1:length(vi.DIRs)) { vi.RSTs[i] <- brick(vi.DIRs[i]) }

# create a reference raster
vi.RSTs$fun <- mean
vi.mosaic <- do.call(mosaic, vi.RSTs)

# stack vegetation index images
vi.resamp.RSTs <- list()
for (i in 1:(length(vi.RSTs)-1))
{
  viRST.resampled <- resample(vi.RSTs[[i]], vi.mosaic, method = "bilinear")
  vi.resamp.RSTs[i] <- viRST.resampled
}
vi.stack <- do.call(stack, vi.resamp.RSTs)

outname <- paste0(out.dir, '/', veg.index, '_time_series.tif')
writeRaster(vi.stack, outname, format="GTiff", overwrite=TRUE)

###### export day of year list for layers in the vi stack
doy.list <- extr_doy(vi.DIRs)
doy.list <- data.frame(doy.list)
colnames(doy.list) <- c('doy')
outname <- paste0(out.dir, '/', 'doy.csv')
write.csv(doy.list, outname)
#*****************************************************************************************#