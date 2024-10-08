###########################################################################################
#
#  this script creates time series vegetation index from planetscope. Note: this script 
#  provides the flexibility to use the mean vegation index of nearby years if no observation'
#  were found during a defined peroid of current year and observation were acquired for
#  other years
#
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
list.of.packages <- c("ggplot2","raster", 'TDPanalysis', 'plyr')  
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
# create a temporarly directory to store compiled VI across years
temp.dir <- file.path(out.dir, 'complied_VI')
if (! file.exists(temp.dir)) dir.create(temp.dir,recursive=TRUE)

# define the years that will be used to borrow data from
years <- c(2017, 2018, 2019, 2020, 2021, 2022, 2023)
# define the year to produce vegetation index time series
year.target <- 2022

# define the date range to build vegetation index time series
day.start <- "2022/04/01"
doy.start <- TDPanalysis::date.to.DOY(day.start, format = 'yyyy/mm/dd')
day.end <- "2022/10/31"
doy.end <- TDPanalysis::date.to.DOY(day.end, format = 'yyyy/mm/dd')
# define max bin size that an observation is needed to build vegetation time series
binsize <- 5

# define the type of vegetation index to use
veg.index <- 'ndgi'

# define site name, this need to be the same with that used in file names
site.name <- 'kfc'

# define output project system
out.trs <- '+proj=utm +zone=03 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'
#*****************************************************************************************#
#*
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
# define file path that contains multiple year vegetation indices. The folder structure
# should be year/day_moasic/..
data.dir <- file.path("")
# search all data contained in the years define in "years"
data.list <- list()
for (year in years)
{
  year.data.dir <- paste0(data.dir, "/", year, "/", 'day_mosaic') # revised folder name
  year.data.lst <- list.dirs(year.data.dir, full.names = TRUE, recursive = FALSE)
  data.list <- c(data.list, list(files = year.data.lst)) 
}
names(data.list) <- as.character(years)
#*****************************************************************************************#

#************************** compile time series VI data  *********************************#
# pull out data list for the target year and create a list of data collection doy
year.target.data.dirs <-ldply(data.list[as.character(year.target)], data.frame)[,2]
### extract day of year that data were collected in the target year
doy.year.target <- extr_doy(year.target.data.dirs)

bin.doy.start <- doy.start
while(bin.doy.start < doy.end)
{
  # detemine the ending date of bin
  bin.doy.end <- bin.doy.start + binsize
  # find file that where collect within the start and end of current bin
  bin.doy.match <- which(doy.year.target > bin.doy.start & 
                           doy.year.target <= bin.doy.end)
  # if data were collected during the bin time range, then pull out the data
  if (length(bin.doy.match) > 0) 
  {
    bin.dir.list <- year.target.data.dirs[bin.doy.match]
    for (dir in bin.dir.list)
    {
      vi.dir <- list.files(dir, pattern = paste0(veg.index, '.tif$'), full.names = TRUE)
      file.copy(vi.dir, temp.dir, overwrite = TRUE)
      print(vi.dir)
    }
  }
  # if no data were collected during the bin time range, then borrow data from
  # other years
  if (length(bin.doy.match) == 0)
  {
    # search data collected in the same doy range from other years
    doy.match.list <- c()
    for (year in years)
    {
      year.data.dirs <-ldply(data.list[as.character(year)], data.frame)[,2]
      year.doy.list <- extr_doy(year.data.dirs)
      year.doy.match <- which(year.doy.list > bin.doy.start & 
                                year.doy.list <= bin.doy.end)
      year.doy.dir <- year.data.dirs[year.doy.match]
      doy.match.list <- c(doy.match.list, year.doy.dir)
    }
    
    if (length(doy.match.list) == 1)
    {
      match.dir <- list.files(doy.match.list, pattern = paste0(veg.index, '.tif$'), 
                              full.names = TRUE)
      vi.rst <- brick(match.dir)
      outdate <- as.Date((bin.doy.start+bin.doy.end)/2, 
                         origin = paste0(year.target, "-01-01"))
      date <- gsub("-", '', outdate)
      out.name <- paste0(temp.dir, '/', site.name, '_', date, '_', veg.index, '.tif')
      writeRaster(vi.rst, out.name, format = 'GTiff', overwrite = TRUE)
    }
    
    if (length(doy.match.list) > 1)
    {
      match.dirs <- c()
      for (dir in doy.match.list)
      {
        match.dir <- list.files(dir, pattern = paste0(veg.index, '.tif$'), 
                                full.names = TRUE)
        match.dirs <- c(match.dirs, match.dir)
      }
      vi.RSTs <- list()
      for(i in 1:length(match.dirs)) { vi.RSTs[i] <- brick(match.dirs[i]) }
      # stack all image
      # create a reference raster
      vi.RSTs$fun <- mean
      vi.mosaic <- do.call(mosaic,vi.RSTs)
      
      # stack vegetation index images
      vi.resamp.RSTs <- list()
      for (i in 1:(length(vi.RSTs)-1))
      {
        viRST.resampled <- resample(vi.RSTs[[i]], vi.mosaic, method = "bilinear")
        vi.resamp.RSTs[i] <- viRST.resampled
      }
      vi.stack <- do.call(stack,vi.resamp.RSTs)
      
      # calculate mean vi of all raster layers
      doyMeanRST <- try(calc(vi.stack, fun = mean, na.rm = T, trim = 0.1), silent = TRUE)
      
      outdate <- as.Date((bin.doy.start+bin.doy.end)/2, 
                         origin = paste0(year.target, "-01-01"))
      date <- gsub("-", '', outdate)
      
      out.name <- paste0(temp.dir, '/', site.name, '_', date, '_', veg.index, '.tif')
      writeRaster(doyMeanRST, out.name, format = 'GTiff', overwrite = TRUE)
    }
  }
  bin.doy.start <- bin.doy.start + binsize
}
#*****************************************************************************************#

#******************** create vegetation index time series stack **************************#
###### calculate ndvi time series
print("generating vegetation index time series")
vi.DIRs <- list.files(temp.dir, pattern = paste0(veg.index, '.tif$'), 
                        full.names = TRUE, recursive = TRUE)
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

###### remove the temporary folder
unlink(temp.dir, recursive = T, force = T)
#*****************************************************************************************#






























