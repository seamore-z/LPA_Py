#---------------------------------------------------------------------
#Load packages
#---------------------------------------------------------------------
library(sf)
library(terra)
library(imager)   #needed for efficient distance to snow calculate
library(ncdf4)

library(iterators)
library(foreach)
library(doMC)

library(matrixStats)   
library(WGCNA)
library(zoo)

library(RcppRoll)
library(Rcpp)

library(rjson)
library(XML)
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#Sort out image chunk boundaries, edited to name list items chunkStart
#and chunkEnd
#Douglas Bolton, edited by Seamore Zhu
#---------------------------------------------------------------------
chunkBoundaries <-function(numPix, numChunks) {
  
  lenChunk <- ceiling(numPix / numChunks) #Length of each chunk
  
  chunkStart <- matrix(0,numChunks,1)
  chunkEnd <- matrix(0,numChunks,1)
  for (n in 1:numChunks) {
    chunkStart[n] <- lenChunk*(n-1) + 1
    chunkEnd[n] <- lenChunk*n}
  
  chunkEnd[numChunks] <- numPix
  boundaries <- list('chunkStart' = chunkStart, 'chunkEnd' = chunkEnd)
  
  return(boundaries)
}


#---------------------------------------------------------------------
#Reconstruct output images from chunks (.Rds) in order to make netCDF
#composite image for each given phenometric
#Seamore Zhu
#---------------------------------------------------------------------
readChunks <- function(numChunks, numPix, tempDir, phenometric) {
  
  boundaries <- chunkBoundaries(numPix, numChunks)  #Get boundaries
  
  #Build matrix for the full image from matrices for each chunk
  mat <- matrix(NA,numPix,6)
  for (n in 1:numChunks) {
    matSub <- matrix(NA,boundaries$chunkEnd[n]-boundaries$chunkStart[n]+1,6)
    fileName <- paste0(tempDir,'c',n,'_p',phenometric,'.Rds')
    matSub <- try(readRDS(fileName),silent=T)
    if (inherits(matSub, 'try-error')) {print('Oh no')}
    mat[boundaries$chunkStart[n]:boundaries$chunkEnd[n],] <- matSub
  }
  
  return(mat)
}


#---------------------------------------------------------------------
#Get spatial info from HLS file
#Will need in order to write spatial attribute info to netcdf file
#Written by Douglas Bolton
#---------------------------------------------------------------------
getNetCDF_projection_info <- function(baseImage) {
  
  #Get extent, and then define pixel centers in the x and y direction
  ext = ext(baseImage)
  res = res(baseImage)[1]
  if (res == 60) {res <-  10}              #Error in metadata for S10, where resolutions is written as 60m. Change to 10m.  
  
  x = seq(ext[1]+res/2,ext[1]+ncol(baseImage)*res, res)
  y = seq(ext[3]+res/2,ext[3]+nrow(baseImage)*res, res)
  
  #Define dimensions for netCDF file
  dimx = ncdim_def(name = 'x', longname = 'x coordinate', units='m', vals = as.double(x))
  dimy = ncdim_def(name = 'y', longname = 'y coordinate', units='m', vals = rev(as.double(y)))
  
  #Get projection in wkt format
  # wkt <- showWKT(projection(baseImage)) 
  wkt <- st_as_text(st_crs(baseImage))
  
  #Need to pull the central meridian from the wkt 
  spt <- unlist(strsplit(gsub(']','',wkt),','))
  central_meridian <- as.numeric(spt[which(spt == "PARAMETER[\"central_meridian\"")+1])
  
  #Need "geoTransform" data for netcdf (top left corner coordinates and resolution)
  geoTransform <- paste(ext[1],res,0,ext[4],0,-res)
  
  return(list(dimx=dimx,dimy=dimy, wkt=wkt, central_meridian=central_meridian, geoTransform=geoTransform))
}


#---------------------------------------------------------------------
#For an image chunk, generate composites for each phenometric (OGI, 
#50PGCI, OGMx, Peak, OGD, 50PCGD, and OGMn) from the mean of images 
#within +/- an interval of the mean phenometric for the entire image
#Seamore Zhu
#---------------------------------------------------------------------
DoNonvegComp <- function(b2, b3, b4, b5, b6, b7, dates, year, params){
  
  #Initialize dates from an image chunk
  yrdoy_char = as.character(format(dates, "%Y%j"))   #Returning to yrdoy from dates
  yrs = substring(yrdoy_char, first = 1, last = 4)   #Get years
  doys = substring(yrdoy_char, first = 5, last = 7)  #Get doys
  
  phen_composites <- vector(mode = "list", length = 7)  #Instantiate list of phenometric composites
  
  #Generate mean phenometrics from the full image for the given year
  phenoPath <- paste0(params$dirs$phenDir,'MSLSP_',tile,'_',as.character(year),'.nc')
  phenoImg <- nc_open(phenoPath)
  #Get 7 DOY phenometrics
  OGI <- ncvar_get(phenoImg, 'OGI'); OGI50 <- ncvar_get(phenoImg, '50PCGI')
  OGMx <- ncvar_get(phenoImg, 'OGMx'); Peak <- ncvar_get(phenoImg, 'Peak')
  OGD <- ncvar_get(phenoImg, 'OGD'); OGD50 <- ncvar_get(phenoImg, '50PCGD')
  OGMn <- ncvar_get(phenoImg, 'OGMn')
  nc_close(phenoImg)  #Close nc file
  #Mean of phenometrics
  OGImean <- round(mean(OGI,na.rm=T)); OGI50mean <- round(mean(OGI50,na.rm=T))
  OGMxmean <- round(mean(OGMx,na.rm=T)); Peakmean <- round(mean(Peak,na.rm=T))
  OGDmean <- round(mean(OGD,na.rm=T)); OGD50mean <- round(mean(OGD50,na.rm=T))
  OGMnmean <- round(mean(OGMn,na.rm=T))
  meanpheno <- c(OGImean,OGI50mean,OGMxmean,Peakmean,OGDmean,OGD50mean,OGMnmean)
  
  #Filter chunked images to year of interest
  phen_yr_i = which(yrs %in% as.character(year))
  phen_doys = doys[phen_yr_i]
  b2_y <- b2[,phen_yr_i]; b3_y <- b3[,phen_yr_i]; b4_y <- b4[,phen_yr_i]
  b5_y <- b5[,phen_yr_i]; b6_y <- b6[,phen_yr_i]; b7_y <- b7[,phen_yr_i]
  
  # Add composites for each phenometric
  i<-1
  for (phenometric in meanpheno) {
    #Filter chunked images to dates of interest
    comp_doys_i = which(as.integer(phen_doys) %in% c((phenometric-7):(phenometric+7)))  #+/-7 days
    b2_r <- b2_y[,comp_doys_i]; b3_r <- b3_y[,comp_doys_i]; b4_r <- b4_y[,comp_doys_i]
    b5_r <- b5_y[,comp_doys_i]; b6_r <- b6_y[,comp_doys_i]; b7_r <- b7_y[,comp_doys_i]
    b2_r[b2_r==32767]<-NA;b3_r[b3_r==32767]<-NA;b4_r[b4_r==32767]<-NA
    b5_r[b5_r==32767]<-NA;b6_r[b6_r==32767]<-NA;b7_r[b7_r==32767]<-NA
    
    #Composites for each band
    b2_c <- round(rowMeans(b2_r,na.rm=TRUE)); b3_c <- round(rowMeans(b3_r,na.rm=TRUE)); b4_c <- round(rowMeans(b4_r,na.rm=TRUE))
    b5_c <- round(rowMeans(b5_r,na.rm=TRUE)); b6_c <- round(rowMeans(b6_r,na.rm=TRUE)); b7_c <- round(rowMeans(b7_r,na.rm=TRUE))
    composite <- cbind(b2_c, b3_c, b4_c, b5_c, b6_c, b7_c)
    phen_composites[[i]] <- composite  #Add composite to list
    i<-i+1
  }
  
  return(phen_composites)
}


#---------------------------------------------------------------------
#Run composite code for an image chunk
#Write composite results for each chunk to disk
#Seamore Zhu
#---------------------------------------------------------------------
runNonvegComposite <- function(chunk, numPix, imgYrs, phenYrs, errorLog, params) {
  
  #pheno_pars <- params$phenology_parameters
  
  #Get all images to process
  ######################
  chunkFold <- paste0(params$dirs$chunkDir,'c',chunk,'/')  # THIS MAY HAVE TO BE UPDATED IF USING NEW CHUNKS BASED ON PREPROCESSING THAT DOES NOT REMOVE WATER
  outFold <-   paste0(params$dirs$tempDir,'outputs/')  
  
  imgList <- list.files(path=chunkFold, pattern=glob2rx("HLS_*.Rds"), full.names=F)
  numImgs = length(imgList)
  
  yrdoy = as.numeric(matrix(NA,numImgs,1))
  for (i in 1:numImgs) {
    imName <- gsub('.Rds','',imgList[i])
    yrdoy[i]= as.numeric(unlist(strsplit(imName,'_',fixed = T))[4])}
  
  ord = order(yrdoy)  #Determine image order
  
  #Read in all imagery for chunk
  b2 <- matrix(NA,numPix,numImgs);  b3 <- matrix(NA,numPix,numImgs);  b4 <- matrix(NA,numPix,numImgs)
  b5 <- matrix(NA,numPix,numImgs);  b6 <- matrix(NA,numPix,numImgs);  b7 <- matrix(NA,numPix,numImgs)
  
  for (i in 1:length(ord)) {
    img <- imgList[ord[i]]
    imgData <- try(matrix(readRDS(paste0(chunkFold,img)),nrow=numPix),silent = TRUE)
    if (inherits(imgData, 'try-error')) {cat(paste('runPhenoChunk: Error for chunk',chunk,img), file=errorLog, append=T);next} 
    
    b2[,i] <- imgData[,1]; b3[,i] <- imgData[,2]; b4[,i] <- imgData[,3]
    b5[,i] <- imgData[,4]; b6[,i] <- imgData[,5]; b7[,i] <- imgData[,6]
    
    remove(imgData)
  }
  
  dates = as.Date(strptime(yrdoy[ord], format="%Y%j"))  #Format as dates
  
  #For each year of imagery/phenology, create phenometric composites
  all_phenocomps <- vector(mode = "list", length = 7)
  for (year in phenYrs) {
    all_phenocomps <- DoNonvegComp(b2,  b3,  b4,  b5,  b6, b7, dates, year)
    
    #Output chunk phenometric images as c[chunk number]_p[phenometric number].Rds
    j<-1
    for (phenocomp in all_phenocomps) {
      yearcheck <- paste0(outFold,'y',year)
      if(!dir.exists(yearcheck)){dir.create(yearcheck)}
      saveRDS(phenocomp, paste0(yearcheck,'/c',chunk,'_p',j,'.Rds'))
      j<-j+1
    }
  }
  
}


#---------------------------------------------------------------------
#Function to update netCDF file for the product by updating 
#phenometric composites to include areas with no phenology
#The layers to be included in the product are defined in MuSLI_LSP_V1_Layers.csv
#Seamore Zhu
#---------------------------------------------------------------------
createComposite <- function(tile, numChunks, numPix, baseImage, phenYrs, params) {
  outFold <-   paste0(params$dirs$tempDir,'outputs/')  #where image chunks are stored
  prj_info <- getNetCDF_projection_info(baseImage)     #getting projection info from baseImage
  
  for (year in phenYrs) {
    phenoPath <- paste0(params$dirs$phenDir,'MSLSP_',tile,'_',as.character(year),'.nc')
    phenoImg <- nc_open(phenoPath, write=TRUE, verbose=FALSE)
    #For each phenometric, reconstruct full composite image, separate out bands, and write to netCDF product file
    for (i in c(1:7)) {
      full_pheno_mat <- readChunks(numChunks, numPix, paste0(outFold,'y',as.character(year),'/'), i)
      for (j in c(1:6)) {
        full_pheno_band <- matrix(full_pheno_mat[,j], dim(baseImage)[1],dim(baseImage)[2])
        full_pheno_band[full_pheno_band < -32767 | full_pheno_band > 32767 | is.na(full_pheno_band)] <- 32767  #Fix incorrect values
        new_var = ncvar_def(paste0('p',i,'_b',j), 'Reflectance', list(prj_info$dimx,prj_info$dimy), 32767, 'Composite band', prec='short', compression=2)
        try(phenoImg <- ncvar_add(phenoImg, new_var), silent=TRUE)  #Add new layer, MAY UPDATE EVENTUALLY TO ADD LINES THAT MERGE COMPOSITE WITH SYNTHETIC
        if (inherits(phenoImg, 'try-error')) {nc_close(phenoImg)}
        ncvar_put(phenoImg, new_var, full_pheno_band)
      }
      remove(full_pheno_mat)
    }
    #Close nc file
    nc_close(phenoImg)
  }
}


################################################################################
# MSLSP_Script additions
######################

#If we are making synthetic images, also make composites 
if (params$phenology_parameters$doComposites) {
  boundaries <- chunkBoundaries(3660*3660, 196)
  i<-1
  for (item in boundaries$chunkEnd) {
    numPix <- item - boundaries$chunkStart[i]+1
    runNonvegComposite(i, numPix, imgYrs, phenYrs, errorLog, params)
    i<-i+1
  }
  
  createComposite(tile, numChunks, numPix, baseImage, phenYrs, params)
}


