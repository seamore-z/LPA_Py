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

# Inputs: list of years doing phenology for (script must loop and create composites for each yr)

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

readChunks <- function(numChunks, numPix, tempDir, phenometric) {
  boundaries <- chunkBoundaries(numPix, numChunks)
  
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

# THIS FUNCTION IS FOR INDIVIDUAL IMAGE CHUNK
DoNonvegComp <- function(b2, b3, b4, b5, b6, b7, dates, year){
  
  # Initialize dates from an image chunk
  yrdoy_char = as.character(format(dates, "%Y%j")) # Returning to yrdoy from dates
  yrs = substring(yrdoy_char, first = 1, last = 4) # Get years
  doys = substring(yrdoy_char, first = 5, last = 7) # Get doys
  
  phen_composites <- vector(mode = "list", length = 7)
  
  # Generate mean phenometrics for the full image
  phenoPath <- paste0('/projectnb/modislc/users/seamorez/Classes/DIP/output/05WPS/phenoMetrics/MSLSP_05WPS_',as.character(year),'.nc')
  phenoImg <- nc_open(phenoPath)
  #Get 7 DOY phenometrics
  OGI <- ncvar_get(phenoImg, 'OGI'); OGI50 <- ncvar_get(phenoImg, '50PCGI')
  OGMx <- ncvar_get(phenoImg, 'OGMx'); Peak <- ncvar_get(phenoImg, 'Peak')
  OGD <- ncvar_get(phenoImg, 'OGD'); OGD50 <- ncvar_get(phenoImg, '50PCGD')
  OGMn <- ncvar_get(phenoImg, 'OGMn')
  nc_close(phenoImg) # Close nc file
  #Mean of phenometrics
  OGImean <- round(mean(OGI,na.rm=T)); OGI50mean <- round(mean(OGI50,na.rm=T))
  OGMxmean <- round(mean(OGMx,na.rm=T)); Peakmean <- round(mean(Peak,na.rm=T))
  OGDmean <- round(mean(OGD,na.rm=T)); OGD50mean <- round(mean(OGD50,na.rm=T))
  OGMnmean <- round(mean(OGMn,na.rm=T))
  meanpheno <- c(OGImean,OGI50mean,OGMxmean,Peakmean,OGDmean,OGD50mean,OGMnmean)
  
  phen_yr_i = which(yrs %in% as.character(year))
  phen_doys = doys[phen_yr_i]
  b2_y <- b2[,phen_yr_i]; b3_y <- b3[,phen_yr_i]; b4_y <- b4[,phen_yr_i]
  b5_y <- b5[,phen_yr_i]; b6_y <- b6[,phen_yr_i]; b7_y <- b7[,phen_yr_i]
  
  i<-1
  for (phenometric in meanpheno) {
    comp_doys_i = which(as.integer(phen_doys) %in% c((phenometric-7):(phenometric+7)))
    b2_r <- b2_y[,comp_doys_i]; b3_r <- b3_y[,comp_doys_i]; b4_r <- b4_y[,comp_doys_i]
    b5_r <- b5_y[,comp_doys_i]; b6_r <- b6_y[,comp_doys_i]; b7_r <- b7_y[,comp_doys_i]
    b2_r[b2_r==32767]<-NA;b3_r[b3_r==32767]<-NA;b4_r[b4_r==32767]<-NA;b5_r[b5_r==32767]<-NA;b6_r[b6_r==32767]<-NA;b7_r[b7_r==32767]<-NA
    
    #Composites for each band
    b2_c <- round(rowMeans(b2_r,na.rm=TRUE)); b3_c <- round(rowMeans(b3_r,na.rm=TRUE)); b4_c <- round(rowMeans(b4_r,na.rm=TRUE))
    b5_c <- round(rowMeans(b5_r,na.rm=TRUE)); b6_c <- round(rowMeans(b6_r,na.rm=TRUE)); b7_c <- round(rowMeans(b7_r,na.rm=TRUE))
    composite <- cbind(b2_c, b3_c, b4_c, b5_c, b6_c, b7_c)
    phen_composites[[i]] <- composite
    i<-i+1
  }
  
  return(phen_composites)
}

#Output as another layer in the temporary folder (For this function, you will need the output nc
#file already created! You need to call this file to get average phenometrics)
runNonvegComposite <- function(chunk, numPix, imgYrs, phenYrs) {#, errorLog, params) {
  
  #pheno_pars <- params$phenology_parameters
  
  #Get all images to process
  ######################
  chunkFold <- paste0('/projectnb/modislc/users/seamorez/Classes/DIP/output/05WPS/imageChunks/','c',chunk,'/')
  outFold <- '/projectnb/modislc/users/seamorez/Classes/DIP/output/05WPS/test_outputs/'
  #chunkFold <- paste0(params$dirs$chunkDir,'c',chunk,'/') 
  #outFold <-   paste0(params$dirs$tempDir,'outputs/')  
  
  
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
    #if (inherits(imgData, 'try-error')) {cat(paste('runPhenoChunk: Error for chunk',chunk,img), file=errorLog, append=T);next} 
    
    b2[,i] <- imgData[,1]; b3[,i] <- imgData[,2]; b4[,i] <- imgData[,3]
    b5[,i] <- imgData[,4]; b6[,i] <- imgData[,5]; b7[,i] <- imgData[,6]
    
    remove(imgData)
  }
  
  dates = as.Date(strptime(yrdoy[ord], format="%Y%j"))  #Format as dates
  
  all_phenocomps <- vector(mode = "list", length = 7)
  for (year in phenYrs) {
    # comp_mat AND WRITING RESULTS PART NEED TO BE EDITED BASED ON UNIQUE PROPERTIES OF MULTI-BAND
    # IMAGE CHUNKS
    all_phenocomps <- DoNonvegComp(b2,  b3,  b4,  b5,  b6, b7, dates, year)
    
    j<-1
    for (phenocomp in all_phenocomps) {
      yearcheck <- paste0(outFold,'y',year)
      if(!dir.exists(yearcheck)){dir.create(yearcheck)}
      #foldercheck <- paste0(outFold,'y',year,'/c',chunk,'/')
      #if(!dir.exists(foldercheck)){dir.create(foldercheck)}
      saveRDS(phenocomp, paste0(yearcheck,'/c',chunk,'_p',j,'.Rds'))
      j<-j+1
    }
  }
  
}

createComposite <- function(numChunks, numPix, baseImage, phenYrs, outFold) {
  prj_info <- getNetCDF_projection_info(baseImage) 
  for (year in phenYrs) {
    phenoPath <- paste0('/projectnb/modislc/users/seamorez/Classes/DIP/output/05WPS/phenoMetrics/MSLSP_05WPS_',as.character(2021),'.nc')
    phenoImg <- nc_open(phenoPath, write=TRUE, verbose=FALSE)
    for (i in c(1:7)) {
      full_pheno_mat <- readChunks(numChunks, numPix, paste0(outFold,'y',as.character(year),'/'), i)
      for (j in c(1:6)) {
        print(dim(full_pheno_mat))
        full_pheno_band <- matrix(full_pheno_mat[,j], dim(baseImage)[1],dim(baseImage)[2])
        print(dim(full_pheno_band))
        print(full_pheno_band[1,1])
        full_pheno_band[full_pheno_band < -32767 | full_pheno_band > 32767 | is.na(full_pheno_band)] <- 32767
        new_var = ncvar_def(paste0('p',i,'_b',j), 'Reflectance', list(prj_info$dimx,prj_info$dimy), 32767, 'Composite band', prec='short', compression=2)
        try(phenoImg <- ncvar_add(phenoImg, new_var), silent=TRUE)
        if (inherits(phenoImg, 'try-error')) {nc_close(phenoImg)}
        ncvar_put(phenoImg, new_var, full_pheno_band)
      }
    }
    # Close nc file
    nc_close(phenoImg)
  }
}

################################################################################PRACTICE
# Initialize bands from an image chunk
#Get all images to process
######################
imgFold <- '/projectnb/modislc/users/seamorez/Classes/DIP/input/HLS30/05WPS/images/'
chunkFold <- '/projectnb/modislc/users/seamorez/Classes/DIP/output/05WPS/imageChunks/c100/'
outFold <- '/projectnb/modislc/users/seamorez/Classes/DIP/output/05WPS/test_outputs/'
#outFold <-   paste0(params$dirs$tempDir,'outputs/')  

chunkList <- list.files(path=chunkFold, pattern=glob2rx("HLS_*.Rds"), full.names=F)

numImgs = length(chunkList)
numChunks = 196

imgList <- list.dirs(path=imgFold, full.names=T)[-1]

yrdoy = as.numeric(matrix(NA,numImgs,1))
for (i in 1:numImgs) {
  imName <- gsub('.Rds','',imgList[i])
  yrdoy[i]= as.numeric(unlist(strsplit(imName,'_',fixed = T))[4])}

ord = order(yrdoy)  #Determine image order

imgName_strip  <- matrix(NA,length(imgList),1)
for (i in 1:length(imgList)) {
  imgName_strip[i] = tail(unlist(strsplit(imgList[i],'/')),n = 1)}

baseImage  <-  rast(paste0(imgList[1],'/',imgName_strip[1],'.Fmask.tif')) #Set up base image that we'll use for outputs
numPix  <-  round(ncol(baseImage)*nrow(baseImage)/numChunks)   #Get number of pixels per chunk

baseChunk  <-  try(matrix(readRDS(paste0(chunkFold,chunkList[1])),nrow=numPix),silent = TRUE) #Set up base image that we'll use for outputs

# START HERE! TRYING TO WORK OUT HOW TO TRACK THIS PROGRESS #
boundaries <- chunkBoundaries(3660*3660, 196)
i<-1
for (item in boundaries$chunkEnd) {
  numPix <- item - boundaries$chunkStart[i]+1
  runNonvegComposite(i, numPix, c(2021), c(2021)) #trying for chunk 100 for 2021
  i<-i+1
}
 #runNonvegComposite('100', numPix, c(2021), c(2021)) #trying for chunk 100 for 2021
createComposite(196, 13395600, baseImage, c(2021), outFold)

#Read in all imagery for chunk
b2 <- matrix(NA,numPix,numImgs);  b3 <- matrix(NA,numPix,numImgs);  b4 <- matrix(NA,numPix,numImgs)
b5 <- matrix(NA,numPix,numImgs);  b6 <- matrix(NA,numPix,numImgs);  b7 <- matrix(NA,numPix,numImgs)

for (i in 1:length(ord)) {
  img <- chunkList[ord[i]]
  imgData <- try(matrix(readRDS(paste0(chunkFold,img)),nrow=numPix),silent = TRUE)
  if (inherits(imgData, 'try-error')) {cat(paste('runPhenoChunk: Error for chunk',chunk,img), file=errorLog, append=T);next} 
  
  b2[,i] <- imgData[,1]; b3[,i] <- imgData[,2]; b4[,i] <- imgData[,3]
  b5[,i] <- imgData[,4]; b6[,i] <- imgData[,5]; b7[,i] <- imgData[,6]
  
  remove(imgData)
}

dates = as.Date(strptime(yrdoy[ord], format="%Y%j"))  #Format as dates

# Initialize dates from an image chunk
yrdoy_char = as.character(format(dates, "%Y%j")) # Returning to yrdoy from dates
yrs = substring(yrdoy_char, first = 1, last = 4) # Get years
doys = substring(yrdoy_char, first = 5, last = 7) # Get doys

phen_yr_i = which(yrs %in% '2021')
phen_imgList = imgList[ord[phen_yr_i]] # Should not have to use ord for bands as they should be ordered by the above code
phen_doys = doys[phen_yr_i]

comp_doys_i = which(as.integer(phen_doys) %in% c(143:157)) # doy 150 +/-7 days
comp_imgList = phen_imgList[ord[comp_doys_i]] # Should not have to use ord for bands as they should be ordered by the above code

# Test compositing
if (params$setup$runNonvegComposite) { 
  #Run compositing code for each image chunk
  imgLog <- foreach(j=1:numChunks) %dopar% {
    log <- try({runNonvegComposite(j, numPixPerChunk[j], imgYrs, phenYrs, errorLog, params)},silent=T)
    if (inherits(log, 'try-error')) {cat(paste('runNonvegComposite: Error for chunk', j,'\n'), file=errorLog, append=T)}
  }
  
  
  #Loop through the years processed and output netcdf files
  yrs <- phenStartYr:phenEndYr  
  log <- foreach(yr = yrs) %dopar% { 
    productFile  <- paste0(params$dirs$phenDir,'MSLSP_',tile,'_',yr,'.nc') 
    qaFile  <- paste0(params$dirs$phenDir,'MSLSP_',tile,'_',yr,'_Extended_QA.nc') 
    
    CreateExtendedQA(yr,qaFile, productTable, baseImage, params)
    CreateProduct(yr,productFile, qaFile, productTable, baseImage, waterMask, params)
  }
}