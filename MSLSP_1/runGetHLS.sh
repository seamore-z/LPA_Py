#!/bin/bash -l
#$ -j y


tile=$1
baseDir=$2
imgStartYr=$3
imgEndYr=$4

imgDir="/projectnb/modislc/users/seamorez/Classes/DIP/input/HLS30/05WPS/images/"


#Get Images from web
##########
getHLS.sh $tile $imgStartYr $imgEndYr $imgDir 

