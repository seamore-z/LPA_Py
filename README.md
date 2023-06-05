# LPA_Py
Landsat Phenology Algorithm

The scripts provided here allow users to apply the Landsat Phenology Algorithm across Landsat ARD tiles. This is a Python version of LPA originally created by Eli Melaas. It has been updated for ARD Collection 2 including Alaska and Hawaii, and for greater parameterization ability to do analysis on different preprocessing corrections, VIs, averaging periods, and peak greenness trends. Here is a brief description of each script:

untar_ARD.py - Untars Landsat ARD files and generate directories for further processing<br>
calc_vi.ipynb - Generates VI GeoTIFF for each image in the ARD stack and applies topographic correction
