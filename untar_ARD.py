## -----------------------------------------------------------------------------------------------------------------##
## Written by Eli Melaas, adapted for Python by Seamore Zhu 
## -----------------------------------------------------------------------------------------------------------------##

from joblib import Parallel, delayed
import sys
import os
import glob
import tarfile

# Register parallel processing with joblib
Parallel(n_jobs=16)

#args = sys.argv[1:]  # Get command-line arguments
#tile = args[0]       # Extract the first argument to tile name
#data_dir = args[1]   # Extract the second argument to data directory

tile = 'h15v02'  # Will be parameterized in json file and passed as argument using shell script (see above code)
data_dir = '/projectnb/modislc/users/seamorez/Landsat_Pheno/ARD/'  # Will be parameterized in json file

#os.chdir(data_dir+tile)
os.chdir(data_dir+tile) 

os.makedirs(os.getcwd()+'/DEM/', exist_ok=True)
os.makedirs(os.getcwd()+'/IMG/', exist_ok=True)
os.makedirs(os.getcwd()+'/MAPS/', exist_ok=True)
os.makedirs(os.getcwd()+'/NLCD/', exist_ok=True)
os.makedirs(os.getcwd()+'/PHENO/', exist_ok=True)
os.makedirs(os.getcwd()+'/SHP/', exist_ok=True)

# Untar all inputs and sort into folders
in_dirs_SR = glob.glob("*SR*", recursive=True)
in_dirs_TA = glob.glob("*TOA*", recursive=True)

dir_names = [os.path.basename(dir_SR)[0:32] for dir_SR in in_dirs_SR]
for dir_name in dir_names:
    os.makedirs(os.path.join(os.getcwd(), 'IMG', dir_name), exist_ok=True)

for i in range(len(in_dirs_TA)):
    print(i)

    dir_name_TA = os.path.basename(in_dirs_TA[i])[0:32]
    tar = tarfile.open(in_dirs_TA[i])
    tar.extractall(os.path.join(os.getcwd(), 'IMG', dir_name_TA))
    tar.close()
    os.remove(in_dirs_TA[i])
    in_dirs_TAB = glob.glob(os.path.join(os.getcwd(), 'IMG', dir_name_TA, '*TOA_B*'), recursive=True)
    for file in in_dirs_TAB:
        os.remove(file)

for i in range(len(in_dirs_SR)):
    print(i)

    dir_name_SR = os.path.basename(in_dirs_SR[i])[0:32]
    tar = tarfile.open(in_dirs_SR[i])
    tar.extractall(os.path.join(os.getcwd(), 'IMG', dir_name_SR))
    tar.close()
    os.remove(in_dirs_SR[i])
