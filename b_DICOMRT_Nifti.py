from dcmrtstruct2nii import dcmrtstruct2nii, list_rt_structs
from pydicom import dcmread
import rpy2.robjects as robjects
import numpy as np
from rpy2.robjects import pandas2ri
pandas2ri.activate()

readRDS = robjects.r['readRDS']

#Reading in the string names from R
filenames = readRDS("filenames.rds")

#turning it into a nested array object
filenames_arr = np.array(filenames)

#Getting the Patient ID, which will be the nifti filename
#Also Including SCANID
#Collecting ID array
ID = []
#For the number of patient files
for x in range(len(filenames_arr)):
  #Get the image of the RTSt object
  image = dcmread(filenames_arr[x][2])
  #Obtain Patient ID
  label = image['PatientID']
  #Obtain SCAN ID
  scanID = filenames_arr[x][0]
  #Combine into One String, gives a way of identifying patient by both
  col_lab =  label.value + scanID
  #Add combined label to the ID list
  ID.append

#Replacing the s/ with an _
#For the sake of cleaness
for x in range(len(ID)):
  ID[x] = ID[x].replace("s/", "_")

#Creating a Destination Filename String
destfiles = []  
for x in range(len(ID)):
  path = "/Users/somasue/Desktop/NIFTI/" + str(ID[x]) + ".nii"
  destfiles.append(path)

 
#Take all the file names and create Nifti
#As a reminder, filenames_arr[x][2] is the RTst file
#filenames_arr[x][1] is the CT DICOM Directory
#destfiles[x] is the destination directory
for x in range(len(filenames_arr)):
  dcmrtstruct2nii(filenames_arr[x][2], filenames_arr[x][1], destfiles[x])


