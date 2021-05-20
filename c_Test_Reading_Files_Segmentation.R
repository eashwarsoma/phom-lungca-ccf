library(dplyr)
library(oro.dicom)
library(oro.nifti)
library(RNifti)
library(RadOnc)
library(tkrplot)
library(stringr)
#Reading in the large NIFTII folders
#Each element represents a patient
niftii.filenames <- list.files("/Users/somasue/Desktop/NIFTI", 
                                 full.names = TRUE) 

#Each list item is a patient
#Within each item are all the associated scans
nscans.filenames <- lapply(niftii.filenames, list.files, full.names = TRUE)

nscans.filenames[[755]]

#Function to get the image and GTV indices
#Selects all GTV files
target.index <- lapply(nscans.filenames,
                       function(x) {
                         c(grep("image.nii", x), grep("GTV", x))
                       }
)

#Which files have no GTV. These may have a CTV or PTV or ITV, but exlcuding to be safe
inds.to.cut <- lapply(target.index, function(x) return(length(x) < 2)) %>%
  unlist %>%
  which
target.index <- target.index[-inds.to.cut]
nscans.filenames <- nscans.filenames[-inds.to.cut]

#Process to select the largest GTV file indicesif there are multiple
for (i in seq_len(length(target.index))) {
  #Only bother with the ones > 2 (2 representing the image index + GTV index)
  if (length(nscans.filenames[[i]][target.index[[i]]]) > 2) {
    #Gets all the GTV file indices
    GTV.inds <- grep("GTV", nscans.filenames[[i]])
    
    #Gets the path names of the GTV segment files to Read in
    paths.to.check <- nscans.filenames[[i]][GTV.inds]
    
    #Reads in Nifti files
    Segs <- lapply(paths.to.check, readNIfTI)
    
    #Gets the matrix of segment coordinates
    inds.with.vertex <- lapply(Segs, function (x)  which(x > 0, arr.ind = TRUE))
    
    #Selects the index of the largest one
    ind.to.take <- which(sapply(inds.with.vertex, dim)[1, ] == max(sapply(inds.with.vertex, dim)[1, ]))
    
    #If there is a tie, keep the first
    ind.to.take <- ind.to.take[1]
    
    #Correlates the index of pathname that produced largest segment with the index from the filenames
    GTV.ind.to.keep <- which(nscans.filenames[[i]] == paths.to.check[ind.to.take])
    
    #Takes that index along with the main image index and replaces the target indices 
    target.index[[i]] <- c(grep("image.nii", nscans.filenames[[i]]), GTV.ind.to.keep)
  }
  
  print(paste("Done with File Index", i))
}


#Creating a new list of files with just the CT scan and Segment scan
true.list.scans <- list()
for (i in seq_len(length(target.index))) {
  true.list.scans[[i]] <- nscans.filenames[[i]][target.index[[i]]]
}


#Function that crops scans using the true list scans with CT first than segment
#Again, Eashwar, need to have list where each element is character[2]
#element 1 is the path name to the image.nii or CT scan
#element 2 is the path name to the GTV.nii or segment scan
scan.cropper <- function(paths) {
  CT <- readNIfTI(paths[1])
  Seg <- readNIfTI(paths[2])
  
  #Selecting the coordinates marked as a segmentation
  #Getting the boundary coordinates of the segment
  inds.with.vertex <- which(Seg > 0, arr.ind = TRUE)
  min.x <- min(inds.with.vertex[, 1])
  max.x <- max(inds.with.vertex[, 1])
  min.y <- min(inds.with.vertex[, 2])
  max.y <- max(inds.with.vertex[, 2])
  min.z <- min(inds.with.vertex[, 3])
  max.z <- max(inds.with.vertex[, 3])
  
  cropped.CT <- CT[min.x:max.x, min.y:max.y, min.z:max.z]
  return(cropped.CT)
}

#Cropping all the scans in the list
list.cropped.scans <- lapply(true.list.scans, scan.cropper)

#naming the cropped scans with patient ID
names(list.cropped.scans) <- str_remove(list.files("/Users/somasue/Desktop/NIFTI",
                                                   full.names = FALSE)[-inds.to.cut],
                                        ".nii")
  
saveRDS(list.cropped.scans, "list.cropped.scans.rds")

#Testing Random Scans
list.cropped.scans <- readRDS("list.cropped.scans.rds")

slices3d(list.cropped.scans[[round(runif(1, 1, length(list.cropped.scans)))]])
