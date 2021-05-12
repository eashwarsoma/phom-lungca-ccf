library(stringr)
library(tidyr)
library(dplyr)
library(oro.dicom)
library(oro.nifti)

#The paths in this represent folders of patients/studies. 
#Within each folder are multiple CT scans, dose plans, and segmentations
patients.filenames <- list.files("/Volumes/Cancr/physics/Peng/Lung_SBRT/Radiomics/Data", 
                                 full.names = TRUE) 


#Excluding 2011_02__Studies to do no RT ST
#Excluding 2020-08__Studies due to no RTST
#Excluding 2020-09__Studies due to no CT
inds.to.drop <- which(patients.filenames %in% c(
                        "/Volumes/Cancr/physics/Peng/Lung_SBRT/Radiomics/Data/2011-02__Studies",
                        "/Volumes/Cancr/physics/Peng/Lung_SBRT/Radiomics/Data/2020-08__Studies",
                        "/Volumes/Cancr/physics/Peng/Lung_SBRT/Radiomics/Data/2020-09__Studies"))
patients.filenames <- patients.filenames[-inds.to.drop]





#Each list item are patient(s)
#There are variable number of patients that can be told apart 
#By the preceeding number
#Within each item are all the associated scans
scans.filenames <- lapply(patients.filenames, list.files, full.names = TRUE)

#Function to extract exact file names from dates folders
filepaths.from.date.folder <- function(date.folder.file) {

  #Gets the basic pattern
  RTst <- str_match(date.folder.file, "s\\/\\d.*_RTst") %>% na.omit()
  CTsc <- str_match(date.folder.file, "s\\/\\d.*_CT") %>% na.omit()
  
  #Omitting the CT and RT effectively creating an ID variable
  #Retaining the CT and RT for later use in getting the correct file path
  CTdf <- cbind.data.frame(ID = str_remove(CTsc, "_CT"), CT = CTsc)
  RTdf <- cbind.data.frame(ID = str_remove(RTst, "_RTst"), RT = RTst)
  
  #Gets all the relevant scans per date folder
  scansdf <- merge(CTdf, RTdf, by="ID")
  
  #Drop duplicate rows
  scansdf <- unique(scansdf[ , 1:3])
  
  #Obtaining the indices of the CT and RT scan file name from the original string
  inds.CT <- lapply(scansdf$CT, grep, date.folder.file)
  inds.RT <- lapply(scansdf$RT, grep, date.folder.file)
  
  #List to collect file names
  list.obs <- list()
  for (i in 1:dim(scansdf)[1]) {
    list.obs[[i]] <- c(date.folder.file[inds.CT[[i]]], date.folder.file[inds.RT[[i]]])
  }
  
  names(list.obs) <- scansdf[,1]
  
  return(list.obs)

}

#Getting the Exact CT and RT St Filenames
true.file.names <- lapply(scans.filenames, filepaths.from.date.folder)

#Outersect function to help select corrected file
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}

#For those studies with multiple RTst, trying to select the corrected file
for (i in 1:length(true.file.names)) {
  for (j in 1:length(true.file.names[[i]]))
    if (length(true.file.names[[i]][[j]]) > 2) {
    rt.ind <- grep("RTst", true.file.names[[i]][[j]])
    cor.ind <- grep("COR", true.file.names[[i]][[j]], ignore.case = TRUE)
    
    #Only drop stuff if the drop value has stuff
    drop <- outersect(cor.ind, rt.ind)
      if (length(drop) != 0) {
      true.file.names[[i]][[j]] <- true.file.names[[i]][[j]][-drop]
      }
    }

}

#Revising the list structure so that each patient is a list item
#Entry 1 is scan ID, E2 is CT, E3 is RT
list.filenames <- list()
counter <- 1 #Here to get good indexing
for (i in 1:length(true.file.names)) {
  for (j in 1:length(true.file.names[[i]])) {
      list.filenames[[counter]] <- c(names(true.file.names[[i]][j]), true.file.names[[i]][[j]])
      counter <- counter + 1
  }
}

#Now Manually Editing the list.file names to select the good RTst File
#For these cases, there are 2 Corrected RTst files, just choosing the more recent one
#IF FILE STRUCTURE EVER CHANGES< CHECK THIS
which(sapply(list.filenames, length) != 3)

list.filenames[[473]] <- list.filenames[[473]][c(1, 2, 4)]
list.filenames[[493]] <- list.filenames[[493]][c(1, 2, 4)]
list.filenames[[500]] <- list.filenames[[500]][c(1, 2, 4)]
list.filenames[[575]] <- list.filenames[[575]][c(1, 2, 4)]
list.filenames[[718]] <- list.filenames[[718]][c(1, 2, 4)]


#This gives us a dcm file of RTst
rtst.dcmlist <- lapply(list.filenames, function (x) list.files(x[3], full.names = TRUE))

#This only keeps the dcm files 
rtst.dcmlist.edit <- lapply(rtst.dcmlist, function (x) x[grep(".dcm", x)])

#Updating the elements within needed.scans.filenames
for (i in 1:length(list.filenames)) {
  list.filenames[[i]][3] <- rtst.dcmlist.edit[[i]][1]
}


saveRDS(list.filenames, "filenames.rds")


#Code to Removing redundant files
#CONTINUE

