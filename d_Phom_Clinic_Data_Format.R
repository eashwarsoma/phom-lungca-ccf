library(ripserr)
library(TDAstats)
library(reshape2)
library(ggplot2)
library(miscset)
library(tidyverse)
library(readxl)
library(survival)
library(survminer)
library(stringr)

#Reading in the cropped scans
#cropped.scans <- readRDS("list.cropped.scans.rds")

#Doing the Cubical Complex and getting raw Phom
#list.raw.phom <- lapply(cropped.scans, cubical)

#Saving Raw PHom as an RDS file in Case
#saveRDS(list.raw.phom, "list.raw.phom.rds")

#De-comment the above to run cubical complexes


#Reading in raw persistent homology 
list.raw.phom <- readRDS("list.raw.phom.rds")

#Creating a list file structure for each patient
list.all.data <- list()

#First adding IDs
for (i in seq_len(length(list.raw.phom))) {
  get_ID <- list.raw.phom[i] %>%
    names %>%
    strsplit(split = "_", n = 2) %>%
    unlist
  
  list.all.data[[i]] <- c(Patient.ID = get_ID[1],
                          Scan.ID = get_ID[2])
}

#Function to count total features that exist at each threshold value
cubical.total.feat.counter <- function(phom, minimum, maximum, res) {
  by <- (maximum - minimum) / res
  seq.to.check <- seq(minimum, maximum, by)
  
  mat <- seq.to.check %>%
    length %>%
    rep(x = 0) %>%
    cbind(seq.to.check)
  
  for (i in seq_len(length(seq.to.check))) {
    tally <- sum(seq.to.check[i] >= phom[,2] & 
                 seq.to.check[i] <= phom[,3])
    
    mat[i, 1] <- tally
  }
  
  return(mat[, c(2, 1)])
}

#Creating a function to count based on individual dim feat
featcounter.vector.cubical.act <- function(phom, minimum, maximum, res, feat) {
  
  phom <- phom[which(phom[, 1] == feat), ]
  
  by <- (maximum - minimum) / res
  seq.to.check <- seq(minimum, maximum, by)
  
  mat <- seq.to.check %>%
    length %>%
    rep(x = 0) %>%
    cbind(seq.to.check)
  
  for (i in seq_len(length(seq.to.check))) {
    tally <- sum(seq.to.check[i] >= phom[, 2] & 
                 seq.to.check[i] <= phom[, 3])
    
    mat[i, 1] <- tally
  }
  
  return(mat[,c(2,1)])
}

#Adding TDA Data
for (i in seq_len(length(list.raw.phom))) {
  #Removing the -1 dimension and the death at 9999
  phom.test.crr <- list.raw.phom[[i]][-which(list.raw.phom[[i]]["dimension"] == -1), ]
  
  #Creating the Topological Feature Curves
  test.0 <- featcounter.vector.cubical.act(phom.test.crr, min = -1024, max = 3071,
                                                          res = 1000, feat = 0)
  test.1 <- featcounter.vector.cubical.act(phom.test.crr, min = -1024, max = 3071,
                                                          res = 1000, feat = 1)
  test.2 <- featcounter.vector.cubical.act(phom.test.crr, min = -1024, max = 3071,
                                                          res = 1000, feat = 2)
  test.tot <- cubical.total.feat.counter(phom.test.crr, min = -1024, max = 3071,
                                                        res = 1000)
  
  all.feat.curves <- cbind(test.0[, 1],
                           test.0[, 2],
                           test.1[, 2],
                           test.2[, 2],
                           test.tot[, 2])
  
  all.feat.curves <- as.data.frame(all.feat.curves)
  
  colnames(all.feat.curves) <- c("HU", "dim0", "dim1", "dim2", "dimtot")
  
  list.all.data[[i]] <- c(list.all.data[[i]], feat_curves = list(all.feat.curves))
  list.all.data[[i]] <- c(list.all.data[[i]], raw_phom = list(list.raw.phom[[i]])) 
}

#Visualizing Random Curves
topfeatcurv.melt <- melt(list.all.data[[167]][["feat_curves"]],
                         id.vars = c("HU"), 
                         measure.vars = c("dim0", "dim1", "dim2", "dimtot"))
colnames(topfeatcurv.melt) <- c("filtration", "feature.type", "feature.count")

ggplot(topfeatcurv.melt, aes(x = filtration, y = feature.count, color = feature.type)) +
  geom_path(size = 1) +
  theme_bw() +
  labs(x = "Filtration (HU)", 
       y = "Number of Features",
       color = "Feature Dimension") 

#Clinic Data Retrieve
key <- read_excel("/Volumes/Cancr/physics/Peng/Lung_SBRT/Radiomics/Key.xlsx")
clinical.data <- read_excel("clinic_data.xlsx")
date_IDs <- readRDS("Scan_Date_ID_Key.rds")

#Incorporate Clinical Variables
for (i in seq_len(length(list.all.data))) {

  #Getting the CCF MRN Number
  pid <- list.all.data[[i]]$Patient.ID
  
  #Getting the Scan Date
  #All analysis time intervals grounded against scan Date
  scan.date <- as.Date(date_IDs[which(date_IDs$PatID == pid), "date"])
  
  #Getting the MRN
  MRN <- as.numeric(key[which(key$anon_id == pid), ]["ccf_number"])
  
  #Subsetting clinic data to the appropriate row
  clinic.MRN <- subset(clinical.data, `CCF #` == MRN)
  
  #Only getting data if an MRN is present
  if(nrow(clinic.MRN) == 1) {
  
    #Getting the clinical Data into a data frame
    clin.dat <- cbind.data.frame(
    #Basic initial ID Information (MRN, Scan Date, Sex, Ethnicity, DOB, Age, Diagnosis Date, Rad Start Date)
      MRN = MRN, 
      Scan.Date = scan.date, 
      Sex = clinic.MRN$Gender, 
      Ethnicity = clinic.MRN$Ethnicity, 
      DOB = clinic.MRN$`Date of Birth`, 
      Age.Diag = clinic.MRN$`Age at Diagnosis`, 
      DODiag = clinic.MRN$`Date of Diagnosis...170`, 
      Date.Rad.Start = clinic.MRN$`RT Start Date...199`,
      
      #Basic Patient Health Information (BMI, Years SMoking, KPS, HGB, CCI comorbidity, 
      #Previous Cancer Diagnosis, #Number of Previous Cancer Diagnoses
      BMI = clinic.MRN$`BMI`, 
      Years.Smoke = clinic.MRN$`Pack Years Smoking`,
      KPS = clinic.MRN$`KPS`, 
      HGB = clinic.MRN$`HGB Pre`, 
      CCI = clinic.MRN$`Total points:`,
      Prev.Canc.Diag = clinic.MRN$`Previous Cancer Diagnosis?`,
      Num.Prev.Canc.Diag = clinic.MRN$`Number of Previous Cancer Diagnoses`, 
      
      #Treatment Info; prior surgery, pre or post radiation chemo, intent of radiation, dose of radiation
      Prior.Surg = clinic.MRN$`Prior Surgery?`, 
      Post.Rad.Chemo = clinic.MRN$`Post-SBRT Chemo/Systemic RX?...657`,
      Pre.Rad.Chemo = clinic.MRN$`Pre-SBRT Chemo?...191`, 
      Rad.Intent = clinic.MRN$`Intent...194`, 
      Tot.Dose = clinic.MRN$`Total Dose (Gy)...201`, 
      
      #Cancer Characteristics including tumor size, TNM, overall stage, Histology, Pathology
      CT.Size = clinic.MRN$`CT (Cm)...177`, 
      T.Stage = clinic.MRN$`T Stage...179`, 
      N.Stage = clinic.MRN$`N Stage...180`,
      M.Stage = clinic.MRN$`M Stage...181`, 
      Overall.Stage = clinic.MRN$Stage...182, 
      Histo = clinic.MRN$Histology...189,
      Path = clinic.MRN$Pathology,
      
      #Survival Information
      Vital.Status = clinic.MRN$`Survival Status`, 
      OS.Length = as.Date(clinic.MRN$`Date of Death`) - scan.date,
      Date.Death = clinic.MRN$`Date of Death`,
      
      #Progression Information
      Progression.Status = ifelse(is.na(clinic.MRN$`1st Failure Date...678`), "Stable", "Progressed"), 
      PFS.Length = as.Date(clinic.MRN$`1st Failure Date...678`) - scan.date,
      PFS.Date = clinic.MRN$`1st Failure Date...678`,
      PFS.Failure.Type = clinic.MRN$`Failure Type...677`,
      
      #Local Control Information
      Local.Status = ifelse(clinic.MRN$Local...666 == "Yes", "Failed", "Controlled"), 
      Local.Length = as.Date(clinic.MRN$`If Local Failure, Date...671`) - scan.date, 
      Date.Local.Failure = clinic.MRN$`If Local Failure, Date...671`,
      
      #Other Progression (Lobar, Nodal, Distal) Information, first Met Site
      Lobar.Status = ifelse(clinic.MRN$Lobar...667 == "No", "Lobar.Controlled", "Lobar.Failed"),
      Date.Lobar = clinic.MRN$`If Lobar Failure, Date`,
      
      Nodal.Status = ifelse(clinic.MRN$Nodal...668 == "No", "Nodal.Controlled", "Nodal.Failed"),
      Date.Nodal = clinic.MRN$`If Nodal Failure, Date...673`,
      
      Distant.Status = ifelse(clinic.MRN$Distant...669 == "No", "Distant.Controlled", "Distant.Failed"),
      Date.Distant = clinic.MRN$`If Distant Failure, Date...674`,
      
      #First Metastatic Site
      First.Met = clinic.MRN$`1st Distant Metastasis Failure Site...679`,
      
      #Follow Up Status
      Date.Last.Update = clinic.MRN$`Last Updated On`,
      Date.Last.Followup = clinic.MRN$`Last Follow-up On`,
      Date.Last.Known.Alive = clinic.MRN$`Last Known Alive`,
      
      Update.Time = as.Date(clinic.MRN$`Last Updated On`) - scan.date,
      Followup.Time = as.Date(clinic.MRN$`Last Follow-up On`) - scan.date,
      Alive.Time = as.Date(clinic.MRN$`Last Known Alive`) - scan.date
    )
  
    list.all.data[[i]] <- c(list.all.data[[i]], clinic.data = list(clin.dat))
  }
}

#Seeing How many entries dont have clinical data
#Follow up with Peng on This
no.MRN.found <- which(sapply(list.all.data, length) == 4)

#Removing the entries with No Clinical Data
list.all.data <- list.all.data[-no.MRN.found]

#Incorporate Moments of the Distribution
#Need Moments 1 to 4 for each Feat Curve (16 possible variables)
#To be super thorough, calculates all 4 flavors of moments
#As a reminder: Mean is the 1st raw moment; variance is the 2nd central moment
#Skewness and Kurtosis are the standardized centralized 3rd and 4th moments respectively

#Calculates raw moment
raw.moment <- function(vec, degree) {
  len <- length(vec)
  val <- sum(vec ^ degree)
  
  return(val / len)
}

#Calculates centralized moment (by subtracting mean)
cent.moment <- function(vec, degree) {
  avg <- mean(vec)
  len <- length(vec)
  val <- sum((vec - avg) ^ degree)
  
  return(val / len)
}

#Calculates standardized moment(by dividing sd^n)
stand.moment <- function(vec, degree) {
  sd <- sd(vec)
  len <- length(vec)
  val <- sum(vec ^ degree) / (sd ^ degree)
  
  return(val / len)
}

#Calculates standardized and centralized moment
standcent.moment <- function (vec, degree) {
  avg <- mean(vec)
  sd <- sd(vec)
  len <- length(vec)
  val <- sum(((vec - avg) / sd) ^ degree)
  return(val / len)
}

for (i in seq_len(length(list.all.data))) {
  dim0 <- list.all.data[[i]]$feat_curves$dim0
  dim1 <- list.all.data[[i]]$feat_curves$dim1
  dim2 <- list.all.data[[i]]$feat_curves$dim2
  dimtot <-list.all.data[[i]]$feat_curves$dimtot
  
  #Raw moments For Dim 0 Curve (Mean is first raw moment)
  dim0.raw.moms <- c(raw.1st.dim0.mean = raw.moment(dim0, 1), 
                     raw.2nd.dim0 = raw.moment(dim0, 2), 
                     raw.3rd.dim0 = raw.moment(dim0, 3), 
                     raw.4th.dim0 = raw.moment(dim0, 4))
  
  #Centralized moments For Dim 0 Curve (Variance is second Centralized moment)
  dim0.cent.moms <- c(cent.1st.dim0 = cent.moment(dim0, 1), 
                      cent.2nd.dim0.variance = cent.moment(dim0, 2), 
                      cent.3rd.dim0 = cent.moment(dim0, 3), 
                      cent.4th.dim0 = cent.moment(dim0, 4))
  
  #standardized moments For Dim 0 Curve
  dim0.stand.moms <- c(stand.1st.dim0 = stand.moment(dim0, 1), 
                       stand.2nd.dim0 = stand.moment(dim0, 2), 
                       stand.3rd.dim0 = stand.moment(dim0, 3), 
                       stand.4th.dim0 = stand.moment(dim0, 4))
  
  #standardized and centralized moments For Dim 0 Curve
  dim0.standcent.moms <- c(standcent.1st.dim0 = standcent.moment(dim0, 1), 
                           standcent.2nd.dim0 = standcent.moment(dim0, 2), 
                           standcent.3rd.dim0.skew = standcent.moment(dim0, 3), 
                           standcent.4th.dim0.kurt = standcent.moment(dim0, 4))
  
  #Collecting all the 0dim moments into one list
  list.dim0.mom <- list(dim0.raw.moms = dim0.raw.moms, dim0.cent.moms = dim0.cent.moms, 
                        dim0.stand.moms = dim0.stand.moms, dim0.standcent.moms = dim0.standcent.moms)
  
  #Doing the Same for Dim 1
  dim1.raw.moms <- c(raw.1st.dim1.mean = raw.moment(dim1, 1), 
                     raw.2nd.dim1 = raw.moment(dim1, 2), 
                     raw.3rd.dim1 = raw.moment(dim1, 3), 
                     raw.4th.dim1 = raw.moment(dim1, 4))
  
  dim1.cent.moms <- c(cent.1st.dim1 = cent.moment(dim1, 1), 
                      cent.2nd.dim1.variance = cent.moment(dim1, 2), 
                      cent.3rd.dim1 = cent.moment(dim1, 3), 
                      cent.4th.dim1 = cent.moment(dim1, 4))
  
  dim1.stand.moms <- c(stand.1st.dim1 = stand.moment(dim1, 1), 
                       stand.2nd.dim1 = stand.moment(dim1, 2), 
                       stand.3rd.dim1 = stand.moment(dim1, 3), 
                       stand.4th.dim1 = stand.moment(dim1, 4))
  
  dim1.standcent.moms <- c(standcent.1st.dim1 = standcent.moment(dim1, 1), 
                           standcent.2nd.dim1 = standcent.moment(dim1, 2), 
                           standcent.3rd.dim1.skew = standcent.moment(dim1, 3), 
                           standcent.4th.dim1.kurt = standcent.moment(dim1, 4))
  
  list.dim1.mom <- list(dim1.raw.moms = dim1.raw.moms, dim1.cent.moms = dim1.cent.moms, 
                        dim1.stand.moms = dim1.stand.moms, dim1.standcent.moms = dim1.standcent.moms)
  
  #Doing the Same for Dim 2
  dim2.raw.moms <- c(raw.1st.dim2.mean = raw.moment(dim2, 1), 
                     raw.2nd.dim2 = raw.moment(dim2, 2), 
                     raw.3rd.dim2 = raw.moment(dim2, 3), 
                     raw.4th.dim2 = raw.moment(dim2, 4))
  
  dim2.cent.moms <- c(cent.1st.dim2 = cent.moment(dim2, 1), 
                      cent.2nd.dim2.variance = cent.moment(dim2, 2), 
                      cent.3rd.dim2 = cent.moment(dim2, 3), 
                      cent.4th.dim2 = cent.moment(dim2, 4))
  
  dim2.stand.moms <- c(stand.1st.dim2 = stand.moment(dim2, 1), 
                       stand.2nd.dim2 = stand.moment(dim2, 2), 
                       stand.3rd.dim2 = stand.moment(dim2, 3), 
                       stand.4th.dim2 = stand.moment(dim2, 4))
  
  dim2.standcent.moms <- c(standcent.1st.dim2 = standcent.moment(dim2, 1), 
                           standcent.2nd.dim2 = standcent.moment(dim2, 2), 
                           standcent.3rd.dim2.skew = standcent.moment(dim2, 3), 
                           standcent.4th.dim2.kurt = standcent.moment(dim2, 4))
  
  list.dim2.mom <- list(dim2.raw.moms = dim2.raw.moms, dim2.cent.moms = dim2.cent.moms, 
                        dim2.stand.moms = dim2.stand.moms, dim2.standcent.moms = dim2.standcent.moms)
  
  #Doing the Same for Dim tot
  dimtot.raw.moms <- c(raw.1st.dimtot.mean = raw.moment(dimtot, 1), 
                       raw.2nd.dimtot = raw.moment(dimtot, 2), 
                       raw.3rd.dimtot = raw.moment(dimtot, 3), 
                       raw.4th.dimtot = raw.moment(dimtot, 4))
  
  dimtot.cent.moms <- c(cent.1st.dimtot = cent.moment(dimtot, 1), 
                        cent.2nd.dimtot.variance = cent.moment(dimtot, 2), 
                        cent.3rd.dimtot = cent.moment(dimtot, 3), 
                        cent.4th.dimtot = cent.moment(dimtot, 4))
  
  dimtot.stand.moms <- c(stand.1st.dimtot = stand.moment(dimtot, 1), 
                         stand.2nd.dimtot = stand.moment(dimtot, 2), 
                         stand.3rd.dimtot = stand.moment(dimtot, 3), 
                         stand.4th.dimtot = stand.moment(dimtot, 4))
  
  dimtot.standcent.moms <- c(standcent.1st.dimtot = standcent.moment(dimtot, 1), 
                             standcent.2nd.dimtot = standcent.moment(dimtot, 2), 
                             standcent.3rd.dimtot.skew = standcent.moment(dimtot, 3), 
                             standcent.4th.dimtot.kurt = standcent.moment(dimtot, 4))
  
  list.dimtot.mom <- list(dimtot.raw.moms = dimtot.raw.moms,
                          dimtot.cent.moms = dimtot.cent.moms, 
                          dimtot.stand.moms = dimtot.stand.moms,
                          dimtot.standcent.moms = dimtot.standcent.moms)
  
  #Creating a super list of all moments
  #This list -> 4 feature curves -> 4 moment flavorts -> 4 degrees (64 possible variables)
  #The special moments (mean, variance, kurtosis, or skewed) are labeled
  feat.curve.moms = list(dim0moms = list.dim0.mom,
                         dim1moms = list.dim1.mom,
                         dim2moms = list.dim2.mom,
                         dimtotmoms = list.dimtot.mom)
  
  #Putting into main list
  list.all.data[[i]] <- c(list.all.data[[i]], feat.curve.moms = list(feat.curve.moms))
}

#Removing DUplicates and Keeping  Most uptodate record
dups <- list.all.data %>%
  sapply(FUN = function(x) x[["Patient.ID"]]) %>%
  duplicated %>%
  which
list.all.data <- list.all.data[-dups]

which(sapply(list.all.data, function (x) dim(x[["clinic.data"]])[1]) == 2)
#42 221 435 468 477 528 650 725 726 763 have multiple clinical records, keeping earliest scan
list.all.data[[42]]$clinic.data <- list.all.data[[42]]$clinic.data[1, ]
list.all.data[[221]]$clinic.data <- list.all.data[[221]]$clinic.data[1, ]
list.all.data[[435]]$clinic.data <- list.all.data[[435]]$clinic.data[2, ]
list.all.data[[468]]$clinic.data <- list.all.data[[468]]$clinic.data[2, ]
list.all.data[[477]]$clinic.data <- list.all.data[[477]]$clinic.data[1, ]
list.all.data[[528]]$clinic.data <- list.all.data[[528]]$clinic.data[1, ]
list.all.data[[650]]$clinic.data <- list.all.data[[650]]$clinic.data[2, ]
list.all.data[[725]]$clinic.data <- list.all.data[[725]]$clinic.data[2, ]
list.all.data[[726]]$clinic.data <- list.all.data[[726]]$clinic.data[2, ]
list.all.data[[763]]$clinic.data <- list.all.data[[763]]$clinic.data[1, ]


#Saving the RDS File
#This file has the patient ID, scan ID, raw PHom, Feature Curve, Moments, and Clinical Data Based on Patient
#Each list item represents a single patient
#This object requires extracter functions to use properly
saveRDS(list.all.data, "formatted_clinical_phom_data.rds")
