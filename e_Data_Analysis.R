library(tidyverse)
library(plyr)
library(survival)
library(survminer)
library(reshape)
library(moments)
library(gridExtra)
library(rstatix)
library(gt)
library(paletteer)
library(ggpubr)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
library(tableone)
library(TDAstats)
library(grid)
library(gridExtra)
library(ggplotify)
library(png)
library(rstatix)
library(ggfortify)
library(stringr)
library(hablar)
library(extrafont)
library(rms)
library(cmprsk)


#Function to linearly scale something to 0 to 1
normalize <- function(col) {
  return ((col - min(col)) / (max(col) - min(col)))
}

list.all.data <- readRDS("formatted_clinical_phom_data.rds")

#Function to extract moments from single list element
#Possible variables here; need to use lapply to get full etraction
mom.dim.vars <- c("dim0", "dim1", "dim2", "dimtot")
mom.types.vars <- c("raw", "cent", "stand", "standcent")

moment.extracter <- function(list.data, mom.dims, mom.types) {
  #Gets the patient ID
  PID <- list.data$Patient.ID
  
  #The dimension feature curve moments to extract
  dim.to.choose <- mom.dims
  
  #The moment types to extract
  moms.to.choose <- mom.types
  
  #Adding string characters to it to make str match unambigous 
  moms.to.choose.mod <- paste("\\.", moms.to.choose, "\\.", sep = "")
  
  #Logical Indices of dimension curves
  dim.inds <- str_detect(names(list.data$feat.curve.moms), paste(dim.to.choose, collapse = "|"))
  
  #Logical Indices of moments to choose
  moms.dims.ind  <- str_detect(names(unlist(list.data$feat.curve.moms[dim.inds], recursive = FALSE)), paste(moms.to.choose.mod, collapse = "|"))
  
  #Lists of the moments we want
  cap.list <- unlist(list.data$feat.curve.moms[dim.inds], recursive = FALSE)[moms.dims.ind]
  
  #Extracting into a vector
  vec.col <- numeric()
  for (i in 1:length(cap.list)) {
    vec.col <- c(vec.col, cap.list[[i]])
  }
  
  fin.vec <- c(PID = PID, vec.col)
  
  #Ensuring Numerica
  fin.vec[-1] <- as.numeric(fin.vec[-1])
  
  return(fin.vec)
}

#Function to extract clinical variables from list element
#Possible variables here; need to use lapply to get full etraction
patient.vars <- c("MRN", "Sex", "Ethnicity", "Age.Diag", "BMI",
                  "Years.Smoke", "KPS", "HGB", "CCI",
                  "Prev.Canc.Diag", "Num.Prev.Canc.Diag", "Scan.Date", 
                  "DODiag", "Date.Rad.Start", "Death.Cause")
treatment.vars <- c("Prior.Surg", "Post.Rad.Chemo", "Pre.Rad.Chemo",
                    "Rad.Intent", "Tot.Dose")
tumor.vars <- c("CT.Size", "T.Stage", "N.Stage", "M.Stage",
                "Overall.Stage", "Histo", "Path", "ALK", "BRAF", "EGFR", "KRAS")
event.vars <- c("Vital.Status", "OS.Length",
                "Progression.Status", "PFS.Length", "PFS.Failure.Type",
                "Local.Status", "Local.Length", 
                "Lobar.Status", "Nodal.Status", "Distant.Status", 
                "First.Met", "Update.Time", "Followup.Time", "Alive.Time",
                "Date.Death", "PFS.Date", "Date.Local.Failure", "Date.Lobar", "Date.Nodal", 
                "Date.Last.Update", "Date.Last.Followup", "Date.Last.Known.Alive")

clinic.extracter <- function(list.data, patient.vars, treatment.vars, tumor.vars, event.vars) {
  
  PID <- list.data$Patient.ID
  
  #Placeholder Dataframe
  plac <- list.data$clinic.data[c(event.vars, patient.vars, treatment.vars, tumor.vars)]
  
  #1 0 convention for events
  # Distant.Failed = Distant Control
  # Lobar.Failed = Lobar Control
  # Nodal.Failed = Nodal Control
  # Failed = Local Control
  # Dead = Vital Status
  # Progressed = Progression Free Survival
  fails.ind <- plac %in% c("Distant.Failed", "Lobar.Failed", "Nodal.Failed", "Failed", "Dead", "Progressed")
  controlled.ind <- plac %in% c("Distant.Controlled", "Lobar.Controlled", "Nodal.Controlled", "Controlled", "Alive", "Stable")
  
  plac[fails.ind] <- 1
  
  plac[controlled.ind] <- 0
  
  
  fin.vec <- cbind.data.frame(PID = PID, plac)
  
  return(fin.vec)
  
}


####Creating Clinical and Moment Data Frame####
moms.list <- lapply(list.all.data, moment.extracter, mom.dims = mom.dim.vars, mom.types = mom.types.vars)
moms.df <- as.data.frame(do.call(rbind, moms.list))

#Converting to numeric
moms.df[, -1] <- moms.df %>% select(-1) %>% sapply(as.character)
moms.df[, -1] <- moms.df %>% select(-1) %>% sapply(as.numeric) 

#Logging, Note errors produced from logging zeroes
moms.df[, -1] <- moms.df %>% select(-1) %>% sapply(log)


clinic.list <- lapply(list.all.data, clinic.extracter, patient.vars = patient.vars, treatment.vars = treatment.vars,
                      tumor.vars = tumor.vars, event.vars = event.vars)
clinic.df <- do.call(rbind.data.frame, clinic.list)


#Combining moments and clinical vairable data frames
clinic.mom.df <- join(clinic.df, moms.df)

sapply(clinic.mom.df, function(x) paste(sum(!is.na(x)), "/", length(x)))


#Creating a corrected Interval column that uses survival or follow up time
#Making other correct.intervals as well (Continue)
clinic.mom.df$Survival.Correct.Interval <- as.numeric(ifelse(clinic.mom.df$Vital.Status == 1, 
                                                             clinic.mom.df$Date.Death - clinic.mom.df$DODiag, 
                                                             clinic.mom.df$Date.Last.Followup - clinic.mom.df$DODiag))

clinic.mom.df$PFS.Correct.Interval <- as.numeric(ifelse(clinic.mom.df$Progression.Status == 1, 
                                                        clinic.mom.df$PFS.Date - clinic.mom.df$Date.Rad.Start, 
                                                        difftime(clinic.mom.df$Date.Last.Update, clinic.mom.df$Date.Rad.Start, units = "days")))

clinic.mom.df$Local.Correct.Interval <- as.numeric(ifelse(clinic.mom.df$Local.Status == 1, 
                                                          clinic.mom.df$Date.Local.Failure - clinic.mom.df$Date.Rad.Start, 
                                                          difftime(clinic.mom.df$Date.Last.Update, clinic.mom.df$Date.Rad.Start, units = "days")))

#Standardizing Variables, Giving these variable slots NA to reflect actual lack of value
clinic.mom.df$T.Stage[clinic.mom.df$T.Stage == "N/A"] <- NA
clinic.mom.df$T.Stage[clinic.mom.df$T.Stage == "T1"] <- NA
clinic.mom.df$T.Stage[clinic.mom.df$T.Stage == "T2"] <- NA
clinic.mom.df$T.Stage[clinic.mom.df$T.Stage == "X"] <- NA

#Giving KPS and CCI standard scores
clinic.mom.df$KPS_adj <- ifelse(clinic.mom.df$KPS <= 70, "KPS \u2264 70", "KPS > 70")

clinic.mom.df$CCI_adj <- ifelse(clinic.mom.df$CCI <= 4, "0-4",
                                ifelse(clinic.mom.df$CCI <= 8, "5-8", "> 8"))

clinic.mom.df$CCI_adj <- factor(clinic.mom.df$CCI_adj, levels = c("0-4", "5-8", "> 8"))

vars.to.check <- c("raw.1st.dim0.mean", "Age.Diag", "Sex", "KPS_adj", "CCI_adj", "T.Stage", "Prior.Surg", "Pre.Rad.Chemo", "Post.Rad.Chemo")


#Consider Using M0 Patients only, and look into what prior surgery means
clinic.mom.df_MO <- subset(clinic.mom.df, M.Stage == "M0")

####Running Cox Model for OS, PFS, Local Control, and Cancer Specific Survival####


#Basic Cox Model Function
univ.Cox <- function(stat, int, vars, dat) {
  #List object to collect data
  list.UR <- list()
  #For loop to run through all the vars
  for (i in 1:length(vars)) {
    #Creating the relevant formula
    form <- as.formula(paste("Surv(", int, ",", stat, ")", "~", vars[i], sep = ""))
    
    #Computing the Cox Model
    res.cox <- coxph(form, data = dat)
    
    # CI
    CI <- exp(confint(res.cox))
    
    #Original HR
    HR <- exp(coef(res.cox))
    
    #Pvalue
    pval <- summary(res.cox)$coefficients[,5]
    
    #Collecting the data and sticking it into the list
    list.UR[[i]] <- cbind(HR, CI, pval)
    
  }
  
  #With the list object, turn it into a data frame
  df.UR <- do.call(rbind, list.UR)
  
  return(df.UR)
  
}
multv.Cox <- function(stat, int, vars, dat) {
  #Creating the relevant formula
  form <- as.formula(paste(paste("Surv(", int, ",", stat, ")"), 
                           paste(vars, collapse=" + "), sep=" ~ "))
  
  #Computing the Cox Model
  res.cox <- coxph(form, data = dat)
  
  #lower CI
  CI <- exp(confint(res.cox))
  
  
  
  #Original HR
  HR <- exp(coef(res.cox))
  
  #Pvalue
  pval <- summary(res.cox)$coefficients[,5]
  
  df.mr <- cbind(HR, CI, pval)
  
  return(df.mr)
}


#UV and MV Models for Survival
univ.Cox("Vital.Status", "Survival.Correct.Interval", vars.to.check, clinic.mom.df_MO) %>% signif(2) %>% as.data.frame %>%
  plyr::mutate(signif = ifelse(pval < .05, "*", "ns"))

multv.Cox("Vital.Status", "Survival.Correct.Interval", vars.to.check, clinic.mom.df_MO) %>% signif(3) %>% as.data.frame %>%
  plyr::mutate(signif = ifelse(pval < .05, "*", "ns"))

#UV and MV Models for PFS
univ.Cox("Progression.Status", "PFS.Correct.Interval", vars.to.check, clinic.mom.df) %>% signif(2) %>% as.data.frame %>%
  mutate(signif = ifelse(pval < .05, "*", "ns"))

multv.Cox("Progression.Status", "PFS.Correct.Interval", vars.to.check, clinic.mom.df) %>% signif(2) %>% as.data.frame %>%
  mutate(signif = ifelse(pval < .05, "*", "ns"))

#UV and MV Models for Local Control
univ.Cox("Local.Status", "Local.Correct.Interval", vars.to.check, clinic.mom.df) %>% signif(2) %>% as.data.frame %>%
  mutate(signif = ifelse(pval < .05, "*", "ns"))

multv.Cox("Local.Status", "Local.Correct.Interval", vars.to.check, clinic.mom.df) %>% signif(2) %>% as.data.frame %>%
  mutate(signif = ifelse(pval < .05, "*", "ns"))

#UV and MV Models for Cancer Specific Survival
#Getting the indices indicating cancer death
#Manually Coded From Reading the reasons
cancer_death <- read.csv("cancer_death.csv")

clinic.mom.df <- join(clinic.mom.df, cbind.data.frame(PID = cancer_death$PID, Canc.Spec = cancer_death$Cancer.Specific))

inds.for.sure <- which(clinic.mom.df$Canc.Spec == "no")

#Creating a Separate Dataframe for this stuff
clinic.mom.df.cancer.spec <- clinic.mom.df

#Zeroing Vital Status for non cancer related deaths
clinic.mom.df.cancer.spec$Vital.Status[inds.for.sure] <- 0

#Recalculating Correct Interval, using data from charts
clinic.mom.df.cancer.spec$Survival.Correct.Interval <- as.numeric(ifelse(clinic.mom.df.cancer.spec$Vital.Status == 1, 
                                                                         clinic.mom.df.cancer.spec$Date.Death - clinic.mom.df.cancer.spec$DODiag, 
                                                                         clinic.mom.df.cancer.spec$Date.Last.Followup - clinic.mom.df.cancer.spec$DODiag))

clinic.mom.df.cancer.spec_M0 <- subset(clinic.mom.df.cancer.spec, M.Stage == "M0")

#Putting a 2 for Vital Status for non cancer related deaths for competing risks
clinic.mom.df.cancer.spec_CR <- clinic.mom.df
clinic.mom.df.cancer.spec_CR$Vital.Status[inds.for.sure] <- 2

clinic.mom.df.cancer.spec_CR$Survival.Correct.Interval <- as.numeric(ifelse(clinic.mom.df.cancer.spec_CR$Vital.Status != 0, 
                                                                            clinic.mom.df.cancer.spec_CR$Date.Death - clinic.mom.df.cancer.spec_CR$DODiag, 
                                                                            clinic.mom.df.cancer.spec_CR$Date.Last.Followup - clinic.mom.df.cancer.spec_CR$DODiag))

clinic.mom.df.cancer.spec_CR_M0 <- subset(clinic.mom.df.cancer.spec_CR, M.Stage == "M0")


univ.Cox("Vital.Status", "Survival.Correct.Interval", vars.to.check, clinic.mom.df.cancer.spec_M0) %>% signif(3) %>% as.data.frame %>%
  plyr::mutate(signif = ifelse(pval < .05, "*", "ns"))

multv.Cox("Vital.Status", "Survival.Correct.Interval", vars.to.check, clinic.mom.df.cancer.spec_M0) %>% signif(3) %>% as.data.frame %>%
  plyr::mutate(signif = ifelse(pval < .05, "*", "ns"))


####Testing Proportional Hazard Assumptions Plot COME BACK####
fit <- coxph(Surv(Survival.Correct.Interval, Vital.Status) ~ raw.1st.dim0.mean+Age.Diag+Sex+KPS_adj+CCI_adj+T.Stage+
               Prior.Surg+Pre.Rad.Chemo+Post.Rad.Chemo,  data=clinic.mom.df.cancer.spec_M0) 
temp <- cox.zph(fit) 
print(temp)                 
plot(temp)


####Kaplan Meier Curves for OS and Cancer Specific Survival####
#Overall Survival
head(clinic.mom.df_MO)

#Cancer Specific Survival
head(clinic.mom.df.cancer.spec_M0)

#CSC for Competing risks
head(clinic.mom.df.cancer.spec_CR_M0)

#Start with OS
data_valid <- subset(clinic.mom.df_MO, 
                     Overall.Stage %in% 
                       c("IA", "IB", "IIA", "IIB") & 
                       Prior.Surg == "No" & 
                       Pre.Rad.Chemo == "No") %>% 
  select(Survival.Correct.Interval, Vital.Status, 
         raw.1st.dim0.mean, Age.Diag, Sex, KPS_adj, 
         CCI_adj, T.Stage, Prior.Surg, Pre.Rad.Chemo, Post.Rad.Chemo, Overall.Stage)

quart.mom <- quantile(data_valid$raw.1st.dim0.mean)

#PRINT FOR MANUSCRIPT
print(quart.mom)


data_valid$cat <- ifelse(data_valid$raw.1st.dim0.mean <= quart.mom[2], "Mom25",
                         ifelse(data_valid$raw.1st.dim0.mean <= quart.mom[3], "Mom50",
                                ifelse(data_valid$raw.1st.dim0.mean <= quart.mom[4], "Mom75",
                                       ifelse(data_valid$raw.1st.dim0.mean <= quart.mom[5], "Mom100", NA))))

data_valid$cat <- factor(data_valid$cat,
                         levels = c("Mom25", "Mom50", "Mom75", "Mom100"))

fit <- survfit(Surv(Survival.Correct.Interval, Vital.Status) ~ cat,
               data = data_valid)

#Organizing the factor levels
data_valid$cat <- revalue(data_valid$cat, c("Mom100" = "Highest PHOM Score Group", 
                                            "Mom25" = "Lowest PHOM Score Group", 
                                            "Mom50" = "Low PHOM Score Group", 
                                            "Mom75" = "High PHOM Score Group"))

data_valid$cat <- factor(data_valid$cat, levels = c("Lowest PHOM Score Group", 
                                                    "Low PHOM Score Group", 
                                                    "High PHOM Score Group",
                                                    "Highest PHOM Score Group"))

#Getting the fit
fit <- survfit(Surv(Survival.Correct.Interval, Vital.Status) ~ cat, data = data_valid)

#Getting the median survival
med.surv <- surv_median(fit)
med.surv$strata <- revalue(med.surv$strata, c("cat=Lowest PHOM Score Group" = "Lowest PHOM Score Group", 
                                              "cat=Low PHOM Score Group" = "Low PHOM Score Group", 
                                              "cat=High PHOM Score Group" = "High PHOM Score Group", 
                                              "cat=Highest PHOM Score Group" = "Highest PHOM Score Group"))
colnames(med.surv) <- c("strata", "median", "lowerlim", "upperlim")

#This log rank compares all groups
log.rank.all <- survdiff(Surv(Survival.Correct.Interval, Vital.Status) ~ cat,
                         data = data_valid)

pval.tot <- 1 - pchisq(log.rank.all$chisq, length(log.rank.all$n) - 1)

#This log rank post hoc compares pairwise groups
statkm <- pairwise_survdiff(Surv(Survival.Correct.Interval, Vital.Status) ~ cat, 
                            data = data_valid,
                            p.adjust.method = "bonferroni")

#Creating a stat object DF for ggplot
stat.km.df <- data.frame(group1 = rep(NA, 6),
                         group2 = rep(NA, 6), 
                         p.adj = rep(NA, 6),
                         p.adj.signif = rep(NA, 6))

stat.km.df$group1 = c(rep("Lowest PHOM Score Group", 3),
                      rep("Low PHOM Score Group", 2),
                      rep("High PHOM Score Group", 1))
stat.km.df$group2 = c("Low PHOM Score Group",
                      "High PHOM Score Group",
                      "Highest PHOM Score Group",
                      "High PHOM Score Group",
                      "Highest PHOM Score Group",
                      "Highest PHOM Score Group")

stat.km.df$p.adj[1:3] <- statkm$p.value[, "Lowest PHOM Score Group"]
stat.km.df$p.adj[4:5] <- na.omit(statkm$p.value[, "Low PHOM Score Group"])
stat.km.df$p.adj[6] <- na.omit(statkm$p.value[, "High PHOM Score Group"])
stat.km.df$p.adj.signif <- ifelse(stat.km.df$p.adj < .05, "*", "ns")

#Converting group names to x positions
stat.km.df$group1 <- case_when(stat.km.df$group1 == med.surv$strata[1] ~ med.surv$median[1],
                               stat.km.df$group1 == med.surv$strata[2] ~ med.surv$median[2],
                               stat.km.df$group1 == med.surv$strata[3] ~ med.surv$median[3],
                               stat.km.df$group1 == med.surv$strata[4] ~ med.surv$median[4])
stat.km.df$group2 <- case_when(stat.km.df$group2 == med.surv$strata[1] ~ med.surv$median[1],
                               stat.km.df$group2 == med.surv$strata[2] ~ med.surv$median[2],
                               stat.km.df$group2 == med.surv$strata[3] ~ med.surv$median[3],
                               stat.km.df$group2 == med.surv$strata[4] ~ med.surv$median[4])

stat.km.df$p.adj <- signif(stat.km.df$p.adj, 2)

#KM Curves  graph
kmcurve <- autoplot(fit,
                    data = data_valid,
                    conf.int = F,
                    censor.shape = "|", 
                    censor.size = 2,
                    censor.colour = "#393D3F",
                    surv.size = 1) + 
  geom_vline(aes(xintercept = median, color = strata), alpha = .5, data = med.surv, size = 1) +
  scale_y_continuous(expand = c(0, .03), breaks = c(0, .2, .4, .6, .8, 1)) + 
  scale_color_manual(values = c("#85BBAC", "#55776F", "#DB5461", "#8A4950")) + 
  theme_classic() + 
  theme(legend.justification=c(1,1), legend.position=c(1,.88),
        plot.title = element_text(hjust = 0.5, size = 20, family = "Arial", color = "#393D3F"),
        axis.line = element_line(size = 1, color = "#393D3F"),
        axis.ticks = element_line(color = "#393D3F"),
        legend.text = element_text(size = 12, family = "Arial", color = "#393D3F"),
        legend.title = element_text(size = 14, family = "Arial", color = "#393D3F"),
        axis.text.x = element_text(size = 16, family = "Arial", color = "#393D3F"),
        axis.text.y = element_text(size = 16, family = "Arial", color = "#393D3F"),
        axis.title = element_text(size = 16, family = "Arial", color = "#393D3F"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)) + 
  labs(title = "Overall Survival By Persistent Homology (PHOM) Score", x = "Survival Time (days)", 
       y = "Overall Survival Probability", color = "Persistent Homology (PHOM) Score Groups") + 
  guides(colour = guide_legend(override.aes = list(shape = 15)))


#Using the ***, **, * Standard
stat.km.df$p.adj.signif <- ifelse(stat.km.df$p.adj < .001, "***", 
                                  ifelse(stat.km.df$p.adj < .01, "**",
                                         ifelse(stat.km.df$p.adj < .05, "*", "ns")))

kmcurve.stat <- kmcurve +
  stat_pvalue_manual(stat.km.df, label = "p.adj.signif",
                     y.position = c(.8, .7, .15),
                     label.size = 8, hide.ns = TRUE, bracket.size = .75,
                     bracket.nudge.y = .1, tip.length = 0, color = "#393D3F") +
  annotate("text", x = 140, y = .1, size = 4,
           #label = paste("p = ", format(signif(pval.tot, 2), scientific = FALSE)))
           label = paste("p <<< .001 "), color = "#393D3F")


ggsave("./Figures/KM_Curves_OS.png", plot = kmcurve.stat,
       scale = 1, width = 10, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)

##Now doing Cancer Specific Survival
data_valid <- subset(clinic.mom.df.cancer.spec_CR_M0, 
                     Overall.Stage %in% 
                       c("IA", "IB", "IIA", "IIB") & 
                       Prior.Surg == "No" & 
                       Pre.Rad.Chemo == "No") %>% 
  select(Survival.Correct.Interval, Vital.Status, 
         raw.1st.dim0.mean, Age.Diag, Sex, KPS_adj, 
         CCI_adj, T.Stage, Prior.Surg, Pre.Rad.Chemo, 
         Post.Rad.Chemo, Overall.Stage)

quart.mom <- quantile(data_valid$raw.1st.dim0.mean)

#PRINT FOR MANUSCRIPT
print(quart.mom)

#Separating MOM scores into 75% and below, then 75% and above
data_valid$cat <- ifelse(data_valid$raw.1st.dim0.mean <= quart.mom[4], "Momlow", "Momhi")

data_valid$cat <- factor(data_valid$cat,
                         levels = c("Momlow", "Momhi"))


#Competing Risks Cumulative Incidence, test statistics Gray's Test
fit3 <- cuminc(ftime = data_valid$Survival.Correct.Interval, 
               fstatus = data_valid$Vital.Status,
               group = data_valid$cat)

names(fit3) <- c("Low PHOM Score_Cancer Death", "Higest PHOM Score_Cancer Death", 
                 "Low PHOM Score_Other Death", "Higest PHOM Score_Other Death", "Tests")

CI_curve <- ggcompetingrisks(fit3, palette = "Dark2",
                             legend = "right",
                             gsep = "_",
                             ggtheme = theme_classic(), multiple_panels = FALSE) +
  scale_color_manual(values = c("#C72006", "#7000B9")) +
  scale_linetype_manual(values = c("solid", "dashed")) + 
  theme(legend.justification=c(1,1), legend.position=c(1,.88),
        plot.title = element_text(hjust = 0.5, size = 20, family = "Arial", color = "#393D3F"),
        axis.line = element_line(size = 1, color = "#393D3F"),
        axis.ticks = element_line(color = "#393D3F"),
        legend.text = element_text(size = 12, family = "Arial", color = "#393D3F"),
        legend.title = element_text(size = 14, family = "Arial", color = "#393D3F"),
        axis.text.x = element_text(size = 16, family = "Arial", color = "#393D3F"),
        axis.text.y = element_text(size = 16, family = "Arial", color = "#393D3F"),
        axis.title = element_text(size = 16, family = "Arial", color = "#393D3F"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)) + 
  labs(title = "Cumulative Incidence by Persistent Homology (PHOM) Score", x = "Survival Time (days)", 
       y = "Probability of Event", color = "Cause of Death", linetype = "PHOM Score Group") + 
  guides(colour = guide_legend(override.aes = list(shape = 15)))


CI_curve <- CI_curve +
  annotate("text", x = 2000, y = .33, size = 5,
           label = paste("p = ", format(signif(fit3$Tests[1,"pv"], 2), 
                                        scientific = FALSE)), 
           color = "#C72006") + 
  annotate("text", x = 2000, y = .65, size = 5,
           label = paste("p = ", format(signif(fit3$Tests[2,"pv"], 2), 
                                        scientific = FALSE)), 
           color = "#7000B9")

ggsave("./Figures/CI_curve_CSS.png", plot = CI_curve,
       scale = 1, width = 10, height = 6, units = "in",
       dpi = 400, limitsize = TRUE)

#####Cancer Specific and Overall Survival Cox Models#####
#Overall Survival
clinic.mom.df_MO

#Cancer Specific Survival
clinic.mom.df.cancer.spec_M0

#Starting with OS
data_valid <- subset(clinic.mom.df_MO, 
                     Overall.Stage %in% 
                       c("IA", "IB", "IIA", "IIB") & 
                       Prior.Surg == "No" & 
                       Pre.Rad.Chemo == "No") %>% 
  select(Survival.Correct.Interval, Vital.Status, 
         raw.1st.dim0.mean, Age.Diag, Sex, KPS_adj, 
         CCI_adj, T.Stage, Prior.Surg, Pre.Rad.Chemo, Post.Rad.Chemo, Overall.Stage)

data_valid <- na.omit(data_valid)
d <- datadist(data_valid)
options(datadist = "d")

##Cox Model RF
MV.Cox <- multv.Cox("Vital.Status", "Survival.Correct.Interval", 
                    c("raw.1st.dim0.mean", "Age.Diag", 
                      "Sex", "KPS_adj", "CCI_adj", 
                      "Overall.Stage", "Post.Rad.Chemo"), 
                    data_valid) %>% signif(3) %>% as.data.frame() %>%
  plyr::mutate(., model = "MV") %>% tibble::rownames_to_column("label")

UV.Cox <- univ.Cox("Vital.Status", "Survival.Correct.Interval", 
                   c("raw.1st.dim0.mean", "Age.Diag", 
                     "Sex", "KPS_adj", "CCI_adj", 
                     "Overall.Stage", "Post.Rad.Chemo"), 
                   data_valid) %>% signif(3) %>% as.data.frame() %>%
  plyr::mutate(., model = "UV")  %>% tibble::rownames_to_column("label")

Cox.tot <- rbind(UV.Cox, MV.Cox)

colnames(Cox.tot) <- c("label", "Hazard Ratio", "Lower Bound", "Upper Bound", 
                       "p-value", "model")

Cox.tot$label <- Cox.tot$label %>%
  revalue(., c("raw.1st.dim0.mean"="PHOM Score", 
               "Age.Diag"="Age at Diagnosis",
               "SexMale"="Male Sex",
               "KPS_adjKPS ≥ 70" = "KPS ≥ 70",
               "CCI_adj5-8" = "CCI = 5-8", 
               "CCI_adj> 8" = "CCI > 8", 
               "Overall.StageIB" = "Stage 1b", 
               "Overall.StageIIA" = "Stage 2a", 
               "Overall.StageIIB" = "Stage 2b", 
               "Post.Rad.ChemoYes" = "Had Post-Radiation Chemo"))

Cox.tot$model <- Cox.tot$model %>%
  revalue(., c("UV"="Univariate Model", 
               "MV"="Multivariate Model"))


Cox.tot$label <- factor(Cox.tot$label, levels = rev(c("PHOM Score", "Age at Diagnosis", "Male Sex", 
                                                      "KPS ≥ 70", "CCI = 5-8", "CCI > 8", 
                                                      "Had Post-Radiation Chemo", 
                                                      "Stage 1b",  "Stage 2a", "Stage 2b")))


Cox.tot$model <- factor(Cox.tot$model, levels = c("Univariate Model", "Multivariate Model"))

coxforest_OS <- ggplot(data=Cox.tot, aes(x=label, y=`Hazard Ratio`, ymin=`Lower Bound`, ymax=`Upper Bound`)) +
  # geom_pointrange(color='black', shape=19, size = .4, fatten = .01) + 
  geom_errorbar(width=0.5, size=1, color = "#393D3F") + #393D3F
  geom_point(size=5, shape=18, color = "#393D3F") +
  facet_wrap(~model) +
  scale_y_continuous(trans='log10', limits = c(.1, 9),
                     breaks = c(.15, .4, 1, 2.5, 6)) + 
  geom_hline(yintercept=1, lty=2, size = .2, color = "#393D3F") +  
  coord_flip() +
  labs(title = "Overall Survival", x = "Cox Variable", 
       y = "Hazard Ratio (log scaled) (95% CI)") + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, family = "Arial", color = "#393D3F"),
        axis.text.x = element_text(size = 13, family = "Arial", color = "#393D3F"),
        axis.text.y = element_text(size = 13, family = "Arial", color = "#393D3F"),
        strip.text.x = element_text(size = 13, family = "Arial", color = "#393D3F"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title=element_text(size=13, family = "Arial", color = "#393D3F"),
        axis.title.y = element_blank(),
        panel.spacing = unit(6, "mm"),
        strip.background = element_rect(size = 1, color = "#393D3F", fill = NA),
        panel.border = element_rect(size = 1, color = "#393D3F", fill = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA))

ggsave("./Figures/coxforest_OS.png", plot = coxforest_OS,
       scale = 1, width = 8, height = 8, units = "in",
       dpi = 400, limitsize = TRUE)

#No dowing Cancer Specific Survival
data_valid <- subset(clinic.mom.df.cancer.spec_M0, 
                     Overall.Stage %in% 
                       c("IA", "IB", "IIA", "IIB") & 
                       Prior.Surg == "No" & 
                       Pre.Rad.Chemo == "No") %>% 
  select(Survival.Correct.Interval, Vital.Status, 
         raw.1st.dim0.mean, Age.Diag, Sex, KPS_adj, 
         CCI_adj, T.Stage, Prior.Surg, Pre.Rad.Chemo, Post.Rad.Chemo, Overall.Stage)

data_valid <- na.omit(data_valid)
d <- datadist(data_valid)
options(datadist = "d")

##Cox Model RF
MV.Cox <- multv.Cox("Vital.Status", "Survival.Correct.Interval", 
                    c("raw.1st.dim0.mean", "Age.Diag", 
                      "Sex", "KPS_adj", "CCI_adj", 
                      "Overall.Stage", "Post.Rad.Chemo"), 
                    data_valid) %>% signif(3) %>% as.data.frame() %>%
  plyr::mutate(., model = "MV") %>% tibble::rownames_to_column("label")

UV.Cox <- univ.Cox("Vital.Status", "Survival.Correct.Interval", 
                   c("raw.1st.dim0.mean", "Age.Diag", 
                     "Sex", "KPS_adj", "CCI_adj", 
                     "Overall.Stage", "Post.Rad.Chemo"), 
                   data_valid) %>% signif(3) %>% as.data.frame() %>%
  plyr::mutate(., model = "UV")  %>% tibble::rownames_to_column("label")

Cox.tot <- rbind(UV.Cox, MV.Cox)

colnames(Cox.tot) <- c("label", "Hazard Ratio", "Lower Bound", "Upper Bound", 
                       "p-value", "model")

Cox.tot$label <- Cox.tot$label %>%
  revalue(., c("raw.1st.dim0.mean"="PHOM Score", 
               "Age.Diag"="Age at Diagnosis",
               "SexMale"="Male Sex",
               "KPS_adjKPS ≥ 70" = "KPS ≥ 70",
               "CCI_adj5-8" = "CCI = 5-8", 
               "CCI_adj> 8" = "CCI > 8", 
               "Overall.StageIB" = "Stage 1b", 
               "Overall.StageIIA" = "Stage 2a", 
               "Overall.StageIIB" = "Stage 2b", 
               "Post.Rad.ChemoYes" = "Had Post-Radiation Chemo"))

Cox.tot$model <- Cox.tot$model %>%
  revalue(., c("UV"="Univariate Model", 
               "MV"="Multivariate Model"))


Cox.tot$label <- factor(Cox.tot$label, levels = rev(c("PHOM Score", "Age at Diagnosis", "Male Sex", 
                                                      "KPS ≥ 70", "CCI = 5-8", "CCI > 8", 
                                                      "Had Post-Radiation Chemo", 
                                                      "Stage 1b",  "Stage 2a", "Stage 2b")))


Cox.tot$model <- factor(Cox.tot$model, levels = c("Univariate Model", "Multivariate Model"))

coxforest_CSC <- ggplot(data=Cox.tot, aes(x=label, y=`Hazard Ratio`, ymin=`Lower Bound`, ymax=`Upper Bound`)) +
  # geom_pointrange(color='black', shape=19, size = .4, fatten = .01) + 
  geom_errorbar(width=0.5, size=1, color = "#393D3F") + #393D3F
  geom_point(size=5, shape=18, color = "#393D3F") +
  facet_wrap(~model) +
  scale_y_continuous(trans='log10', limits = c(.1, 9),
                     breaks = c(.15, .4, 1, 2.5, 6)) + 
  geom_hline(yintercept=1, lty=2, size = .2, color = "#393D3F") +  
  coord_flip() +
  labs(title = "Cancer Specific Survival", x = "Cox Variable", 
       y = "Hazard Ratio (log scaled) (95% CI)") + 
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, family = "Arial", color = "#393D3F"),
        axis.text.x = element_text(size = 13, family = "Arial", color = "#393D3F"),
        axis.text.y = element_text(size = 13, family = "Arial", color = "#393D3F"),
        strip.text.x = element_text(size = 13, family = "Arial", color = "#393D3F"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title=element_text(size=13, family = "Arial", color = "#393D3F"),
        axis.title.y = element_blank(),
        panel.spacing = unit(6, "mm"),
        strip.background = element_rect(size = 1, color = "#393D3F", fill = NA),
        panel.border = element_rect(size = 1, color = "#393D3F", fill = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA))

ggsave("./Figures/coxforest_CSC.png", plot = coxforest_CSC,
       scale = 1, width = 8, height = 8, units = "in",
       dpi = 400, limitsize = TRUE)

comb_For <- egg::ggarrange(coxforest_OS + 
                             theme(axis.title.x = element_blank()), 
                           coxforest_CSC + theme(axis.text.y = element_blank(),
                                                 axis.ticks.y = element_blank(),
                                                 axis.title.y = element_blank(),
                                                 axis.title.x = element_blank()), 
                           nrow = 1)

# Annotate the figure by adding a common labels
comb_For <- annotate_figure(comb_For,
                            bottom = text_grob("Hazard Ratio (log scaled)", 
                                               color = "#393D3F", size = 13,
                                               hjust = 0, vjust = 0)
)


ggexport(comb_For, filename = "./Figures/comb_for.png", 
         width = 4800, height = 2400, res = 400)



#####Prediction Cox Model, Assessing Importance of PHOM####
#Only Doing OS Because 
#CSC Looked Wonky
data_valid <- subset(clinic.mom.df_MO, 
                     Overall.Stage %in% 
                       c("IA", "IB", "IIA", "IIB") & 
                       Prior.Surg == "No" & 
                       Pre.Rad.Chemo == "No") %>% 
  select(Survival.Correct.Interval, Vital.Status, 
         raw.1st.dim0.mean, Age.Diag, Sex, KPS_adj, 
         CCI_adj, T.Stage, Prior.Surg, Pre.Rad.Chemo, Post.Rad.Chemo, Overall.Stage)



data_valid <- na.omit(data_valid)
d <- datadist(data_valid)
options(datadist = "d")

lungsurv <- Surv(data_valid$Survival.Correct.Interval, data_valid$Vital.Status)

#Somer's D and C index calculation table
#Seeing how much this index drops per dropping variable
#And looking at proprtion of overall xi squared
#Basically says how much each variable reduces the sum of squared errors
model_full <- cph(formula = lungsurv ~ 
                    raw.1st.dim0.mean + 
                    Age.Diag + 
                    Sex + 
                    KPS_adj + 
                    CCI_adj + 
                    Overall.Stage +
                    Post.Rad.Chemo, 
                  data=data_valid, x = TRUE, y = TRUE, surv = TRUE, time.inc = 500) 
anova_tab <- anova(model_full, what='proportion chisq')
model_full_SD <- rms::validate(model_full, method = "boot", B=100)


#Knocking PHOM
model_phom <- cph(formula = lungsurv ~ 
                    #raw.1st.dim0.mean + 
                    Age.Diag + 
                    Sex + 
                    KPS_adj + 
                    CCI_adj + 
                    Overall.Stage +
                    Post.Rad.Chemo, 
                  data=data_valid, x = TRUE, y = TRUE, surv = TRUE)
model_phom_SD <- rms::validate(model_phom, method = "boot", B=1000)

#Knocking Age.Diag
model_age <- cph(formula = lungsurv ~ 
                   raw.1st.dim0.mean + 
                   #Age.Diag + 
                   Sex + 
                   KPS_adj + 
                   CCI_adj + 
                   Overall.Stage +
                   Post.Rad.Chemo, 
                 data=data_valid, x = TRUE, y = TRUE, surv = TRUE)
model_age_SD <- rms::validate(model_age, method = "boot", B=1000)


#Knocking Sex
model_sex <- cph(formula = lungsurv ~ 
                   raw.1st.dim0.mean + 
                   Age.Diag + 
                   #Sex + 
                   KPS_adj + 
                   CCI_adj + 
                   Overall.Stage +
                   Post.Rad.Chemo, 
                 data=data_valid, x = TRUE, y = TRUE, surv = TRUE)
model_sex_SD <- rms::validate(model_sex, method = "boot", B=1000)


#Knocking KPS
model_kps <- cph(formula = lungsurv ~ 
                   raw.1st.dim0.mean + 
                   Age.Diag + 
                   Sex + 
                   #KPS_adj + 
                   CCI_adj + 
                   Overall.Stage +
                   Post.Rad.Chemo, 
                 data=data_valid, x = TRUE, y = TRUE, surv = TRUE)
model_kps_SD <- rms::validate(model_kps, method = "boot", B=1000)

#Knocking CCI
model_cci <- cph(formula = lungsurv ~ 
                   raw.1st.dim0.mean + 
                   Age.Diag + 
                   Sex + 
                   KPS_adj + 
                   #CCI_adj + 
                   Overall.Stage +
                   Post.Rad.Chemo, 
                 data=data_valid, x = TRUE, y = TRUE, surv = TRUE)
model_cci_SD <- rms::validate(model_cci, method = "boot", B=1000)

#Knocking Stage
model_stage <- cph(formula = lungsurv ~ 
                     raw.1st.dim0.mean + 
                     Age.Diag + 
                     Sex + 
                     KPS_adj + 
                     CCI_adj + 
                     #Overall.Stage +
                     Post.Rad.Chemo, 
                   data=data_valid, x = TRUE, y = TRUE, surv = TRUE)
model_stage_SD <- rms::validate(model_stage, method = "boot", B=1000)

#Knocking Chemo
model_chemo <- cph(formula = lungsurv ~ 
                     raw.1st.dim0.mean + 
                     Age.Diag + 
                     Sex + 
                     KPS_adj + 
                     CCI_adj + 
                     Overall.Stage, # +
                   #Post.Rad.Chemo, 
                   data=data_valid, x = TRUE, y = TRUE, surv = TRUE)
model_chemo_SD <- rms::validate(model_chemo, method = "boot", B=1000)

#Arranging all data into table
names(anova_tab[, "Chi-Square"]) <- c("PHOM Score", "Age", "Sex", 
                                      "KPS", "CCI", "Overall Stage", 
                                      "Post Radiation Chemo", "Full Model")
SD_tab <- c(
  model_phom_SD["Dxy", "index.corrected"],
  model_age_SD["Dxy", "index.corrected"],
  model_sex_SD["Dxy", "index.corrected"],
  model_kps_SD["Dxy", "index.corrected"],
  model_cci_SD["Dxy", "index.corrected"],
  model_stage_SD["Dxy", "index.corrected"],
  model_chemo_SD["Dxy", "index.corrected"],
  model_full_SD["Dxy", "index.corrected"]
)

names(SD_tab) <- c("PHOM Score", "Age", "Sex", 
                   "KPS", "CCI", "Overall Stage", 
                   "Post Radiation Chemo", "Full Model")

disc_tab <- cbind.data.frame(`C index` = SD_tab, `Chi-Square` = anova_tab[,"Chi-Square"])

#Actually Converting to C index
disc_tab$`C index` <- (disc_tab$`C index`/2) + .5

#Calculating Drop in C index
disc_tab$`Drop in C index` <- disc_tab[["Full Model", "C index"]] - disc_tab$`C index`

disc_tab <- disc_tab %>% signif(4)

#Saving the Table Output
#Again C index reflection of discriminatory power, model is weak, but best we have
#KPS, PHOM Score, and CCI have relatively greatest contribution to discriminatory pwoer
#And greatest contribution in reducing mean square error in overall survival (Chi square)
#After a 1000 boot straps
print(xtable(disc_tab, type = "latex", digits = c(4,4,4,4)), 
      file = "./Figures/disc_tab.tex")


####Prediction Models Cox Calibration Curves####
#Only Doing OS Because 
#CSC Looked Wonky
data_valid <- subset(clinic.mom.df_MO, 
                     Overall.Stage %in% 
                       c("IA", "IB", "IIA", "IIB") & 
                       Prior.Surg == "No" & 
                       Pre.Rad.Chemo == "No") %>% 
  select(Survival.Correct.Interval, Vital.Status, 
         raw.1st.dim0.mean, Age.Diag, Sex, KPS_adj, 
         CCI_adj, T.Stage, Prior.Surg, Pre.Rad.Chemo, Post.Rad.Chemo, Overall.Stage)



data_valid <- na.omit(data_valid)
d <- datadist(data_valid)
options(datadist = "d")

lungsurv <- Surv(data_valid$Survival.Correct.Interval, data_valid$Vital.Status)


#1 year
model_1 <- cph(formula = lungsurv ~ 
                 raw.1st.dim0.mean + 
                 Age.Diag + 
                 Sex + 
                 KPS_adj + 
                 CCI_adj + 
                 Overall.Stage +
                 Post.Rad.Chemo, 
               data=data_valid, x = TRUE, y = TRUE, surv = TRUE, time.inc = 365.25)

#2 year
model_2 <- cph(formula = lungsurv ~ raw.1st.dim0.mean + 
                 Age.Diag + 
                 Sex + 
                 KPS_adj + 
                 CCI_adj + 
                 Overall.Stage +
                 Post.Rad.Chemo, 
               data=data_valid, x = TRUE, y = TRUE, surv = TRUE, time.inc = 730.5)

#5 year
model_5 <- cph(formula = lungsurv ~ raw.1st.dim0.mean + 
                 Age.Diag + 
                 Sex + 
                 KPS_adj + 
                 CCI_adj + 
                 Overall.Stage +
                 Post.Rad.Chemo, 
               data=data_valid, x = TRUE, y = TRUE, surv = TRUE, time.inc = 1826.25)

#8 year
model_8 <- cph(formula = lungsurv ~ raw.1st.dim0.mean + 
                 Age.Diag + 
                 Sex + 
                 KPS_adj + 
                 CCI_adj + 
                 Overall.Stage +
                 Post.Rad.Chemo, 
               data=data_valid, x = TRUE, y = TRUE, surv = TRUE, time.inc = 2922)

#Creating all the calibration objects
dats_1 <- calibrate(model_1, B=1000, u=365.25)
dats_2 <- calibrate(model_2, B=1000, u=730.5)
dats_5 <- calibrate(model_5, B=1000, u=1826.25)
dats_8 <- calibrate(model_8, B=1000, u=2922)

#Saving as png
png("./Figures/calplots.png", width = 4800, height = 3600, res = 400)
#Seeting up margins
par(mar = c(5, 5, 3, 1), mfrow=c(2,2))
plot(dats_1, lwd=1, lty=1, 
     cex.axis=1, cex.main=1.2, cex.sub=0.6, 
     xlab="",
     ylab = "", subtitles = FALSE)
legend(x=.815, y=.77, legend=c("Apparent", "Bias-corrected", "Ideal"), 
       lty=c(1, 1, 1), bty="n", col=c("#202020", "#193fff", "#e9e9e9"))
title(ylab="Actual death (proportion)", line=2.3, cex.lab=1.2)
title(xlab="Predicted death (proportion)", line=2.3, cex.lab=1.2)
title(main="Predicted 1 year survival", line=.9, cex.lab=1.2)

plot(dats_2, lwd=1, lty=1, 
     cex.axis=1, cex.main=1.2, cex.sub=0.6, 
     xlab="",
     ylab = "", subtitles = FALSE)
legend(x=.575, y=.467, legend=c("Apparent", "Bias-corrected", "Ideal"), 
       lty=c(1, 1, 1), bty="n", col=c("#202020", "#193fff", "#e9e9e9"))
title(ylab="Actual death (proportion)", line=2.3, cex.lab=1.2)
title(xlab="Predicted death (proportion)", line=2.3, cex.lab=1.2)
title(main="Predicted 2 year survival", line=.9, cex.lab=1.2)

plot(dats_5, lwd=1, lty=1, 
     cex.axis=1, cex.main=1.2, cex.sub=0.6, 
     xlab="",
     ylab = "", subtitles = FALSE)
legend(x=.28, y=.17, legend=c("Apparent", "Bias-corrected", "Ideal"), 
       lty=c(1, 1, 1), bty="n", col=c("#202020", "#193fff", "#e9e9e9"))
title(ylab="Actual death (proportion)", line=2.3, cex.lab=1.2)
title(xlab="Predicted death (proportion)", line=2.3, cex.lab=1.2)
title(main="Predicted 5 year survival", line=.9, cex.lab=1.2)

plot(dats_8, lwd=1, lty=1, 
     cex.axis=1, cex.main=1.2, cex.sub=0.6, 
     xlab="",
     ylab = "", subtitles = FALSE)
legend(x=.15, y=.1, legend=c("Apparent", "Bias-corrected", "Ideal"), 
       lty=c(1, 1, 1), bty="n", col=c("#202020", "#193fff", "#e9e9e9"))
title(ylab="Actual death (proportion)", line=2.3, cex.lab=1.2)
title(xlab="Predicted death (proportion)", line=2.3, cex.lab=1.2)
title(main="Predicted 8 year survival", line=.9, cex.lab=1.2)

dev.off()

####Creating Nomogram for survival####
phmodel <- cph(formula = lungsurv ~ raw.1st.dim0.mean + 
                 Age.Diag + 
                 Sex + 
                 KPS_adj + 
                 CCI_adj + 
                 Overall.Stage +
                 Post.Rad.Chemo, 
               data=data_valid, x = TRUE, y = TRUE, surv = TRUE, 
               time.inc = 8*365.25)

surv<-Survival(phmodel)

surv1<-function(x) surv(1*365.25,lp=x)
surv2<-function(x) surv(2*365.25,lp=x)
surv5<-function(x) surv(5*365.25,lp=x)
surv8<-function(x) surv(8*365.25,lp=x)

quan<-Quantile(phmodel)
med<-function(x) quan(lp=x)/365.25
ss<-c(.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95)

nom <- nomogram(phmodel,fun=list(surv1, surv2, surv5, surv8, med), 
                lp = FALSE,
                funlabel=c('1-year Survival','2-year Survival',
                           '5-year Survival', '8-year Survival', 
                           'Median Survival (Years)'),
                fun.at=list(ss, ss, ss, ss, c(0, .5 ,1:8)))
names(nom) <- c("PHOM Score", "Age", "Sex", "KPS", 
                "CCI", "Overall Stage", "Post Radiation Chemo", 
                "total.points", '1-year Survival','2-year Survival',
                '5-year Survival', '8-year Survival', 
                'Median Survival (Years)')

#Saving as png
png("./Figures/nomogram.png", width = 4800, height = 3600, res = 400)
plot(nom, cex.var=1, cex.axis=.75, lmgp=.25)
dev.off()

###More Detailed Nomogram to Turn Into Excel Sheet##
ss_d<- seq(0, 1, .01) 
ssm_d<-seq(0, 8, .01) 
nom_d <- nomogram(phmodel,fun=list(surv1, surv2, surv5, surv8, med), 
                  lp = FALSE,
                  funlabel=c('1-year Survival','2-year Survival',
                             '5-year Survival', '8-year Survival', 
                             'Median Survival'),
                  fun.at=list(ss_d, ss_d, ss_d, ss_d, ssm_d))
names(nom_d) <- c("PHOM Score", "Age", "Sex", "KPS", 
                  "CCI", "Overall Stage", "Post Radiation Chemo", 
                  "total.points", '1-year Survival','2-year Survival',
                  '5-year Survival', '8-year Survival', 
                  'Median Survival')

#Gathering all the Points
point_risks <- rbind(
  cbind.data.frame(predictor = "PHOM Score", 
                   value = nom_d$`PHOM Score`$raw.1st.dim0.mean,
                   points = nom_d$`PHOM Score`$points),
  
  cbind.data.frame(predictor = "Age", 
                   value = nom_d$Age$Age.Diag,
                   points = nom_d$Age$points),
  
  cbind.data.frame(predictor = "Sex", 
                   value = nom_d$Sex$Sex,
                   points = nom_d$Sex$points),
  
  cbind.data.frame(predictor = "KPS", 
                   value = nom_d$KPS$KPS_adj,
                   points = nom_d$KPS$points),
  
  cbind.data.frame(predictor = "CCI", 
                   value = nom_d$CCI$CCI_adj,
                   points = nom_d$CCI$points),
  
  cbind.data.frame(predictor = "Overall Stage", 
                   value = nom_d$`Overall Stage`$Overall.Stage,
                   points = nom_d$`Overall Stage`$points),
  
  cbind.data.frame(predictor = "Post Radiation Chemo", 
                   value = nom_d$`Post Radiation Chemo`$Post.Rad.Chemo,
                   points = nom_d$`Post Radiation Chemo`$points)
  
)

one_year_dat <- cbind.data.frame(
  length = "1 year",
  points = nom_d$`1-year Survival`$x,
  prob = nom_d$`1-year Survival`$x.real
)

two_year_dat <- cbind.data.frame(
  length = "2 year",
  points = nom_d$`2-year Survival`$x,
  prob = nom_d$`2-year Survival`$x.real
)

five_year_dat <- cbind.data.frame(
  length = "5 year",
  points = nom_d$`5-year Survival`$x,
  prob = nom_d$`5-year Survival`$x.real
)

eight_year_dat <- cbind.data.frame(
  length = "8 year",
  points = nom_d$`8-year Survival`$x,
  prob = nom_d$`8-year Survival`$x.real
)

med_surv <- cbind.data.frame(
  length = "med",
  points = nom_d$`Median Survival`$x,
  prob = nom_d$`Median Survival`$x.real
)

#Rounding Point_Risk to single digits for simplicity
point_risks$points <- point_risks$points %>% round(0)

#The max number of points is 231, so we need to make a dataframe with 232 rows
max_points <- point_risks %>% group_by(predictor) %>% 
  dplyr::summarise(max = max(points)) %>% 
  select(max) %>% sum()

#Creating a final tally matrix to make a csv out of
mat_fin <- matrix(nrow = max_points + 1, ncol = 6)
colnames(mat_fin) <- c("Points", "1 Year Survival", "2 Year Survival",
                       "5 Year Survival", "8 Year Survival", "Median Survival")

#Adding the Points to this matrix
mat_fin[, "Points"] <- 0:max_points

#Adding 1 year survival to this table
for (i in 1:nrow(mat_fin)) {
  #Gets the point of interest
  point <- mat_fin[i, "Points"]
  
  #Gets the index of which point value is closes
  ind <- which(abs(one_year_dat$points-point)==min(abs(one_year_dat$points-point)))
  
  #Gets the probability of the point and adds it to the table
  mat_fin[i, "1 Year Survival"] <- one_year_dat$prob[ind]
}

#Adding 2 year survival to this table
for (i in 1:nrow(mat_fin)) {
  #Gets the point of interest
  point <- mat_fin[i, "Points"]
  
  #Gets the index of which point value is closes
  ind <- which(abs(two_year_dat$points-point)==min(abs(two_year_dat$points-point)))
  
  #Gets the probability of the point and adds it to the table
  mat_fin[i, "2 Year Survival"] <- two_year_dat$prob[ind]
}

#Adding 5 year survival to this table
for (i in 1:nrow(mat_fin)) {
  #Gets the point of interest
  point <- mat_fin[i, "Points"]
  
  #Gets the index of which point value is closes
  ind <- which(abs(five_year_dat$points-point)==min(abs(five_year_dat$points-point)))
  
  #Gets the probability of the point and adds it to the table
  mat_fin[i, "5 Year Survival"] <- five_year_dat$prob[ind]
}

#Adding 8 year survival to this table
for (i in 1:nrow(mat_fin)) {
  #Gets the point of interest
  point <- mat_fin[i, "Points"]
  
  #Gets the index of which point value is closes
  ind <- which(abs(eight_year_dat$points-point)==min(abs(eight_year_dat$points-point)))
  
  #Gets the probability of the point and adds it to the table
  mat_fin[i, "8 Year Survival"] <- eight_year_dat$prob[ind]
}

#Adding median survival to this table
for (i in 1:nrow(mat_fin)) {
  #Gets the point of interest
  point <- mat_fin[i, "Points"]
  
  #Gets the index of which point value is closes
  ind <- which(abs(med_surv$points-point)==min(abs(med_surv$points-point)))
  
  #Gets the probability of the point and adds it to the table
  mat_fin[i, "Median Survival"] <- med_surv$prob[ind]
}

#Writing these as csvs
write.csv(mat_fin, "./NomApp/mat_fin.csv")
write.csv(point_risks, "./NomApp/point_risks.csv")


#####Older Code, Retained for posterity####
#Random Survival Forests
library(randomForestSRC)
library(ggRandomForests)
library(pec)


#Overall Survival
clinic.mom.df_MO

#Cancer Specific Survival
clinic.mom.df.cancer.spec_M0

data_valid <- subset(clinic.mom.df.cancer.spec_M0, 
                     Overall.Stage %in% 
                       c("IA", "IB", "IIA", "IIB") & 
                       Prior.Surg == "No" & 
                       Pre.Rad.Chemo == "No") %>% 
  select(Survival.Correct.Interval, Vital.Status, 
         raw.1st.dim0.mean, Age.Diag, Sex, KPS_adj, 
         CCI_adj, T.Stage, Post.Rad.Chemo)

data_valid$Sex <- as.factor(data_valid$Sex)
data_valid$KPS_adj <- as.factor(data_valid$KPS_adj)
data_valid$CCI_adj <- as.factor(data_valid$CCI_adj)
data_valid$T.Stage <- as.factor(data_valid$T.Stage)
data_valid$Post.Rad.Chemo <- as.factor(data_valid$Post.Rad.Chemo)

names(data_valid) <- c("Survival.Correct.Interval", "Vital.Status", 
                       "PHOM_Score", 
                       "Age_Diagnosis", "Sex", "KPS_adj", "CCI_adj", 
                       "T_Stage", "Post_Rad_Chemo")

rfsrc_model <- rfsrc(Surv(Survival.Correct.Interval, Vital.Status) ~ 
                       PHOM_Score +
                       Age_Diagnosis + 
                       Sex +
                       KPS_adj + 
                       CCI_adj + 
                       T_Stage + 
                       Post_Rad_Chemo, 
                     data = data_valid, 
                     nsplit = 10, 
                     na.action = "na.impute", 
                     nodesize = 15,
                     tree.err = TRUE, importance = TRUE)
print(rfsrc_model)
plot(rfsrc_model)

#10 fold Cross Validation
# Sample the data and create a training subset.
train <- sample(1:nrow(data_valid), round(nrow(data_valid) * 0.9))
# Train the model.
data_valid.grow <- rfsrc(Surv(Survival.Correct.Interval, Vital.Status) ~ 
                           PHOM_Score +
                           Age_Diagnosis + 
                           Sex +
                           KPS_adj + 
                           CCI_adj + 
                           T_Stage + 
                           Post_Rad_Chemo, 
                         data = data_valid[train, ], 
                         nsplit = 10, 
                         na.action = "na.impute", 
                         nodesize = 15,
                         tree.err = TRUE, importance = TRUE)

# Test the model.
pred <- predict(data_valid.grow, data_valid[-train , ])

# Compare the results.
print(data_valid.grow)
print(pred)

var.select(rfsrc_model)



####Cox Forest PLots###



vars.to.check <- c("raw.1st.dim0.mean", "Age.Diag", "Sex", "KPS_adj", "CCI_adj", "T.Stage", "Prior.Surg", "Pre.Rad.Chemo", "Post.Rad.Chemo")



#UV and MV Models for Survival
UV.init <- univ.Cox("Vital.Status", "Survival.Correct.Interval", vars.to.check, clinic.mom.df) %>% signif(3) %>% as.data.frame 

MV.init <- multv.Cox("Vital.Status", "Survival.Correct.Interval", vars.to.check, clinic.mom.df) %>% signif(3) %>% as.data.frame 

UV.Cox <- cbind(UV.init, "UV")

MV.Cox <- cbind(MV.init, "MV")

UV.Cox <- tibble::rownames_to_column(UV.Cox, "row_names")
MV.Cox <- tibble::rownames_to_column(MV.Cox, "row_names")


colnames(UV.Cox) <- c("label", "HR", "LL", "UL", "pval", "model")

colnames(MV.Cox) <- c("label", "HR", "LL", "UL", "pval", "model")

Cox.tot <- rbind(UV.Cox, MV.Cox)

colnames(Cox.tot) <- c("label", "Hazard Ratio", "Lower Bound", "Upper Bound", 
                       "p-value", "model")

Cox.tot$label <- Cox.tot$label %>%
  revalue(., c("raw.1st.dim0.mean"="PHOM Score", 
               "Age.Diag"="Age at Diagnosis",
               "SexMale"="Male Sex",
               "KPS_adjKPS > 70" = "KPS > 70",
               "CCI_adj5-8" = "CCI = 5-8", 
               "CCI_adj> 8" = "CCI > 8", 
               "T.StageT1b" = "T Stage 1b", 
               "T.StageT2a" = "T Stage 2a", 
               "T.StageT2b" = "T Stage 2b", 
               "T.StageT3" = "T Stage T3", 
               "T.StageT4" = "T Stage T4", 
               "Prior.SurgYes" = "Had Prior Surgery", 
               "Pre.Rad.ChemoYes" = "Had Pre-radiation Chemo", 
               "Post.Rad.ChemoYes" = "Had Post-radiation Chemo"))

Cox.tot$model <- Cox.tot$model %>%
  revalue(., c("UV"="Univariate Model", 
               "MV"="Multivariate Model"))


Cox.tot$label <- factor(Cox.tot$label, levels = rev(c("PHOM Score", "Age at Diagnosis", "Male Sex", "KPS > 70", "CCI = 5-8", "CCI > 8", 
                                                      "Had Prior Surgery", "Had Pre-radiation Chemo", "Had Post-radiation Chemo", 
                                                      "T Stage 1b",  "T Stage 2a", "T Stage 2b", "T Stage T3", "T Stage T4")))

#Significant Truth Table
Truth_Labels <- join(
  subset(Cox.tot, model == "Univariate Model") %>% mutate(sigUV = ifelse(`p-value` < .05, "sig", "nsig")) %>% select(label, sigUV),
  subset(Cox.tot, model == "Multivariate Model") %>% mutate(sigMV = ifelse(`p-value` < .05, "sig", "nsig")) %>% select(label, sigMV)
)

Cox.tot$model <- factor(Cox.tot$model, levels = c("Univariate Model", "Multivariate Model"))

coxforest <- ggplot(data=Cox.tot, aes(x=label, y=`Hazard Ratio`, ymin=`Lower Bound`, ymax=`Upper Bound`)) +
  # geom_pointrange(color='black', shape=19, size = .4, fatten = .01) + 
  geom_errorbar(width=0.5, size=1, color = "#393D3F") + #393D3F
  geom_point(size=5, shape=18, color = "#393D3F") +
  facet_wrap(~model) +
  scale_y_continuous(trans='log10', limits = c(.2,4)) + 
  geom_hline(yintercept=1, lty=2, size = .2, color = "#393D3F") +  
  coord_flip() +
  labs(title = "Survival Forest Plot", x = "Cox Variable", 
       y = "log(Hazard Ratio) (95% CI)") + 
  theme_classic() +
  scale_x_discrete(labels = rev(c(expression(bold("PHOM Score")),
                                  expression(bold("Age at Diagnosis")),
                                  expression(bold("Male Sex")),
                                  expression(bold("KPS ≥ 70")),
                                  expression(bold("CCI = 5-8")),
                                  expression(bold("CCI > 8")),
                                  "Had Prior Surgery",
                                  "Had Pre-radiation Chemo",
                                  "Had Post-radiation Chemo",
                                  expression(italic("T Stage 1b")),
                                  expression(italic("T Stage 2a")),
                                  expression(italic("T Stage 2b")),
                                  "T Stage T3",
                                  "T Stage T4"))) +
  theme(plot.title = element_text(hjust = 0.5, size = 16, family = "Arial", color = "#393D3F"),
        axis.text.x = element_text(size = 13, family = "Arial", color = "#393D3F"),
        axis.text.y = element_text(size = 13, family = "Arial", color = "#393D3F"),
        strip.text.x = element_text(size = 13, family = "Arial", color = "#393D3F"),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title=element_text(size=13, family = "Arial", color = "#393D3F"),
        axis.title.y = element_blank(),
        panel.spacing = unit(6, "mm"),
        strip.background = element_rect(size = 1, color = "#393D3F", fill = NA),
        panel.border = element_rect(size = 1, color = "#393D3F", fill = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA))


coxforest


ggsave("./Figures/coxforest.png", plot = coxforest,
       scale = 1, width = 8, height = 8, units = "in",
       dpi = 400, limitsize = TRUE, bg = "transparent")




###KM Curve###
#Adding A Moment 1 quartile Variable
# quart.mom <- quantile(clinic.mom.df$raw.mom.mean.log)
quart.mom <- quantile(clinic.mom.df$raw.1st.dim0.mean)

clinic.mom.df$cat <- ifelse(clinic.mom.df$raw.1st.dim0.mean <= quart.mom[2], "Mom25",
                            ifelse(clinic.mom.df$raw.1st.dim0.mean <= quart.mom[3], "Mom50",
                                   ifelse(clinic.mom.df$raw.1st.dim0.mean <= quart.mom[4], "Mom75",
                                          ifelse(clinic.mom.df$raw.1st.dim0.mean <= quart.mom[5], "Mom100", NA))))

clinic.mom.df$cat <- factor(clinic.mom.df$cat,
                            levels = c("Mom25", "Mom50", "Mom75", "Mom100"))

fit <- survfit(Surv(Survival.Correct.Interval, Vital.Status) ~ cat,
               data = clinic.mom.df)

#Organizing the factor levels
clinic.mom.df$cat <- revalue(clinic.mom.df$cat, c("Mom100" = "Highest PHOM Score Group (0.90 to 6.05)", 
                                                  "Mom25" = "Lowest PHOM Score Group (-3.33 to -0.91)", 
                                                  "Mom50" = "Low PHOM Score Group (-0.91 to -0.10)", 
                                                  "Mom75" = "High PHOM Score Group (-0.10 to 0.90)"))

clinic.mom.df$cat <- factor(clinic.mom.df$cat, levels = c("Lowest PHOM Score Group (-3.33 to -0.91)", 
                                                          "Low PHOM Score Group (-0.91 to -0.10)", 
                                                          "High PHOM Score Group (-0.10 to 0.90)",
                                                          "Highest PHOM Score Group (0.90 to 6.05)"))

#Getting the fit
fit <- survfit(Surv(Survival.Correct.Interval, Vital.Status) ~ cat, data = clinic.mom.df)

#Getting the median survival
med.surv <- surv_median(fit)
med.surv$strata <- revalue(med.surv$strata, c("cat=Lowest PHOM Score Group (-3.33 to -0.91)" = "Lowest PHOM Score Group (-3.33 to -0.91)", 
                                              "cat=Low PHOM Score Group (-0.91 to -0.10)" = "Low PHOM Score Group (-0.91 to -0.10)", 
                                              "cat=High PHOM Score Group (-0.10 to 0.90)" = "High PHOM Score Group (-0.10 to 0.90)", 
                                              "cat=Highest PHOM Score Group (0.90 to 6.05)" = "Highest PHOM Score Group (0.90 to 6.05)"))
colnames(med.surv) <- c("strata", "median", "lowerlim", "upperlim")

#This log rank compares all groups
log.rank.all <- survdiff(Surv(Survival.Correct.Interval, Vital.Status) ~ cat,
                         data = clinic.mom.df)

pval.tot <- 1 - pchisq(log.rank.all$chisq, length(log.rank.all$n) - 1)

#This log rank post hoc compares pairwise groups
statkm <- pairwise_survdiff(Surv(Survival.Correct.Interval, Vital.Status) ~ cat, 
                            data = clinic.mom.df,
                            p.adjust.method = "bonferroni")

#Creating a stat object DF for ggplot
stat.km.df <- data.frame(group1 = rep(NA, 6),
                         group2 = rep(NA, 6), 
                         p.adj = rep(NA, 6),
                         p.adj.signif = rep(NA, 6))

stat.km.df$group1 = c(rep("Lowest PHOM Score Group (-3.33 to -0.91)", 3),
                      rep("Low PHOM Score Group (-0.91 to -0.10)", 2),
                      rep("High PHOM Score Group (-0.10 to 0.90)", 1))
stat.km.df$group2 = c("Low PHOM Score Group (-0.91 to -0.10)",
                      "High PHOM Score Group (-0.10 to 0.90)",
                      "Highest PHOM Score Group (0.90 to 6.05)",
                      "High PHOM Score Group (-0.10 to 0.90)",
                      "Highest PHOM Score Group (0.90 to 6.05)",
                      "Highest PHOM Score Group (0.90 to 6.05)")

stat.km.df$p.adj[1:3] <- statkm$p.value[, "Lowest PHOM Score Group (-3.33 to -0.91)"]
stat.km.df$p.adj[4:5] <- na.omit(statkm$p.value[, "Low PHOM Score Group (-0.91 to -0.10)"])
stat.km.df$p.adj[6] <- na.omit(statkm$p.value[, "High PHOM Score Group (-0.10 to 0.90)"])
stat.km.df$p.adj.signif <- ifelse(stat.km.df$p.adj < .05, "*", "ns")

#Converting group names to x positions
stat.km.df$group1 <- case_when(stat.km.df$group1 == med.surv$strata[1] ~ med.surv$median[1],
                               stat.km.df$group1 == med.surv$strata[2] ~ med.surv$median[2],
                               stat.km.df$group1 == med.surv$strata[3] ~ med.surv$median[3],
                               stat.km.df$group1 == med.surv$strata[4] ~ med.surv$median[4])
stat.km.df$group2 <- case_when(stat.km.df$group2 == med.surv$strata[1] ~ med.surv$median[1],
                               stat.km.df$group2 == med.surv$strata[2] ~ med.surv$median[2],
                               stat.km.df$group2 == med.surv$strata[3] ~ med.surv$median[3],
                               stat.km.df$group2 == med.surv$strata[4] ~ med.surv$median[4])

stat.km.df$p.adj <- signif(stat.km.df$p.adj, 2)

#KM Curves  graph
kmcurve <- autoplot(fit,
                    data = tab.KM,
                    conf.int = F,
                    censor.shape = "|", 
                    censor.size = 2,
                    censor.colour = "#393D3F",
                    surv.size = 1) + 
  geom_vline(aes(xintercept = median, color = strata), alpha = .5, data = med.surv, size = 1) +
  scale_y_continuous(expand = c(0, .03), breaks = c(0, .2, .4, .6, .8, 1)) + 
  scale_color_manual(values = c("#85BBAC", "#55776F", "#DB5461", "#8A4950")) + 
  theme_classic() + 
  theme(legend.justification=c(1,1), legend.position=c(1,.88),
        plot.title = element_text(hjust = 0.5, size = 20, family = "Arial", color = "#393D3F"),
        axis.line = element_line(size = 1, color = "#393D3F"),
        axis.ticks = element_line(color = "#393D3F"),
        legend.text = element_text(size = 12, family = "Arial", color = "#393D3F"),
        legend.title = element_text(size = 14, family = "Arial", color = "#393D3F"),
        axis.text.x = element_text(size = 16, family = "Arial", color = "#393D3F"),
        axis.text.y = element_text(size = 16, family = "Arial", color = "#393D3F"),
        axis.title = element_text(size = 16, family = "Arial", color = "#393D3F"),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA)) + 
  labs(title = "Survival By Persistent Homology (PHOM) Score", x = "Survival Time (days)", 
       y = "Survival Probability", color = "Persistent Homology (PHOM) Score Groups") + 
  guides(colour = guide_legend(override.aes = list(shape = 15)))


#Using the ***, **, * Standard
stat.km.df$p.adj.signif <- ifelse(stat.km.df$p.adj < .001, "***", 
                                  ifelse(stat.km.df$p.adj < .01, "**",
                                         ifelse(stat.km.df$p.adj < .05, "*", "ns")))

kmcurve.stat <- kmcurve +
  stat_pvalue_manual(stat.km.df, label = "p.adj.signif",
                     y.position = c(.8, .7, .15),
                     label.size = 8, hide.ns = TRUE, bracket.size = .75,
                     bracket.nudge.y = .1, tip.length = 0, color = "#393D3F") +
  annotate("text", x = 140, y = .1, size = 4,
           #label = paste("p = ", format(signif(pval.tot, 2), scientific = FALSE)))
           label = paste("p <<< .001 "), color = "#393D3F")


ggsave("./Figures/KM_Curves.png", plot = kmcurve.stat,
       scale = 1, width = 10, height = 6, units = "in",
       dpi = 400, limitsize = TRUE, bg = "transparent")


####Competing Risks Model Continue TRYING####
library(riskRegression)
library(prodlim)
library(survival)

#Using Competing Risks Data frame
data_valid <- subset(clinic.mom.df.cancer.spec_CR, 
                     Overall.Stage %in% 
                       c("IA", "IB", "IIA", "IIB") & 
                       Prior.Surg == "No" & 
                       Pre.Rad.Chemo == "No") %>% 
  select(Survival.Correct.Interval, Vital.Status, 
         raw.1st.dim0.mean, Age.Diag, Sex, KPS_adj, 
         CCI_adj, T.Stage, Prior.Surg, Pre.Rad.Chemo, 
         Post.Rad.Chemo, Overall.Stage) %>% na.omit()

data_valid <- na.omit(data_valid)
d <- datadist(data_valid)
options(datadist = "d")

data_valid$Sex <- as.factor(data_valid$Sex)
data_valid$KPS_adj <- as.factor(data_valid$KPS_adj)
data_valid$CCI_adj <- as.factor(data_valid$CCI_adj)
data_valid$T.Stage <- as.factor(data_valid$T.Stage)
data_valid$Post.Rad.Chemo <- as.factor(data_valid$Post.Rad.Chemo)

csc <- CSC(Hist(Survival.Correct.Interval, Vital.Status) ~ raw.1st.dim0.mean + 
             Age.Diag + 
             Sex + 
             KPS_adj + 
             CCI_adj + 
             Overall.Stage +
             Post.Rad.Chemo,
           data = data_valid) 

fgr <- FGR(Hist(Survival.Correct.Interval, Vital.Status) ~ raw.1st.dim0.mean + 
             Age.Diag + 
             Sex + 
             KPS_adj + 
             CCI_adj + 
             Overall.Stage +
             Post.Rad.Chemo,
           data = data_valid, cause = 1)


fgr_lim <- FGR(Hist(Survival.Correct.Interval, Vital.Status) ~ raw.1st.dim0.mean + 
                 Age.Diag + 
                 Sex + 
                 Overall.Stage,
               data = data_valid, cause = 1)


score <- Score(list("csc" = csc,
                    "fgr" = fgr,
                    "fgr_limg" = fgr_lim ),
               formula = Hist(Survival.Correct.Interval, Vital.Status) ~ 1, 
               data = data_valid, 
               times = seq(0, 3000, 100),
               plots = "calibration",
               split.method = "bootcv",
               B = 10,
               summary = "risks",
               contrasts=TRUE)


plotCalibration(score, times = 1000, cens.method="local", method = "nne")



# simulated data to test 
library(cmprsk)

ftime <- data_valid$Survival.Correct.Interval

fstatus <- data_valid$Vital.Status

cov <- data_valid[ , 3]

z <- crr(data_valid$Survival.Correct.Interval, data_valid$Vital.Status, 
         data_valid[ , 3:4], failcode=1)

print(z)
summary(z)



data_valid$Vital.Status <- factor(data_valid$Vital.Status, 0:2, 
                                  labels=c("censor", "cancer", "other"))

data <- finegray(Surv(Survival.Correct.Interval, Vital.Status) ~ raw.1st.dim0.mean + 
                   Age.Diag + 
                   Sex + 
                   KPS_adj + 
                   CCI_adj + 
                   T.Stage +
                   Post.Rad.Chemo, data=data_valid)
summary(data)







