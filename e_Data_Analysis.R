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
<<<<<<< HEAD
library(stringr)
library(hablar)
=======
>>>>>>> 72af4f36942260c26257abb83c725a1433290221

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
  dim.inds <- list.data$feat.curve.moms %>%
    names %>%
    str_detect(paste(dim.to.choose, collapse = "|"))
  
  #Logical Indices of moments to choose
  moms.dims.ind <- list.data$feat.curve.moms[dim.inds] %>%
    unlist(recursive = FALSE) %>%
    names %>%
    str_detect(paste(moms.to.choose.mod, collapse = "|"))
  
  #Lists of the moments we want
  cap.list <- list.data$feat.curve.moms[dim.inds] %>%
    unlist(recursive = FALSE) %>%
    extract(moms.dims.ind)
  
  #Extracting into a vector
  vec.col <- numeric()
  for (i in seq_len(length(cap.list))) {
    vec.col <- c(vec.col, cap.list[[i]])
  }
  
  fin.vec <- c(PID = PID, vec.col)
  
  return(fin.vec)
}


#Function to extract clinical variables from list element
#Possible variables here; need to use lapply to get full etraction
patient.vars <- c("Sex", "Ethnicity", "Age.Diag", "BMI",
                  "Years.Smoke", "KPS", "HGB", "CCI",
                  "Prev.Canc.Diag", "Num.Prev.Canc.Diag")
treatment.vars <- c("Prior.Surg", "Post.Rad.Chemo", "Pre.Rad.Chemo",
                    "Rad.Intent", "Tot.Dose")
tumor.vars <- c("CT.Size", "T.Stage", "N.Stage", "M.Stage",
                "Overall.Stage", "Histo", "Path")
event.vars <- c("Vital.Status", "OS.Length",
                "Progression.Status", "PFS.Length", "PFS.Failure.Type",
                "Local.Status", "Local.Length", 
                "Lobar.Status", "Nodal.Status", "Distant.Status", 
                "First.Met", "Update.Time", "Followup.Time", "Alive.Time")

clinic.extracter <- function(list.data, patient.vars, treatment.vars, tumor.vars, event.vars) {
  
  #Keeping patient ID
  PID <- list.data$Patient.ID
  
  #Placeholder Dataframe
  plac <- list.data$clinic.data[c(event.vars, patient.vars, treatment.vars, tumor.vars)]
  
  #Storing Names
  names <- names(plac)
  
  #1 0 convention for events
  # Distant.Failed = Distant Control
  # Lobar.Failed = Lobar Control
  # Nodal.Failed = Nodal Control
  # Failed = Local Control
  # Dead = Vital Status
  # Progressed = Progression Free Survival
  new.vec <- case_when(
    plac %in% c("Distant.Failed", "Lobar.Failed", "Nodal.Failed", "Failed", "Dead", "Progressed") ~ "1",
    plac %in% c("Distant.Controlled", "Lobar.Controlled", "Nodal.Controlled", "Controlled", "Alive", "Stable") ~ "0",
    TRUE ~ as.character(plac)
  )
  
  names(new.vec) <- names
  
  fin.vec <- c(PID = PID, new.vec)
  
  return(fin.vec)
  
}
#COMEBACK WITH VARIABLE SELECTION

#Replicating Old Paper Results
<<<<<<< HEAD
moms.list <- lapply(list.all.data, moment.extracter, mom.dims = mom.dim.vars, mom.types = mom.types.vars)
moms.df <- as.data.frame(do.call(rbind, moms.list))

#Rounding and filtering out 1 and 0 columns
moms.df[, -1] <- moms.df[, -1] %>% sapply(as.numeric) %>% round(5) %>% as.data.frame()

vecs <- sapply(moms.df, function (x) all(x == 0))
vecs2 <- sapply(moms.df, function (x) all(x == .999))

#Due to NAs, need to remove standcent dim 2 stuff manually
moms.df.corr <- moms.df[, -which(vecs == TRUE | vecs2 == TRUE)] %>% select(-c("standcent.1st.dim2", "standcent.2nd.dim2"))

#Logging...note produces expected NAs
moms.df.corr[moms.df.corr == 0] <- NA
moms.df.corr[, -1] <- moms.df.corr[, -1] %>% log()


clinic.list <- lapply(list.all.data, clinic.extracter, patient.vars = patient.vars, treatment.vars = treatment.vars,
                      tumor.vars = tumor.vars, event.vars = event.vars)
clinic.df <- as.data.frame(do.call(rbind, clinic.list))
=======
moms.list <- lapply(list.all.data, moment.extracter,
                    mom.dims = "dim0", mom.types = "raw")
moms.df <- do.call(rbind, moms.list) %>% as.data.frame

clinic.list <- lapply(list.all.data, clinic.extracter, patient.vars = c("Sex", "Age.Diag"), treatment.vars = NULL,
                      tumor.vars = "Overall.Stage", event.vars = c("Vital.Status", "OS.Length", "Followup.Time"))
clinic.df <- do.call(rbind, clinic.list) %>% as.data.frame
>>>>>>> 72af4f36942260c26257abb83c725a1433290221

#making proper numeric
clinic.df[, c("Vital.Status", "OS.Length", "Progression.Status", "PFS.Length", 
              "Local.Status", "Local.Length", "Lobar.Status", 
              "Nodal.Status", "Distant.Status", "Update.Time", 
              "Followup.Time", "Alive.Time", "Age.Diag", "CCI",
              "BMI", "Years.Smoke", "HGB", "Num.Prev.Canc.Diag", "Tot.Dose", "CT.Size")] <- 
  sapply(clinic.df[, c("Vital.Status", "OS.Length", "Progression.Status", "PFS.Length", 
                "Local.Status", "Local.Length", "Lobar.Status", 
                "Nodal.Status", "Distant.Status", "Update.Time", 
                "Followup.Time", "Alive.Time", "Age.Diag", "CCI",
                "BMI", "Years.Smoke", "HGB", "Num.Prev.Canc.Diag", "Tot.Dose", "CT.Size")], as.numeric)


#Combining moments and clinical vairable data frames
clinic.mom.df <- join(clinic.df, moms.df.corr)

sapply(clinic.mom.df, function(x) paste(sum(!is.na(x)), "/", length(x)))

#Creating a corrected Interval column that uses survival or follow up time
<<<<<<< HEAD
clinic.mom.df$Correct.Interval <- as.numeric(ifelse(clinic.mom.df$Vital.Status == "1", clinic.mom.df$OS.Length, clinic.mom.df$Followup.Time))

#Basic Cox Model
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


univ.Cox("Vital.Status", "Correct.Interval", 
         c("raw.1st.dim0.mean", "raw.1st.dim1.mean", 
           "cent.2nd.dim0.variance", "cent.2nd.dim1.variance", 
           "Age.Diag", "Sex", 
           "KPS", "CCI", 
           "Overall.Stage"), clinic.mom.df) %>% signif(3) %>% as.data.frame %>%
  mutate(signif = ifelse(pval < .05, "*", "ns"))

multv.Cox("Vital.Status", "Correct.Interval", 
         c("raw.1st.dim0.mean", "raw.1st.dim1.mean", 
           "cent.2nd.dim0.variance", "cent.2nd.dim1.variance", 
           "Age.Diag", "Sex", 
           "KPS", "CCI", 
           "Overall.Stage"), clinic.mom.df) %>% signif(3) %>% as.data.frame %>%
  mutate(signif = ifelse(pval < .05, "*", "ns"))





#Lasso model



=======
clinic.mom.df <- clinic.mom.df %>%
  mutate(Correct.Interval = ifelse(Vital.Status == "1", OS.Length, Followup.Time),
         Correct.Interval = as.numeric(Correct.Interval),
         Age.Diag = as.numeric(Age.Diag),
         raw.1st.dim0.mean = as.numeric(raw.1st.dim0.mean),
         raw.2nd.dim0 = as.numeric(raw.2nd.dim0),
         raw.3rd.dim0 = as.numeric(raw.3rd.dim0),
         raw.4th.dim0 = as.numeric(raw.4th.dim0),
         Vital.Status = as.numeric(Vital.Status))

coxph(Surv(Correct.Interval, Vital.Status) ~ log(raw.1st.dim0.mean),
      data = clinic.mom.df) %>%
  summary

coxph(Surv(Correct.Interval, Vital.Status) ~ log(raw.1st.dim0.mean) + Age.Diag + Overall.Stage + Sex,
      data = clinic.mom.df) %>%
  summary
>>>>>>> 72af4f36942260c26257abb83c725a1433290221

plot(log(clinic.mom.df$raw.1st.dim0.mean),
     clinic.mom.df$Correct.Interval)

clinic.mom.df <- mutate(clinic.mom.df, raw.mom.mean.log = log(raw.1st.dim0.mean))

clinic.mom.df$cat <- factor(clinic.mom.df$cat,
                            levels = c("Highest", "High", "Low", "Lowest"))

ggsurvplot(
  fit = survfit(Surv(Correct.Interval, Vital.Status) ~ cat, data = clinic.mom.df), 
  xlab = "Days", 
  legend.title = "Log of Raw Moment 1 of 0 Dimension Feature Curve",
  ylab = "Overall survival probability"
)

###PROPER PLOT HERE####
#Adding A Moment 1 quartile Variable
# quart.mom <- quantile(clinic.mom.df$raw.mom.mean.log)
quart.mom <- quantile(clinic.mom.df$raw.1st.dim0.mean)

cbind(log.transform = quantile(clinic.mom.df$raw.mom.mean.log),
      original = quantile(clinic.mom.df$raw.1st.dim0.mean))


clinic.mom.df$cat <- ifelse(clinic.mom.df$raw.mom.mean.log <= quart.mom[2], "Mom25",
                     ifelse(clinic.mom.df$raw.mom.mean.log <= quart.mom[3], "Mom50",
                            ifelse(clinic.mom.df$raw.mom.mean.log <= quart.mom[4], "Mom75",
                                   ifelse(clinic.mom.df$raw.mom.mean.log <= quart.mom[5], "Mom100", NA))))

clinic.mom.df$cat <- factor(clinic.mom.df$cat,
                            levels = c("Mom25", "Mom50", "Mom75", "Mom100"))

fit <- survfit(Surv(Correct.Interval, Vital.Status) ~ cat,
               data = clinic.mom.df)

#Organizing the factor levels
clinic.mom.df$cat <- revalue(clinic.mom.df$cat, c("Mom100" = "First Moment 75-100 percentile", 
                                                  "Mom25" = "First Moment 0-25 percentile", 
                                                  "Mom50" = "First Moment 25-50 percentile", 
                                                  "Mom75" = "First Moment 50-75 percentile"))

clinic.mom.df$cat <- factor(clinic.mom.df$cat, levels = c("First Moment 0-25 percentile", 
                                                          "First Moment 25-50 percentile", 
                                                          "First Moment 50-75 percentile",
                                                          "First Moment 75-100 percentile"))

#Getting the fit
fit <- survfit(Surv(Correct.Interval, Vital.Status) ~ cat, data = clinic.mom.df)

#Getting the median survival
med.surv <- surv_median(fit)
med.surv$strata <- revalue(med.surv$strata, c("cat=First Moment 0-25 percentile" = "First Moment 0-25 percentile", 
                                              "cat=First Moment 25-50 percentile" = "First Moment 25-50 percentile", 
                                              "cat=First Moment 50-75 percentile" = "First Moment 50-75 percentile", 
                                              "cat=First Moment 75-100 percentile" = "First Moment 75-100 percentile"))
colnames(med.surv) <- c("strata", "median", "lowerlim", "upperlim")

#This log rank compares all groups
log.rank.all <- survdiff(Surv(Correct.Interval, Vital.Status) ~ cat,
                         data = clinic.mom.df)
pval.tot <- 1 - pchisq(log.rank.all$chisq, length(log.rank.all$n) - 1)

#This log rank post hoc compares pairwise groups
statkm <- pairwise_survdiff(Surv(Correct.Interval, Vital.Status) ~ cat, 
                            data = clinic.mom.df,
                            p.adjust.method = "bonferroni")

#Creating a stat object DF for ggplot
stat.km.df <- data.frame(group1 = rep(NA, 6),
                         group2 = rep(NA, 6), 
                         p.adj = rep(NA, 6),
                         p.adj.signif = rep(NA, 6))

stat.km.df$group1 = c(rep("First Moment 0-25 percentile", 3),
                      rep("First Moment 25-50 percentile", 2),
                      rep("First Moment 50-75 percentile", 1))
stat.km.df$group2 = c("First Moment 25-50 percentile",
                      "First Moment 50-75 percentile",
                      "First Moment 75-100 percentile",
                      "First Moment 50-75 percentile",
                      "First Moment 75-100 percentile",
                      "First Moment 75-100 percentile")

stat.km.df$p.adj[1:3] <- statkm$p.value[, "First Moment 0-25 percentile"]
stat.km.df$p.adj[4:5] <- na.omit(statkm$p.value[, "First Moment 25-50 percentile"])
stat.km.df$p.adj[6] <- na.omit(statkm$p.value[, "First Moment 50-75 percentile"])
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
                    censor.size = 2) + 
  geom_vline(aes(xintercept = median, color = strata), alpha = .5, data = med.surv) +
  scale_y_continuous(expand = c(0, .03), breaks = c(0, .2, .4, .6, .8, 1)) + 
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07", "#9932CC")) + 
  theme_classic() + 
  theme(legend.justification=c(1,1), legend.position=c(1,.88),
        plot.title = element_text(hjust = 0.5, size = 18), 
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        axis.title = element_text(size = 18)) + 
  labs(title = "KM Curves for 0D Feature Curve Scaled Moment 1 Groups (Log Transformed)", x = "Survival Time (days)", 
       y = "Survival Probability", color = "Zero Feature Curve\nScaled Moment 1 Quartiles") + 
  guides(colour = guide_legend(override.aes = list(shape = 15)))
kmcurve


kmcurve.stat <- kmcurve +
  stat_pvalue_manual(stat.km.df, label = "p.adj",
                     y.position = c(.9, .7, .6, .3, .15),
                     size = 4, hide.ns = TRUE,
                     bracket.nudge.y = .1, tip.length = 0) +
  annotate("text", x = 140, y = .1, size = 4,
           label = paste("p = ", format(signif(pval.tot, 2), scientific = FALSE)))

kmcurve.stat


#####DIAGNOSIS BOX####
library(glmnet)


<<<<<<< HEAD

clinic.mom.df.lim <- clinic.mom.df %>% select(-c("Num.Prev.Canc.Diag", "HGB", "Path", "Histo", 
                                                 "OS.Length", "Progression.Status", "PFS.Length", 
                                                 "PFS.Failure.Type", "Local.Status", "Local.Length", "Lobar.Status", 
                                                 "Nodal.Status", "Distant.Status", "First.Met", "Update.Time", 
                                                 "Followup.Time", "Alive.Time"))

clinic.mom.df.lim <- clinic.mom.df.lim %>% na.omit()



#Surv Variables
y <- clinic.mom.df.lim[, c("Correct.Interval", "Vital.Status")]
y.surv <- Surv(time = y$Correct.Interval, event = y$Vital.Status)

#Predictor Vars
x <- clinic.mom.df.lim[, c("Sex", "Ethnicity", "Age.Diag", "BMI", 
                       "Years.Smoke", "KPS", "CCI", "Prev.Canc.Diag", "Prior.Surg", 
                       "Post.Rad.Chemo", "Pre.Rad.Chemo", "Rad.Intent", "Tot.Dose", 
                       "CT.Size", "T.Stage", "N.Stage", "M.Stage", "Overall.Stage", 
                       "raw.1st.dim0.mean", "raw.2nd.dim0", "raw.3rd.dim0", "raw.4th.dim0", 
                       "cent.2nd.dim0.variance", "cent.3rd.dim0", "cent.4th.dim0", "stand.1st.dim0", 
                       "stand.2nd.dim0", "stand.3rd.dim0", "stand.4th.dim0", "standcent.3rd.dim0.skew", 
                       "standcent.4th.dim0.kurt", "raw.1st.dim1.mean", "raw.2nd.dim1", 
                       "raw.3rd.dim1", "raw.4th.dim1", "cent.2nd.dim1.variance", "cent.3rd.dim1", 
                       "cent.4th.dim1", "stand.1st.dim1", "stand.2nd.dim1", "stand.3rd.dim1", 
                       "stand.4th.dim1", "standcent.3rd.dim1.skew", "standcent.4th.dim1.kurt", 
                       "raw.1st.dim2.mean", "raw.2nd.dim2", "raw.3rd.dim2", "raw.4th.dim2", 
                       "cent.2nd.dim2.variance", "cent.3rd.dim2", "cent.4th.dim2", "stand.1st.dim2", 
                       "stand.2nd.dim2", "stand.3rd.dim2", "stand.4th.dim2", "standcent.3rd.dim2.skew", 
                       "standcent.4th.dim2.kurt", "raw.1st.dimtot.mean", "raw.2nd.dimtot", 
                       "raw.3rd.dimtot", "raw.4th.dimtot", "cent.2nd.dimtot.variance", 
                       "cent.3rd.dimtot", "cent.4th.dimtot", "stand.1st.dimtot", "stand.2nd.dimtot", 
                       "stand.3rd.dimtot", "stand.4th.dimtot", "standcent.3rd.dimtot.skew", 
                       "standcent.4th.dimtot.kurt")]

fit <- glmnet(x, y.surv, family = "cox")

coef(fit, s = .05)

cvfit <- cv.glmnet(as.matrix(x), y.surv, family = "cox", type.measure = "C")


coxph(y.surv ~ x.mat)

x.mat <- as.matrix.data.frame(x)

#####DIAGNOSIS BOX####





=======
#####DIAGNOSIS BOX####
>>>>>>> 72af4f36942260c26257abb83c725a1433290221
