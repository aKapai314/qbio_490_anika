# set working directory to analysis data
setwd("/Users/anikakapai/desktop/qbio_490_anika/week4_clinical_anika/analysis_data")

# install packages needed
if(!require("survival")){
  install.packages("survival")
}
if(!require("survminer")){
  install.packages("survminer")
}
if(!require("ggplot2")){
  install.packages("ggplot2")
}
# load packages
library(BiocManager)
library(TCGAbiolinks)
library(survival)
library(survminer)
library(ggplot2)

# read in clinical data
clinical <- read.csv("/Users/anikakapai/desktop/qbio_490_anika/week4_clinical_anika/analysis_data/brca_clinical_data.csv")

# query data
clinical_query <- GDCquery(project = "TCGA-BRCA", data.category = "Clinical", file.type = "xml")
#prepare drug and radiation data
clinical.drug <- GDCprepare_clinic(query = clinical_query, clinical.info = "drug")
clinical.rad <- GDCprepare_clinic(query = clinical_query, clinical.info = "radiation")

# first variable: histological_type = 0 NAs
sum(is.na(clinical$histological_type)) 

# second variable: radiation_dosage = 139 NAs
sum(is.na(clinical$lymph_node_examined_count)) 

#  Activity #1:  boxplot that compares distribution of lymph nodes with different histological subtypes
jpeg("/Users/anikakapai/desktop/qbio_490_anika/week4_clinical_anika/boxplot_histological_type.jpeg")
boxplot_hist = ggboxplot(clinical, x = "histological_type", y = "lymph_node_examined_count") +
  theme(axis.text.x = element_text(angle = 90))

boxplot_hist

# Activity #2:

# create survival time column by days to death or days to last follow up
clinical$survival_time <- ifelse(is.na(clinical$days_to_death), clinical$days_to_last_followup, clinical$days_to_death)
#create death event column
clinical$death_event <- ifelse(clinical$vital_status == "Dead", TRUE, FALSE)

# Initialize a 'survival' object for histology
surv_object_histology <- Surv(time = clinical$survival_time,
                              event = clinical$death_event)

# Create a fit object for histology
histology_fit <- surv_fit( surv_object_histology ~ clinical$histological_type,
                           data = clinical )

survplot_histology = ggsurvplot(histology_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")

jpeg("/Users/anikakapai/desktop/qbio_490_anika/week4_clinical_anika/KM_histological_subtype.jpeg")
KM_plot_histology = survplot_histology$plot + 
  theme_bw() + # visual
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

KM_plot_histology


# Looks at the distribution of lymph nodes
summary(clinical$lymph_node_examined_count)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.00    3.00    9.00   10.45   15.00   44.00     139

# stratify patients by the # of lymph nodes they had examined using nested ifelse structure
# categories: 0-9, 10-19, 20-29, 30-39, 40+
clinical$ln_category <- ifelse(clinical$lymph_node_examined_count < 10, "0-9", ifelse(clinical$lymph_node_examined_count < 20, "10-19", ifelse(clinical$lymph_node_examined_count < 30, "20-29", ifelse(clinical$lymph_node_examined_count < 40, "30-39", "40+"))))

# Initialize a 'survival' object for lymph nodes
surv_object_ln <- Surv(time = clinical$survival_time,
                       event = clinical$death_event)

# Create a fit object for lymph nodes
ln_fit <- surv_fit( surv_object_ln ~ clinical$ln_category, data = clinical )

survplot_ln = ggsurvplot(ln_fit, pval=TRUE, ggtheme = theme(plot.margin = unit(c(1,1,1,1), "cm")), legend = "right")

jpeg("/Users/anikakapai/desktop/qbio_490_anika/week4_clinical_anika/KM_ln.jpeg")
KM_plot_ln = survplot_ln$plot + 
  theme_bw() + # visuals
  theme(axis.title = element_text(size=20), 
        axis.text = element_text(size=16),
        legend.title = element_text(size=14),
        legend.text = element_text(size=12))

KM_plot_ln
