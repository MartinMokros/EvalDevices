library(readr)
library(dplyr)
library(Metrics)
library(ggplot2)
library(purrr)
library(ggpubr)
library(tidyr)
library(broom)

###### Setting main folder
folder_path = "D:/Devices_Bench/R_Results/"

#############################################################################################
############################### Importing reference data ####################################
#############################################################################################

setwd(paste0(folder_path,"Ref"))
temp = list.files(pattern="*.csv")
Ref_files = lapply(temp, read_csv)
Ref_files = lapply(Ref_files, select, tree, diameter, TLS_new, HMLS, ipad, multicam)
Ref_files = lapply(Ref_files, rename, DBH_ref = diameter,
                   HMLS_ID = HMLS,
                   TLS_ID = TLS_new,
                   ipad_ID = ipad,
                   MC_ID = multicam)

#Example of plot A
head(Ref_files[1])


#Reference data only with trees over 7 cm of DBH
Trees_over_7 = lapply(Ref_files, function(x) filter(x, DBH_ref > 0.07))

print(paste("On all eight plots",sum(sapply(Ref_files, nrow)),"trees were measured.", 
            "And",sum(sapply(Trees_over_7, nrow)),"trees have DBH >=7cm"))

#############################################################################################
###################### Processing Terrestrial laser scanning data ###########################
#############################################################################################

#Importing data from TLS and joining them with referene data
setwd(paste0(folder_path,"TLS"))
temp = list.files(pattern="*.csv")
TLS_files = lapply(temp, read_csv, skip = 1, col_names= c("TLS_ID", "X", "Y", "TLS_DBH"))
TLS_files = lapply(TLS_files, select, TLS_ID, TLS_DBH)
TLS_join = map2(Ref_files, TLS_files, inner_join, by = "TLS_ID")

#Example of plot A before and after join with reference data
head(TLS_files[1])
head(TLS_join[1])

#Calculate number of detected trees
Tree_detected_TLS = 0
for(i in 1:8) {
  TLS_PLOT = as.data.frame(TLS_files[i])
  Tree_detected_TLS = Tree_detected_TLS+ nrow(TLS_PLOT)
}
#Number of all detected trees
print(paste("Numebr of detected trees across all plots by TLS is:", sum(Tree_detected_TLS)))

#Calculating falsly detected trees
False_TLS_Trees = c()
for(i in 1:8) {
  TLS_Plot = as.data.frame(TLS_files[i])
  TLS_join_Plot = as.data.frame(TLS_join[i])
  False_TLS_Tree_new = nrow(TLS_Plot) - (nrow(TLS_join_Plot) - sum(is.na(TLS_join_Plot$TLS_ID)))
  False_TLS_Trees = c(False_TLS_Trees, False_TLS_Tree_new)
}

print(paste("Numebr of falsly detected trees across all plots by TLS is:", sum(False_TLS_Trees)))

#Filtering out reference trees under 7 cm in DBH
TLS_join = lapply(TLS_join, function(x) filter(x, DBH_ref >= 0.07))

#Pre-Calculating error
TLS_join = lapply(TLS_join, function(x) mutate(x, Error = TLS_DBH - DBH_ref))

#Excluding matched trees with error over 100%
TLS_join = lapply(TLS_join, function(x) mutate(x, twenty=  abs(Error / (DBH_ref/100))))

TLS_join = lapply(TLS_join, function(x) filter(x, twenty <= 100))

#Calculating tree detection rate (TDR)
TDR_TLS = c()
nrow(as.data.frame(TLS_join[2]))
for(i in 1:8) {
  detected_trees = nrow(as.data.frame(TLS_join[i]))
  TLS_TDR_new = ((detected_trees)/(nrow(as.data.frame(Trees_over_7[i]))/100))
  TDR_TLS = c(TDR_TLS, TLS_TDR_new)
}
plots = c("A" ,"B","C","D","E","F","G","H")

print(paste("Percentage of tree detection rate for plot", plots, "is:", as.vector(round(TDR_TLS))))

#Calculating RMSE, rRMSE, Bias, rBias
RMSE_TLS = c()
rRMSE_TLS = c()
Bias_TLS = c()
rBias_TLS = c()

for(i in 1:8) {
  TLS = as.data.frame(TLS_join[i])
  RMSE_new = sqrt(mean(TLS$Error^2, na.rm=TRUE))
  rRMSE_new = (RMSE_new / mean(TLS$DBH_ref)) * 100
  Bias_new = mean(TLS$Error, na.rm=TRUE)
  rBias_new = abs(Bias_new / mean(TLS$DBH_ref)) * 100
  
  RMSE_TLS = c(RMSE_TLS, RMSE_new)
  rRMSE_TLS = c(rRMSE_TLS, rRMSE_new)
  Bias_TLS = c(Bias_TLS, Bias_new)
  rBias_TLS = c(rBias_TLS, rBias_new)
}

TLS_Eval = data.frame(RMSE_TLS*100, rRMSE_TLS, Bias_TLS*100, rBias_TLS)
row.names(TLS_Eval) = c("A", "B", "C", "D", "E", "F", "G", "H")
colnames(TLS_Eval) = c("RMSE (cm)", "rRMSE (%)", "Bias (cm)", "rBias (%)")
(TLS_Eval = round(TLS_Eval, digits = 1))

#Visualization by scatter plot - Reference DBH versus Estimated DBH
TLS_scatter = list()

for(i in 1:8) {
  TLS_data =as.data.frame(TLS_join[i])
  TLS_scatter_new =  ggplot(TLS_data, aes(DBH_ref, TLS_DBH))+
    geom_point(size=3)+
    geom_smooth(method = lm, se = FALSE)+
    theme_classic()+
    xlab("Reference DBH (m)")+ylab("TLS DBH (m)")
  
  TLS_scatter[[i]] = TLS_scatter_new
}

ggarrange(plotlist=TLS_scatter,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          ncol = 4, nrow = 2)

#Visualizations of errors by box plots
Plot_labels = c("A", "B", "C", "D", "E", "F", "G", "H")
TLS_errors = lapply(TLS_join, function(x) select(x, DBH_ref, TLS_DBH))
TLS_errors = lapply(TLS_errors, function(x) rename(x, est_DBH = TLS_DBH))

TLS_errors_ = list()
for(i in 1:8){
  TLS_error = as.data.frame(TLS_errors[i])
  TLS_error = mutate(TLS_error, Error = est_DBH - DBH_ref, Device = "TLS", Plot = Plot_labels[i])
  TLS_errors_[[i]] = TLS_error
}

TLS_errors_join = do.call(rbind, TLS_errors_)

ggplot(TLS_errors_join, aes(x=Plot, y=Error*100))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.25)+
  theme_bw()+
  xlab("Plots")+ylab("TLS errors (cm)")+
  geom_hline(aes(yintercept = 0))

#Calculating ovela statistics across all plots
RMSE_TLS_all = sqrt(mean(TLS_errors_join$Error^2, na.rm=TRUE))
rRMSE_TLS_all = (RMSE_TLS_all / mean(TLS_errors_join$DBH_ref)) * 100
Bias_TLS_all = mean(TLS_errors_join$Error, na.rm=TRUE)
rBias_TLS_all = abs(Bias_TLS_all / mean(TLS_errors_join$DBH_ref)) * 100
TDR_TLS_all = nrow(TLS_errors_join)/sum(sapply(Trees_over_7, nrow))*100

TLS_Eval_overal = data.frame(RMSE_TLS_all*100, rRMSE_TLS_all, Bias_TLS_all*100, rBias_TLS_all, TDR_TLS_all)
row.names(TLS_Eval_overal) = c("All plots")
colnames(TLS_Eval_overal) = c("RMSE (cm)", "rRMSE (%)", "Bias (cm)", "rBias (%)", "TDR (%)")
(TLS_Eval_overal = round(TLS_Eval_overal, digits = 1))

#############################################################################################
###################### Processing hand-held personal laser scanning data ####################
#############################################################################################

#Importing data from HMLS and joining them with referene data
setwd(paste0(folder_path,"HMLS"))
temp = list.files(pattern="*.csv")
HMLS_files = lapply(temp, read_csv, skip = 1, col_names= c("HMLS_ID", "X", "Y", "HMLS_DBH"))
HMLS_files = lapply(HMLS_files, select, HMLS_ID, HMLS_DBH)
HMLS_join = map2(Ref_files, HMLS_files, inner_join, by = "HMLS_ID")

#Example of plot A before and after join with reference data
head(HMLS_files[1])
head(HMLS_join[1])

#Calculate number of detected trees
Tree_detected_HMLS = 0
for(i in 1:8) {
  HMLS_PLOT = as.data.frame(HMLS_files[i])
  Tree_detected_HMLS = Tree_detected_HMLS+ nrow(HMLS_PLOT)
}
#Number of all detected trees
print(paste("Numebr of detected trees across all plots by PLShh is:", sum(Tree_detected_HMLS)))

#Calculating falsly detected trees
False_HMLS_Trees = c()
for(i in 1:8) {
  HMLS_Plot = as.data.frame(HMLS_files[i])
  HMLS_join_Plot = as.data.frame(HMLS_join[i])
  False_HMLS_Tree_new = nrow(HMLS_Plot) - (nrow(HMLS_join_Plot) - sum(is.na(HMLS_join_Plot$HMLS_ID)))
  False_HMLS_Trees = c(False_HMLS_Trees, False_HMLS_Tree_new)
}

print(paste("Numebr of falsly detected trees across all plots by PLShh is:", sum(False_HMLS_Trees)))

#Filtering out reference trees under 7 cm in DBH
HMLS_join = lapply(HMLS_join, function(x) filter(x, DBH_ref >= 0.07))

#Pre-Calculating error
HMLS_join = lapply(HMLS_join, function(x) mutate(x, Error = HMLS_DBH - DBH_ref))

#Excluding matched trees with error over 100%
HMLS_join = lapply(HMLS_join, function(x) mutate(x, twenty=  abs(Error / (DBH_ref/100))))

HMLS_join = lapply(HMLS_join, function(x) filter(x, twenty <= 100))

#Calculating tree detection rate (TDR)
TDR_HMLS = c()

for(i in 1:8) {
  detected_trees = nrow(as.data.frame(HMLS_join[i]))
  HMLS_TDR_new = ((detected_trees)/(nrow(as.data.frame(Trees_over_7[i]))/100))
  TDR_HMLS = c(TDR_HMLS, HMLS_TDR_new)
}
plots = c("A" ,"B","C","D","E","F","G","H")

print(paste("Percentage of tree detection rate for plot", plots, "is:", as.vector(round(TDR_HMLS))))

#Calculating RMSE, rRMSE, Bias, rBias
RMSE_HMLS = c()
rRMSE_HMLS = c()
Bias_HMLS = c()
rBias_HMLS = c()

for(i in 1:8) {
  HMLS = as.data.frame(HMLS_join[i])
  RMSE_new = sqrt(mean(HMLS$Error^2, na.rm=TRUE))
  rRMSE_new = (RMSE_new / mean(HMLS$DBH_ref)) * 100
  Bias_new = mean(HMLS$Error, na.rm=TRUE)
  rBias_new = abs(Bias_new / mean(HMLS$DBH_ref)) * 100
  
  RMSE_HMLS = c(RMSE_HMLS, RMSE_new)
  rRMSE_HMLS = c(rRMSE_HMLS, rRMSE_new)
  Bias_HMLS = c(Bias_HMLS, Bias_new)
  rBias_HMLS = c(rBias_HMLS, rBias_new)
}

HMLS_Eval = data.frame(RMSE_HMLS*100, rRMSE_HMLS, Bias_HMLS*100, rBias_HMLS)
row.names(HMLS_Eval) = c("A", "B", "C", "D", "E", "F", "G", "H")
colnames(HMLS_Eval) = c("RMSE (cm)", "rRMSE (%)", "Bias (cm)", "rBias (%)")
(HMLS_Eval = round(HMLS_Eval, digits = 1))

#Visualization by scatter plot - Reference DBH versus Estimated DBH
HMLS_scatter = list()

for(i in 1:8) {
  HMLS_data =as.data.frame(HMLS_join[i])
  HMLS_scatter_new =  ggplot(HMLS_data, aes(DBH_ref, HMLS_DBH))+
    geom_point(size=3)+
    geom_smooth(method = lm, se = FALSE)+
    theme_classic()+
    xlab("Reference DBH (m)")+ylab("PLShh DBH (m)")
  
  HMLS_scatter[[i]] = HMLS_scatter_new
}

ggarrange(plotlist=HMLS_scatter,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          ncol = 4, nrow = 2)

#Visualizations of errors by box plots
Plot_labels = c("A", "B", "C", "D", "E", "F", "G", "H")
HMLS_errors = lapply(HMLS_join, function(x) select(x, DBH_ref, HMLS_DBH))
HMLS_errors = lapply(HMLS_errors, function(x) rename(x, est_DBH = HMLS_DBH))

HMLS_errors_ = list()
for(i in 1:8){
  HMLS_error = as.data.frame(HMLS_errors[i])
  HMLS_error = mutate(HMLS_error, Error = est_DBH - DBH_ref, Device = "HMLS", Plot = Plot_labels[i])
  HMLS_errors_[[i]] = HMLS_error
}

HMLS_errors_join = do.call(rbind, HMLS_errors_)

ggplot(HMLS_errors_join, aes(x=Plot, y=Error*100))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.25)+
  theme_bw()+
  xlab("Plots")+ylab("PLShh errors (cm)")+
  geom_hline(aes(yintercept = 0))

#Calculating ovela statistics across all plots
RMSE_HMLS_all = sqrt(mean(HMLS_errors_join$Error^2, na.rm=TRUE))
rRMSE_HMLS_all = (RMSE_HMLS_all / mean(HMLS_errors_join$DBH_ref)) * 100
Bias_HMLS_all = mean(HMLS_errors_join$Error, na.rm=TRUE)
rBias_HMLS_all = abs(Bias_HMLS_all / mean(HMLS_errors_join$DBH_ref)) * 100
TDR_HMLS_all = nrow(HMLS_errors_join)/sum(sapply(Trees_over_7, nrow))*100

HMLS_Eval_overal = data.frame(RMSE_HMLS_all*100, rRMSE_HMLS_all, Bias_HMLS_all*100, rBias_HMLS_all, TDR_HMLS_all)
row.names(HMLS_Eval_overal) = c("All plots")
colnames(HMLS_Eval_overal) = c("RMSE (cm)", "rRMSE (%)", "Bias (cm)", "rBias (%)", "TDR (%)")
(HMLS_Eval_overal = round(HMLS_Eval_overal, digits = 1))

#############################################################################################
############################## Processing iPad Pro data #####################################
#############################################################################################

#Importing data from ipad and joining them with referene data
setwd(paste0(folder_path,"ipad"))
temp = list.files(pattern="*.csv")
ipad_files = lapply(temp, read_csv, skip = 1, col_names= c("ipad_ID", "X", "Y", "ipad_DBH"))
ipad_files = lapply(ipad_files, select, ipad_ID, ipad_DBH)
ipad_join = map2(Ref_files, ipad_files, inner_join, by = "ipad_ID")

#Example of plot A before and after join with reference data
head(ipad_files[1])
head(ipad_join[1])

#Calculate number of detected trees
Tree_detected_ipad = 0
for(i in 1:8) {
  ipad_PLOT = as.data.frame(ipad_files[i])
  Tree_detected_ipad = Tree_detected_ipad+ nrow(ipad_PLOT)
}
#Number of all detected trees
print(paste("Numebr of detected trees across all plots by ipad is:", sum(Tree_detected_ipad)))

#Calculating falsly detected trees
False_ipad_Trees = c()
for(i in 1:8) {
  ipad_Plot = as.data.frame(ipad_files[i])
  ipad_join_Plot = as.data.frame(ipad_join[i])
  False_ipad_Tree_new = nrow(ipad_Plot) - (nrow(ipad_join_Plot) - sum(is.na(ipad_join_Plot$ipad_ID)))
  False_ipad_Trees = c(False_ipad_Trees, False_ipad_Tree_new)
}

print(paste("Numebr of falsly detected trees across all plots by ipad is:", sum(False_ipad_Trees)))

#Filtering out reference trees under 7 cm in DBH
ipad_join = lapply(ipad_join, function(x) filter(x, DBH_ref >= 0.07))

#Pre-Calculating error
ipad_join = lapply(ipad_join, function(x) mutate(x, Error = ipad_DBH - DBH_ref))

#Excluding matched trees with error over 100%
ipad_join = lapply(ipad_join, function(x) mutate(x, twenty=  abs(Error / (DBH_ref/100))))

ipad_join = lapply(ipad_join, function(x) filter(x, twenty <= 100))

#Calculating tree detection rate (TDR)
TDR_ipad = c()

for(i in 1:8) {
  detected_trees = nrow(as.data.frame(ipad_join[i]))
  ipad_TDR_new = ((detected_trees)/(nrow(as.data.frame(Trees_over_7[i]))/100))
  TDR_ipad = c(TDR_ipad, ipad_TDR_new)
}
plots = c("A" ,"B","C","D","E","F","G","H")

print(paste("Percentage of tree detection rate for plot", plots, "is:", as.vector(round(TDR_ipad))))

#Calculating RMSE, rRMSE, Bias, rBias
RMSE_ipad = c()
rRMSE_ipad = c()
Bias_ipad = c()
rBias_ipad = c()

for(i in 1:8) {
  ipad = as.data.frame(ipad_join[i])
  RMSE_new = sqrt(mean(ipad$Error^2, na.rm=TRUE))
  rRMSE_new = (RMSE_new / mean(ipad$DBH_ref)) * 100
  Bias_new = mean(ipad$Error, na.rm=TRUE)
  rBias_new = abs(Bias_new / mean(ipad$DBH_ref)) * 100
  
  RMSE_ipad = c(RMSE_ipad, RMSE_new)
  rRMSE_ipad = c(rRMSE_ipad, rRMSE_new)
  Bias_ipad = c(Bias_ipad, Bias_new)
  rBias_ipad = c(rBias_ipad, rBias_new)
}

ipad_Eval = data.frame(RMSE_ipad*100, rRMSE_ipad, Bias_ipad*100, rBias_ipad)
row.names(ipad_Eval) = c("A", "B", "C", "D", "E", "F", "G", "H")
colnames(ipad_Eval) = c("RMSE (cm)", "rRMSE (%)", "Bias (cm)", "rBias (%)")
(ipad_Eval = round(ipad_Eval, digits = 1))

#Visualization by scatter plot - Reference DBH versus Estimated DBH
ipad_scatter = list()

for(i in 1:8) {
  ipad_data =as.data.frame(ipad_join[i])
  ipad_scatter_new =  ggplot(ipad_data, aes(DBH_ref, ipad_DBH))+
    geom_point(size=3)+
    geom_smooth(method = lm, se = FALSE)+
    theme_classic()+
    xlab("Reference DBH (m)")+ylab("iPad DBH (m)")
  
  ipad_scatter[[i]] = ipad_scatter_new
}

ggarrange(plotlist=ipad_scatter,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          ncol = 4, nrow = 2)

#Visualizations of errors by box plots
Plot_labels = c("A", "B", "C", "D", "E", "F", "G", "H")
ipad_errors = lapply(ipad_join, function(x) select(x, DBH_ref, ipad_DBH))
ipad_errors = lapply(ipad_errors, function(x) rename(x, est_DBH = ipad_DBH))

ipad_errors_ = list()
for(i in 1:8){
  ipad_error = as.data.frame(ipad_errors[i])
  ipad_error = mutate(ipad_error, Error = est_DBH - DBH_ref, Device = "ipad", Plot = Plot_labels[i])
  ipad_errors_[[i]] = ipad_error
}

ipad_errors_join = do.call(rbind, ipad_errors_)

ggplot(ipad_errors_join, aes(x=Plot, y=Error*100))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.25)+
  theme_bw()+
  xlab("Plots")+ylab("iPad errors (cm)")+
  geom_hline(aes(yintercept = 0))

#Calculating ovela statistics across all plots
RMSE_ipad_all = sqrt(mean(ipad_errors_join$Error^2, na.rm=TRUE))
rRMSE_ipad_all = (RMSE_ipad_all / mean(ipad_errors_join$DBH_ref)) * 100
Bias_ipad_all = mean(ipad_errors_join$Error, na.rm=TRUE)
rBias_ipad_all = abs(Bias_ipad_all / mean(ipad_errors_join$DBH_ref)) * 100
TDR_ipad_all = nrow(ipad_errors_join)/sum(sapply(Trees_over_7, nrow))*100

ipad_Eval_overal = data.frame(RMSE_ipad_all*100, rRMSE_ipad_all, Bias_ipad_all*100, rBias_ipad_all, TDR_ipad_all)
row.names(ipad_Eval_overal) = c("All plots")
colnames(ipad_Eval_overal) = c("RMSE (cm)", "rRMSE (%)", "Bias (cm)", "rBias (%)", "TDR (%)")
(ipad_Eval_overal = round(ipad_Eval_overal, digits = 1))

#############################################################################################
###################### Processing Multi-camera system data ##################################
#############################################################################################

#Importing data from MC and joining them with referene data
setwd(paste0(folder_path,"MC"))
temp = list.files(pattern="*.csv")
MC_files = lapply(temp, read_csv, skip = 1, col_names= c("MC_ID", "X", "Y", "MC_DBH"))
MC_files = lapply(MC_files, select, MC_ID, MC_DBH)
MC_join = map2(Ref_files, MC_files, inner_join, by = "MC_ID")

#Example of plot A before and after join with reference data
head(MC_files[1])
head(MC_join[1])

#Calculate number of detected trees
Tree_detected_MC = 0
for(i in 1:8) {
  MC_PLOT = as.data.frame(MC_files[i])
  Tree_detected_MC = Tree_detected_MC+ nrow(MC_PLOT)
}
#Number of all detected trees
print(paste("Numebr of detected trees across all plots by MC is:", sum(Tree_detected_MC)))

#Calculating falsly detected trees
False_MC_Trees = c()
for(i in 1:8) {
  MC_Plot = as.data.frame(MC_files[i])
  MC_join_Plot = as.data.frame(MC_join[i])
  False_MC_Tree_new = nrow(MC_Plot) - (nrow(MC_join_Plot) - sum(is.na(MC_join_Plot$MC_ID)))
  False_MC_Trees = c(False_MC_Trees, False_MC_Tree_new)
}

print(paste("Numebr of falsly detected trees across all plots by MC is:", sum(False_MC_Trees)))

#Filtering out reference trees under 7 cm in DBH
MC_join = lapply(MC_join, function(x) filter(x, DBH_ref >= 0.07))

#Pre-Calculating error
MC_join = lapply(MC_join, function(x) mutate(x, Error = MC_DBH - DBH_ref))

#Excluding matched trees with error over 100%
MC_join = lapply(MC_join, function(x) mutate(x, twenty=  abs(Error / (DBH_ref/100))))

MC_join = lapply(MC_join, function(x) filter(x, twenty <= 100))

#Calculating tree detection rate (TDR)
TDR_MC = c()

for(i in 1:8) {
  detected_trees = nrow(as.data.frame(MC_join[i]))
  MC_TDR_new = ((detected_trees)/(nrow(as.data.frame(Trees_over_7[i]))/100))
  TDR_MC = c(TDR_MC, MC_TDR_new)
}
plots = c("A" ,"B","C","D","E","F","G","H")

print(paste("Percentage of tree detection rate for plot", plots, "is:", as.vector(round(TDR_MC))))

#Calculating RMSE, rRMSE, Bias, rBias
RMSE_MC = c()
rRMSE_MC = c()
Bias_MC = c()
rBias_MC = c()

for(i in 1:8) {
  MC = as.data.frame(MC_join[i])
  RMSE_new = sqrt(mean(MC$Error^2, na.rm=TRUE))
  rRMSE_new = (RMSE_new / mean(MC$DBH_ref)) * 100
  Bias_new = mean(MC$Error, na.rm=TRUE)
  rBias_new = abs(Bias_new / mean(MC$DBH_ref)) * 100
  
  RMSE_MC = c(RMSE_MC, RMSE_new)
  rRMSE_MC = c(rRMSE_MC, rRMSE_new)
  Bias_MC = c(Bias_MC, Bias_new)
  rBias_MC = c(rBias_MC, rBias_new)
}

MC_Eval = data.frame(RMSE_MC*100, rRMSE_MC, Bias_MC*100, rBias_MC)
row.names(MC_Eval) = c("A", "B", "C", "D", "E", "F", "G", "H")
colnames(MC_Eval) = c("RMSE (cm)", "rRMSE (%)", "Bias (cm)", "rBias (%)")
(MC_Eval = round(MC_Eval, digits = 1))


#Visualization by scatter plot - Reference DBH versus Estimated DBH
MC_scatter = list()

for(i in 1:8) {
  MC_data =as.data.frame(MC_join[i])
  MC_scatter_new =  ggplot(MC_data, aes(DBH_ref, MC_DBH))+
    geom_point(size=3)+
    geom_smooth(method = lm, se = FALSE)+
    theme_classic()+
    xlab("Reference DBH (m)")+ylab("MultiCam DBH (m)")
  
  MC_scatter[[i]] = MC_scatter_new
}

ggarrange(plotlist=MC_scatter,
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
          ncol = 4, nrow = 2)

#Visualizations of errors by box plots
Plot_labels = c("A", "B", "C", "D", "E", "F", "G", "H")
MC_errors = lapply(MC_join, function(x) select(x, DBH_ref, MC_DBH))
MC_errors = lapply(MC_errors, function(x) rename(x, est_DBH = MC_DBH))

MC_errors_ = list()
for(i in 1:8){
  MC_error = as.data.frame(MC_errors[i])
  MC_error = mutate(MC_error, Error = est_DBH - DBH_ref, Device = "MC", Plot = Plot_labels[i])
  MC_errors_[[i]] = MC_error
}

MC_errors_join = do.call(rbind, MC_errors_)

ggplot(MC_errors_join, aes(x=Plot, y=Error*100))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.25)+
  theme_bw()+
  xlab("Plots")+ylab("MultiCam errors (cm)")+
  geom_hline(aes(yintercept = 0))

#Calculating ovela statistics across all plots
RMSE_MC_all = sqrt(mean(MC_errors_join$Error^2, na.rm=TRUE))
rRMSE_MC_all = (RMSE_MC_all / mean(MC_errors_join$DBH_ref)) * 100
Bias_MC_all = mean(MC_errors_join$Error, na.rm=TRUE)
rBias_MC_all = abs(Bias_MC_all / mean(MC_errors_join$DBH_ref)) * 100
TDR_MC_all = nrow(MC_errors_join)/sum(sapply(Trees_over_7, nrow))*100

MC_Eval_overal = data.frame(RMSE_MC_all*100, rRMSE_MC_all, Bias_MC_all*100, rBias_MC_all, TDR_MC_all)
row.names(MC_Eval_overal) = c("All plots")
colnames(MC_Eval_overal) = c("RMSE (cm)", "rRMSE (%)", "Bias (cm)", "rBias (%)", "TDR (%)")
(MC_Eval_overal = round(MC_Eval_overal, digits = 1))

#############################################################################################
########################### Comparing results from devices ##################################
#############################################################################################

#Merging statistics of all devices for each plot to tables
dir.create(file.path(folder_path, "Stats_Figures"))
setwd(paste0(folder_path, "Stats_Figures"))

rRMSE = data.frame(rRMSE_TLS, rRMSE_HMLS, rRMSE_ipad, rRMSE_MC)
row.names(rRMSE) = c("A", "B", "C", "D", "E", "F", "G", "H")
colnames(rRMSE) = c("TLS", "PLShh", "iPad", "MultiCam")
(rrRMSE = round(rRMSE, digits = 1))
write.csv(x=rrRMSE, file="rRMSE.csv")

RMSE = data.frame(RMSE_TLS, RMSE_HMLS, RMSE_ipad, RMSE_MC)
row.names(RMSE) = c("A", "B", "C", "D", "E","F",  "G", "H")
colnames(RMSE) = c("TLS", "PLShh", "iPad", "MultiCam")
(roundRMSE = round(RMSE*100, digits = 1))
write.csv(x=roundRMSE, file="RMSE.csv")

Bias = data.frame(Bias_TLS, Bias_HMLS, Bias_ipad, Bias_MC)
row.names(Bias) = c("A", "B", "C", "D", "E", "F", "G", "H")
colnames(Bias) = c("TLS", "PLShh", "iPad", "MultiCam")
(roundBias = round(Bias*100, digits = 2))
write.csv(x=roundBias, file="Bias.csv")

rBias = data.frame(rBias_TLS, rBias_HMLS, rBias_ipad, rBias_MC)
row.names(rBias) = c("A", "B", "C", "D", "E", "F", "G", "H")
colnames(rBias) = c("TLS", "PLShh", "iPad", "MultiCam")
(roundBias = round(rBias, digits = 2))
write.csv(x=roundBias, file="rBias.csv")

TDR = data.frame(TDR_TLS, TDR_HMLS, TDR_ipad, TDR_MC)
row.names(TDR) = c("A", "B", "C", "D", "E", "F", "G", "H")
colnames(TDR) = c("TLS", "PLShh", "iPad", "MultiCam")
(roundTDR = round(TDR, digits = 1))
write.csv(x=roundTDR, file="TDR.csv")

FTDT = data.frame(False_TLS_Trees, False_HMLS_Trees, 
                  False_ipad_Trees, False_MC_Trees)
row.names(FTDT) = c("A", "B", "C", "D", "E", "F", "G", "H")
colnames(FTDT) = c("TLS", "PLShh", "iPad", "MultiCam")
FTDT
write.csv(x=FTDT, file="FTDT.csv")

#Merging overal statistics of all devices

RMSE_all = c(RMSE_TLS_all, RMSE_HMLS_all, RMSE_ipad_all, RMSE_MC_all)
rRMSE_all = c(rRMSE_TLS_all, rRMSE_HMLS_all, rRMSE_ipad_all, rRMSE_MC_all)
Bias_all = c(Bias_TLS_all, Bias_HMLS_all, Bias_ipad_all, Bias_MC_all)
rBias_all = c(rBias_TLS_all, rBias_HMLS_all, rBias_ipad_all, rBias_MC_all)
TDR_all = c(TDR_TLS_all, TDR_HMLS_all, TDR_ipad_all, TDR_MC_all)
FTDT_all = c(sum(False_TLS_Trees), sum(False_HMLS_Trees), 
              sum(False_ipad_Trees), sum(False_MC_Trees))

stat_all = data.frame(RMSE_all*100, rRMSE_all, Bias_all*100, rBias_all, TDR_all, FTDT_all)
colnames(stat_all) = c("RMSE (cm)", "rRMSE (%)", "Bias (cm)", "rBias (%)", "TDR (%)", "FTD (n)")
row.names(stat_all) = c("TLS", "PLShh", "iPad", "MultiCam")
(round_stat_all = round(stat_all, digits = 2))
write.csv(x=round_stat_all, file="stat_all.csv")

#Calculating two way ANOVA and post-hoc Tukey test
Errors_join = rbind(TLS_errors_join, HMLS_errors_join, ipad_errors_join, MC_errors_join)


res.aov <- aov(Error ~ Device * Plot, data = Errors_join)
summary(res.aov)
(res.tuk <- TukeyHSD(res.aov))

write.csv( tidy(res.aov) , "ANOVA.csv" )
write.csv( tidy(res.tuk) , "TUKEY.csv" )

#Testing Bias agianst 0

TLS_Ttest = filter(Errors_join, Device == "TLS")
t.test(TLS_Ttest$Error,mu=0)

HMLS_Ttest = filter(Errors_join, Device == "HMLS")
t.test(HMLS_Ttest$Error,mu=0)

ipad_Ttest = filter(Errors_join, Device == "ipad")
t.test(ipad_Ttest$Error,mu=0)

MC_Ttest = filter(Errors_join, Device == "MC")
t.test(MC_Ttest$Error,mu=0)

############################# Creating Figure 5 #############################

rRMSE$Plots = row.names(rRMSE)
rRMSE = gather(rRMSE, "Device", "rRMSE", -Plots)
TDR = gather(TDR, "Device", "TDR")
TDR_rRmse = cbind(rRMSE, TDR[2])

TDR_rRmse$Device = factor(TDR_rRmse$Device, 
                          levels=c("TLS", "PLShh",
                                   "iPad", "MultiCam"))

png("Figure_5.png", width=3000, height=1500, res=300)
ggplot(TDR_rRmse, aes(x=Plots, y=TDR, group=Device))+
  geom_line(aes(color=Device), size=2)+
  geom_point(aes(color=Device), size=5)+
  theme_bw()+
  ylab("Tree detection rate (%)")+xlab("Test Sites")+
  scale_color_manual(values = 
                       c("#09015f", "#af0069", 
                         "#55b3b1", "#f6c065"), 
                     labels = c("TLS", expression(PLS['hh']),
                                "iPad", "MultiCam"))+
  theme(legend.title=element_blank(),
        legend.text = element_text(size=20),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        strip.text = element_text(size=20),
        legend.position = "bottom")
dev.off()

############################# Creating Figure 6 #############################

dat_text <- data.frame(
  Device = c("TLS", "HMLS", "ipad", "MC"))

TLS_lm = Errors_join %>% 
  filter(Device=="TLS")
TLS_lm = lm(DBH_ref ~ est_DBH , data = TLS_lm) 
TLS_lm = summary(TLS_lm)$r.squared

HMLS_lm = Errors_join %>% 
  filter(Device=="HMLS")
HMLS_lm = lm(DBH_ref ~ est_DBH , data = HMLS_lm) 
HMLS_lm = summary(HMLS_lm)$r.squared

ipad_lm = Errors_join %>% 
  filter(Device=="ipad")
ipad_lm = lm(DBH_ref ~ est_DBH , data = ipad_lm) 
ipad_lm = summary(ipad_lm)$r.squared

MC_lm = Errors_join %>% 
  filter(Device=="MC")
MC_lm = lm(DBH_ref ~ est_DBH , data = MC_lm) 
MC_lm = summary(MC_lm)$r.squared

Errors_join = rbind(TLS_errors_join, HMLS_errors_join, ipad_errors_join, MC_errors_join)
Errors_join$Device = factor(Errors_join$Device, levels = c("TLS","HMLS","ipad","MC"))

my_lab = as_labeller(c(TLS="TLS", HMLS="PLS[hh]", ipad="iPad", MC="MultiCam"), default = label_parsed)

png("Figure_6.png", width=3000, height=2000, res=300)
Scat = ggplot(Errors_join, aes(DBH_ref*100, est_DBH*100))+
  geom_point(size=3, alpha=0.75, aes(color=Device))+
  geom_smooth(method = lm, se = FALSE, color = "black")+
  theme_bw()+
  theme(legend.position = "none")+
  xlab("Reference DBH (cm)")+ylab("Estimated DBH (cm)")+
  scale_color_manual(values = 
                       c("#09015f", "#af0069", 
                         "#55b3b1", "#f6c065"))+
  theme(legend.title=element_blank(),
        legend.text = element_text(size=20),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        strip.text = element_text(size=18),
        strip.background = element_rect(colour="black",
                                        fill="white"))+
  facet_wrap(~Device, labeller = my_lab)
Scat + geom_text(
  data    = dat_text,
  mapping = aes(x = 16, y = 68, label = 
                  c(paste("r^2 ==", round(TLS_lm, 3)), paste("r^2 ==", round(HMLS_lm, 3)), 
                    paste("r^2 ==", round(ipad_lm, 3)),paste("r^2 ==", round(MC_lm, 3)))),
  parse=T, size=6)
dev.off()

############################# Creating Figure 7 #############################

Errors_join = rbind(TLS_errors_join, HMLS_errors_join, ipad_errors_join, MC_errors_join)

Errors_join$Device = factor(Errors_join$Device, levels=c("MC", 'ipad', 'HMLS','TLS'))

png("Figure_7.png", width=5000, height=1500, res=350)
ggplot(Errors_join, aes(x=Device, y=Error*100))+
  geom_boxplot(aes(fill=Device), outlier.size = 4,
               outlier.alpha = 0.75, outlier.shape = 21)+
  theme_bw()+
  xlab("")+ylab("Errors (cm)")+
  scale_x_discrete(labels=c("HMLS"=expression(PLS['hh']), 
                            "ipad" = "iPad", "MC" = "MultiCam"))+
  scale_y_continuous(breaks=seq(-100,100,10))+
  scale_fill_manual(values = 
                      c("#f6c065","#55b3b1","#af0069","#09015f"))+
  geom_hline(aes(yintercept = 0), color="red", size=1)+
  theme(legend.title=element_blank(),
        legend.text = element_text(size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=20),
        strip.text = element_text(size=20),
        legend.position = "none")+
  coord_flip()
dev.off()

############################# Creating Figure 8 #############################

Errors_join = mutate(Errors_join, rError = Error/(DBH_ref/100))

png("Figure_8.png", width=4000, height=2000, res=300)
ggplot(Errors_join, aes(x=Device, y=Error*100, fill = Plot))+
  #geom_point(pch=21, position = position_jitterdodge(), show.legend=FALSE)+
  geom_boxplot(size=0.5, outlier.shape = NA, alpha=0.6)+
  #geom_jitter(aes(color = Plot), width = 0.15, alpha=0.5, size=2)+
  theme_bw()+
  xlab("")+ylab("Errors (cm)")+
  #labs(fill="Plots")+
  scale_x_discrete(labels=c("HMLS"=expression(PLS['hh']),
                            "ipad" = "iPad", "MC" = "MultiCam"))+
  scale_y_continuous(breaks=seq(-1000,1000,5))+
  scale_fill_discrete(name="Plots")+
  geom_hline(aes(yintercept = 0), color="black", size=0.5)+
  theme(legend.position = "bottom",
        legend.text = element_text(size=20),
        legend.title = element_text(size=20),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size=20),
        strip.text = element_text(size=20))+
  scale_color_discrete(guide = FALSE)+
  guides(fill = guide_legend(nrow = 1))+
  #facet_wrap(~Device, scales="free", labeller = labeller(Device = device.labs))
  #coord_cartesian(ylim=c(-100,100))+
  coord_flip(ylim=c(-20,20))
dev.off()

### Calculating overal statistics for trees with DBH >10 and >20 cm 
############################# Creating Figure 9 #############################

Errors_join = rbind(TLS_errors_join, HMLS_errors_join, ipad_errors_join, MC_errors_join)
head(Errors_join)

Trees_over_20 = lapply(Ref_files, function(x) filter(x, DBH_ref >= 0.2))
Trees_over_10 = lapply(Ref_files, function(x) filter(x, DBH_ref >= 0.1))

print(paste("On all eight plots",sum(sapply(Ref_files, nrow)),"trees were measured.", 
            "And",sum(sapply(Trees_over_7, nrow)),"trees have DBH >7cm,",
            sum(sapply(Trees_over_10, nrow)), "trees have DBH >10cm,",
            sum(sapply(Trees_over_20, nrow)), "trees have DBH >20cm,"))


Errors_over_20cm = filter(Errors_join, DBH_ref>0.2)
Errors_over_20cm %>%
  group_by(Device) %>%
  summarise(n=n())
    
(Errors_over_20cm = Errors_over_20cm %>%
  group_by(Device) %>%
  summarise(RMSE = sqrt(mean(Error^2, na.rm=TRUE)), rRMSE = (RMSE / mean(DBH_ref)) * 100,
            Bias = mean(Error, na.rm=TRUE), 
            rBias = abs(Bias / mean(DBH_ref)) * 100,
            n=n(),
            TDR = (n/153)*100,
            DBH = 20))

Errors_over_10cm = filter(Errors_join, DBH_ref>0.1)
Errors_over_10cm %>%
  group_by(Device) %>%
  summarise(n=n())

(Errors_over_10cm = Errors_over_10cm %>%
    group_by(Device) %>%
    summarise(RMSE = sqrt(mean(Error^2, na.rm=TRUE)), rRMSE = (RMSE / mean(DBH_ref)) * 100,
              Bias = mean(Error, na.rm=TRUE), 
              rBias = abs(Bias / mean(DBH_ref)) * 100,
              n=n(),
              TDR = (n/229)*100,
              DBH = 10))

Errors_over_7cm = filter(Errors_join, DBH_ref>0.07)
Errors_over_7cm %>%
  group_by(Device) %>%
  summarise(n=n())

(Errors_over_7cm = Errors_over_7cm %>%
    group_by(Device) %>%
    summarise(RMSE = sqrt(mean(Error^2, na.rm=TRUE)), rRMSE = (RMSE / mean(DBH_ref)) * 100,
              Bias = mean(Error, na.rm=TRUE), 
              rBias = abs(Bias / mean(DBH_ref)) * 100,
              n=n(),
              TDR = (n/268)*100,
              DBH = 7))

Errors_ove = rbind(Errors_over_7cm, Errors_over_10cm, Errors_over_20cm)

Errors_ove$Device = factor(Errors_ove$Device, 
                                 levels=c("TLS", "HMLS",
                                          "ipad", "MC"))

png("Figure_9_a.png", width=3000, height=2000, res=300)
ggplot(Errors_ove, aes(x=Device, y=TDR, fill=factor(DBH),label =round(TDR)))+
  geom_col(position = "dodge", alpha=0.80, color="black")+
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 6)+
  theme_bw()+
  ylab("Tree detection rate (%)")+
  scale_fill_manual(values = 
                      c("#e8e8e8", "#bbbfca", 
                        "#495464"), 
                    labels = c(">7", ">10", ">20"))+
  scale_x_discrete(labels = c("TLS", expression(PLS['hh']),
                     "iPad", "MultiCam"))+
  theme(legend.title=element_blank(),
        legend.text = element_text(size=20),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        strip.text = element_text(size=20),
        axis.title.x = element_blank(),
        legend.position = "bottom")
dev.off()

png("Figure_9_b.png", width=3000, height=2000, res=300)
ggplot(Errors_ove, aes(x=Device, y=rRMSE, fill=factor(DBH),label =round(rRMSE, 1)))+
  geom_col(position = "dodge", alpha=0.80, color="black")+
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 6)+
  theme_bw()+
  ylab("Relative RMSE (%)")+
  scale_fill_manual(values = 
                      c("#e8e8e8", "#bbbfca", 
                        "#495464"), 
                    labels = c(">7", ">10", ">20"))+
  scale_x_discrete(labels = c("TLS", expression(PLS['hh']),
                              "iPad", "MultiCam"))+
  theme(legend.title=element_blank(),
        legend.text = element_text(size=20),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        strip.text = element_text(size=20),
        axis.title.x = element_blank(),
        legend.position = "bottom")
dev.off()

############################# Creating Figure 10 #############################

stat_all$Device = row.names(stat_all)
stat_all$Device = factor(stat_all$Device, 
                         levels=c("TLS", "PLShh",
                                  "iPad", "MultiCam"))


png("Figure_10.png", width=8000, height=5000, res=600)
ggplot(TDR_rRmse, aes(TDR, rRMSE, fill=Device))+
  stat_ellipse(aes(color=Device), type = "t", size=1.75, alpha=0.75, show.legend=FALSE)+
  geom_point(size=8, alpha=0.75, shape=21)+
  geom_point(data=stat_all, aes(`TDR (%)`, rRMSE, color=Device), shape=13, size = 12, show.legend=FALSE)+
  theme_bw()+
  xlab("Tree detection rate (%)")+ylab("Relative RMSE (%)")+
  scale_x_continuous(breaks=seq(0,100,10))+
  scale_y_continuous(breaks=seq(0,100,10))+
  scale_fill_manual(values = 
                      c("#09015f", "#af0069", 
                        "#55b3b1", "#f6c065"), 
                    labels = c("TLS", expression(PLS['hh']),
                               "iPad", "MultiCam"))+
  scale_color_manual(values = 
                       c("#09015f", "#af0069", 
                         "#55b3b1", "#f6c065"))+
  theme(legend.title=element_blank(),
        legend.text = element_text(size=22),
        axis.text = element_text(size=18),
        axis.title = element_text(size=22),
        strip.text = element_text(size=22),
        legend.position = c(0.75, 0.845))+
  guides(fill = guide_legend(nrow = 1))
#geom_text(aes(label=Plots), color="black", hjust=0, vjust=-1)
dev.off()

############################# Creating Figure S1 #############################
Errors_join = rbind(TLS_errors_join, HMLS_errors_join, 
                    ipad_errors_join, MC_errors_join)

Errors_join$Device = factor(Errors_join$Device, 
                            levels=c('TLS', 'HMLS','ipad',"MC"))

Plot_labels = c("Site A", "Site B", "Site C", "Site D",
                "Site E", "Site F", "Site G", "Site H")
names(Plot_labels) = c("A", "B", "C", "D", "E", "F", "G", "H")

png("Figure_S1.png", width=6000, height=2500, res=300)
ggplot(Errors_join, aes(DBH_ref*100, est_DBH*100))+
  geom_point(aes(fill = Device), size=5, alpha=0.4, shape=21)+
  geom_smooth(aes(color=Device),method = lm, se = FALSE, size=1.2, alpha=1)+
  scale_x_continuous(breaks = seq(0, 100, by = 10))+
  scale_y_continuous(breaks = seq(0, 100, by = 10))+
  #scale_shape_manual(values=c(3, 16, 17, 20))+
  scale_color_manual(values = 
                       c("#09015f","#af0069", "#55b3b1","#f6c065"), labels=c("HMLS"=expression(PLS['hh']),
                                                                             "ipad" = "iPad", "MC" = "MultiCam"))+
  scale_fill_manual(values = 
                      c("#09015f","#af0069", "#55b3b1","#f6c065"), labels=c("HMLS"=expression(PLS['hh']),
                                                                            "ipad" = "iPad", "MC" = "MultiCam"))+
  theme_bw()+
  theme(legend.position="bottom") +
  theme(legend.title=element_blank(),
        legend.text = element_text(size=24),
        axis.text = element_text(size=16),
        axis.title = element_text(size=20),
        strip.text = element_text(size=20))+
  #scale_color_manual(values=c("#003f5c", "#7a5195", "#ef5675", "#ffa600"))+
  xlab("Reference DBH (cm)")+ylab("Estimated DBH (cm)")+
  facet_wrap(~ Plot, nrow = 2, labeller = labeller(Plot = Plot_labels))
dev.off()
