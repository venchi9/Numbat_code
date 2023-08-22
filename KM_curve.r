#Do KM curves for Ali
# https://bioconnector.github.io/workshops/r-survival.html#cox_regression

#This is also helpful
#https://www.datacamp.com/tutorial/survival-analysis-R

#Actually used this one for my data:
#https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html#Kaplan-Meier_plots

until qrsh -l mem_requested=20G -pe smp 4; do sleep 2; done
conda activate r_4
cd /directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/figures
R

#Load directories
library(survival)
library(survminer)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggsurvfit)

#Read in file
res<-read_tsv(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Survival_stats.txt", col_names = TRUE)

#Do KM curve

surv_object<-Surv(time = res$Time_to_recurrence, event = res$Recurrence)
surv_object

survfit2(Surv(Time_to_recurrence, Recurrence) ~ 1, data = res) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Time to recurrence"
  ) + add_risktable()
  dev.off()


#Worked


#Now do death
survfit2(Surv(Time_to_death, Death) ~ 1, data = res) %>% 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Time to death"
  ) + add_risktable()
  dev.off()
