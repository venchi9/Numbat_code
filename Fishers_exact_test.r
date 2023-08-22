#I will try to do Fisher's exact test on the clone percentages.
#This assumes that if there was no difference between the primary and LN, that the percentage would be 50%
#https://www.reneshbedre.com/blog/fisher-exact-test.html


fpga
until qrsh -l mem_requested=50G -pe smp 4; do sleep 2; done
conda activate R_4_2
R

library(tidyverse)
library(glue)
library(gtools)
library(ggplot2)
library(ggpubr)
library(stats)
library(rstatix)


data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/clones/data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/clones/figures"

#First read in the table
sample_name<-"HN200519A"
sample_name<-"HN120819A"
sample_name<-"HN021219A"
sample_name<-"HN230620A"

#HN200519A
x<-read.csv(file = str_glue("{data_dir}/{sample_name}_clone_percentages.csv"))
#Enter in the table manually because I'm stupid
df<-data.frame("Clone1" = c(1.71,0.94),  "Null"=c(50,50), row.names =c("Primary", "LN"))
fisher<-pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")
  group1  group2     n     p p.adj p.adj.signif
* <chr>   <chr>  <dbl> <dbl> <dbl> <chr>       
1 Primary LN      103.     1     1 ns         

df<-data.frame("Clone2" =c(3.09, 0.313),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")
  group1  group2     n     p p.adj p.adj.signif
* <chr>   <chr>  <dbl> <dbl> <dbl> <chr>       
1 Primary LN      103. 0.243 0.243 ns 

df<-data.frame("Clone3" =c(1.22, 71.2),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n        p    p.adj p.adj.signif
* <chr>   <chr>  <dbl>    <dbl>    <dbl> <chr>       
1 Primary LN      172. 5.08e-14 5.08e-14 ****   

df<-data.frame("Clone4" =c(31.94, 1.17),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n           p       p.adj p.adj.signif
* <chr>   <chr>  <dbl>       <dbl>       <dbl> <chr>       
1 Primary LN      133. 0.000000226 0.000000226 ****  

df<-data.frame("Clone5" =c(0.64, 26.31),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n         p     p.adj p.adj.signif
* <chr>   <chr>  <dbl>     <dbl>     <dbl> <chr>       
1 Primary LN      127. 0.0000038 0.0000038 ****     

df<-data.frame("Clone6"=c(61.36, 0.039),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")
  group1  group2     n        p    p.adj p.adj.signif
* <chr>   <chr>  <dbl>    <dbl>    <dbl> <chr>       
1 Primary LN      161. 8.01e-14 8.01e-14 ****        

#Need to divide the p-values by 6 before putting into the paper
#HN120819A
x<-read.csv(file = str_glue("{data_dir}/{sample_name}_clone_percentages.csv"))
#Clone1
df<-data.frame("Clone1"=c(1.04, 0.697),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n     p p.adj p.adj.signif
* <chr>   <chr>  <dbl> <dbl> <dbl> <chr>       
1 Primary LN      102.     1     1 ns   

#CLone2
df<-data.frame("Clone1"=c(35.16, 13.163),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n      p  p.adj p.adj.signif
* <chr>   <chr>  <dbl>  <dbl>  <dbl> <chr>       
1 Primary LN      148. 0.0125 0.0125 *     

#Clone3
df<-data.frame("Clone1"=c(63.795, 86.138),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n     p p.adj p.adj.signif
* <chr>   <chr>  <dbl> <dbl> <dbl> <chr>       
1 Primary LN      250.   0.3   0.3 ns      

#Divide by three before using

#HN021219A
sample_name<-"HN021219A"
x<-read.csv(file = str_glue("{data_dir}/{sample_name}_clone_percentages.csv"))

#Clone1
df<-data.frame("Clone1"=c(0.924, 1.97),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n     p p.adj p.adj.signif
* <chr>   <chr>  <dbl> <dbl> <dbl> <chr>       
1 Primary LN      103.     1     1 ns  

#Clone2
df<-data.frame("Clone1"=c(0.206, 0.383),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n     p p.adj p.adj.signif
* <chr>   <chr>  <dbl> <dbl> <dbl> <chr>       
1 Primary LN      101.     1     1 ns   

#Clone3
df<-data.frame("Clone1"=c(18.969, 65.024),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n        p    p.adj p.adj.signif
* <chr>   <chr>  <dbl>    <dbl>    <dbl> <chr>       
1 Primary LN      184. 0.000136 0.000136 ***       

#Clone4
df<-data.frame("Clone1"=c(0.618, 23.371),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n         p     p.adj p.adj.signif
* <chr>   <chr>  <dbl>     <dbl>     <dbl> <chr>       
1 Primary LN      124. 0.0000151 0.0000151 ****  

#Clone5
df<-data.frame("Clone1"=c(21.237, 6.951),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n      p  p.adj p.adj.signif
* <chr>   <chr>  <dbl>  <dbl>  <dbl> <chr>       
1 Primary LN      128. 0.0303 0.0303 *  

#Clone6
df<-data.frame("Clone1"=c(6.804, 2.29),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")
  group1  group2     n     p p.adj p.adj.signif
* <chr>   <chr>  <dbl> <dbl> <dbl> <chr>       
1 Primary LN      109. 0.165 0.165 ns          

#Clone7
df<-data.frame("Clone1"=c(51.340, 0),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n        p    p.adj p.adj.signif
* <chr>   <chr>  <dbl>    <dbl>    <dbl> <chr>       
1 Primary LN      151. 3.68e-12 3.68e-12 ****   

#Divide by 7 before using p values

#HN230620A
sample_name<-"HN230620A"
x<-read.csv(file = str_glue("{data_dir}/{sample_name}_clone_percentages.csv"))

#Clone1
df<-data.frame("Clone1"=c(16.915, 1.326),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n        p    p.adj p.adj.signif
* <chr>   <chr>  <dbl>    <dbl>    <dbl> <chr>       
1 Primary LN      118. 0.000448 0.000448 ***    

#Clone2
df<-data.frame("Clone1"=c(17.910, 0.515),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n       p   p.adj p.adj.signif
* <chr>   <chr>  <dbl>   <dbl>   <dbl> <chr>       
1 Primary LN      118. 0.00022 0.00022 ***   

#Clone3
df<-data.frame("Clone1"=c(10.447, 0.884),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n      p  p.adj p.adj.signif
* <chr>   <chr>  <dbl>  <dbl>  <dbl> <chr>       
1 Primary LN      111. 0.0106 0.0106 * 

#Clone4
df<-data.frame("Clone1"=c(6.965, 28.002),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n       p   p.adj p.adj.signif
* <chr>   <chr>  <dbl>   <dbl>   <dbl> <chr>       
1 Primary LN      135. 0.00258 0.00258 **      

#Clone5
df<-data.frame("Clone1"=c(9.950, 68.828),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n            p        p.adj p.adj.signif
* <chr>   <chr>  <dbl>        <dbl>        <dbl> <chr>       
1 Primary LN      179. 0.0000000913 0.0000000913 ****   

#Clone6
df<-data.frame("Clone1"=c(37.810, 0.44),  "Null"=c(50,50), row.names =c("Primary", "LN"))
pairwise_fisher_test(as.matrix(df), p.adjust.method = "fdr")

  group1  group2     n             p         p.adj p.adj.signif
* <chr>   <chr>  <dbl>         <dbl>         <dbl> <chr>       
1 Primary LN      138. 0.00000000105 0.00000000105 **** 

#Divide by 6 before putting in the paper

