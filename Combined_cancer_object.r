#Make sure we have a good combined cancer only combined object.  This is going to be quite painful.

#Data is here:
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Freemux_MT"

a<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Freemux_MT/RZ726_LN2_freemux_MT.rds")
b<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Freemux_MT/pri34_freemux_MT.rds")
c<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Freemux_MT/RZ726_pri_freemux_MT.rds")
d<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Freemux_MT/RZ726_LN_freemux_MT.rds")
e<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Freemux_MT/pri234_freemux_MT.rds")
f<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Freemux_MT/LN23_rep2_freemux_MT.rds")
g<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Freemux_MT/LN23_rep1_freemux_MT.rds")
h<-readRDS(file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/Freemux_MT/LN14_freemux_MT.rds")

#Put in location
a$Location<-"LN"
a$MajoritySinglet_Individual_Assignment<-"HN021219A"
b$Location<-"Primary"
c$Location<-"Primary"
d$Location<-"LN"
e$Location<-"Primary"
f$Location<-"LN"
g$Location<-"LN"

#Filter out the patients we don't need.

Idents(a)<-"MajoritySinglet_Individual_Assignment"
a<-subset(a, idents = c("HN021219A"))

Idents(b)<-"MajoritySinglet_Individual_Assignment"
b<-subset(b, idents = c("HN120819A", "HN200519A"))

Idents(c)<-"MajoritySinglet_Individual_Assignment"
c<-subset(c, idents = c("HN021219A", "HN230620A"))

Idents(d)<-"MajoritySinglet_Individual_Assignment"
d<-subset(d, idents = c("HN230620A"))

Idents(e)<-"MajoritySinglet_Individual_Assignment"
e<-subset(e, idents = c("HN120819A", "HN200519A"))

Idents(f)<-"MajoritySinglet_Individual_Assignment"
f<-subset(f, idents = c("HN120819A"))

Idents(g)<-"MajoritySinglet_Individual_Assignment"
g<-subset(g, idents = c("HN120819A"))


#LN14
Idents(h)<-"MajoritySinglet_Individual_Assignment"
h<-subset(h, idents =c("HN200519A"))

#Now perform a "Merge"

combined<-merge(a, y=c(b,c,d,e,f,g,h), project = "Combined")

saveRDS(combined, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/combined/merged_object_paper_notransform.rds")

#Now work out what seurat has done with the cell ids...
Idents(combined)<-"MajoritySinglet_Individual_Assignment"
test<-subset(combined, idents = "HN200519A")
AAACCCAAGAAATTGC-1_2

test<-subset(combined, idents = "HN120819A")
AAACCCAAGCATTGTC-1_2

test<-subset(combined, idents = "HN021219A")
AAACCCAAGCGTTCAT-1_1

test<-subset(combined, idents = "HN230620A")
AAACCCAAGGGTGGGA-1_3

#Now put in the curated names from the individual objects
data_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data"
fig_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/clones/figures"
list_dir<-"/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/reactomeGeneLists"

sample_name<-"HN200519A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))
Idents(x)<-"DEGcluster"
x<-subset(x, idents = c("Lymph_expanding", "Primary_4", "Primary_6"))
#Collapse into Primary and LN
x$combined<-ifelse(((x$DEGcluster == "Primary_4") | (x$DEGcluster == "Primary_6")), "Primary", "LN")
#Extract the cell ids
new<-x$combined
#Create data frame
df<-data.frame(new)
#Paste the suffix
rownames(df)<-gsub('.{4}$', '-1', rownames(df))
pt1<-df

#Add to invidiual groups
htest<-AddMetaData(object = h, metadata = pt1)
worked

btest<-AddMetaData(object = b, metadata = pt1)

etest<-AddMetaData(object = e, metadata = pt1)

#Worked!!!

#Hn120819A
sample_name<-"HN120819A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))
Idents(x)<-"DEGcluster"
x$combined<-ifelse(((x$DEGcluster =="Lymph_2") | (x$DEGcluster =="Lymph_3")), "LN", as.character(x$DEGcluster))
x$combined<-ifelse(((x$DEGcluster =="primary_2") | (x$DEGcluster =="primary_3")), "Primary", as.character(x$combined))
Idents(x)<-"combined"
x<-subset(x, idents = c("Primary", "LN"))
#Extract the cell ids
new1<-x$combined
#Create data frame
df<-data.frame(new1)
rownames(df)<-gsub('.{4}$', '-1', rownames(df))
pt2<-df

#Put into individual objects
btest<-AddMetaData(object = btest, metadata = pt2)
etest<-AddMetaData(object = etest, metadata = pt2)
ftest<-AddMetaData(object = f, metadata = pt2)
gtest<-AddMetaData(object = g, metadata = pt2)

#HN021219A
sample_name<-"HN021219A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))
Idents(x)<-"DEGcluster"
x$combined<-ifelse(((x$DEGcluster =="Lymph_3") | (x$DEGcluster =="Lymph_4")), "LN", as.character(x$DEGcluster))
x$combined<-ifelse(((x$DEGcluster =="primary")), "Primary", as.character(x$combined))
Idents(x)<-"combined"
x<-subset(x, idents = c("Primary", "LN"))
#Extract the cell ids
new2<-x$combined
#Create data frame
df<-data.frame(new2)
rownames(df)<-gsub('.{2}$', '-1-1', rownames(df))
pt3<-df

#Add to objects
ctest<-AddMetaData(object = c, metadata = pt3)
atest<-AddMetaData(object = a, metadata = pt3)

#HN230620A
sample_name<-"HN230620A"
x<-readRDS(file = str_glue("{data_dir}/{sample_name}/{sample_name}_cancer_object_reclustered.rds"))
Idents(x)<-"DEGcluster"
x$combined<-ifelse(((x$DEGcluster =="Lymph_node")), "LN", as.character(x$DEGcluster))
x$combined<-ifelse(((x$DEGcluster =="primary")), "Primary", as.character(x$combined))
Idents(x)<-"combined"
x<-subset(x, idents = c("Primary", "LN"))
#Extract the cell ids
new3<-x$combined
#Create data frame
df<-data.frame(new3)
rownames(df)<-gsub('.{2}$', '-1-3', rownames(df))
pt4<-df

#Add into the objects
ctest<-AddMetaData(object = ctest, metadata = pt4)
dtest<-AddMetaData(object = d, metadata = pt4)

#Okay, that has all worked.
#This next part is going to SUCK
#atest - only has one new column
Idents(atest)<-"new2"
atest<-subset(atest, idents =c("LN"))
atest$final<-atest$new2

#btest
#Has new1 and new2
btest$final<-ifelse(((btest$new == "Primary") | (btest$new1 =="Primary")), "Primary", "Other")
Idents(btest)<-"final"
btest<-subset(btest,idents =c("Primary"))

#ctest has new 2 and new3
ctest$final<-ifelse(((ctest$new2 == "Primary") | (ctest$new3 =="Primary")), "Primary", "Other")
Idents(ctest)<-"final"
ctest<-subset(ctest,idents =c("Primary"))

#dtest, only has "new3"
Idents(dtest)<-"new3"
dtest<-subset(dtest,idents =c("LN"))
dtest$final<-dtest$new3

#etest - has new and new1
etest$final<-ifelse(((etest$new == "Primary") | (etest$new1 =="Primary")), "Primary", "Other")
Idents(etest)<-"final"
etest<-subset(etest,idents =c("Primary"))

#ftest just new1
Idents(ftest)<-"new1"
ftest<-subset(ftest,idents =c("LN"))
ftest$final<-ftest$new1

#gtest just new1
Idents(gtest)<-"new1"
gtest<-subset(gtest,idents =c("LN"))
gtest$final<-gtest$new1

#htest new and new2
htest$final<-ifelse(((htest$new == "LN") | (htest$new2 =="LN")), "LN", "Other")
Idents(htest)<-"final"
htest<-subset(htest,idents =c("LN"))

#All worked, now merge and then sctransform

merged<-merge(atest, y =c(btest,ctest,dtest,etest,ftest,gtest,htest), project = "combined")
saveRDS(merged, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/combined/merged_object_paper_notransform.rds")

merged <- SCTransform(merged, vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(merged, file = "/directflow/SCCGGroupShare/projects/venchi/HN/HN_2022_new_demultiplexed_data/Data/combined/merged_object_paper_SCT.rds")



