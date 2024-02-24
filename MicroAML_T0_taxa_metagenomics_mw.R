### Laure Bindels & Sarah Pötgens
### MicroAML metagenomics dataset
### This script concerns the T0 dataset

rm(list=ls())


##########################################################
###            IMPORT OF THE MAPPING FILE              ###
##########################################################

mapping = read.csv("mapping_file_MicroAML.csv", header = TRUE, sep = ",", row.names = 1)
mapping = mapping[,5:12]  #to select for specific variables
mapping = mapping[ order(row.names(mapping)), ]
a <- which(rownames(mapping) == 'aml109bis')
b <- which(rownames(mapping) == 'aml113bis')
c <- which(rownames(mapping) == 'aml510bis')
mapping <- mapping[-c(a,b,c),] # check line positions to remove
mapping_T0 = subset(mapping, time == 'T0')

## check the structure
str(mapping)
str(mapping_T0)


##########################################################
###            IMPORT OF THE TAXA TABLE                ###
##########################################################

## load the taxa table after taxa renaming (script-rename-metaphlan-output.R)
## Need rearrangement from the output of the bioinformatic pipeline
# open the file in openoffice for transposition, removal of the NCBI ID and renaming of the first column
# csv-delimited txt file, with samples as rows and taxa as columns
# percentage values

relabund_p = read.table("metaphlan_output_MicroAML_taxa_renamed.csv", header = TRUE, sep = ",")
relabund_p = as.data.frame(relabund_p) 
relabund_p = relabund_p[order(relabund_p$sample_ID), ]

## check the structure
str(relabund_p)
# check also visually, is the order correct? 

database = cbind(mapping, relabund_p)
a = database$sample_ID
b = database$patient
compare = cbind.data.frame(a, b) # comparison ok, validation of the merging

database_T0 = subset(database, time == 'T0')
relabund_p_T0 = database_T0[,-c(1:9)]
relabund_p_T0 = relabund_p_T0[ ,which(colSums(relabund_p_T0) > 0)] # remove the taxa with only 0 value, SP double checked the function in OpenOffice Calc
#92 taxa removed


#########################################################
###               FILTER LIST FOR TAXA                 ###
##########################################################

## establish the filter list to remove taxa with a mean below than 0.01% in the full dataset (T0)
data = relabund_p_T0
percent = 0.01
keep.taxa = which(colMeans(data) > percent)
length(keep.taxa)
#[1] 425
data.f1 = data[,keep.taxa]

# remove the taxa that are present in less than 25% group AML and CT
data.AML = data.f1[1:30,]
data.CT = data.f1 [31:60,]
presence_p = 25 
number_of_samples = 30 #TBU
threshold = presence_p/100 * number_of_samples
keep.taxa2 = which(colSums(data.AML !=0) > threshold)
keep.taxa3 = which(colSums(data.CT !=0) > threshold)
keep.taxa25p = c(keep.taxa2, keep.taxa3)
keep.taxa25p = unique(keep.taxa25p) #to keep only once each taxa column
length(keep.taxa25p)
#[1] 320

relabund_p_T0.filter.1 <- relabund_p_T0[,keep.taxa]
relabund_p_T0.filter.2 <- relabund_p_T0.filter.1[, keep.taxa25p]


##########################################################
###              UNIVARIATE ANALYSES                   ###
##########################################################

data3 = cbind(mapping_T0, relabund_p_T0.filter.2) 
relabund_p_ct = subset(data3, Group =="CT")
relabund_p_ct = relabund_p_ct[,9:ncol(data3)]
relabund_p_aml = subset(data3, Group =="AML")
relabund_p_aml = relabund_p_aml[,9:ncol(data3)]


## MW tests on percentages
pvalue_mw = vector('numeric') 
qvalue_mw = vector('numeric') 


for (i in 1:ncol(relabund_p_ct))
{ y1 = relabund_p_ct[,i]
y2 = relabund_p_aml[,i]
mw = wilcox.test(y1,y2)   # https://www.statmethods.net/stats/nonparametric.html
p = mw$p.value
pvalue_mw = append(pvalue_mw,p)
}

# qvalue Mann-Whitney baseline
# To correct the pvalue for the multiple testing, FDR correction. http://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
qvalue_mw = p.adjust(pvalue_mw, method = "fdr")	


## final stats table
data_stats = data.frame(cbind(pvalue_mw, qvalue_mw))                  
rownames(data_stats) = colnames(relabund_p_T0.filter.2)


write.table(data_stats, file = 'MicroAML_metagenomics_taxa_with_filters_stats_mw.txt', sep = "\t", quote = FALSE, append = FALSE)

