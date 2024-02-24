### Laure Bindels 
### MicroAML metagenomics dataset
### This script concerns the T0 dataset
### Taxa analyses at the genus level

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
# open the file in open office for transposition, removal of the NCBI ID and renaming of the first column
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
compare = cbind.data.frame(a, b) 
# comparison ok, validation of the merging

database_T0 = subset(database, time == 'T0')
relabund_p_T0 = database_T0[,-c(1:9)]
relabund_p_T0 = relabund_p_T0[ ,which(colSums(relabund_p_T0) > 0)] 
# remove the taxa with only 0 value, SP double checked the function in OpenOffice Calc
#92 taxa removed, ok


##########################################################
###           SELECTION OF THE GENUS LEVEL             ###
##########################################################

## select only the data at the genus level 
relabund_p_T0_genus =relabund_p_T0[,grep("g_",colnames(relabund_p_T0))]
# selection of 156 genera


#########################################################
###               FILTER LIST FOR TAXA                 ###
##########################################################

## establish the filter list to remove genus with a mean below than 0.01% in the full dataset (T0)
data = relabund_p_T0_genus
percent = 0.01
keep.taxa = which(colMeans(data) > percent)
length(keep.taxa)
#[1] 92
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
#[1] 78

relabund_p_T0.filter.1 <- relabund_p_T0_genus[,keep.taxa]
relabund_p_T0.filter.2 <- relabund_p_T0.filter.1[, keep.taxa25p]


##########################################################
###              CLR TRANSFORMATION                    ###
##########################################################

library(mixOmics)
library(vegan)

#relabund_p_T0 is in %
X = data.matrix(relabund_p_T0.filter.2)
# data.matrix to transform the data frame in numeric matrix, needed for the PCA


min(X[X > 0])
# [1] 0.00006
# I fix the offset as half the minimal value => 0.00003

offset.data = X + 0.00003
# offset.data = t(offset.data)

clr.data = data.matrix(logratio.transfo(offset.data, logratio = "CLR"))
# to log transform the data


##########################################################
###          PCA ANALYSIS AT THE GENUS LEVEL           ###
##########################################################

## select only the data at the genus level on the filtered data
clr.data_genus = clr.data

# to define the group
Group <- mapping_T0$Group
Group <- as.factor(Group)
Group <- relevel(Group, ref = "CT")
library(plyr)
count(mapping_T0, 'Group') # to know how many samples in each group
#Group freq
#1   AML   30
#2    CT   30

res.pca.1 <- pca(clr.data_genus, ncomp = 10, scale = TRUE, center = TRUE)
plot(res.pca.1) # to see the eighenvalues of each component, we can select 3 components

res.pca.1 <- pca(clr.data_genus, ncomp = 3, scale = TRUE, center = TRUE)

plotVar(res.pca.1)  # to see the correlation circle plot
plotVar(res.pca.1, cutoff = 0.5)  # to see the correlation circle plot

plot(sort(abs(res.pca.1$loadings$X[,1]), decreasing = TRUE), type = 'h', ylab = 'Loading vector, abs value in decreasing order')
# plot the abs value of the loading vector associated to component 1

plotIndiv(res.pca.1, ind.names = FALSE, group = Group, col.per.group = c("grey", "orange"), legend = TRUE, ellipse = TRUE,
          ellipse.level = 0.95, pch = 16, cex = 1, title = 'PCA at the genus level')
# see help "points" to select for a shape for each dot.
plotIndiv(res.pca.1, ind.names = TRUE, group = Group, col.per.group = c("grey", "orange"), legend = TRUE, ellipse = TRUE,
          ellipse.level = 0.95, title = 'PCA at the genus level')

#save image for publication
jpeg(filename = 'pca_genus_metagenomics_CLR.jpeg', width = 6*300, height = 5*300, res = 300)
print(plotIndiv(res.pca.1, ind.names = FALSE, group = Group, col.per.group = c("grey", "orange"), pch = c(18,20), style = "ggplot2", ellipse = TRUE, ellipse.level = 0.95, 
                title = "PCA at the genus level", size.title = rel(2), size.xlabel = rel(1.5), size.ylabel = rel(1.5), size.axis = rel(1), legend = FALSE ))
dev.off()


##########################################################
###              UNIVARIATE ANALYSES                   ###
##########################################################

## MW tests on clr transformed data, I start with the clr transformed data
clr.data_matrix = clr.data_genus[,grep("",colnames(relabund_p_T0.filter.2))] 
# an artifice to transform the clr object in matrix
data4 = cbind(mapping_T0, clr.data_matrix)
clr.data_ct = subset(data4, Group =="CT")
clr.data_ct = clr.data_ct[,9:ncol(data4)]
clr.data_aml = subset(data4, Group =="AML")
clr.data_aml = clr.data_aml[,9:ncol(data4)]

pvalue_mw_clr = vector('numeric') 
qvalue_mw_clr = vector('numeric') 

for (i in 1:ncol(clr.data_matrix))
{ y1 = clr.data_ct[,i]
y2 = clr.data_aml[,i]
mw = wilcox.test(y1,y2)   # https://www.statmethods.net/stats/nonparametric.html  
p = mw$p.value
pvalue_mw_clr = append(pvalue_mw_clr,p)
}

# qvalue Mann-Whitney baseline
# To correct the pvalue for the multiple testing, FDR correction. http://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
qvalue_mw_clr = p.adjust(pvalue_mw_clr, method = "fdr")	


## Aldex2-MW on count data

library("ALDEx2")
# browseVignettes("ALDEx2") # for help
# to define the conds vector if not already done

aldex_data = 100*(relabund_p_T0.filter.2)
aldex_data = t(aldex_data) # data should be with taxa in rows and samples in columns
aldex_data_rounded = round(aldex_data,0) # count should be rounded (aldex expects integers)
x.all <- aldex(aldex_data_rounded, Group, mc.samples=128, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE)

aldex.plot(x.all, type="MA", test="welch", xlab="Log-ratio abundance", ylab="Difference")

pvalue_aldex = x.all$wi.ep 
qvalue_aldex = x.all$wi.eBH 
aldex_effect = x.all$effect

## final stats table
data_stats = data.frame(cbind(pvalue_mw_clr, qvalue_mw_clr, aldex_effect, pvalue_aldex, qvalue_aldex))                  
rownames(data_stats) = colnames(relabund_p_T0.filter.2)


write.table(data_stats, file = 'MicroAML_metagenomics_genus_CLR-MW_aldex2.txt', sep = "\t", quote = FALSE, append = FALSE)
