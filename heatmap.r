#This File was used to generate the Heatmap for DOE Grant and for Publication of RPT4a RPT4b and PIPs Paper

#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }


#########################################################
### B) Reading in data and transform it into matrix format
#########################################################


#Modify Data Read-in
setwd("C:\\Users\\Workstation\\Documents\\Revolution\\Project0\\datasets")
data <- read.csv("test_csv_for_r_heatmap.csv", comment.char="#")
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-5 into a matrix
rownames(mat_data) <- rnames                  # assign row names 


#########################################################
### C) Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("black", "yellow", "orange", "red"))(n = 110)

# (optional) defines the color breaks manually for a "skewed" color transition
col_breaks = c(seq(0,0,length=1),  # for black
  seq(0.0001,0.0010,length=10),              # for yellow
  seq(0.0011,0.02573,length=50),			#for orange
  seq(0.02574,0.08577,length=50))              # for red

# creates a 5 x 5 inch image
pdf("heatmaps_in_r.pdf")

hclust2 <- function(x, method="centroid", ...)
  hclust(x, method=method, ...)
dist2 <- function(x, ...)
  as.dist(1-cor(t(x), method="pearson"))

heatmap.2(mat_data,
  #cellnote = mat_data,  # same data set for cell labels
  main = "LFQ MS/MS Native MG132 PAG1", # heat map title
  hclustfun = hclust2,
  distfun = dist2,
  notecol="black",      # change font color of cell labels to black
  #density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns off trace lines inside the heat map
  margins =c(8,22),     # widens margins around plot
  col=my_palette,       # use on color palette defined earlier 
  breaks=col_breaks,    # enable color transition at specified limits
  key = TRUE,
  keysize = 1.5,
  symkey = FALSE,
  key.title = "Color Key",
  dendrogram="row",     # only draw a row dendrogram
  Colv="NA")            # turn off column clustering

dev.off()               # close the PNG device
