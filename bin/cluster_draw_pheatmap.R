#! /usr/bin/env Rscript

#
# Draw a heatmap by clustering with Cluster 3.0 and then plotting with pheatmap
# Author: George Bell, Bioinformatics and Research Computing
#
# Version 1.0
# Version 1.1   -- Keep special characters when reading files; Handle cases without gene and sample clustering [12 Nov 2020 - GB]
# Version 1.1.1 -- Add message about 0s getting turned into NAs [4 Dec 2020 - GB]
# Version 1.2   -- Fix color range so 0 is in the center.  Users can still modify 'breaks'.  Add label fontsize variables [4 Dec 2020 - GB]
# Version 1.3   -- Add case so clustering can be skipped (with '-g 0 -e 0') and the heatmap will still be created. [21 Dec 2021 - GB] 
# Version 2.0   -- Use ComplexHeatmap: https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
#                  so that sample labels can be printed on top and feature labels on left [22 Aug 2025 - KdS] 

####################################  beginning of user-adjustable variables  ####################################

# Maximum number of genes to print (or names will be left out from the figure)
max.genes.to.print = 100
# Maximum number of samples to print (or names will be left out from the figure)
max.samples.to.print = 100
# Heatmap colors
neg.color = "blue"		# For samples with gene level less than gene average
zero.color = "white"		# For samples with gene level equal to gene average
pos.color = "red"		# For samples with gene level greater than gene average
#NA.color = "#B3B3B3"		# Color for missing data
NA.color = "yellow"		# Color for missing data
na_col = "yellow" 	# Color for missing data
relative.label.size = 1		# For heatmap sample and gene labels; default is 1
color.range.length = 255	# How many colors in our color range?
display_numbers = FALSE		# Should we print log2 ratios on heatmap?  (yes=>TRUE; no=>FALSE)
fontsize_row = 6		# Font size for row (gene) labels
fontsize_col = 10 		# Font size for column (sample) labels
breaks = seq(-3, 3, length.out = color.range.length+1)	# Set the max and min color range creates 256 break points
#breaks = seq(-3, 3, length.out = color.range.length)	# Set the max and min color range creates 255 colors You need n-1 colors for n breaks (since colors fill the intervals between breaks)
# breaks = NA			# Let the program choose the color range ('breaks = NA') or define the color range such as
legend.title = "Expression Level"	# Title that appears on the color scale legend

# Cluster 3.0 clustering details
# -l    => log-transform matrix
# -cg a => mean center (or '-cg m' => median center)
# -g 1  => cluster by gene (row) with Uncentered correlation
# -e 1  => cluster by sample (column) with Uncentered correlation
# -m a  => merge clusters by Pairwise average-linkage

cluster.bin = Sys.which("cluster")
if (nchar(cluster.bin) == 0) cluster.bin = "/usr/local/bin/cluster"
cluster.base.command = paste(cluster.bin, "-l -cg a -g 1 -e 1 -m a -f")	# followed by input file

####################################  end of user-adjustable variables  ####################################


# Make this backwards-compatible with "R --vanilla < foo.R" syntax
offset = 0
if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
  message("Running as R --vanilla ....")
  offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./cluster_draw_pheatmap.R" }

if (length(commandArgs()) < (7 - offset))
{
	message("\nCluster a matrix with Cluster 3.0 and then create the heatmap with R's pheatmap()\n")
	message(paste("USAGE:   ", this.script, "matrix outputPDF [sample_annotations feature_annotations]"))
	message(paste("Example: ", this.script, "Expr_values.txt MyHeatmap.pdf\n"))
	message("Notes:\n  The input matrix should not be log-transformed.")
	message("  Cluster 3.0 command: ", cluster.base.command)
	message("  Output heatmap file is 2 pages, with slightly different color schemes.")
	message("  To change clustering or display details, modify the top of this script.\n")
	q()
}

# Tab-delimited input file
input.file = commandArgs()[6 - offset]
# Name of output file
output.file = commandArgs()[7 - offset]

# Sample annotations (optional)
sample.anno.file = commandArgs()[8 - offset]

# Feature annotations (optional)
feature.anno.file = commandArgs()[9 - offset]

library(pheatmap)
library(ComplexHeatmap)


###
###  Cluster with Cluster 3.0
###

# -l    => log-transform matrix
# -cg a => mean center (or '-cg m' => median center)
# -g 1  => cluster by gene (row) with Uncentered correlation
# -e 1  => cluster by sample (column) with Uncentered correlation
# -m a  => merge clusters by Pairwise average-linkage

cluster.command = paste(cluster.base.command, input.file)
message("\nClustering matrix with Cluster 3.0, using the command")
message(cluster.command)
message("\nNote that log-transformation will cause 0 values to appear in the heatmap as NA (indicated by gray, the NA color).  If you want to avoid these NA entries, add pseudocount(s) to your matrix before running this code.\n")


system(cluster.command)
message("Clustering completed.  Now assembling results to create image ...\n")

if (! is.na(sample.anno.file) && sample.anno.file != "NA")
{
	message("Using ", sample.anno.file, " as the sample annotation file.")
	sample.anno = data.frame(read.delim(sample.anno.file, row.names=1, check.names=FALSE), check.names=FALSE)
	annotation_col = sample.anno
} else {
	message("No sample annotation file found.")
	annotation_col = NA
}

if (! is.na(feature.anno.file) && feature.anno.file != "NA")
{
	message("Using ", feature.anno.file, " as the feature annotation file.")
	feature.anno = data.frame(read.delim(feature.anno.file, row.names=1, check.names=FALSE), check.names=FALSE)
	annotation_row = feature.anno
} else {
	message("No feature annotation file found.")
	annotation_row = NA
}


######################  Code from https://rdrr.io/github/uc-bd2k/gimmR/man/importCdt.html (and slightly modified to work with Cluster 3.0)  ######################

importAtr <-
function(atrFile){
	if (!file.exists(atrFile)) {
		warning(paste("File not found: ", atrFile, ". Trying ", atrFile, ".atr", sep=""))
		atrFile <- paste(atrFile, "atr", sep = ".")
	}
	colc<-c("character","character","character","numeric")
	# Keep special characters when reading file (GB)
	# atr<-read.table(atrFile,colClasses=colc)
	atr<-read.table(atrFile,colClasses=colc, check.names=FALSE)
	height<-1-atr[,4]
	clusterNum<-dim(atr)[[1]]
	merge<-matrix(0,clusterNum,2)
	clusterSize<-integer(clusterNum)
	for(i in 1:clusterNum) {
		for( j in 1:2) {
			tempStr<-atr[i,(j+1)]
			positive<-substr(tempStr,1,4) == "NODE"
			id<-substr(tempStr,5,(nchar(tempStr)-1))
			if(positive) merge[i,j]<-as.integer(id)
			else  merge[i,j]<- -as.integer(id)
			if(merge[i,j] > 0) {
				clusterSize[i]<-clusterSize[i]+clusterSize[merge[i,j]]
			}else {
				clusterSize[i]<-clusterSize[i]+1
			}
		}
	}
	# Subtract 1 from all ARRAY IDs (because it's expected that there's no ARRY0X) ==> GWB
	merge = ifelse(merge<=0, merge-1, merge)	
	perm<-clusterNum
	allLeafs <- FALSE
	while(!allLeafs) {
		allLeafs <- all(perm < 0)
		if(!allLeafs) {
			index <- which(perm>=0)[1]
			newPerm <- as.vector(merge[perm[index], ])
			if(index > 1) newPerm <- c(perm[1:(index-1)], newPerm)
			if(index < length(perm)) newPerm <- c(newPerm, perm[(index+1):length(perm)])
			perm <- newPerm
		}
	}
	res <- list(merge=merge, height=height, order=-perm, labels = 1:length(perm), method = "average")
	class(res) <- "hclust"
	res
}

importGtr <-
function(gtrFile){
	if (!file.exists(gtrFile)) {
		warning(paste("File not found: ", gtrFile, ". Trying ", gtrFile, ".gtr", sep=""))
		gtrFile <- paste(gtrFile, "gtr", sep = ".")
	}
	colc<-c("character","character","character","numeric")
	# Keep special characters when reading file (GB)
	# gtr<-read.table(gtrFile,colClasses=colc)
	gtr<-read.table(gtrFile,colClasses=colc, check.names=FALSE)
	height<-1-gtr[,4]
	clusterNum<-dim(gtr)[[1]]
	merge<-matrix(0,clusterNum,2)
	clusterSize<-integer(clusterNum)
	for(i in 1:clusterNum) {
		for( j in 1:2) {
			tempStr<-gtr[i,(j+1)]
			positive<-substr(tempStr,1,4) == "NODE"
			id<-substr(tempStr,5,(nchar(tempStr)-1))
			if(positive) merge[i,j]<-as.integer(id)
			else  merge[i,j]<- -as.integer(id)
			if(merge[i,j] > 0) {
				clusterSize[i]<-clusterSize[i]+clusterSize[merge[i,j]]
			}else {
				clusterSize[i]<-clusterSize[i]+1
			}
		}
	}
	# Subtract 1 from all GENE IDs (because it's expected that there's no GENE0)
	merge = ifelse(merge<=0, merge-1, merge)
	perm<-clusterNum
	allLeafs <- FALSE
	while(!allLeafs) {
		allLeafs <- all(perm < 0)
		if(!allLeafs) {
			index <- which(perm>=0)[1]
			newPerm <- as.vector(merge[perm[index], ])
			if(index > 1) newPerm <- c(perm[1:(index-1)], newPerm)
			if(index < length(perm)) newPerm <- c(newPerm, perm[(index+1):length(perm)])
			perm <- newPerm
		}
	}
	res <- list(merge=merge, height=height, order=-perm, labels = 1:length(perm), method = "average")
	class(res) <- "hclust"
	res
}

importCdt <-
function(cdtFile) {
	if (!file.exists(cdtFile)) {
		warning(paste("File not found: ", cdtFile, ". Trying ", cdtFile, ".cdt", sep=""))
		cdtFile <- paste(cdtFile, "cdt", sep = ".")
	}
	# Keep special characters when reading file (GB)
	# cdt <- read.table(cdtFile, sep = "\t", stringsAsFactors=FALSE, header = TRUE, quote = "")
	cdt <- read.table(cdtFile, sep = "\t", stringsAsFactors=FALSE, header = TRUE, quote = "", check.names=FALSE)
	data.start.row <- which(cdt[,1] == "EWEIGHT") + 1
	data.start.col <- which(colnames(cdt) == "GWEIGHT") + 1
	if(length(data.start.row) < 1) { #there was no "EWEIGHT" row
		data.start.row <- 1
	}
	if(length(data.start.col) < 1) { #there was no "GWEIGHT" column
		if(colnames(cdt)[1] == "GID") data.start.col <- 4
		else data.start.col <- 3
	}
	if(colnames(cdt)[1] == "GID"){
		ID.col <- 2
		NAME.col <- 3
	} else {
		ID.col <- 1
		NAME.col <- 2
	}
	data.last.row <- dim(cdt)[1]
	data.last.col <- dim(cdt)[2]
	cdt <- cdt[data.start.row:data.last.row, c(ID.col,NAME.col,data.start.col:data.last.col)]
	for(i in 3:dim(cdt)[2]) cdt[, i]<-as.numeric(as.character(cdt[, i]))
	# Remove feature IDs from matrix and set to row names ==> GWB
	rownames(cdt) = cdt[,1]
	cdt = cdt[,3:ncol(cdt)]
	cdt
}

# CDT clustered matrix file
cdtFile = gsub(".txt", ".cdt", input.file)
# Gene (column) tree file
gtrFile = gsub(".txt", ".gtr", input.file)
# Sample (row) tree file
atrFile = gsub(".txt", ".atr", input.file)

message("Reading and processing Cluster 3.0 gene and sample tree files ...")

rawExpr = read.delim(input.file, row.names=1, check.names=FALSE)
num.genes = nrow(rawExpr)
num.samples = ncol(rawExpr)

if (file.exists(gtrFile))
{
	gclust = importGtr(gtrFile)
	gclust$labels = rownames(rawExpr)

	feature.pdf = gsub(".pdf", ".feature_dendrogram.pdf", output.file)
	pdf(feature.pdf, w=max(11, num.genes/8), h=8.5, useDingbats=FALSE)
	# plot(gclust, axes=FALSE, ylab="", main="Clustered features")
	# Try plot.dendrogram
	gclust.d = as.dendrogram(gclust)
	nodePar = list(lab.cex = 0.6, pch = c(NA, 19), cex = 0.7, col = "blue")
	plot(gclust.d, type = "rectangle", ylab = "", axes=F, main="Clustered features", nodePar = nodePar)
	foo = dev.off()
} else {
	message("No gene dendrogram file found, so assuming that genes weren't clustered.")
	gclust = FALSE
}

if (file.exists(atrFile))
{
	sclust = importAtr(atrFile)
	sclust$labels = colnames(rawExpr)

	sample.pdf = gsub(".pdf", ".sample_dendrogram.pdf", output.file)
	pdf(sample.pdf, w=max(11, num.samples/8), h=8.5, useDingbats=FALSE)
	# plot(sclust, axes=FALSE, ylab="", main="Clustered samples")
	# Try plot.dendrogram
	sclust.d = as.dendrogram(sclust)
	nodePar = list(lab.cex = 0.6, pch = c(NA, 19), cex = 0.7, col = "blue")
	plot(sclust.d, type = "rectangle", ylab = "", axes=F, main="Clustered samples", nodePar = nodePar)
	foo = dev.off()
} else {
	message("No sample dendrogram file found, so assuming that samples weren't clustered.")
	sclust = FALSE
}

# If we cluster
if ( file.exists(gtrFile) || file.exists(atrFile) )
{
	cdt = importCdt(cdtFile)
	cdt = cdt[rownames(rawExpr), colnames(rawExpr)]
} else {
	# We didn't cluster
	message("It looks like you aren't doing any clustering.  That's OK.  Reading nrm file ....")
	
	# Get nrm filename by replacing .txt with .nrm
	nrmFile = gsub(".txt", ".nrm", input.file)
	
	cdt = importCdt(nrmFile)
	cdt = cdt[rownames(rawExpr), colnames(rawExpr)]
}

show_rownames = TRUE
show_colnames = TRUE
if (num.genes > max.genes.to.print)
{
	message("Since number of genes (", num.genes, ") exceeds the display max (", max.genes.to.print, ") they will not be printed.")
	show_rownames = FALSE
}
if (num.samples > max.samples.to.print)
{
	message("Since number of samples (", num.samples, ") exceeds the display max (", max.samples.to.print, ") they will not be printed.")
	show_colnames = FALSE
}

# For ideas about drawing nice dendrograms, see
# http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning

pdf(output.file, w=11, h=8.5, useDingbats=FALSE)
title = paste(gsub(".txt", "", input.file), "heatmap")
# First heatmap with linear interpolation
ComplexHeatmap::pheatmap(as.matrix(cdt), cluster_rows=gclust, cluster_cols=sclust, show_rownames=show_rownames, show_colnames=show_colnames, col=colorRampPalette(c(neg.color, zero.color, pos.color))(color.range.length), main=title, na_col=NA.color, annotation_col=annotation_col, annotation_row=annotation_row, display_numbers=display_numbers, breaks=breaks, fontsize_row=fontsize_row, fontsize_col=fontsize_col, row_names_side="left", row_dend_side="right", column_names_side="top", column_dend_side="bottom", heatmap_legend_param = list(title = legend.title))
# Second heatmap with spline interpolation
ComplexHeatmap::pheatmap(as.matrix(cdt), cluster_rows=gclust, cluster_cols=sclust, show_rownames=show_rownames, show_colnames=show_colnames, col=colorRampPalette(c(neg.color, zero.color, pos.color),interpolate="spline")(color.range.length), main=title, na_col=NA.color, annotation_col=annotation_col, annotation_row=annotation_row, display_numbers=display_numbers, breaks=breaks, fontsize_row=fontsize_row, fontsize_col=fontsize_col, row_names_side="left", row_dend_side="right", column_names_side="top", column_dend_side="bottom", heatmap_legend_param = list(title = legend.title))

foo = dev.off()

message(paste("\nAll done!  See final figures of\n  complete heatmap (", output.file, ")", sep=""))
if (file.exists(gtrFile))
{
	message(paste("  gene tree (", feature.pdf, ")", sep=""))
}
if (file.exists(atrFile))
{
	message(paste("  sample tree (", sample.pdf, ")", sep=""))
}
message()
