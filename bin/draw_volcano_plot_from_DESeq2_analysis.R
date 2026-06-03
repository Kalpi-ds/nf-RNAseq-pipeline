#!/usr/bin/env Rscript

###
###  Draw a MA plot for a DESeq2 analysis
###
###  George Bell - BaRC
###  Version 1.0: 2 Nov 2018
###  Input file: Output file from DESeq2's results() object
###    with gene symbols in the first column
###  Draw MA plot, coloring genes that meet threshold, 
###    and labeling up to num.genes.to.label (set to 100) genes.
###  Add geom_text_repel(max.overlaps=100) to prevent much skipping of labels for nearby points
###

# version = "Version 2.0 -- 2 November 2018"
# version = "Version 2.1 -- 8 February 2019"
# version = "Version 2.2 -- 4 June 2021"
# version = "Version 2.3 -- 8 December 2021"	# Fix bug that labeled too many genes [See lines 131-136]
version = "Version 2.4 -- 2 May 2022"	# Keep gene IDs as is Change print ordering of points so labeled points are on top 

suppressMessages(library(ggplot2))
library(ggrepel)

# Maximum number of genes to label with gene names
num.genes.to.label = 100
# Color for up genes
up.color = "#FFAAAA"
# Color for down genes
down.color = "#AAAAFF"
# Offset to add before calculating log2(baseMean)
counts.mean.offset = 1
# Y-axis limits
# y.axis.min = -10
# y.axis.max = 10

offset = 0
if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
	message("Running as R --vanilla ....")
	offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./draw_volcano_plot_from_DESeq_analysis.R" }

input.filename = commandArgs()[6 - offset]
sampleA.name = commandArgs()[7 - offset]
sampleB.name = commandArgs()[8 - offset]
fdr.threshold = as.numeric(commandArgs()[9 - offset])
log2FC.threshold = as.numeric(commandArgs()[10 - offset])
# gene.symbols.filename = commandArgs()[11 - offset]

if ( ! is.na(log2FC.threshold) ) {
	message(paste("\nUsing DESeq2 output file *", input.filename, "* for input", sep=""))
	message(paste("Group A name is *", sampleA.name, "*", sep=""))
	message(paste("Group B name is *", sampleB.name, "*", sep=""))
	message(paste("FDR threshold is *", fdr.threshold, "*", sep=""))
	message(paste("Log2 fold change threshold is *", log2FC.threshold, "*", sep=""))
} else {
	message("\nDraw a volcano plot from DESeq2 output,")
	message("  coloring points that meet FDR and fold change thresholds,")
	message("  and labeling top 100 points based on p-value.")
	message("USAGE: ", this.script, " inputFile Sample.A.name Sample.B.name FDR.threshold log2FC.threshold")
	message(paste(version, "\n"))
	quit()
}

# Read DESeq output file, and make the first column row names
# DESeq2.stats = read.delim(input.filename, row.names=1)
###  Modified to permit gene matrices with duplicated gene symbols [GB, 7 Sep 2018]
# Input gene matrix (no row names)
DESeq2.stats = read.delim(input.filename, header=T, sep="\t", check.names=F)
# Order matrix by decreasing baseMean, adding a numeric suffix to each symbol we've already seen.
gene.symbols.original = as.character(DESeq2.stats[,1])
# gene.symbols.unique = .Internal(make.names(gene.symbols.original, unique=TRUE))
gene.symbols.unique = make.unique(gene.symbols.original)
o = order(gene.symbols.original != gene.symbols.unique, DESeq2.stats$baseMean, decreasing=TRUE)
gene.symbols.unique[o] = make.unique(gene.symbols.unique[o])
# Make row names from unique gene symbols.
rownames(DESeq2.stats) = gene.symbols.unique
# Drop the gene symbol (first column) since we have them as row names
DESeq2.stats = DESeq2.stats[,-1]

# head(DESeq2.stats)

# Drop rows with a FDR of NA
DESeq2.stats = DESeq2.stats[! is.na(DESeq2.stats$padj),]

# Order genes by increasing p-values
DESeq2.stats = DESeq2.stats[order(DESeq2.stats$pvalue, decreasing=FALSE), ]

# For plot labels
comparison = paste("log2 ( ", sampleA.name, " / ", sampleB.name, " )", sep="")
figure.title = paste("Volcano plot:", sampleA.name, "vs.", sampleB.name)

# Get IDs
gene.symbols = rownames(DESeq2.stats)
message("Getting p-value from ", colnames(DESeq2.stats)[4])
p.value = DESeq2.stats$pvalue
# Get the FDR and then do -log10
message("Getting FDR from ", colnames(DESeq2.stats)[5])
fdr = DESeq2.stats$padj
minus.log10.fdr = -log(fdr, 10)
# Get the log2 fold change
log2FC = DESeq2.stats[,2]
# log2 mean
log2.counts.mean = log(DESeq2.stats[,1] + counts.mean.offset, 2)
# 0 - Log10 FDR
min.non0.FDR = min(fdr[fdr > 0], na.rm=T)

# Set FDRs of 0 to min.non0.FDR
fdr.no0 = sapply(fdr, function(x) ifelse(x==0,min.non0.FDR,x))
min(fdr.no0)
negLog10FDR = 0 - log10(fdr.no0)

# Filter by FDR and fold change (for coloring)
filtered.rows =      ((fdr < fdr.threshold) & (abs(log2FC) > log2FC.threshold))
filtered.rows.up =   ((fdr < fdr.threshold) & (abs(log2FC) > log2FC.threshold) & log2FC > 0)
filtered.rows.down = ((fdr < fdr.threshold) & (abs(log2FC) > log2FC.threshold) & log2FC < 0)

# Use p-value threshold of 100th gene in this list
DESeq2.stats.filtered = DESeq2.stats[filtered.rows,]
# Select the top num.genes.to.label genes (or all of them, if DESeq2.stats.filtered is smaller)
# Fixed 2 Feb 2019
# pvalue.threshold.labels = DESeq2.stats.filtered[num.genes.to.label,5]
pvalue.threshold.labels = DESeq2.stats.filtered[min(num.genes.to.label,nrow(DESeq2.stats.filtered)),5]

# Testing
message(paste("Raw p-value threshold for labels (top", num.genes.to.label, "genes by p-value) = ", pvalue.threshold.labels))

# filtered.rows.labels =      ((p.value <= pvalue.threshold.labels) & (abs(log2FC) > log2FC.threshold))
# filtered.rows.up.labels =   ((p.value <= pvalue.threshold.labels) & (abs(log2FC) > log2FC.threshold) & log2FC > 0)
# filtered.rows.down.labels = ((p.value <= pvalue.threshold.labels) & (abs(log2FC) > log2FC.threshold) & log2FC < 0)
filtered.rows.labels =      ((p.value <= pvalue.threshold.labels) & (abs(log2FC) > log2FC.threshold) & (fdr < fdr.threshold))
filtered.rows.up.labels =   ((p.value <= pvalue.threshold.labels) & (abs(log2FC) > log2FC.threshold) & log2FC > 0 & (fdr < fdr.threshold))
filtered.rows.down.labels = ((p.value <= pvalue.threshold.labels) & (abs(log2FC) > log2FC.threshold) & log2FC < 0 & (fdr < fdr.threshold))

num.points.to.label = length(log2FC[filtered.rows])
message(paste(num.points.to.label, "points meet FDR and FC thresholds."))
message(paste(length(log2FC[filtered.rows.labels]), "points will be colored (top ~", num.genes.to.label, "by p-value, after filtering by FDR and FC)."))
gene.symbols[filtered.rows.labels]

comparison = paste(sampleB.name, "vs",  sampleA.name)
log2FC.label = paste("log2 ( ", sampleA.name, " / ",  sampleB.name, " )", sep="")

# Create plot with ggplot2
input = gsub(".txt", "",   input.filename)
pdf.filename = paste("Volcano_plot", input,  sampleA.name, "vs",  sampleB.name, fdr.threshold, log2FC.threshold, "pdf", sep=".")
pdf(pdf.filename, w=11, h=8.5, useDingbats=FALSE)

if (num.points.to.label < num.genes.to.label)
{
	ggplot() +
	# xlim(0, 20) + 
	# ylim(y.axis.min, y.axis.max) + 
	geom_point(aes(log2FC, negLog10FDR), color="darkgray", size=1) + 
	geom_point(aes(log2FC[filtered.rows.up], negLog10FDR[filtered.rows.up]), color=up.color, size=1) + 
	geom_point(aes(log2FC[filtered.rows.down], negLog10FDR[filtered.rows.down]), color=down.color, size=1) + 
	geom_point(aes(log2FC[filtered.rows.up.labels], negLog10FDR[filtered.rows.up.labels]), color="red", size=3) + 
	geom_point(aes(log2FC[filtered.rows.down.labels], negLog10FDR[filtered.rows.down.labels]), color="blue", size=3) + 
	geom_vline(xintercept = 0, color="grey80") + 
	labs(x = log2FC.label) +
	labs(y = "-log10 (FDR)") + 
	theme_classic(base_size = 16) +
	ggtitle(figure.title) +
	theme(plot.title = element_text(hjust = 0.5)) +
	geom_text_repel(aes(log2FC[filtered.rows.labels], negLog10FDR[filtered.rows.labels]), label = gene.symbols[filtered.rows.labels], size=4, max.overlaps=100)
} else {
	ggplot() +
	# xlim(0, 20) + 
	# ylim(y.axis.min, y.axis.max) + 
	geom_point(aes(log2FC, negLog10FDR), color="darkgray", size=1) + 
	geom_point(aes(log2FC[filtered.rows.up], negLog10FDR[filtered.rows.up]), color=up.color, size=1) + 
	geom_point(aes(log2FC[filtered.rows.down], negLog10FDR[filtered.rows.down]), color=down.color, size=1) + 
	geom_point(aes(log2FC[filtered.rows.up.labels], negLog10FDR[filtered.rows.up.labels]), color="red", size=3) + 
	geom_point(aes(log2FC[filtered.rows.down.labels], negLog10FDR[filtered.rows.down.labels]), color="blue", size=3) + 
	geom_vline(xintercept = 0, color="grey80") + 
	labs(x = log2FC.label) +
	labs(y = "-log10 (FDR)") + 
	theme_classic(base_size = 16) +
	ggtitle(figure.title) +
	theme(plot.title = element_text(hjust = 0.5)) +
	geom_text_repel(aes(log2FC[filtered.rows.labels], negLog10FDR[filtered.rows.labels]), label = gene.symbols[filtered.rows.labels], size=4, max.overlaps=100)
}

foo = dev.off()

message(paste("Created plot!  See", pdf.filename, "\n"))

# Save details of your session
sink("sessionInfo.txt")
date()
print(sessionInfo(), locale=FALSE)
sink()
