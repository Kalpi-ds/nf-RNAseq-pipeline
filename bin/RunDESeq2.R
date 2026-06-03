#!/usr/bin/env Rscript

###
###  Run DESeq2 from the command line
###  George Bell - Bioinformatics and Research Computing, Whitehead Institute
###
###  USAGE: R --vanilla < RunDESeq2.R inputCounts OutputFile [group for each sample]
###  EX: R --vanilla < RunDESeq2.R InputDEseq2.txt DESeq2_output.txt UHR UHR brain brain
###
###  Note that comparison is based on alphabetical order of samples, 
###  so with samples A and B output ratios will be for B/A
###  If you want the other direction, modify the names of the groups
###
###  Updated 20 Nov 2012:
###    Change estimateVarianceFunctions() to estimateDispersions() [newest command, needed for DESeq 1.8.3]
###    Change handling of sample names to always compare experimental/control
###      and include sample names in output file
###  Updated 19 Mar 2013:
###    Include normalized counts in output file
###  Updated 10 July 2013:
###    Make PCA figure using vst-transformed counts matrix
###  Updated 24 September 2013:
###    Adapt for DESeq2
###  Updated 15 Jun 2017:
###    v1.16 does not use betaPrior by default, need to use lfcShrink function ==> otherwise none of the DE genes will be
###    in common with edgeR (see line 89 below)
###    Also see discussion from M.Love https://support.bioconductor.org/p/95695/
###  Updated 13 December 2017
###    Change results(dds) to results(dds, independentFiltering=FALSE) [to get FDR values for all genes]
###  Updated 26 Feb 2018
###    Round output statistics to significant.figures [initially set to 6]
###  Updated 12 Apr 2018
###    add log2 normalized reads by rlog, and change from vst to rlog for PCA plot 
###  Updated 6 June 2018
###    Add raw p-value histogram for QC
###  Updated 7 Jul 2018
###		Removed "res=res" in lfcShrink.  
###  Updated 21 Aug 2018
###		Added MA plot and volcano plot.  
###  Updated 28 Aug 2018
###		Added MA plot for pre-lfcShrink values.  
###  Updated 4 February 2021
###		Allow for redundant feature (gene) IDs as input.  If present, a unique numerical suffix is added. 
###  Updated 25 May 2022
###     Use normal shrinkage and add a call to draw_MA_plot_from_DESeq2_analysis.R and draw_volcano_plot_from_DESeq2_analysis.R
###  Updated 20 July 2022
###		Force input matrix to be integers (so add a rounding step) [GB]
###  Updated 5 March 2024
###		Change default log2FC threshold (for plots) from 0.7 to 1 and move variables to top [GB]

offset = 0
significant.figures = 6
fdr.threshold    = as.numeric(Sys.getenv("PADJ_CUTOFF",    unset = "0.05"))
log2FC.threshold = as.numeric(Sys.getenv("LFC_THRESHOLD",  unset = "1"))

if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
	message("Running as R --vanilla ....")
	offset = 3
}

# Get the full path to this script
this.script = commandArgs()[4]
this.script = sub("^--file=", "", this.script)
if (is.na(this.script)) { this.script = "./DESeq2_normalize_only.R" }

if (length(commandArgs()) < (7 - offset))
{
	message("\nAnalyze a matrix of raw counts with DESeq2 to identify differentially expressed genes/features.")
	message(paste("USAGE:   ", this.script, "inputCounts OutputFile group1 group2 [...]"))
	message(paste("Example: ", this.script, "Input_DESeq2.txt DESeq2_output.txt UHR UHR brain brain\n"))
	message(paste("Note that the first group is selected as the control, and the last group is selected as \'experimental\'."))
	message(paste("Use for DESeq2 v1.16 or higher.  Otherwise use RunDESeq2_old_v1.10.R\n"))
	q()
}

# Load library (quietly)
message("\nLoading DESeq2 and other required R packages ....")
suppressMessages(library(DESeq2))
library("vsn")


message(paste("Running DESeq2 version", packageDescription("DESeq2")$Version, "\n"))

# First argument is input file
input.filename = commandArgs()[6 - offset]
# Second argument is output file
output.filename = commandArgs()[7 - offset]
# Third and fourth arguments are groups (in same order as input file)
groups = commandArgs()[(8 - offset):length(commandArgs())]

message(paste("Input filename is", input.filename))
message(paste("Output filename is", output.filename))
pca.figure.filename = paste(output.filename, "PCA.pdf", sep=".")
# Remove any ".txt." in the name
pca.figure.filename = gsub(".txt.", ".", pca.figure.filename)

disp.figure.filename = paste(output.filename, "dispersion.pdf", sep=".")

message("Reading counts matrix, making IDs unique (by adding a numerical suffix) if required ....")

# Read matrix of counts
# Change method so duplicate row names will be tolerated.  28 January 2021
# counts.matrix = read.delim(data.matrix.file, row.names=1, check.names=F)
# Permit gene matrices with duplicated gene symbols [GB, 7 Sep 2018]	Added "check.names=FALSE" on 28 Jan 2021
counts.matrix = read.delim(input.filename, header=T, sep="\t", check.names=FALSE)
# Get total counts per gene (for sorting redundant rows)
counts.sum.by.gene = apply(counts.matrix[,2:ncol(counts.matrix)], 1, sum)
# Order matrix by decreasing counts.sum.by.gene, adding a numeric suffix to each symbol we've already seen.
gene.symbols.original = as.character(counts.matrix[,1])
# gene.symbols.unique = .Internal(make.names(gene.symbols.original, unique=TRUE))
gene.symbols.unique = make.unique(gene.symbols.original)
o = order(gene.symbols.original != gene.symbols.unique, counts.sum.by.gene, decreasing=TRUE)
gene.symbols.unique[o] = make.unique(gene.symbols.unique[o])
# Make row names from unique gene symbols.
rownames(counts.matrix) = gene.symbols.unique
# Drop the gene symbol (first column) since we have them as row names
counts = counts.matrix[,-1]

# Force counts values into integers -- 20 Jul 2022 [GB]
counts = round(counts, 0)

message("Groups are")
groups
control.group = groups[1]
exp.group = groups[length(groups)]
message(paste("Control group is", control.group))
message(paste("Experimental group is", exp.group))

# Make a DESeqDataSet
# Need to relevel because of inconsistent strangeness
# dds = DESeqDataSetFromMatrix(countData = counts, colData = DataFrame(condition=factor(groups)), design = ~ condition)
dds = DESeqDataSetFromMatrix(countData = counts, colData = DataFrame(condition=relevel(factor(groups), ref=control.group)), design = ~ condition)
# factor(groups)

###  MAJOR COMMAND:
# Estimate size factors and dispersions and fit generalized linear model
dds = DESeq(dds)

# Print size factors
sizeFactors(dds)


# check dispersion distribution
# If fitted curve is off, you can try with local dispersion in stead of parametric dispersion.
# You can estimated the dispersion with residue as mentioned in the post by Michael Love:
# https://support.bioconductor.org/p/81094/
disp.figure.filename = gsub(".txt.", ".", disp.figure.filename)
message(paste("Create dispersion distribution figure: ", disp.figure.filename, "\n"))
pdf(disp.figure.filename, useDingbats=FALSE)
plotDispEsts(dds)
foo = dev.off()


# Do stats
#NOTE: by default, results function has independentFiltering=TRUE => in general, genes with low counts (based on mean normalized count) are removed (ie. adj p-val is NA for these)
#use results(dds, independentFiltering=FALSE) to get adj p-val for all genes, however, you may want to do extra filtering (eg. remove genes with low counts).
#see: https://bioconductor.org/packages/3.7/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#indfilt
# res = results(dds)
res = results(dds, independentFiltering=FALSE)
#add for v1.16, see https://support.bioconductor.org/p/95695/, including res in lfcShrink will slightly affect adj p-val e.g. lfcShrink(dds,coef=2,res=res).
# resultsNames(dds)
# res.shrunk = lfcShrink(dds, coef=2)	# Change from 'normal' to 'apeglm' (11 Feb 2021)
# 05 May 2022 #change back to normal shrink
message("Running lfcShrink() using the normal method ...")
res.shrunk = lfcShrink(dds, type="normal", coef=2)

# change header for log2FC column
colnames(res.shrunk)[2] = paste("log2(", paste(exp.group, control.group, sep="/"), ")", sep="")

# Print a histogram of raw p-values for QC
# For interpretation, see http://varianceexplained.org/statistics/interpreting-pvalue-histogram/
pval.figure.filename = paste(output.filename, "pvalue_histogram.pdf", sep=".")
# Remove any ".txt." in the name
pval.figure.filename = gsub(".txt.", ".", pval.figure.filename)
message(paste("Creating a histogram of p-values for QC:", pval.figure.filename, "\n"))
pdf(pval.figure.filename, useDingbats=FALSE)
hist(res.shrunk$pvalue, breaks=20, main="Histogram of raw p-values (for QC)", xlab="raw p-value", col="wheat")
foo = dev.off()

# based on 
dds.rlog = rlog(dds, blind=FALSE )


#for multiple conditions or specific contrasts use
#res_brain_UHR<-results(dds,contrast=c("condition","brain","UHR"),independentFiltering=FALSE)
#res_brain_UHR<-lfcShrink(dds,coef="condition_brain_UHR")
#use res_brain_UHR
#Note: lfcShrink function has changed since 1.16, https://support.bioconductor.org/p/98833/, https://support.bioconductor.org/p/95695/


# Add normalized counts for this version (GB - 19 Mar 2013)
counts.normalized = round(t(t(counts(dds))/sizeFactors(dds)), 2)
colnames(counts.normalized) = paste(colnames(counts), "norm", sep=".")

# Add normalized counts for this version (BY - 12 Mar 2018)
counts.rlog = round(assay(dds.rlog),2)
colnames(counts.rlog ) = paste(colnames(counts), "rlog", sep=".")

# Print output (including norm counts)
message(paste("Rounding output statistics to", significant.figures, "significant figures....\n"))
output.table = cbind(rownames(dds), signif(as.matrix(res.shrunk), significant.figures), counts, counts.normalized, counts.rlog)
colnames(output.table)[1] = "Feature.ID"
write.table(output.table, file=output.filename, sep="\t", quote=F, row.names=F)

# We may need to modify the function plotPCA() because it makes a max of 12 colors (so >12 samples cause problems)
plotPCA_BaRC = function (x, intgroup = "condition", ntop = 500)
{
	library(genefilter)
	library(RColorBrewer)
	rv = rowVars(exprs(x), useNames=TRUE)
	select = order(rv, decreasing = TRUE)[seq_len(ntop)]
	pca = prcomp(t(exprs(x)[select, ]))
	fac = factor(apply(pData(x)[, intgroup, drop = FALSE], 1,
		paste, collapse = " : "))
	if (nlevels(fac) >= 3) {
		colours = brewer.pal(nlevels(fac), "Paired")
	} else { colours = c("green", "blue") }
	xyplot(PC2 ~ PC1, groups = fac, data = as.data.frame(pca$x), 
	  pch = 16, las=1, cex = 2, aspect = "iso", col = colours, main = draw.key(key = list(rect = list(col = colours), 
	  text = list(levels(fac)), rep = FALSE)))
}

DESeq2_plotPCA_BaRC = 
function (object, intgroup = "condition", ntop = 500, returnData = FALSE) 
{
	library(ggplot2)
    rv <- rowVars(assay(object), useNames=TRUE)
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
        length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
    percentVar <- pca$sdev^2/sum(pca$sdev^2)
    if (!all(intgroup %in% names(colData(object)))) {
        stop("the argument 'intgroup' should specify columns of colData(dds)")
    }
    intgroup.df <- as.data.frame(colData(object)[, intgroup, 
        drop = FALSE])
    group <- if (length(intgroup) > 1) {
        factor(apply(intgroup.df, 1, paste, collapse = " : "))
    }
    else {
        colData(object)[[intgroup]]
    }
    d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
        intgroup.df, name = colnames(object))
    if (returnData) {
        attr(d, "percentVar") <- percentVar[1:2]
        return(d)
    }
    print(ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
        geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
        100), "% variance")) + coord_fixed())
    # Add variation of this plot that fills the page
    print(ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group")) + 
        geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 
        100), "% variance")) + ylab(paste0("PC2: ", round(percentVar[2] * 
        100), "% variance")))
}

### Do VST transformation (recommended for PCA plot)
#dds.vst = varianceStabilizingTransformation(dds, blind=TRUE)
#rld = rlogTransformation(dds, blind=TRUE)
### Run Principal Components Analysis
#pdf(pca.figure.filename, useDingbats=FALSE)
#plotPCA(rld, intgroup=c("condition"))
#plotPCA_BaRC(rld)
#dev.off()


pdf(pca.figure.filename, useDingbats=FALSE)
# PCA plot
DESeq2_plotPCA_BaRC(dds.rlog, intgroup=c("condition"))
# sd vs ranked mean
meanSdPlot(assay(dds.rlog))
boxplot(assay(dds.rlog), las=1)
foo = dev.off()

# MA plot (with unshrunken log fold changes)
counts.mean = res$baseMean
log2.ratio = res[,2]
# Add 1 pseudocount before taking logs
log2.counts.mean = log2(counts.mean + 1)
comparison = paste("log2 (", paste(exp.group, control.group, sep=" / "), ")", sep="")
MA.plot.noShrink.filename = paste(output.filename, "MA_plot.no_lfcShrink.pdf", sep=".")
# Remove any ".txt." in the name
MA.plot.noShrink.filename = gsub(".txt.", ".", MA.plot.noShrink.filename)
pdf(MA.plot.noShrink.filename, w=11, h=8.5, useDingbats=FALSE)
plot(log2.counts.mean, log2.ratio, main=paste("MA plot (no lfcShrink):", exp.group, "vs.", control.group),
  pch=20, ylab=comparison, xlab="log2 (mean expression level)", cex=0.5, xlim=c(0, max(log2.counts.mean)), las=1)
abline(h=0, col="gray")
foo = dev.off()

# MA and Volcano plots
# Try BaRC NFS scripts first; fall back to built-in ggplot2 versions so the
# pipeline works outside the Whitehead cluster without any NFS access.
message("For MA and volcano plots, FDR threshold = ", fdr.threshold, ", log2FC threshold = ", log2FC.threshold)

ma.script      = Sys.which("draw_MA_plot_from_DESeq2_analysis.R")
volcano.script = Sys.which("draw_volcano_plot_from_DESeq2_analysis.R")

if (nchar(ma.script) > 0 && nchar(volcano.script) > 0) {
    MA.plot.command      = paste(ma.script,      output.filename, exp.group, control.group, fdr.threshold, log2FC.threshold)
    volcano.plot.command = paste(volcano.script, output.filename, exp.group, control.group, fdr.threshold, log2FC.threshold)
    message(MA.plot.command)
    system(MA.plot.command)
    message(volcano.plot.command)
    system(volcano.plot.command)
} else {
    message("BaRC NFS scripts not found — generating built-in ggplot2 MA and Volcano plots.")
    suppressMessages(library(ggplot2))
    suppressMessages(library(ggrepel))

    # Build a shared data frame from shrunk results
    plot.df = data.frame(
        baseMean = as.numeric(res.shrunk$baseMean),
        log2FC   = as.numeric(res.shrunk[, 2]),
        padj     = as.numeric(res.shrunk$padj),
        gene     = rownames(dds),
        stringsAsFactors = FALSE
    )
    plot.df$significant = !is.na(plot.df$padj) &
                           plot.df$padj  <= fdr.threshold &
                           abs(plot.df$log2FC) >= log2FC.threshold
    plot.df$color = ifelse(plot.df$significant, "Significant", "Not significant")

    # ── MA plot ────────────────────────────────────────────────────────────────
    ma.filename = paste(output.filename, "MA_plot.pdf", sep=".")
    ma.filename = gsub(".txt.", ".", ma.filename)
    message("Writing built-in MA plot: ", ma.filename)
    pdf(ma.filename, width=8, height=6, useDingbats=FALSE)
    print(
        ggplot(plot.df, aes(x = log10(baseMean + 1), y = log2FC, color = color)) +
        geom_point(size = 0.6, alpha = 0.6) +
        scale_color_manual(values = c("Significant" = "red", "Not significant" = "grey60")) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +
        geom_hline(yintercept = c(-log2FC.threshold, log2FC.threshold),
                   linetype = "dotted", color = "darkblue") +
        labs(
            title  = paste("MA plot:", exp.group, "vs.", control.group),
            x      = "log10 (mean expression + 1)",
            y      = paste0("log2(", exp.group, "/", control.group, ")"),
            color  = ""
        ) +
        theme_bw(base_size = 12) +
        theme(legend.position = "bottom")
    )
    foo = dev.off()

    # ── Volcano plot ───────────────────────────────────────────────────────────
    plot.df$neglog10padj = -log10(ifelse(is.na(plot.df$padj), 1, plot.df$padj) + 1e-300)
    volcano.filename = paste(output.filename, "Volcano_plot.pdf", sep=".")
    volcano.filename = gsub(".txt.", ".", volcano.filename)
    message("Writing built-in Volcano plot: ", volcano.filename)
    pdf(volcano.filename, width=8, height=6, useDingbats=FALSE)
    print(
        ggplot(plot.df, aes(x = log2FC, y = neglog10padj, color = color)) +
        geom_point(size = 0.6, alpha = 0.6) +
        scale_color_manual(values = c("Significant" = "red", "Not significant" = "grey60")) +
        geom_vline(xintercept = c(-log2FC.threshold, log2FC.threshold),
                   linetype = "dashed", color = "blue") +
        geom_hline(yintercept = -log10(fdr.threshold),
                   linetype = "dashed", color = "blue") +
        labs(
            title  = paste("Volcano plot:", exp.group, "vs.", control.group),
            x      = paste0("log2(", exp.group, "/", control.group, ")"),
            y      = "-log10(adjusted p-value)",
            color  = ""
        ) +
        theme_bw(base_size = 12) +
        theme(legend.position = "bottom")
    )
    foo = dev.off()
}

############

# Save details of your session
sink(paste("sessionInfo", Sys.Date(), "txt", sep="."))
date()
print(sessionInfo(), locale=FALSE)
sink()

###################
