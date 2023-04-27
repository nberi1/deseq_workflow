##### -------- MUTU I 4-condition Differential Expression Analysis ------------
###
### Visualize data quality, identify differentially expressed genes, and ######
### graph important characteristics in any RNA-seq dataset               ######

# 1 Call libraries ----
# 1.1   call libraries ---------------------------------------------------- --------------------------------------------------------------------------------------------------------------------

myLibraries <- c("DESeq2",
                 "tximport",
                 "ensembldb",
                 "EnsDb.Hsapiens.v75",
                 "readr",
                 "tximportData",
                 "GenomicFeatures",
                 suppressPackageStartupMessages("GenomicFeatures"),
                 suppressPackageStartupMessages("SummarizedExperiment"),
                 "pheatmap",
                 "RColorBrewer",
                 "biomaRt",
                 "ggplot2",
                 "DEGreport",
                 "ggrepel",
                 "tidyverse",
                 "AnnotationDbi",
                 "AnnotationHub",
                 "RMariaDB",
                 "dendextend",
                 "vsn",
                 "hexbin",
                 "tximeta",
                 "rtracklayer",
                 "stringr",
                 "jsonlite",
                 "GenomeInfoDb",
                 "magrittr",
                 "genbankr",
                 "rentrez",
                 "BiocIO",
                 "gsubfn",
                 "reshape2",
                 "stringr")

invisible(lapply(myLibraries, library, ch = TRUE))

# ----
# 2 Declare functions ----
# 2.n declarations matches with 3.n eg 2.2 fxns are used in    section 
# 2.1   list required input variables ----
# directory, working (dir)                                       3.1
# contains cDNA-alignments (store all your input files here)
# contains output-files (go here to find your results)
# akataGtfDir (where to find the GTF for akata)                  3.1
# sample metadata - should always be samples.txt                 3.2
# reference condition number (refNum)                            3.2
# assembly, human (current functions only let you use Ensembl)   3.3
# assembly, EBV                                                  4.2
### List variables that will change by experiment
# inputDir                                                       3.1
# project name (projName)/sampleAttributes                       3.2
# filenames for pairwise comparisons                             3.9a
### Functions that save a file                                   2.4, 2.7, 2.8, 2.9
# 2.2   declare functions to get sample metadata -------------------------- ----------------------------------

# get sample metadata
# input is path to working directory
# returns a table of sample metadata
# recommend to call the output "samps"
getSampleMetadata <- function(workingDir) {
  # this fxn organizes the sample metadata from a premade table
  # input a table listing each file in the set of samples
  # one column must be named "trt" (may chagne later on)
  # because right now my treatments are all trts
  # and one column must be named "cond" as in condition
  # get sample metadata
  samples <- read.table(file.path(workingDir,"samples.txt"), header=TRUE)
  samples$cond <- as.factor(samples$cond)
  return(samples)
  
}

# get sample attributes
# note, condition 1 being WT is still the best functioning way to do this
# but I should work on changing that.
# inputs: project name (arbitrary), human genome build, EBV genome build,
# samples dataframe, condition # of reference, 

# get attributes by inputting the samples.txt data table
getAttributes <- function(proj, human, ebv, samples, refNum) {
  
  statement <- paste0("Your project is named ", proj)
  print(statement)
  
  # information about the reference condition
  refName <- unique(samples$trt[which(samps$cond == refNum)])
  
  # get the metadata for all samples that are in the reference condition
  refMeta <- samps[which(samps$trt == refName),]
  refMeta
  numRefCond <- unique(refMeta$cond)
  numRefCond
  
  # information about the experimental condition(s)
  expMeta <- samps[which(samps$trt != refName),]

  # count the number of experimental conditions
  numExpCond <- unique(expMeta$cond)

  # get the names of the experimental conditions (as a list)
  expNames <- unique(expMeta$trt)
  
  cat("Your reference treatment is:", refName, "\nYour experimental treatment(s) are:", expNames)

  
  meta <- list("projectName" = proj,
            "refName" = refName, 
            "numRefCond" = numRefCond, 
            "expNames" = expNames, 
            "numExpCond" = numExpCond,
            "humanGenome" = human,
            "ebvGenome" = ebv)
  
  return(meta)
  
}


# fxn for getting file paths based on genome used
# inputs are dir, genome, and nucType (type of nucleotide indexed)
# currently options rae hg19 and akata
# but whatever you have should work too 
# genome is either hg19 or akata
# nucType is either gDNA or cDNA
getFilePaths <- function(dir, genome, nucType) {
  fold <- paste0(genome, "_", nucType, "_alignments/") #, samps$run)
  #  print(fold)
  filename <- paste0(samps$run, "_quant.sf")
  files <- file.path(dir, fold, filename)
  return(files)
  
}
# files <- getFilePaths(dir, "hg19", "cDNA")

# 2.3   declare functions for human genome annotation --------------------- --------------------------------------------------

# human genome build from ensembl
# and annotation shenanigans part one 
fetchFromEnsembl <- function(org, rel) {
  myDB <- makeTxDbFromEnsembl(organism = org, release = rel)
  return(myDB)
}

# tximport [[ FROM ENSEMBL ONLY ]]
preImport <- function(db, files) {
  k <- keys(db, keytype = "TXNAME")
  tx2gene <- AnnotationDbi::select(db, k, "GENEID", "TXNAME")
  txi <- tximport(files, type = "salmon", 
                  tx2gene = tx2gene,
                  countsFromAbundance = "lengthScaledTPM",
                  ignoreTxVersion = TRUE)
  return(txi)
}

# get coldata for txi
getColdata <- function(txi) {
  cts <- txi$counts
  cts <- cts[rowSums(cts) > 0, ]
  coldata <- data.frame(files,
                        names = paste0(samps$trt, "_", samps$rep),
                        cond = samps$cond,
                        trt = samps$trt,
                        stringsAsFactors = FALSE
  )
  return(coldata)
  
}

# 2.4   declare fxn to graph hg19-aligned data quality -------------------- ---------------------------------
# save files with experiment-wide specific name
# ie if the experiment is basically all about bcl6 KO
# then give it some name that makes sense like bcl6KO
# arguments w = datasetName variable; x = "arbitrary naming string", y = dataset
saveMeanSdPlot <- function(datasetName, plotType, dataset) {
  filename <- paste0(outputDir, datasetName, "_", plotType, "_meanSdPlot.png")
  cat(filename)
  
  # save a .png image of the graph
  png(file = filename,
      width = 600,
      height = 350)
  
  meanSdPlot(assay(dataset))

  dev.off()


}
#

# for whatever reason tihs function and htis one only
# is more consistent if I do dev.off() outside of the fxn
savePCAPlot <- function(datasetName, plotType, dataset) {
  
  filename <- paste0(outputDir, datasetName, "_", plotType, "_PCAPlot.png")
  cat(filename)
  
  png(file = filename,
      width = 600,
      height = 350)
  
  plotPCA(dataset,
          intgroup = "trt") +
    theme_classic() +
    scale_colour_discrete(name = "Treatment",
                          labels = c(sampAttributes$refName, sampAttributes$expNames))
  
  # dev.off()
  
}
#

# make sample correlation heatmap
makeSampleCorrelationHeatmap <- function(datasetName, plotType, dataset) {
  vsd_mat <- assay(dataset)
  vsd_mat <- vsd_mat[rowSums(vsd_mat) > 0,]
  
  vsd_cor <- cor(vsd_mat)
  rownames(vsd_cor) <- paste0(samps$trt,
                              "_rep_",
                              samps$rep)
  colnames(vsd_cor) <- rownames(vsd_cor)
  
  filename <- paste0(outputDir, datasetName, "_", plotType, "_correlation.png")
  filename
  
  png(file = filename,
      width = 600,
      height = 350)
  
  pheatmap(vsd_cor)
  dev.off() 
}
#

# make sample distribution heatmap
makeSampleDistributionHeatmap <- function(datasetName, plotType, dataset) {
  
  sampleDists <- dist(t(assay(dataset)))
  
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- vsd$trt
  colnames(sampleDistMatrix) <- colnames(vsd)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  
  filename <- paste0(outputDir, datasetName, "_", plotType, "_distribution.png")
  filename
  
  png(file = filename,
      width = 600,
      height = 350)
  
  pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors)
  
  dev.off()
  
}
#

# 2.5   placeholder ----
# 2.6   placeholder ----
# 2.7   declare functions to get and plot pairwise results ----

# get pairwise log2 fold change for each "exp" vs "wt" trt---
# wt code commented out so i can come back to it later and change
# to pairwise comparison of any two conditions
getPairwiseLFC <- function(exp, ref, dds) { #}, wt) {
  coef <- paste0("cond_", exp, "_vs_", ref) # "_vs_", wt)
  print(coef)
  lfcShrunk <- lfcShrink(dds,
                         coef = coef,
                         #                         coef = paste0("condition_", con1, "_", con2),
                         type = "apeglm")
  return(lfcShrunk)
  
}
#

# ggplot2 to make a volcano plots for LFC output
# takes a few minutes to run 
# inputs: x is the subset (delimited by "$") of the lapply output (a)
# ie, a$x that you want to plot and the name "name" is
# whatever you want to call the plot
plotOutput <- function(x, name, xmax, ymax, genome){
  # Visualize the comparison from the "compare" function
  # Input is the result of the compare function
  # Returns a ggplot2 object
  xmin <- (xmax * -1)
  
  df <- as.data.frame(x)
  df$ensembl_gene_id <- rownames(df)
  rownames(df) <- NULL

  x$ensembl_gene_id <- rownames(x)
  rownames(x) <- NULL

  g <- ggplot(data = df,
              aes(x = log2FoldChange,
                  y = -log10(padj),
                  xmin = xmin, xmax = xmax,
                  ymin = 0, ymax = ymax)) +
    geom_point() +
#    geom_text_repel(aes(label = ensembl_gene_id)) +
    theme_classic()

  ggsave(g, filename = paste0(outputDir, "/", sampAttributes$projectName, "_", genome, "_", name, "_unlabeled_volcano_plot.png"), height = 5, width = 5, device = "png")
  print(paste0("Your volcano plot comparing treatment ", name, " to WT has been saved to ", outputDir))
  

}
#

# 2.8   declare functions to annotate genes via Ensembl db ---------------- --------------------------
# Add annotations as desired
# write to files

# set up biomart with your chosen human genome build from Ensembl
# get Anns meaning get Annotations
getAnns <- function(dataset, biomart, version){
  # return list of genes
  # needs ensembl, biomarts
  marts <- listMarts()
  mart <- useEnsembl(dataset = dataset,
                     biomart = biomart,
                     version = version)
  
  myMart <- useDataset(dataset, mart = mart)
  
  attributes <- c("ensembl_transcript_id",
                  "ensembl_gene_id",
                  "chromosome_name",
                  "start_position",
                  "end_position",
                  "hgnc_symbol",
                  "description",
                  "gene_biotype")
  
  genes <- getBM(attributes = attributes,
                 values = rownames(normalized_counts),
                 mart = myMart,
                 useCache = FALSE)
  #  head(genes, 10)
  
  genes1 <- genes[,c(2:8)]
  return(genes1)
}
#

# declare fxn to annotate genes
annotateGenes <- function(dat, anns) {
  dat$ensembl_gene_id = rownames(dat)
  rownames(dat) = NULL
  annotations <- merge(x = dat, 
             y = anns, 
             by = "ensembl_gene_id", 
             all.x = TRUE, 
             all.y = FALSE)
  return(annotations)
}
#

# declare function to annotate genes and save tables to files
# makeAnnFiles includes the annotateGenes() function nested within it
makeAnnFiles <- function(dat, anns, filename) {
  file = paste0(outputDir, filename, ".tsv")
  
  dat <- drop_na(dat)
  data_ann <- annotateGenes(dat, anns)
  data_ann <- unique(data_ann)
  
  write.table(data_ann,
              file = file,
              quote = FALSE,
              sep = "\t",
              row.names = FALSE)
  
  cat(file)
  return(data_ann)
}
#

# declare function to annotate genes and save tables to files
# makeAnnFiles includes the annotateGenes() function nested within it
# makeAnnFiles_CompleteCases <- function(dat, anns, filename) {
#   
#   cases <- complete.cases(dat)
#   completed <- dat[cases, ]
#   data_ann <- annotateGenes(completed, anns)
#   data_ann <- unique(data_ann)
#   
#   write.table(data_ann,
#               file = paste0(outputDir, filename, ".tsv"),
#               quote = FALSE,
#               sep = "\t",
#               row.names = FALSE)
# #  cat("your file was written to", paste0(outputDir, filename, ".tsv"))
#   
# #  return(data_ann)
# }
#

# 2.9   declare fxns to get and visualize DEGs ---------------------------- ----
# declare function to identify DEGs
getDEGs <- function(x, hgnc_symbol = hgnc_symbol) {
  x$fc <- NA
  x$sig <- NA
  x$isDEG <- NA
  
  x$fc[which(x$log2FoldChange > 1)] <- "up"
  x$fc[which(x$log2FoldChange < -1)] <- "down"
  x$sig[which(x$padj < 0.05)] <- "yes"
  
  x$isDEG[which(x$fc == "up" & x$sig == "yes")] <- "red"
    x$isDEG[which(x$fc == "down" & x$sig == "yes")] <- "green"
      
    x$name <- NA
    x <- x[ ,!names(x) %in% c("ensembl_gene_id")]
    
    x <- x %>% mutate(name = ifelse(isDEG == "green" | isDEG == "red",
                                    hgnc_symbol, NA))
    x <- unique(x)
    return(x)
}
#

# plot DEGs 
 plotDEGs <- function(dat, xmax, ymax, genome) {
  xmin = -1*(xmax)
  name <- deparse(substitute(dat))
  filename = paste0(outputDir, "/", sampAttributes$projectName, name, ".png")

  #
    g <- ggplot(data = dat,
              aes(x = log2FoldChange,
                  y = -log10(padj),
                  label = name,
                  xmin = xmin, xmax = xmax,
                  ymin = 0, ymax = ymax)) + #,
#                  color = isDEG)) +
      geom_point() +
      geom_text_repel() +
      theme_classic() # +
      # scale_shape_manual(name = "Expression change",
      #                    labels = c("Down", "None", "Up"))
    #  
  
  #
  ggsave(filename, device = "png", height = 5, width = 6)
  #
  return(filename)
}
#

# ----
# 3 Set up environment and process data for hg19-aligned reads ----
# 3.1   assign directories ------------------------------------------------ ---------------------------------------------------

# assign working directory to dir
dir <- "C:/Users/NRB37/Dropbox (Partners HealthCare)/Nina/coworkers/Yifei/mutu/data_analysis/"
cat("Your working directory is", dir)

# go one level into the directory for where all the input files are stored
# note this is all cDNA (ie, RNA-seq) files
inputDir <- paste0(dir, "hg19_cDNA_alignments")
list.files(inputDir)
cat("Input files will be taken from", inputDir)

# same for output files, meaning all graphs and plots and DESeq2 results
outputDir <- paste0(dir, "output_files/")
list.files(outputDir)
cat("Output will appear in ", outputDir)

# 3.2   get sample metadata and hg19 alignments --------------------------- -------

# get sample metadata
samps <- getSampleMetadata(dir)
samps

# first argument is project name, it's arbitrary
sampAttributes <- as.list(getAttributes("YL_cytokines_MUTU", "hg19", "akata", samps, 1))
sampAttributes

# extract reference treatment number and treatment name for use in colData 
refNum <- sampAttributes$numRefCond
names(refNum) <- c(sampAttributes$refName)
refNum

# same for the experimental treatment number and name
expNums <- sampAttributes$numExpCond
names(expNums) <- sampAttributes$expNames
expNums

# get hg19 cDNA transcriptome alignment file paths 
files <- getFilePaths(dir, "hg19", "cDNA")
cat("Please check that these are the quant files you want to analyze")
files
file.exists(files)

# assign the project name to datasetName - include which genome
projName <- sampAttributes$projectName
datasetName <- paste0(projName, "_", sampAttributes$humanGenome)
datasetName

# 3.3   get hg19 transcriptome database ----------------------------------- --------------------------------------------

# fetch human genome hg19 build 75 from Ensembl
# if not already loaded
# checker var
is_ensembl_loaded = exists("ensDBHg19")

# if statement to check if ensembl loaded
if (is_ensembl_loaded) {
  print("ENSEMBL database loaded.")
  } else {
  ensDBHg19 <- fetchFromEnsembl("Homo sapiens", 75)
}

# run preImport
txi <- preImport(ensDBHg19, files)
head(txi, 10)

# get coldata for txi
coldata <- getColdata(txi)

# confirm that the directory is correct - should say "hg19" and "cDNA"
# view coldata
coldata

# make DESeqDataSet
myData <- DESeqDataSetFromTximport(txi, 
                                colData = coldata,
                                design = ~ cond)

# 3.4a  visualize hg19-aligned data quality unblinded --------------------- ---------------------------------

# check data quality with non-blinded data!
# non norm transform data
ntd <- normTransform(myData)

# get mean SD plot
saveMeanSdPlot(datasetName, "nonNormedData", ntd)
# dev.off()

# r log transform data, unblinded
rldUnblind <- rlog(myData,
            blind = FALSE)

# make mean SD plot and PCA plot for this transform
saveMeanSdPlot(datasetName, "rLogTransformedData_unblinded", rldUnblind)
savePCAPlot(datasetName, "rLogTransformedData_unblinded", rldUnblind)
# yes this dev.off is necessary
dev.off()

# Next up, variance stabilizing transform
vsdUnblind <- vst(myData,
           fitType = "local",
           blind = FALSE)

saveMeanSdPlot(datasetName, "varStabilizedData_unblinded", vsdUnblind)
# dev.off()
savePCAPlot(datasetName, "varStabilizedData_unblinded", vsdUnblind)
dev.off()

# 3.4b  graph BLINDED data with best transformation for figures ----------- ----

#VSD looks like it fits better
#so make all graphs again with blind = true
vsd <- vst(myData,
           fitType = "local",
           blind = TRUE)

# type of normalization applied
normType <- "blinded_varStabilizedData"

# save all plots needed
saveMeanSdPlot(datasetName, normType, vsd)
savePCAPlot(datasetName, normType, vsd)
dev.off()
makeSampleCorrelationHeatmap(datasetName, normType, vsd)
makeSampleDistributionHeatmap(datasetName, normType, vsd)

plotPCA(vsd,
        intgroup = "cond") +
  theme_classic() +
  scale_colour_discrete(name = "Treatment",
                        labels = c(sampAttributes$refName, sampAttributes$expNames))
  

# 3.5   run DESeq2 -------------------------------------------------------- -----------------------------------------------

# To run DESeq2 on dataset named ddsHg19
dds <- DESeq(myData)

# 3.6   plot normed counts ------------------------------------------------ -----------------------------------------
### runs post-DESeq2
normalized_counts <- counts(dds, normalized = TRUE)
# fill in

# make a filename
filename <- paste0(outputDir, sampAttributes$projectName, "_", sampAttributes$humanGenome, "_normedCounts.png")
filename
#

# save as a .png file
png(file = filename,
    width = 600,
    height = 350)

plotDispEsts(dds)
dev.off()
#

# 3.7   get and plot pairwise results ------------------------------------- -----------------------------------------------

# run getPairwiseLFC on all treatments vs WT-like
pairCompareHuman <- lapply(expNums, getPairwiseLFC, refNum, dds)

# check that df (line two of three for each sample) is labeled correctly
# (pairCompare$[data frame name])

# plot and get dataframe for each individually
plotOutput(pairCompareHuman$`MutuI_IL4-CD40L_Day1`, "_MUTU_IL4-CD40L_Day1", 6, 60, sampAttributes$humanGenome)
cond_2_vs_1_human <- as.data.frame(pairCompareHuman$`MutuI_IL4-CD40L_Day1`)
cond_2_vs_1_human$exon_id <- rownames(cond_2_vs_1_human)

plotOutput(pairCompareHuman$MutuI_IL10_Day1, "_MUTU_IL10_Day1", 6, 60, sampAttributes$humanGenome)
cond_3_vs_1_human <- as.data.frame(pairCompareHuman$MutuI_IL10_Day1)
cond_3_vs_1_human$exon_id <- rownames(cond_3_vs_1_human)

plotOutput(pairCompareHuman$MutuI_IL21_Day1, "_MUTU_IL21_Day1", 6, 60, sampAttributes$humanGenome)
cond_4_vs_1_human <- as.data.frame(pairCompareHuman$MutuI_IL21_Day1)
cond_4_vs_1_human$exon_id <- rownames(cond_4_vs_1_human)



# 3.8   annotate human genes from ensembl --------------------------------- --------------------------------

# get annotations from Ensembl
# dataset, biomart, version 
anns_human <- getAnns("hsapiens_gene_ensembl", "ensembl", 75)
#
# save annotation file
# obj <- makeAnnFiles(as.data.frame(pairCompare$[condition]), anns, paste0(datasetName, "_[condition]", "_anns"))

cond_2_vs_1_human_anns <- makeAnnFiles(as.data.frame(pairCompareHuman$`MutuI_IL4-CD40L_Day1`), anns_human, paste0(datasetName, "_IL4-CD40L_Day1", "_anns"))
cond_3_vs_1_human_anns <- makeAnnFiles(as.data.frame(pairCompareHuman$MutuI_IL10_Day1), anns_human, paste0(datasetName, "_IL10_Day1", "_anns"))
cond_4_vs_1_human_anns <- makeAnnFiles(as.data.frame(pairCompareHuman$MutuI_IL21_Day1), anns_human, paste0(datasetName, "_IL21_Day1", "_anns"))


# 3.9a  get differentially expressed genes -------------------------------- -------------------------------------

# get human DEGs in each condition
MUTU_IL4_CD40L_Day1_DEGs_human <- getDEGs(cond_2_vs_1_human_anns)
MUTU_IL10_Day1_DEGs_human <- getDEGs(cond_3_vs_1_human_anns)
MUTU_IL21_Day1_DEGs_human <- getDEGs(cond_4_vs_1_human_anns)
#

# make a volcano plot with labeled DEGs
plotDEGs(MUTU_IL4_CD40L_Day1_DEGs_human, 4, 6, "hg19")
plotDEGs(MUTU_IL10_Day1_DEGs_human, 4, 6, "hg19")
plotDEGs(MUTU_IL21_Day1_DEGs_human, 4, 6, "hg19")
#

# 4 Process data for Akata-aligned reads ----
# 4.1   placeholder ----
# 4.2   setup env to annotate Akata-aligned transcript reads -------------- ------
# First confirm that section # 3.1 set up env --- has been run
# get akata cDNA transcriptome alignment file paths 
files <- getFilePaths(dir, "akata", "cDNA")
files
file.exists(files)

# 4.3   get Akata transcriptome database ---------------------------------- ----
# copy of ebv index file in both data-analysis and sequencing folders
akataGtfDir <- "C:/Users/NRB37/Dropbox (Partners HealthCare)/Nina/coworkers/"
akata_gtf <- file.path(akataGtfDir, "akata_gtf_Mingxiang_Rui_me.gtf")
file.exists(akata_gtf)

akataDB <- makeTxDbFromGFF(akata_gtf, format = "gtf")

k <- keys(akataDB, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(akataDB, k, "EXONID", "TXNAME")
txi <- tximport(files, type = "salmon", 
                tx2gene = tx2gene)
(txi)
#

# get coldata for txi
coldata <- getColdata(txi)
#

# confirm that the directory is correct - should say "akata" and "cDNA"
# view coldata
coldata

# make DESeqDataSet
myData <- DESeqDataSetFromTximport(txi, 
                                    colData = coldata,
                                    design = ~ cond)
#

# assign the project name to datasetName - include which genome
projName <- sampAttributes$projectName
datasetName <- paste0(projName, "_", sampAttributes$ebvGenome)
datasetName

# 4.4a  visualize akata-aligned data quality unblinded -------------------- ---------------------------------

# check data quality with non-blinded data!
# non norm transform data
ntd <- normTransform(myData)

# get mean SD plot
saveMeanSdPlot(datasetName, "nonNormedData", ntd)
# dev.off()

# r log transform data, unblinded
rldUnblind <- rlog(myData,
                   blind = FALSE,
                   fitType = "local")

# make mean SD plot and PCA plot for this transform
saveMeanSdPlot(datasetName, "rLogTransformedData_unblinded", rldUnblind)
savePCAPlot(datasetName, "rLogTransformedData_unblinded", rldUnblind)
# yes this dev.off is necessary
dev.off()

# Next up, variance stabilizing transform
vsdUnblind <- varianceStabilizingTransformation(myData,
                  fitType = "local",
                  blind = FALSE)

saveMeanSdPlot(datasetName, "varStabilizedData_unblinded", vsdUnblind)
# dev.off()
savePCAPlot(datasetName, "varStabilizedData_unblinded", vsdUnblind)
dev.off()

# 4.4b  graph BLINDED data with best transformation for figures ----------- ----

#VSD looks like it fits better
#so make all graphs again with blind = true
vsd <- varianceStabilizingTransformation(myData,
           fitType = "local",
           blind = TRUE)

# type of normalization applied
normType <- "blinded_varStabilizedData"

# save all plots needed
saveMeanSdPlot(datasetName, normType, vsd)
savePCAPlot(datasetName, normType, vsd)
dev.off()
makeSampleCorrelationHeatmap(datasetName, normType, vsd)
makeSampleDistributionHeatmap(datasetName, normType, vsd)



plotPCA(vsd,
        intgroup = "cond") +
  theme_classic() +
  scale_colour_discrete(name = "Treatment",
                        labels = c(sampAttributes$refName, sampAttributes$expNames))


# 4.5   run DESeq2 -------------------------------------------------------- ------------------------------------------------------------

dds <- DESeq(myData,
             fitType = "local")

# 4.6   get normed counts ----
### runs post-DESeq2
normalized_counts <- counts(dds, normalized = TRUE)
# fill in

# make a filename
filename <- paste0(outputDir, sampAttributes$projectName, "_", sampAttributes$ebvGenome, "_normedCounts.png")
filename
#

# save as a .png file
png(file = filename,
    width = 600,
    height = 350)

plotDispEsts(dds)
dev.off()
#

# 4.7   get pairwise results ---------------------------------------------- --------------------------------------------------

# # make a list of conditions you want to compare
# myList <- c("2","3")
# names(myList) <- c("GM_IL15_Day6", "GM_IL21_Day6")
# #
# run getPairwiseLFC on all treatments vs WT-like
pairCompareEBV <- lapply(expNums, getPairwiseLFC, refNum, dds)

# plot output
plotOutput(pairCompareEBV$`MutuI_IL4-CD40L_Day1`, "_MUTU_IL4_CD40L_Day1", 6, 60, sampAttributes$ebvGenome)
cond_2_vs_1_ebv <- as.data.frame(pairCompareEBV$`MutuI_IL4-CD40L_Day1`)
cond_2_vs_1_ebv$exon_id <- rownames(cond_2_vs_1_ebv)

plotOutput(pairCompareEBV$MutuI_IL10_Day1, "_MUTU_IL10_Day1", 6, 60, sampAttributes$ebvGenome)
cond_3_vs_1_ebv <- as.data.frame(pairCompareEBV$MutuI_IL10_Day1)
cond_3_vs_1_ebv$exon_id <- rownames(cond_3_vs_1_ebv)

plotOutput(pairCompareEBV$MutuI_IL21_Day1, "_MUTU_IL21_Day1", 6, 60, sampAttributes$ebvGenome)
cond_4_vs_1_ebv <- as.data.frame(pairCompareEBV$MutuI_IL21_Day1)
cond_4_vs_1_ebv$exon_id <- rownames(cond_4_vs_1_ebv)

# 4.8a  organize the akata gtf to extract gene IDs ---------------------------- 


# get the gene IDs - must repeat once per condition
cond_2_vs_1_ebv_anns <- merge(x = cond_2_vs_1_ebv, y = tx2gene, by.x = "exon_id", by.y = "EXONID")
#save as .tsv file
write_tsv(cond_2_vs_1_ebv_anns, file = paste0(outputDir, datasetName, "_IL4_CD40L_Day1", "_anns.tsv"))

#
# get the gene IDs - must repeat once per condition
cond_3_vs_1_ebv_anns <- merge(x = cond_3_vs_1_ebv, y = tx2gene, by.x = "exon_id", by.y = "EXONID")
#save as .tsv file
write_tsv(cond_3_vs_1_ebv_anns, file = paste0(outputDir, datasetName, "_IL10_Day1", "_anns.tsv"))

#
# get the gene IDs - must repeat once per condition
cond_4_vs_1_ebv_anns <- merge(x = cond_4_vs_1_ebv, y = tx2gene, by.x = "exon_id", by.y = "EXONID")
#save as .tsv file
write_tsv(cond_4_vs_1_ebv_anns, file = paste0(outputDir, datasetName, "_IL21_Day1", "_anns.tsv"))

# rows below are for annotating, which I'm skippin for now 
#myAnns <- unique(myAnns[which(complete.cases(myAnns)),])
#ebvAnns <- myAnns

#myAnns <- unique(myAnns[which(complete.cases(myAnns)),])
#ebvAnns <- myAnns
#

# get the gene IDs - must repeat once per condition
#anns <- merge(x = p, y = tx2gene, by.x = "exon_id", by.y = "EXONID")

# #save as .tsv file
# myAnns <- write_tsv(anns, file = paste0(outputDir, datasetName, names(myList)[4], "_anns.tsv"))
# myAnns <- unique(myAnns[which(complete.cases(myAnns)),])
# ebvAnns <- myAnns
# #
# 


# 4.9   get and visualize EBV DEGs----
# get list of diff expressed genes
#gm_IL15_Day6_DEG <- getDEGs(myAnns, cond_2_vs_1$n)
#

# make a volcano plot with labeled DEGs
#plotDEGs(sgJunB_2_DEG, 4, 6, "akata")
#
