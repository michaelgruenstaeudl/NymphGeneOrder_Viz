#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD <m.gruenstaeudl@fu-berlin.de>"
#copyright = "Copyright (C) 2017 Michael Gruenstaeudl"
#version = "2017.05.19.1700"

# ANALYSIS TITLE
analysis_title = paste( "Phylogenetic Tree with Highest Likelihood Score\n",
                        "Software used: RAxML v.8.2.9 (for Bayesian inference: Bayesian: MrBayes v.3.2.5)\n",
                        
                        #"Data partitioning: None; the entire dataset is analyzed as a single data partition\n",
                        #"Data partitioning: By gene; branch length optimization is linked across partitions\n",
                        #"Data partitioning: By codon; branch length optimization is linked across partitions\n",
                        "Data partitioning: Best partitioning scheme as identified by PartitionFinder2; branch length optimization is linked across partitions\n",
                        "Nucleotide substitution model: ML analysis: GTRGAMMAI (for Bayesian inference: best-fitting model per partition)\n",
                        "Resampling support: 1000 rapid bootstrap runs (for Bayesian inference: 1000 post-burnin MCMC runs)\n",
                        #"Tree rooting: none / unrooted\n",
                        "Date: ", Sys.time(),
                        sep ='')

#############
# Libraries #
#############
library(ggtree)
library(grid)
library(gridExtra)
library(svglite)    # For improved svg drivers
library(tcltk)      # For dialog boxes
library(tools)      # For function 'file_path_sans_ext'
library(treeio)     # For function 'as.treedata'
#library(random)     # To generate random number


####################################
# STEP 1. Specify in- and outfiles #
####################################
# SPECIFYING INFILES
inFile_plotTree = tk_choose.files(caption='Select .tre file to be plotted')
inFile_suppVals = tk_choose.files(caption='Select .tre file for which support values are calculated')

# SPECIFYING OUTFILES
outDir = dirname(inFile_plotTree)
outFn = paste(file_path_sans_ext(inFile_plotTree), '_', sep='')
outFile_suppVals = paste(outFn, '_suppVals.tre', sep='')
outFile_plotTree = paste(outFn, '_viz.svg', sep='')

Rplot_cladogram = paste(outFn, '_cladogramPlot.Rdata', sep='')
Rtree_cladogram = paste(outFn, '_cladogramTreeObj.Rdata', sep='')
Rplot_phylogram = paste(outFn, '_phylogramPlot.Rdata', sep='')
Rtree_phylogram = paste(outFn, '_phylogramTreeObj.Rdata', sep='')

# SPECIFYING OUTGROUP
my_outgroup = "Amborella_trichopoda_NC_005086"

###################################################################
# STEP SPECIAL. Summarize posterior probabilities on best ML tree #
###################################################################
cmd1 = paste("rm ", outFile_suppVals)
system(cmd1)

cmd2 = paste("python2.7 /home/michael_science/.local/bin/sumtrees.py",
                            "-o", outFile_suppVals,
                            "-t", inFile_plotTree,
                                  inFile_suppVals)
system(cmd2)

###############################
# STEP 2. Load and root trees #
###############################

# LOAD TREE TO BE PLOTTED
#tree <- ggtree::read.raxml(inFile_plotTree)    # ABANDONING FUNCTION read.raxml B/C TOO BUGGY
tree <- ape::read.tree(inFile_plotTree)         # Loading RAxML-file ".bipartitions" via APE-function

# LOAD TREE WITH THE ADDITIONAL SUPPORT VALUES
tree_suppVals <- ape::read.nexus(outFile_suppVals)
#nl <- tree_suppVals$node.label                 # COMMENTED OUT B/C ABANDONING FUNCTION read.raxml
#suppVals = as.treedata(tree_suppVals, boot=nl) # COMMENTED OUT B/C ABANDONING FUNCTION read.raxml

# ROOT TREES
##tree_MCC@phylo <- ape::root(tree_MCC@phylo, node = X, edgelabel=TRUE)             # COMMENTED OUT B/C ABANDONING FUNCTION read.raxml
#tree@phylo <- ape::root(tree@phylo, outgroup=my_outgroup, edgelabel=TRUE)          # COMMENTED OUT B/C ABANDONING FUNCTION read.raxml
#suppVals@phylo <- ape::root(suppVals@phylo, outgroup=my_outgroup, edgelabel=TRUE)  # COMMENTED OUT B/C ABANDONING FUNCTION read.raxml
tree <- ape::root(tree, outgroup=my_outgroup, edgelabel=TRUE)
tree_suppVals <- ape::root(tree_suppVals, outgroup=my_outgroup, edgelabel=TRUE)

# EXTRACT SUPPORT VALS
nl <- tree_suppVals$node.label
suppVals = as.treedata(tree_suppVals, boot=nl)

################################
# STEP 3a. CONSTRUCT CLADOGRAM #
################################
# CONSTRUCT BASE CLADOGRAM
p1 <- ggtree(tree, branch.length="none", size=0.75, ladderize=TRUE)
# RE-FORMAT TIP LABELS
p1_lbls_new <- gsub("_", " ", p1$data$label[which(p1$data$isTip==TRUE)])
# ADD TIP LABELS
p1 <- p1 + geom_tiplab(label = p1_lbls_new, fontface="italic", offset=0.2)
# RE-FORMAT BOOTSTRAP VALUES
d <- p1$data # Get from phylo, not from p1
d <- d[!d$isTip,]

if ((!("bootstrap" %in% colnames(d))) & ("label" %in% colnames(d))) {
    colnames(d)[which("label"==colnames(d))] <- "bootstrap"
} else {
    stop('ERROR: Cannot find column label either `bootstrap` or `label` in tree.')
}

d$bootstrap <- round(as.numeric(d$bootstrap), digits=2)
p1_bootstrap <- d[d$bootstrap > 50,]
# ADD BOOTSTRAP VALUES
p1 <- p1 + geom_text(data=p1_bootstrap, aes(label=bootstrap), hjust=1.25, vjust=-0.5)
# TFL would be more efficient, but would not have any rounding
# p1 + geom_text2(aes(label=label, subset = !is.na(as.numeric(posterior)) & as.numeric(posterior) > 80))

# Convert POSTERIOR PROBABILITIES FROM SUPPVALS TREE
pp <- suppVals@data
pp <- pp[which(pp[,"node"] %in% p1_bootstrap[,"node"]),"bootstrap"]
ppV <- round(as.numeric(as.character(pp)), digits=2)
new_bootstrap <- p1_bootstrap
#new_bootstrap[,"bootstrap"] <- c(NA, ppV)
new_bootstrap[,"bootstrap"] <- c(NA, NA, 1.00, ppV)
# GET PP VALUES > 0.5
new_bootstrap <- new_bootstrap[new_bootstrap$bootstrap > 0.5,]

p1 <- p1 + geom_text(data=new_bootstrap, aes(label=bootstrap), hjust=2.25, vjust=-0.5, fontface="italic")

# FORMAT TREE
p1 <- p1 + theme(plot.margin=unit(c(1,1,1,1),"cm")) +
        geom_rootpoint() +                                                      # Add a dot to indicate the rootpoint
        geom_treescale(x=0, y=0, width=0.01, color='white') +                   # To adjust height compared to phylogram plot
        ggtitle('Best ML tree as cladogram\nBootstrap values above branches') +
        ggplot2::xlim(0, 15)                                                    # Important for long tip labels (i.e., long taxon names)

# SAVING AS PLOT AND TREE
save(p1, file=Rplot_cladogram)
tree_p1 = as.phylo(p1)
save(tree_p1, file=Rtree_cladogram)

################################
# STEP 3b. CONSTRUCT PHYLOGRAM #
################################
# CONSTRUCT BASE PHYLOGRAM
p2 <- ggtree(tree, size=0.75, ladderize=TRUE)
# RE-FORMAT AND ADD TIP LABELS
p2_lbls_new <- gsub("_", " ", p2$data$label[which(p2$data$isTip==TRUE)])
p2 <- p2 + geom_tiplab(label = p2_lbls_new, fontface="italic", offset=0.001)
# ADD SCALEBAR AND FORMAT TREE
p2 <- p2 + theme(plot.margin=unit(c(1,1,1,1),"cm")) +
        geom_rootpoint() +                                                      # Add a dot to indicate the rootpoint
        geom_treescale(x=0, y=0, width=0.01) +
        ggtitle('Best ML tree') +
        #ggplot2::xlim(0, 0.5)                                                   # Important for long tip labels (i.e., long taxon names)
        ggplot2::xlim(0, 0.23)                                                  # Important for long tip labels (i.e., long taxon names)

# SAVING AS PLOT AND TREE
save(p2, file=Rplot_phylogram)
tree_p2 = as.phylo(p2)
save(tree_p2, file=Rtree_phylogram)

##############################################
# STEP 4. Construct final plot, save to file #
##############################################

# CONSTRUCT MULTIPLOT AND SAVE
svglite(outFile_plotTree, width=20, height=10, standalone=TRUE)
    grid.arrange(p1, p2, top=textGrob(analysis_title),
                 layout_matrix = matrix(c(1,2), ncol=2, byrow=TRUE),
                 widths=c(0.3,0.7))
    #multiplot(p1, p2, ncol=2, widths=c(0.3,0.7)) # Multiplot is a cheap wrapper around grid.arrange
dev.off()
