#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD <m.gruenstaeudl@fu-berlin.de>"
#copyright = "Copyright (C) 2017 Michael Gruenstaeudl"
#version = "2017.02.27.1600"

# ANALYSIS TITLE
analysis_title = paste( "Phylogenetic Tree with Highest Likelihood Score\n",
                        "Software used: RAxML v.8.2.9\n",
                        "Data partitioning: 2 distinct data partitions with joint branchlength estimation\n",
                        "Nucleotide substitution model: GTRGAMMAI; ",
                        "Resampling support: 1000 rapid bootstrap runs\n",
                        "Tree rooting: none / unrooted\n",
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

####################################
# STEP 1. Specify in- and outfiles #
####################################
# SPECIFYING INFILE
inFile = tk_choose.files(caption='Select .tre file')

# SPECIFYING OUTFILES
outDir = dirname(inFile)
outFn = paste(file_path_sans_ext(inFile), '_', sep='')
outFile = paste(outFn, '_viz.svg', sep='')

# SPECIFYING OUTGROUP
my_outgroup = "Amborella_trichopoda_NC_005086"

###############################
# STEP 2. Load and root trees #
###############################

# LOAD TREE
tree <- read.raxml(inFile)

# INFER NODE NUMBER OUT OUTGROUP

# ROOT TREE
#tree_MCC@phylo <- ape::root(tree_MCC@phylo, node = X, edgelabel=TRUE)
tree@phylo <- ape::root(tree@phylo, outgroup=my_outgroup, edgelabel=TRUE)

###########################
# STEP 3. Construct Trees #
###########################

# CONSTRUCT BASE CLADOGRAM
p1 <- ggtree(tree, branch.length="none", size=0.75)
# RE-FORMAT TIP LABELS
p1_lbls_new <- gsub("_", " ", p1$data$label[which(p1$data$isTip==TRUE)])
# ADD TIP LABELS
p1 <- p1 + geom_tiplab(label = p1_lbls_new, fontface="italic", offset=0.2)
# RE-FORMAT BOOTSTRAP VALUES
d <- p1$data # Get from phylo, not from p1
d <- d[!d$isTip,]
d$bootstrap <- round(as.numeric(d$bootstrap), digits=2)
p1_bootstrap <- d[d$bootstrap > 50,]
# ADD BOOTSTRAP VALUES
p1 <- p1 + geom_text(data=p1_bootstrap, aes(label=bootstrap), hjust=1.25, vjust=-0.5)
# TFL would be more efficient, but would not have any rounding
# p1 + geom_text2(aes(label=label, subset = !is.na(as.numeric(posterior)) & as.numeric(posterior) > 80))
# FORMAT TREE
p1 <- p1 + theme(plot.margin=unit(c(1,1,1,1),"cm")) +
        geom_rootpoint() +                                                      # Add a dot to indicate the rootpoint
        geom_treescale(x=0, y=0, width=0.01, color='white') +                   # To adjust height compared to phylogram plot
        ggtitle('Best ML tree as cladogram\nBootstrap values above branches') +
        ggplot2::xlim(0, 40)                                                    # Important for long tip labels (i.e., long taxon names)

# CONSTRUCT BASE PHYLOGRAM
p2 <- ggtree(tree, size=0.75)
# RE-FORMAT TIP LABELS
p2_lbls_new <- gsub("_", " ", p2$data$label[which(p2$data$isTip==TRUE)])
# ADD TIP LABELS
p2 <- p2 + geom_tiplab(label = p2_lbls_new, fontface="italic", offset=0.001)
# ADD SCALEBAR AND FORMAT TREE
p2 <- p2 + theme(plot.margin=unit(c(1,1,1,1),"cm")) +
        geom_rootpoint() +                                                      # Add a dot to indicate the rootpoint
        geom_treescale(x=0, y=0, width=0.01) +
        ggtitle('Best ML tree') +
        #ggplot2::xlim(0, 0.5)                                                   # Important for long tip labels (i.e., long taxon names)
        ggplot2::xlim(0, 0.6)                                                   # Important for long tip labels (i.e., long taxon names)

##############################################
# STEP 4. Construct final plot, save to file #
##############################################

# CONSTRUCT MULTIPLOT AND SAVE
svglite(outFile, width=25, height=25, standalone=TRUE)
    grid.arrange(p1, p2, top=textGrob(analysis_title),
                 layout_matrix = matrix(c(1,2), ncol=2, byrow=TRUE),
                 widths=c(0.3,0.7))
    #multiplot(p1, p2, ncol=2, widths=c(0.3,0.7)) # Multiplot is a cheap wrapper around grid.arrange
dev.off()

