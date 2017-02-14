#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD <mi.gruenstaeudl@gmail.com>"
#copyright = "Copyright (C) 2017 Michael Gruenstaeudl"
#version = "2017.02.14.1800"

# ANALYSIS TITLE
analysis_title = paste( "Phylogenetic Tree with Highest Likelihood Score\n",
                        "Software used: RAxML v.8.2.9; ",
                        "Data partitioning: 1 distinct model/data partition with joint branchlength optimization\n",
                        "Nucleotide substitution model: GTRGAMMAI; ",
                        "Resampling support: 1000 rapid bootstrap runs\n",
                        "Tree rooting: none / unrooted\n",
                        sep ='')

#############
# Libraries #
#############
library(ggtree)
library(grid)
library(gridExtra)
library(svglite) # For improved svg drivers
library(tcltk)   # For dialog boxes
library(tools)   # For function 'file_path_sans_ext'

####################################
# STEP 1. Specify in- and outfiles #
####################################
# SPECIFYING INFILES
inFile = tk_choose.files(caption='Select .tre file')

# SPECIFYING OUTFILES
outDir = dirname(inFile)
outFn = paste(file_path_sans_ext(inFile), '_', sep='')
outFile = paste(outFn, '_viz.svg', sep='')

##############################
# STEP 2. Adjust taxon names #
##############################

# LOAD ALIGNMENTS
tree <- read.raxml(inFile)

## FOR FUTURE
#lbls_old = get.tree(tree)$tip.label
#l1 = gsub("_", " ", lbls)
##l2 = lapply(l1, function(x) paste("paste(italic('", x, "'))", sep=""))
##lbls_new = lapply(l2, function(x) paste("'", x, "'", sep=""))
##d = as.data.frame(cbind(label=lbls_old, label2=lbls_new))
#d = as.data.frame(cbind(label=lbls_old, label2=l))
#ggtree(tree) %<+% d + geom_tiplab(aes(label=label2))

###########################
# STEP 3. Construct Trees #
###########################

# CONSTRUCT CLADOGRAM WITH BOOTSTRAP VALUES
p1 <- ggtree(tree, branch.length="none", size=0.75) +
        theme(plot.margin=unit(c(1,1,1,1),"cm")) +
        #geom_tiplab(offset=0.2, fontface="italic") +
        #geom_label(aes(label=bootstrap), hjust=1.5, vjust=-0.5) +
        geom_text(aes(label=bootstrap), hjust=1.25, vjust=-0.5) +
        geom_treescale(x=0, y=0, width=0.01, color='white') +                   # To adjust height compared to phylogram plot
        ggtitle('Cladogram\nBootstrap values above branches') +
        ggplot2::xlim(0, 30)                                                    # Important for long tip labels (i.e., long taxon names)

# Properly formatting labels
p1 <- p1 + geom_tiplab(label = gsub("_", " ", p1$data$label[which(p1$data$isTip==TRUE)]),
                       fontface="italic",
                       offset=0.2)

# CONSTRUCT PHYLOGRAM WITH SCALEBAR
p2 <- ggtree(tree, size=0.75) +
        theme(plot.margin=unit(c(1,1,1,1),"cm")) +
        #geom_tiplab(offset=0.001, fontface="italic") +
        geom_treescale(x=0, y=0, width=0.01) +
        ggtitle('Phylogram') +
        ggplot2::xlim(0, 0.5)                                                    # Important for long tip labels (i.e., long taxon names)

# Properly formatting labels
p2 <- p2 + geom_tiplab(label = gsub("_", " ", p2$data$label[which(p2$data$isTip==TRUE)]),
                       fontface="italic",
                       offset=0.001)

##############################################
# STEP 3. Construct final plot, save to file #
##############################################

# CONSTRUCT MULTIPLOT AND SAVE
svglite(outFile, width=25, height=25, standalone=TRUE)
    grid.arrange(p1, p2, top=textGrob(analysis_title),
                 layout_matrix = matrix(c(1,2), ncol=2, byrow=TRUE),
                 widths=c(0.3,0.7))
    #multiplot(p1, p2, ncol=2, widths=c(0.3,0.7)) # Multiplot is a cheap wrapper around grid.arrange
dev.off()

