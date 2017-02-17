#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD <m.gruenstaeudl@fu-berlin.de>"
#copyright = "Copyright (C) 2017 Michael Gruenstaeudl"
#version = "2017.02.17.1300"

# ANALYSIS TITLE
title_MCC = paste("Maximum Clade Credibility Tree of Posterior Tree Distribution\n",
                        "Software used: MrBayes v.3.2.5; ",
                        "Data partitioning: 1 distinct model/data partition with joint branchlength optimization\n",
                        "Nucleotide substitution model: GTR+I+G; ",
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
inFile_MCC = tk_choose.files(caption='Select file of Maximum Clade Credibility tree')

# SPECIFYING OUTFILES
outDir = dirname(inFile_MCC)
outFn = paste(file_path_sans_ext(inFile_MCC), '_', sep='')
outFile = paste(outFn, '_viz.svg', sep='')

# REMOVING NEGATIVE BRANCH LENGTHS - IMPORTANT B/C TREEANNOTATOR SOMETIMES GENERATES THEM
inFile_MCC_adj = paste(file_path_sans_ext(inFile_MCC), '_noNegBrLens.tre', sep='')
# Actual command looks like: sed -e 's,:-[0-9\.]\+,:0.0,g' in.tree > out.tree
cmd = paste("sed -e 's,:-[0-9\\.]\\+,:0.0,g'", inFile_MCC, ">", inFile_MCC_adj)
system(cmd)

# LOAD ALIGNMENTS
tree_MCC <- ggtree::read.beast(inFile_MCC_adj) # The function read.beast exists in several packages, but often generates different output.

###########################
# STEP 2. Construct Trees #
###########################

# DISPLAY MCC AS CLADOGRAM WITH PP-VALUES
p1 <- ggtree(tree_MCC, branch.length="none", size=0.75)
# RE-FORMAT AND ADD TIP LABELS
p1_tiplabels <- gsub("_", " ", p1$data$label[which(p1$data$isTip==TRUE)])
p1 <- p1 + geom_tiplab(label = p1_tiplabels, fontface="italic", offset=0.2)
# RE-FORMAT AND ADD POSTERIOR VALUES
d <- p1$data # Get from phylo, not from p1
d <- d[!d$isTip,]
d$posterior <- round(as.numeric(d$posterior), digits=2)
p1_posterior <- d[d$posterior > 0.5,]
p1 <- p1 + geom_text(data=p1_posterior, aes(label=posterior), hjust=1.25, vjust=-0.5)
# TFL would be more efficient, but would not have any rounding
# p1 + geom_text2(aes(label=label, subset = !is.na(as.numeric(posterior)) & as.numeric(posterior) > 80))
# FORMAT TREE
p1 <- p1 + theme(plot.margin=unit(c(1,1,1,1),"cm")) +
        geom_treescale(x=0, y=0, width=0.01, color='white') +                   # To adjust height compared to phylogram plot
        ggtitle('Maximum Clade Credibility tree as cladogram\nPosterior probability values greater than 0.5 are displayed above branches') +
        ggplot2::xlim(0, 30)                                                    # Important for long tip labels (i.e., long taxon names)

# DISPLAY MCC AS PHYLOGRAM
p2 <- ggtree(tree_MCC, size=0.75)
# RE-FORMAT AND ADD TIP LABELS
p2_tiplabels <- gsub("_", " ", p2$data$label[which(p2$data$isTip==TRUE)])
p2 <- p2 + geom_tiplab(label = p2_tiplabels, fontface="italic", offset=0.001)
# ADD SCALEBAR AND FORMAT TREE
p2 <- p2 + theme(plot.margin=unit(c(1,1,1,1),"cm")) +
        geom_treescale(x=0, y=0, width=0.01) +
        ggtitle('Maximum Clade Credibility tree\nPosterior probability values above branches') +
        ggplot2::xlim(0, 0.40)                                                  # Important for long tip labels (i.e., long taxon names)

##############################################
# STEP 3. Construct final plot, save to file #
##############################################

# CONSTRUCT MULTIPLOT AND SAVE
svglite(outFile, width=25, height=25, standalone=TRUE)
    grid.arrange(p1, p2, top=textGrob(title_MCC),
                 layout_matrix = matrix(c(1,2), ncol=2, byrow=TRUE),
                 widths=c(0.3,0.7))
    #multiplot(p1, p2, ncol=2, widths=c(0.3,0.7)) # Multiplot is a cheap wrapper around grid.arrange
dev.off()

