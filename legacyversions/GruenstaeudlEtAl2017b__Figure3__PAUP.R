#!/usr/bin/R
#author = "Michael Gruenstaeudl, PhD <m.gruenstaeudl@fu-berlin.de>"
#copyright = "Copyright (C) 2017 Michael Gruenstaeudl"
#version = "2017.02.21.1600"

# ANALYSIS TITLE
analysis_title = paste( "50% Majority Rule Consensus Tree of All Most Parsimonious Trees\n",
                        "Software used: PAUP 4.0b10\n",
                        #"Data partitioning: 2 distinct data partitions with joint branchlength estimation\n",
                        "Resampling support: 100 jackknife replicates\n",
                        "Tree rooting: none / unrooted\n",
                        "Date: ", Sys.time(),
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

# LOAD ALIGNMENTS
library(ape)
tree <- read.nexus(inFile)

###########################
# STEP 2. Construct Trees #
###########################

# CONSTRUCT BASE CLADOGRAM
p1 <- ggtree(tree, branch.length="none", size=0.75)
# RE-FORMAT AND ADD TIP LABELS
p1_lbls_new <- gsub("_", " ", p1$data$label[which(p1$data$isTip==TRUE)])
p1 <- p1 + geom_tiplab(label = p1_lbls_new, fontface="italic", offset=0.2)
# RE-FORMAT AND ADD POSTERIOR VALUES
d <- p1$data # Get from phylo, not from p1
d <- d[!d$isTip,]
d$label <- round(as.numeric(d$label), digits=2)
p1_label <- d[d$label > 50,]
p1 <- p1 + geom_text(data=p1_label, aes(label=label), hjust=1.25, vjust=-0.5)
# TFL would be more efficient, but would not have any rounding
# p1 + geom_text2(aes(label=label, subset = !is.na(as.numeric(posterior)) & as.numeric(posterior) > 80))
# ADD BOOTSTRAP VALUES AND FORMAT TREE
p1 <- p1 + theme(plot.margin=unit(c(1,1,1,1),"cm")) +
        #geom_text(aes(label=bootstrap), hjust=1.25, vjust=-0.5) +
        geom_treescale(x=0, y=0, width=0.01, color='white') +                   # To adjust height compared to phylogram plot
        ggtitle('50% Majority Rule Consensus Tree as cladogram\nJackknife values above branches') +
        ggplot2::xlim(0, 40)                                                    # Important for long tip labels (i.e., long taxon names)


# CONSTRUCT BASE PHYLOGRAM
p2 <- ggtree(tree, size=0.75)
# RE-FORMAT AND ADD TIP LABELS
p2_lbls_new <- gsub("_", " ", p2$data$label[which(p2$data$isTip==TRUE)])
p2 <- p2 + geom_tiplab(label = p2_lbls_new, fontface="italic", offset=0.001)
# ADD SCALEBAR AND FORMAT TREE
p2 <- p2 + theme(plot.margin=unit(c(1,1,1,1),"cm")) +
        #geom_text(aes(label=bootstrap), hjust=1.25, vjust=-0.5) + 
        geom_treescale(x=0, y=0, width=0.01) +
        ggtitle('50% Majority Rule Consensus Tree') +
        ggplot2::xlim(0, 20)                                                  # Important for long tip labels (i.e., long taxon names)

##############################################
# STEP 3. Construct final plot, save to file #
##############################################

# CONSTRUCT MULTIPLOT AND SAVE
svglite(outFile, width=25, height=25, standalone=TRUE)
    grid.arrange(p2, p1, top=textGrob(analysis_title),
                 layout_matrix = matrix(c(1,2), ncol=2, byrow=TRUE),
                 widths=c(0.3,0.7))
    #multiplot(p1, p2, ncol=2, widths=c(0.3,0.7)) # Multiplot is a cheap wrapper around grid.arrange
dev.off()

