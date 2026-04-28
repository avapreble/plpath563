library(ggtree)    
library(ape)         
library(ggplot2)
library(treeio)

# load tree
iqtree <- read.tree("data/alignment_results/MAFFT_aligned.fasta.treefile") 


# tip labels, colors, formating
colorRegions <- data.frame(label = iqtree$tip.label)
colorRegions$Region <- ifelse(grepl("C.MZ", colorRegions$label), "Mozambique", 
                              ifelse(grepl("TH", colorRegions$label), "Thailand", "USA"))

colorRegions$new_label <- ave(colorRegions$Region, colorRegions$Region,
                              FUN = function(x) paste0(x, "_", seq_along(x)))

iqtreePlot <- ggtree(iqtree) %<+% colorRegions +
  geom_tiplab(aes(label = new_label, color = Region), size = 3) +
  geom_nodelab(aes(label = label), size = 2, hjust = -0.5, vjust = -0.3) +
  xlim(0, max(ggtree(iqtree)$data$x) * 1.4)


iqtreePlotx <- iqtreePlot + theme_tree2() + 
  theme(legend.position = "left", text = element_text(size = 10)) + 
  guides(color = guide_legend(override.aes = list(label = "•", size = 10))) 

iqtreePlotx

# save image
ggsave("figures/iqtree1.png", plot = iqtreePlotx, width = 8, height = 6, dpi = 300)





