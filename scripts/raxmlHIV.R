# running RAxML on alignment generated with MAFFT (15 taxa) 

#loading packages
library(ggtree)    
library(ape)         
library(ggplot2)
library(treeio)


# loading RAxML tree
raxmltree <- read.tree("results/mafft_support.raxml.support")
raxmltree <- ladderize(raxmltree)


# tip labels, colors, formating
colorRegionsraxml <- data.frame(label = iqtree$tip.label)
colorRegionsraxml$Region <- ifelse(grepl("C.MZ", colorRegions$label), "Mozambique", 
                              ifelse(grepl("TH", colorRegions$label), "Thailand", "USA"))

colorRegionsraxml$new_label <- ave(colorRegions$Region, colorRegions$Region,
                              FUN = function(x) paste0(x, "_", seq_along(x)))

#plot tree
raxmltreeplot <- ggtree(raxmltree) %<+% colorRegionsraxml + 
  geom_tiplab(aes(label = new_label, color = Region), size = 3) +
  xlim(0, max(ggtree(raxmltree) $data$x) * 1.4) +
  geom_text2(aes(label = label, subset = !isTip), 
           size = 2.5, hjust = -2.5) 


raxmltreex <- raxmltreeplot + 
  theme_tree2() + 
  theme(legend.position = "left", text = element_text(size = 10)) + guides(color = guide_legend(override.aes = list(label = "•", size = 10))) 

raxmltreex

# save image
ggsave("figures/raxmltree1.png", plot = raxmltreex, width = 8, height = 6, dpi = 300)







