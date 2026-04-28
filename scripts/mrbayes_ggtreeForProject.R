# I did this code from scratch originally and am going to use it to build the other ones faster

# load packages
library(ggtree)    
library(ape)         
library(ggplot2)
library(treeio)



# MrBayes - read in tree
MBtree <- read.nexus("data/alignment_results/MAFFT_aligned.nex.con.tre")
MBtreex <- read.mrbayes("data/alignment_results/MAFFT_aligned.nex.con.tre")

# base tree
plot(MBtree)
nodelabels()

# tip labels
colorRegions <- data.frame(label = MBtree$tip.label)
colorRegions$Region <- ifelse(grepl("C.MZ", colorRegions$label), "Mozambique", 
                              ifelse(grepl("TH", colorRegions$label), "Thailand", "USA"))

colorRegions$new_label <- ave(colorRegions$Region, colorRegions$Region,
                              FUN = function(x) paste0(x, "_", seq_along(x)))

# aesthetics

mbTreePlot <- ggtree(MBtreex) + geom_tiplab()
mbTreePlot <- ggtree(MBtreex) %<+% colorRegions + 
  
  geom_tiplab(aes(label = new_label, color = Region), size = 3) + 
  
  geom_text2(aes(label = sprintf("%.2f", as.numeric(prob)), subset = !isTip), 
             size = 2.5, hjust = -2.5) +
  
  xlim(0, max(ggtree(MBtreex)$data$x) * 1.3) + 
  
  
  theme_tree2() + 
  theme(legend.position = "left", text = element_text(size = 10)) + 
  guides(color = guide_legend(override.aes = list(label = "~", size = 10))) 


# save image
ggsave("figures/mbtree1.png", plot = mbTreePlot, width = 8, height = 6, dpi = 300)
  










