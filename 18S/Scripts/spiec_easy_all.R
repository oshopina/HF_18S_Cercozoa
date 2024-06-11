library(SpiecEasi)
library(igraph)
library(stringr)

env = read.csv('18S/Data/env_16S.csv')
rownames(env) = env$Sample
bacteria = read.csv('18S/Data/otu_16S_4500.RA.csv')
fungi = read.csv('18S/Data/otu_ITS_7150.RA.csv')
rownames(bacteria) = bacteria$X
rownames(fungi) = fungi$X
euks = read.csv2('18S/Data/euk_no_fungi_2300.csv', row.names = 1)

bacteria = bacteria[,-1]
fungi = fungi[,-1]
euks = euks[,1:98]

uni_names = str_extract(colnames(bacteria), '[0-9]+') |> as.numeric()
colnames(bacteria) = paste0('HF', uni_names)
uni_names = str_extract(colnames(fungi), '[0-9]+') |> as.numeric()
colnames(fungi) = paste0('HF', uni_names)

euks = euks[,colnames(euks) %in% colnames(bacteria)]
bacteria = bacteria[,colnames(euks)]
fungi = fungi[,colnames(euks)]

rownames(bacteria) = paste0("B.", rownames(bacteria))
rownames(fungi) = paste0("F.", rownames(fungi))
rownames(euks) = paste0('E.', rownames(euks))

################################ SpiecEasi #####################################

# Define a list of groups and their criteria

groups <- list(
  list(clust = 1, name = "g4"),
  list(clust = 2, name = "t45"),
  list(clust = 3, name = "g6"),
  list(clust = 4, name = "t65"),
  list(clust = 5, name = "g7")
)

bacteria_group_tables = list()
fungi_group_tables = list()
euks_group_tables = list()
env_group_tables = list()

#Loop through the groups
for (i in 1:length(groups)) {
  group = groups[[i]]
  
  env_group <- env[env$OG_clustering == group$clust | env$Tipping_points == group$clust,]
  
  bacteria_group <- bacteria[, colnames(bacteria) %in% env_group$Sample]
  bacteria_group <- bacteria_group[which(apply(bacteria_group, 1, max) > 2),] |> t()
  non_zero_counts <- apply(bacteria_group, 2, function(x) sum(x != 0))
  bacteria_group = bacteria_group[, non_zero_counts > 3]
  
  fungi_group <- fungi[, colnames(fungi) %in% env_group$Sample]
  fungi_group <- fungi_group[which(apply(fungi_group, 1, max) > 3),] |> t()
  non_zero_counts <- apply(fungi_group, 2, function(x) sum(x != 0))
  fungi_group = fungi_group[, non_zero_counts > 3]
  
  euks_group <- euks[, colnames(euks) %in% env_group$Sample]
  euks_group <- euks_group[which(apply(euks_group, 1, max) > 3),] |> t()
  non_zero_counts <- apply(euks_group, 2, function(x) sum(x != 0))
  euks_group = euks_group[, non_zero_counts > 3]
  
  bacteria_group_tables[[i]] = bacteria_group
  fungi_group_tables[[i]] = fungi_group
  euks_group_tables[[i]] = euks_group
  env_group_tables[[i]] = env_group
  
  ## Perform the network analysis and save the result (uncomment this section when needed)
  spieceasi <- spiec.easi(list(bacteria_group, fungi_group, euks_group),
                          method = 'mb', lambda.min.ratio = 1e-2, nlambda = 19,
                          icov.select.params = list(rep.num = 50, ncores = 10))

  saveRDS(spieceasi, paste0("network_", group$name, ".rds"))
}


############################ Network preparation ###############################
matrixes = list()
nets =  list()
nets_dist = list()
nets_abs = list()


for (i in 1:length(groups)) {
  group = groups[[i]]
  bacteria_group = bacteria_group_tables[[i]]
  fungi_group = fungi_group_tables[[i]]
  
  spiec = readRDS(paste0('Data/SpiecEasi_3_groups/network_', group$name, '.rds'))
  matrix = symBeta(getOptBeta(spiec), mode='maxabs') |> as.matrix()
  colnames(matrix) = c(colnames(bacteria_group), colnames(fungi_group))
  rownames(matrix) = c(colnames(bacteria_group), colnames(fungi_group))
  
  net = graph.adjacency(matrix,mode = "undirected", weighted = TRUE, diag = FALSE)
  V(net)$name = colnames(matrix)
  
  net_dist <- net
  E(net_dist)$weight <- 1 - abs(E(net_dist)$weight)
  
  net_abs <- net
  E(net_abs)$weight <- abs(E(net_abs)$weight)
  
  E(net)[weight>0]$color <- "#FF37D5"
  E(net)[weight<0]$color <- "#4335C9"
  
  V(net)$color <- ifelse(substr(V(net)$name, 1, 1) == "B", "black", "white")
  
  matrixes[[i]] = matrix
  nets[[i]] =  net
  nets_dist[[i]] = net_dist
  nets_abs[[i]] = net_abs
  
}

############################ Network graph #####################################
library(grid)
library(gridGraphics)
library(ggplotify)
library(patchwork)
library(ggplot2)

l = layout_with_mds(nets_dist[[1]])

plot(nets[[1]], edge.color = E(nets[[1]])$colors, vertex.size = 2.5, edge.curved = 1,
     vertex.color = V(nets[[1]])$colors, vertex.label = "", rescale=F, layout=l*0.55)

grid.echo()
p1 <- grid.grab()
p1 = ggplotify::as.ggplot(p1)

l = layout_with_mds(nets_dist[[2]])

plot(nets[[2]], edge.color = E(nets[[2]])$colors, vertex.size = 2.5, edge.curved = 1,
     vertex.color = V(nets[[2]])$colors, vertex.label = "", rescale=F, layout=l*0.55)

grid.echo()
p2 <- grid.grab()
p2 = ggplotify::as.ggplot(p2)

l = layout_with_mds(nets_dist[[3]])

plot(nets[[3]], edge.color = E(nets[[3]])$colors, vertex.size = 2.5, edge.curved = 1,
     vertex.color = V(nets[[3]])$colors, vertex.label = "", rescale=F, layout=l*0.55)

grid.echo()
p3 <- grid.grab()
p3 = ggplotify::as.ggplot(p3)

########################## Fake legend ########################################

data <- data.frame(
  Year = 2000:2010,
  Value_A = cumsum(runif(11)),
  Value_B = cumsum(runif(11))
)

dummy_plot = ggplot(data, aes(x = Year)) +
  geom_line(aes(y = Value_A, color = 'Positive')) +
  geom_line(aes(y = Value_B, color = 'Negative')) + 
  scale_color_manual(values = c("#FF37D5", "#4335C9"), name = 'Association') +
  theme_bw()

l1 = cowplot::get_legend(dummy_plot)
l1 = as.ggplot(l1)

dummy_plot2 <- ggplot(data, aes(x = Year)) +
  geom_point(aes(y = Value_A, fill = 'Bacteria'), shape = 21, size = 3) +
  geom_point(aes(y = Value_B, fill = 'Fungi'), shape = 21, size = 3) +
  scale_fill_manual(values = c('Bacteria' = 'black', 'Fungi' = 'white'), name = 'Guild') +
  theme_bw()

l2 = cowplot::get_legend(dummy_plot2)
l2 = as.ggplot(l2)

############################### Final plot ####################################

layout = c(
  area(t = 1, b = 10, l = 1, r = 10),
  area(t = 1, b = 10, l = 12, r = 22),
  area(t = 1, b = 10, l = 20, r = 30),
  area(t = 7, b = 9, l = 29, r = 30)
)

l = c(area(t = 1, b = 15, l = 1, r = 5),
      area(t = 10, b = 15, l = 1, r = 5))

legend = l1 + theme(panel.spacing = unit(0, "lines")) + l2 + plot_layout(design = l) 

final_plot = p1 + ggtitle('pH 3.7~4.5') + theme(plot.title = element_text(face = "bold", size = 20)) + 
  p2 + ggtitle('pH 4.5~6.1') + theme(plot.title = element_text(face = "bold", size = 20)) + 
  p3 + ggtitle('pH 6.1~8.0') + theme(plot.title = element_text(face = "bold", size = 20)) +
  legend + plot_layout(design = layout)

# ggsave('Figures/NETWORK.png', 
#        final_plot, 
#        device = 'png', 
#        width = 21,
#        height = 10)

# saveRDS(list(net = nets, dist = nets_dist, abs = nets_abs), 'Data/full_networks_igraph.rds')
# saveRDS(matrixes, 'matrixes_3_groups.rds')
