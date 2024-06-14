library(SpiecEasi)
library(igraph)
library(stringr)
library(dplyr)

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
  bacteria_group <- bacteria_group[which(apply(bacteria_group, 1, max) > 4),] |> t()
  non_zero_counts <- apply(bacteria_group, 2, function(x) sum(x != 0))
  bacteria_group = bacteria_group[, non_zero_counts > 2]
  
  fungi_group <- fungi[, colnames(fungi) %in% env_group$Sample]
  fungi_group <- fungi_group[which(apply(fungi_group, 1, max) > 6),] |> t()
  non_zero_counts <- apply(fungi_group, 2, function(x) sum(x != 0))
  fungi_group = fungi_group[, non_zero_counts > 2]
  
  euks_group <- euks[, colnames(euks) %in% env_group$Sample]
  euks_group <- euks_group[which(apply(euks_group, 1, max) > 6),] |> t()
  non_zero_counts <- apply(euks_group, 2, function(x) sum(x != 0))
  euks_group = euks_group[, non_zero_counts > 2]
  
  bacteria_group_tables[[i]] = bacteria_group
  fungi_group_tables[[i]] = fungi_group
  euks_group_tables[[i]] = euks_group
  env_group_tables[[i]] = env_group
  
  ## Perform the network analysis and save the result (uncomment this section when needed)
  # spieceasi <- spiec.easi(list(bacteria_group, fungi_group, euks_group),
  #                         method = 'mb', lambda.min.ratio = 1e-2, nlambda = 19,
  #                         icov.select.params = list(rep.num = 50, ncores = 10))
  # 
  # saveRDS(spieceasi, paste0("network_", group$name, ".rds"))
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
  euks_group = euks_group_tables[[i]]
  
  spiec = readRDS(paste0('18S/Results/Networks/Networks_all/network_', group$name, '.rds'))
  matrix = symBeta(getOptBeta(spiec), mode='maxabs') |> as.matrix()
  colnames(matrix) = c(colnames(bacteria_group), colnames(fungi_group), colnames(euks_group))
  rownames(matrix) = c(colnames(bacteria_group), colnames(fungi_group), colnames(euks_group))
  
  net = graph.adjacency(matrix,mode = "undirected", weighted = TRUE, diag = FALSE)
  V(net)$name = colnames(matrix)
  
  net_dist <- net
  E(net_dist)$weight <- 1 - abs(E(net_dist)$weight)
  
  net_abs <- net
  E(net_abs)$weight <- abs(E(net_abs)$weight)
  
  E(net)[weight>0]$color <- "#bf812d"
  E(net)[weight<0]$color <- "#4575b4"
  
  V(net)$color <- case_when(
    substr(V(net)$name, 1, 1) == "B" ~ "blue",
    substr(V(net)$name, 1, 1) == "F" ~ "green",
    substr(V(net)$name, 1, 1) == "E" ~ "orange",
    TRUE ~ "white"
  )
  
  matrixes[[i]] = matrix
  nets[[i]] =  net
  nets_dist[[i]] = net_dist
  nets_abs[[i]] = net_abs
  
}

############################ Network graph #####################################
# library(grid)
# library(gridGraphics)
# library(ggplotify)
# library(patchwork)
# library(ggplot2)
# 
# l = layout_with_mds(nets_dist[[1]])
# 
# plot(nets[[1]], edge.color = E(nets[[1]])$colors, vertex.size = 2.5, edge.curved = 1,
#      vertex.color = V(nets[[1]])$colors, vertex.label = "", rescale=F, layout=l*0.50)
# 
# grid.echo()
# p1 <- grid.grab()
# p1 = ggplotify::as.ggplot(p1)
# 
# l = layout_with_mds(nets_dist[[2]])
# 
# plot(nets[[2]], edge.color = E(nets[[2]])$colors, vertex.size = 2.5, edge.curved = 1,
#      vertex.color = V(nets[[2]])$colors, vertex.label = "", rescale=F, layout=l*0.25)
# 
# grid.echo()
# p2 <- grid.grab()
# p2 = ggplotify::as.ggplot(p2)
# 
# l = layout_with_mds(nets_dist[[3]])
# 
# plot(nets[[3]], edge.color = E(nets[[3]])$colors, vertex.size = 2.5, edge.curved = 1,
#      vertex.color = V(nets[[3]])$colors, vertex.label = "", rescale=F, layout=l*0.55)
# 
# grid.echo()
# p3 <- grid.grab()
# p3 = ggplotify::as.ggplot(p3)
# 
# l = layout_with_mds(nets_dist[[4]])
# 
# plot(nets[[4]], edge.color = E(nets[[4]])$colors, vertex.size = 2.5, edge.curved = 1,
#      vertex.color = V(nets[[4]])$colors, vertex.label = "", rescale=F, layout=l*0.40)
# 
# grid.echo()
# p4 <- grid.grab()
# p4 = ggplotify::as.ggplot(p4)
# 
# l = layout_with_mds(nets_dist[[5]])
# 
# plot(nets[[5]], edge.color = E(nets[[5]])$colors, vertex.size = 2.5, edge.curved = 1,
#      vertex.color = V(nets[[5]])$colors, vertex.label = "", rescale=F, layout=l*0.70)
# 
# grid.echo()
# p5 <- grid.grab()
# p5 = ggplotify::as.ggplot(p5)
# 
# ########################## Fake legend ########################################
# 
# data <- data.frame(
#   Year = 2000:2010,
#   Value_A = cumsum(runif(11)),
#   Value_B = cumsum(runif(11)),
#   Value_C = cumsum(runif(11))
# )
# 
# dummy_plot = ggplot(data, aes(x = Year)) +
#   geom_line(aes(y = Value_A, color = 'Positive')) +
#   geom_line(aes(y = Value_B, color = 'Negative')) + 
#   scale_color_manual(values = c("#4575b4", "#bf812d"), name = 'Association') +
#   theme_bw()
# 
# l1 = cowplot::get_legend(dummy_plot)
# l1 = as.ggplot(l1)
# 
# dummy_plot2 <- ggplot(data, aes(x = Year)) +
#   geom_point(aes(y = Value_A, fill = 'Bacteria'), shape = 21, size = 3) +
#   geom_point(aes(y = Value_B, fill = 'Fungi'), shape = 21, size = 3) +
#   geom_point(aes(y = Value_C, fill = 'Eukaryotes'), shape = 21, size = 3) +
#   scale_fill_manual(values = c('Bacteria' = 'blue', 'Fungi' = 'green', 'Eukaryotes' = 'orange'), name = 'Guild') +
#   theme_bw()
# 
# l2 = cowplot::get_legend(dummy_plot2)
# l2 = as.ggplot(l2)
# ############################### Final plot ####################################
# 
# layout = c(
#   area(t = 1, b = 20, l = 1, r = 30),
#   area(t = 20, b = 40, l = 15, r = 35), 
#   area(t = 1, b = 20, l = 25, r = 52),
#   area(t = 20, b = 40, l = 40, r = 70),
#   area(t = 1, b = 20, l = 50, r = 80),
#   area(t = 20, b = 40, l = 70, r = 80),
#   area(t = 25, b = 45, l = 70, r = 80))
# 
# 
# final_plot = p1 + ggtitle('pH 3.7~4.5') + theme(plot.title = element_text(face = "bold", size = 15)) +
#   p2 + ggtitle('pH 4.3~4.7') + theme(plot.title = element_text(face = "bold", size = 15, colour = 'red')) +
#   p3 + ggtitle('pH 4.5~6.1') + theme(plot.title = element_text(face = "bold", size = 15)) +
#   p4 + ggtitle('pH 6.3~6.7') + theme(plot.title = element_text(face = "bold", size = 15, colour = 'red')) +
#   p5 + ggtitle('pH 6.1~8.0') + theme(plot.title = element_text(face = "bold", size = 15)) +
#   l1 + l2 + plot_layout(design = layout)


# ggsave('18S/Figures/NETWORKS_all.svg',
#        final_plot,
#        device = 'svg',
#        width = 21,
#        height = 16)

# saveRDS(list(net = nets, dist = nets_dist, abs = nets_abs), '18S/Data/full_networks_igraph_all.rds')
# saveRDS(matrixes, 'matrixes_3_groups.rds')

