library(SpiecEasi)
library(igraph)

env = read.csv('18S/Data/env_16S.csv')
rownames(env) = env$Sample
otu = read.csv2('18S/Data/euk_6160.csv', row.names = 1)
################################ SpiecEasi #####################################

# Define a list of groups and their criteria

groups <- list(
  list(clust = 1, name = "g4"),
  list(clust = 2, name = "t45"),
  list(clust = 3, name = "g6"),
  list(clust = 4, name = "t65"),
  list(clust = 5, name = "g7")
)

otu_group_tables = list()
env_group_tables = list()

#Loop through the groups
for (i in 1:length(groups)) {
  group = groups[[i]]
  
  env_group <- env[env$OG_clustering == group$clust | env$Tipping_points == group$clust,]
  
  otu_group <- otu[, colnames(otu) %in% env_group$Sample]
  otu_group <- otu_group[which(apply(otu_group, 1, max) > 5),] |> t()
  
  otu_group_tables[[i]] = otu_group
  env_group_tables[[i]] = env_group
  
  ## Perform the network analysis and save the result (uncomment this section when needed)
  # spieceasi <- spiec.easi(list(otu_group, fungi_group),
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
  otu_group = otu_group_tables[[i]]
  
  spiec = readRDS(paste0('18S/Results/Networks/network_', group$name, '.rds'))
  matrix = symBeta(getOptBeta(spiec), mode='maxabs') |> as.matrix()
  colnames(matrix) = colnames(otu_group)
  rownames(matrix) = colnames(otu_group)
  
  net = graph_from_adjacency_matrix(matrix,mode = "undirected", weighted = TRUE, diag = FALSE)
  V(net)$name = colnames(matrix)
  
  net_dist <- net
  E(net_dist)$weight <- 1 - abs(E(net_dist)$weight)
  
  net_abs <- net
  E(net_abs)$weight <- abs(E(net_abs)$weight)
  
  E(net)[weight>0]$color <- "#bf812d"
  E(net)[weight<0]$color <- "#4575b4"
  
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
     vertex.color = V(nets[[1]])$colors, vertex.label = "", rescale=F, layout=l*0.23)

grid.echo()
p1 <- grid.grab()
p1 = ggplotify::as.ggplot(p1)

l = layout_with_mds(nets_dist[[2]])

plot(nets[[2]], edge.color = E(nets[[2]])$colors, vertex.size = 2.5, edge.curved = 1,
     vertex.color = V(nets[[2]])$colors, vertex.label = "", rescale=F, layout=l*0.035)

grid.echo()
p2 <- grid.grab()
p2 = ggplotify::as.ggplot(p2)

l = layout_with_mds(nets_dist[[3]])

plot(nets[[3]], edge.color = E(nets[[3]])$colors, vertex.size = 2.5, edge.curved = 1,
     vertex.color = V(nets[[3]])$colors, vertex.label = "", rescale=F, layout=l*0.45)

grid.echo()
p3 <- grid.grab()
p3 = ggplotify::as.ggplot(p3)

l = layout_with_mds(nets_dist[[4]])

plot(nets[[4]], edge.color = E(nets[[4]])$colors, vertex.size = 2.5, edge.curved = 1,
     vertex.color = V(nets[[4]])$colors, vertex.label = "", rescale=F, layout=l*0.30)

grid.echo()
p4 <- grid.grab()
p4 = ggplotify::as.ggplot(p4)

l = layout_with_mds(nets_dist[[5]])

plot(nets[[5]], edge.color = E(nets[[5]])$colors, vertex.size = 2.5, edge.curved = 1,
     vertex.color = V(nets[[5]])$colors, vertex.label = "", rescale=F, layout=l*0.62)

grid.echo()
p5 <- grid.grab()
p5 = ggplotify::as.ggplot(p5)

########################## Fake legend ########################################

data <- data.frame(
  Year = 2000:2010,
  Value_A = cumsum(runif(11)),
  Value_B = cumsum(runif(11))
)

dummy_plot = ggplot(data, aes(x = Year)) +
  geom_line(aes(y = Value_A, color = 'Positive')) +
  geom_line(aes(y = Value_B, color = 'Negative')) + 
  scale_color_manual(values = c("#bf812d", "#4575b4"), name = 'Association') +
  theme_bw()

l1 = cowplot::get_legend(dummy_plot)
l1 = as.ggplot(l1)


############################### Final plot ####################################

layout = c(
  area(t = 1, b = 20, l = 1, r = 30),
  area(t = 20, b = 40, l = 15, r = 35), 
  area(t = 1, b = 20, l = 22, r = 52),
  area(t = 20, b = 40, l = 40, r = 70),
  area(t = 1, b = 20, l = 50, r = 95),
  area(t = 20, b = 40, l = 70, r = 80))


final_plot = p1 + ggtitle('pH 3.7~4.5') + theme(plot.title = element_text(face = "bold", size = 15)) +
  p2 + ggtitle('pH 4.3~4.7') + theme(plot.title = element_text(face = "bold", size = 15, colour = 'red')) +
  p3 + ggtitle('pH 4.5~6.1') + theme(plot.title = element_text(face = "bold", size = 15)) +
  p4 + ggtitle('pH 6.3~6.7') + theme(plot.title = element_text(face = "bold", size = 15, colour = 'red')) +
  p5 + ggtitle('pH 6.1~8.0') + theme(plot.title = element_text(face = "bold", size = 15)) +
  l1 + plot_layout(design = layout)


# ggsave('18S/Figures/NETWORK.svg',
#        final_plot,
#        device = 'svg',
#        width = 21,
#        height = 15)

# saveRDS(list(net = nets, dist = nets_dist, abs = nets_abs), '18S/Data/full_networks_igraph.rds')
# saveRDS(matrixes, 'matrixes_3_groups.rds')
