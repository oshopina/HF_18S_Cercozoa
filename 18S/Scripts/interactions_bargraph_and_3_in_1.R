library(dplyr)
library(igraph)
library(ggplot2)

all_networks = readRDS('18S/Data/full_networks_igraph.rds')

edge_weights = data.frame()

for (i in 1:5) {
  network = all_networks[['net']][[i]]
  edge_weight = data.frame(weight = E(network)$weight, network = i)
  edge_weights = rbind(edge_weights, edge_weight)
}

rm(edge_weight, network, i)

edge_weights$Association = ifelse(edge_weights$weight >= 0, "Positive", "Negative")

summary_data <- edge_weights %>%
  group_by(network, Association) %>%
  summarise(count = n())

positive = summary_data %>% tidyr::pivot_wider(names_from = Association, values_from = count, values_fill = 0) %>%
mutate(Positive_Percentage = Positive / (Positive + Negative) * 100)

a <-  (5000 - 0) / (70 - 50) ## Max no of interactions - min number / max % positive - min % positive
b <- 0 - a * 50

color <- c('black', 'red', 'black', 'red', 'black')

final_plot = ggplot() +
  geom_bar(data = summary_data, aes(x = as.factor(network), y = count, fill = Association), stat = "identity", position = "stack", color = 'black') +
  geom_line(data = positive, aes(x = network, y = a * Positive_Percentage + b), linewidth = 1.2, color = 'red') +
  scale_fill_manual(values = c("Positive" = "#bf812d", "Negative" = "#4575b4")) +
  scale_y_continuous(expand = c(0,0), sec.axis = sec_axis(~ (. - b) / a, name = "Positive Interactions (%)")) +
  labs(x = "pH", y = 'No. interactions') +
  theme_classic(base_size = 20) +
  scale_x_discrete(labels = c('3.7~4.5', '4.3~4.7', '4.5~6.1', '6.3~6.7', '6.1~8.0')) +
  theme(axis.text.x = element_text(colour = color))

# ggsave('Figures/interactions.png', final_plot, device = 'png')

############################### 3 in 1 graph ###################################

# plot1 = readRDS('Data/Shared_data/3_plots_05.rds')
# plot234 = readRDS('Data/Shared_data/basic_metrics_plots.rds')
# 
# plot2 = plot234[[1]] + ylab('Degree') + plot234[[2]] + ylab('Betweenness') +
#   plot234[[3]] + ylab('Closeness') + plot234[[4]] + ylab('PageRank Walker')
# 
# library(patchwork)
# 
# l = c(
#   area(t = 1, b = 10, l = 1, r = 10),
#   area(t = 1, b = 10, l = 10, r = 20),
#   area(t = 1, b = 10, l = 20, r = 30),
#   area(t = 7, b = 10, l = 27, r = 30),
#   area(t = 11, b = 20, l = 1, r = 20),
#   area(t = 11, b = 20, l = 22, r = 29)
# )
# 
# three_in_one = plot1 + plot2 + final_plot + plot_layout(design = l)
# 
# ggsave('Figures/three_in_one.svg', three_in_one,
#        device = 'svg',
#        width = 15,
#        height = 13)
