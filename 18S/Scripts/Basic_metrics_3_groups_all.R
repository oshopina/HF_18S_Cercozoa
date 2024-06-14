## Script for calculating basic metrics for network analysis 

## Load required libraries
library(SpiecEasi)
library(igraph)
library(centiserve)
library(dplyr)
library(ggplot2)
library(multcomp)

# #################### otu metrics calculation ##################################

# all_networks = readRDS('18S/Data/full_networks_igraph_all.rds')
# 
# ## Calculate centrality metrics for each network
# otu_metrics <- list()
# 
# ## Loop through the files
# for (i in 1:5) {
# 
#   # Get the network objects using the constructed variable names
#   net_g <- all_networks[["net"]][[i]]
#   net_g_dist <- all_networks[["dist"]][[i]]
#   net_g_abs <- all_networks[["abs"]][[i]]
# 
#   # Calculate the centrality metrics
#   result <- data.frame(
#     degree = degree(net_g),
#     #alpha_centrality = alpha_centrality(net_g),
#     #strength = strength(net_g),
#     betweenness = betweenness(net_g_dist, v = V(net_g_dist)),
#     closeness = closeness(net_g_dist),
#     #transitivity = transitivity(net_g, type = 'localundirected'),
#     #eigen_centrality = eigen_centrality(net_g_dist)$vector,
#     #page_rank = page_rank(net_g_dist)$vector,
#     #bottleneck = bottleneck(net_g_dist, v = V(net_g_dist)),
#     #authority_score = authority_score(net_g)$vector,
#     #hub_score = hub_score(net_g)$vector,
#     #centralization = centr_degree(net_g)$res,
#     markov = markovcent(net_g)
#   )
# 
#   # Save the resulting list for the file
#   otu_metrics[[i]] <- result
# }
# rm(net_g, net_g_dist, net_g_abs, result, i)
# 
# ## Merge tables by metric
# tables_by_metric <- list()
# 
# 
# ##Iterate over the column names of the '4' table in otu_metrics
# for (i in colnames(otu_metrics[[1]])) {
#   table_names <- c(1:5)
# 
#   # Extract the column of interest from each table
#   table_list <- lapply(table_names, function(table_num) {
#     otu_metrics[[table_num]] %>% dplyr::select(matches(i))
#   })
# 
#   # Add row names as a column to each table
#   table_list <- lapply(table_list, function(table) {
#     table$row_names <- rownames(table)
#     table
#   })
# 
#   # Set the names of the new list based on the original table names
#   names(table_list) <- table_names
# 
#   # Rename the metric column in each table with the corresponding suffix
#   table_list <- lapply(names(table_list), function(table_name) {
#     suffix <- gsub("[^0-9]", "", table_name)  # Extract the numeric suffix from the table name
#     table <- table_list[[table_name]]
#     colnames(table)[colnames(table) == i] <- paste0(i, "_", table_name)
#     table
#   })
# 
#   # Merge the columns of all tables by row names with custom suffixes for duplicate columns
#   merged_table <- Reduce(function(x, y) merge(x, y, by = "row_names", all = TRUE),
#                          table_list)
# 
#   # Store the merged table in tables_by_metric
#   tables_by_metric[[i]] <- merged_table
# }
# rm(table_names, table_list, merged_table, i)
# 
# saveRDS(tables_by_metric, '18S/Data/tables_by_metric_all.rds')

########################## Violin plots #######################################
tables_by_metric = readRDS('18S/Data/tables_by_metric_all.rds')

violin_plots = list()
for (i in names(tables_by_metric)) {
  df_long <- tidyr::pivot_longer(tables_by_metric[[i]], 
                                 cols = -row_names, names_to = "sample", values_to = "value")
  df_long$sample = as.factor(df_long$sample)
  anova_result = aov(value ~ sample, data = df_long)
  tuk = glht(anova_result, linfct = mcp(sample = "Tukey"))
  tuk_cld = cld(tuk)
  cld_letters = tuk_cld$mcletters$Letters
  labels_level = max(na.omit(df_long$value)) - sd(na.omit(df_long$value))

  # Plot the violin plot
  color <- c('black', 'red', 'black', 'red', 'black')
  plot <- ggplot(df_long, aes(x = sample, y = value, fill = sample, width = 2)) +
    geom_violin(trim = F) +
    xlab("pH") +
    ylab(i) +
    theme_classic() +
    scale_fill_manual(values = c("#B51945",
                                 "#d53e4f",
                                 "#FDCF7D", 
                                 "#68C2A3",
                                 "#3288bd")) +
    scale_x_discrete(labels = c('3.7~4.5', '4.3~4.7', '4.5~6.1', '6.3~6.7', '6.1~8.0')) +
    scale_y_continuous(limits = c(0, max(df_long$value))) +
    guides(fill = 'none') +
    annotate('text', x = names(cld_letters), y = labels_level*1.35, label = cld_letters,
               size = 6) +
    theme(axis.text.x = element_text(colour = color))
  
  # Add the plot to the list
  violin_plots[[i]] <- plot
}

rm(df_long, plot, anova_result, i)
    
# pdf("Figures/violin_plots_3_groups.pdf")
# for (plot in violin_plots) {
#   print(plot)
# }
# dev.off()

# saveRDS(tables_by_metric, 'Data/Shared_data/tables_by_metric.rds')
# saveRDS(list(violin_plots$degree, violin_plots$betweenness, violin_plots$closeness, violin_plots$page_rank),
#         '18S/Data/basic_metrics_plots.rds')

