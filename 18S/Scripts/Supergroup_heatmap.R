library(ComplexHeatmap)
library(gsveasyr)

otu = read.csv2('18S/Data/euk_6160.csv', row.names = 1)
env <- read.csv('18S/Data/env_16S.csv', row.names = 1)
env <- env[order(env$pH),] ## Sort the way you want the samples to be ordered on heatmap
env = env[rownames(env) %in% colnames(otu),]

collapsed_otu = aggregate(otu[,1:98], by = list(Supergroup = otu$Supergroup), FUN = 'sum')
rownames(collapsed_otu) = collapsed_otu$Supergroup
collapsed_otu = collapsed_otu[,-1] |> t() |> as.data.frame()
collapsed_otu$Other = collapsed_otu$Eukaryota_X + collapsed_otu$Haptista + collapsed_otu$Provora
collapsed_otu = collapsed_otu[,c(1:3, 5, 7, 9, 10)]

collapsed_otu <- collapsed_otu[rownames(env),]
anova <- auto_aov_fixed(collapsed_otu, ~ pH, env_df = env)$Results
anova = subset(anova, str_detect(Parameter, 'pH'))[, c('Data', 'F_value', 'p_value', 'Signif')]
rownames(anova) = anova$Data
anova = anova[colnames(collapsed_otu),]
anova$F_value = round(anova$F_value, digits = 2)
anova$p_value = round(anova$p_value, digits = 4)
anova$Signif <- gsub("\\*+", "*", anova$Signif)



df_scaled <- collapsed_otu |> vegan::decostand(method = 'normalize', MARGIN = 2) |> scale() |> t()

my_palette = colorRampPalette(c(
  'white',
  '#eeeeee',
  '#aaaaaa',
  '#444444',
  '#3a3a3a',
  '#2d2d2d',
  'black'
))

col_fun = circlize::colorRamp2(
  c(3.7, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
  c(
    "#9e0142",
    "#d53e4f",
    "#f46d43",
    "#fdae61",
    "#fee08a",
    "#e6f598",
    "#aadda4",
    "#66a2a5",
    "#3288ad",
    "#5e4fa2"
  )
)

# Get unique pH values
unique_pH <- sort(unique(env$pH))
names_for_pH <- rep("", length(env$pH))
approx_positions <- c(3.7, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0)
# Iterate over each approximate pH value
for (pH_value in approx_positions) {
  # Find the index in unique_pH that is closest to the approximate pH value
  closest_index <- which.min(abs(unique(env$pH) - pH_value))
  
  # Find the index in env_temp$pH that corresponds to the first occurrence of the unique pH value
  first_occurrence_index <- which(env$pH == unique(env$pH)[closest_index])[1]
  
  # Set the corresponding position in names_for_pH to the approximate pH value
  names_for_pH[first_occurrence_index] <- pH_value
}

ha <- HeatmapAnnotation(
  pH = env$pH,
  empty = anno_empty(border = FALSE, height = unit(3, "mm")),
  pH_labels = anno_text(names_for_pH, rot = 0),
  col = list(pH = col_fun),
  show_legend = FALSE,
  gp = gpar(col = "black")
)

ha_f <- rowAnnotation(
  anova = anno_text(anova$Signif)
)

Heatmap(
  df_scaled,
  row_order = rownames(df_scaled),
  column_order = env$Hoosfield.ID,
  col = my_palette(100),
  show_column_names = FALSE,
  row_title = NULL,
  bottom_annotation = ha,
  left_annotation = ha_f,
  show_heatmap_legend = FALSE
)
