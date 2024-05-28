library(ComplexHeatmap)
library(changepoint)
library(changepoint.geo)
library(gsveasyr)

################################ Prepare data ##################################

## Read environmental data
env <- read.csv('18S/Data/env_16S.csv', row.names = 1)
env <- env[order(env$pH),] ## Sort the way you want the samples to be ordered on heatmap

## Load community data
otu = read.csv2('18S/Data/metazoa_100.csv', row.names = 1)
tax = otu[,91:99]
otu = otu[1:90] |> t() |> as.data.frame()
env = env[rownames(env) %in% rownames(otu),]
otu <- otu[rownames(env),] ## sort otu table by env table

df <- otu[, apply(otu, 2, max) >= 36]  # Keep columns with max value >= 36 (3%)

############################### Change point analysis ##########################

df_mat = df |> as.matrix() |> vegan::decostand(method = 'normalize', MARGIN = 2) 
cpt = geomcp(df_mat)
plot(cpt)

dist_cpt = cpt.meanvar(
  distance(cpt),
  method = "PELT",
  penalty = 'CROPS',
  pen.value = c(5, 500)
)

pen.value.full(dist_cpt)
dist_var = cpts.full(dist_cpt)
tail(dist_var)
plot(dist_cpt, diagnostic = T)
plot(dist_cpt, ncpts = 1)
######################### WRITE CHOSEN NUMBER OF CHANGEPOINTS FOR ANGLE AND DISTANCE
cp_dist = readline('Number of changepoint for distance data: ') |> as.numeric()

ang_cpt = cpt.meanvar(
  angle(cpt),
  method = "PELT",
  penalty = 'CROPS',
  pen.value = c(5, 500)
)

pen.value.full(ang_cpt)
ang_var = cpts.full(ang_cpt)
tail(ang_var)
plot(ang_cpt, diagnostic = T)
plot(ang_cpt, ncpts = 3)
######################### WRITE CHOSEN NUMBER OF CHANGEPOINTS FOR ANGLE AND DISTANCE
cp_angle = readline('Number of changepoint for angle data: ') |> as.numeric()

######################## Prepare taxonomy labels for heatmap ##################
tax = tax[rownames(tax) %in% colnames(df),]
tax = tax[!is.na(tax$Class),]
tax <- tax[order(tax$Class, tax$Order, tax$Family, tax$Genus),] ## Order taxonomy alphabetically
tax[is.na(tax)] <- 'Unc.'

tax$phylum_label <- paste0('c_',tax$Class)
tax$phylum_label[duplicated(tax$phylum_label)] <- ""
tax$order_label <- paste0('o_', tax$Order)
tax$order_label[duplicated(tax$order_label)] <- ""
tax$f_g_label <- paste0('f_', tax$Family, ';g_', tax$Genus)
tax$f_g_label[duplicated(tax$f_g_label)] <- ""
tax$otu_label <- paste0('[', substring(rownames(tax), 4), ']')

df <- df[, rownames(tax)]

## Perform ANOVA analysis
anova_results = data.frame()
for(i in colnames(df)) {
  anova = anova(aov(df[,i] ~ pH + I(pH^2), data = env))
  anova <- as.data.frame(anova[, c(1, 4, 5)])
  colnames(anova) <- c("Df", "F_value", "p_value")
  anova <- anova %>% mutate(Signif = case_when(p_value < 0.05 ~ "*", TRUE ~ " "))
  filtered_anova <- anova %>%
    filter(Signif == "*" | row_number() <= 2) %>%
    arrange(desc(Signif)) |> 
    slice(1)
  rownames(filtered_anova) = i
  anova_results = rbind(anova_results, filtered_anova)
}

anova_results = anova_results[colnames(df),]
anova_results$F_value = round(anova_results$F_value, digits = 2)
anova_results$p_value = round(anova_results$p_value, digits = 4)
anova = anova_results

## Save ANOVA results to a Word document
# gtsave(gt::gt(anova), 'singleM_community/Results/anova_ra.docx')

############################## Cluster analysis ###############################

set.seed(1)
clust = kmeans(as.dist(df_mat), centers=3, iter.max = 999)
clust = clust$cluster
clust = as.factor(clust)

############################## Heatmap ########################################

## Scale data for heatmap
df_scaled <- df |> vegan::decostand(method = 'normalize', MARGIN = 2) |> scale() |> t()
medians = apply(df, 2, function(x){
  median(x)
})

## Define color palette for pH levels and names for pH levels
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

## Create taxonomy annotations

median_col = circlize::colorRamp2(c(0, 4), c('white', 'aquamarine4'))

ha_c <- rowAnnotation(
  class = anno_text(tax$phylum_label),
  order = anno_text(tax$order_label),
  empty = anno_empty(border = FALSE, height = unit(1, "mm")),
  Median = medians,
  gp = gpar(col = "black"),
  col = list(Median = median_col),
  show_legend = F, 
  show_annotation_name = F
)

ha_f <- rowAnnotation(
  anova = anno_text(anova$Signif),
  otu = anno_text(tax$otu_label),
  family = anno_text(tax$f_g_label),
  empty = anno_empty(border = FALSE, height = unit(1, "mm")))


## Create additional annotation for change point frequency
ha2 = HeatmapAnnotation(
  change_point_dist = anno_lines(
    data.set(dist_cpt),
    axis = F
  ),
  change_point_ang = anno_lines(
    data.set(ang_cpt),
    axis = F
  ),
  height = unit(3, "cm"), 
  cluster = clust,
  show_legend = F,
  show_annotation_name = F
)

## Create Heatmap
my_palette = colorRampPalette(c(
  'white',
  '#eeeeee',
  '#aaaaaa',
  '#444444',
  '#3a3a3a',
  '#2d2d2d',
  'black'
))

draw(Heatmap(
  df_scaled,
  row_order = rownames(df_scaled),
  column_order = env$Hoosfield.ID,
  col = my_palette(100),
  show_column_names = FALSE,
  show_row_names = FALSE,
  row_title = NULL,
  row_split = tax$Class,
  top_annotation = ha2,
  bottom_annotation = ha,
  left_annotation = ha_c,
  right_annotation = ha_f,
  show_heatmap_legend = FALSE
))

change_points_dist = cpts(dist_cpt, cp_dist)[!is.na(cpts(dist_cpt, cp_dist))]

#creating vertical line for each change point in distance
for(i in change_points_dist) {
  decorate_annotation("change_point_dist", {
    grid.lines(
      x = unit(c(i, i), 'native'),
      y = unit(c(min(
        data.set(dist_cpt)
      ), max(
        data.set(dist_cpt)
      )), 'native'),
      gp = gpar(col = "red", lwd = 3)
    )
  })
}

#Adding title
decorate_annotation("change_point_dist", {
  grid.text(
    "Mean\nchange\npoint",
    x = unit(0, "npc"),
    y = unit(0.5, "npc"),
    rot = 90,
    vjust = -0.1, 
    gp = gpar(fontsize = 9)
  )
})

change_points_ang = cpts(ang_cpt, cp_angle)[!is.na(cpts(ang_cpt, cp_angle))]

#creating vertical line for each change point in angle
for(i in change_points_ang) {
  decorate_annotation("change_point_ang", {
    grid.lines(
      x = unit(c(i, i), 'native'),
      y = unit(c(min(
        data.set(ang_cpt)
      ), max(
        data.set(ang_cpt)
      )), 'native'),
      gp = gpar(col = "red", lwd = 3)
    )
  })
}

#Adding title
decorate_annotation("change_point_ang", {
  grid.text(
    "Variance\nchange\npoint",
    x = unit(0, "npc"),
    y = unit(0.5, "npc"),
    rot = 90,
    vjust = -0.1,
    gp = gpar(fontsize = 9)
  )
})

lgd = Legend(title = 'Median abundance', col_fun = median_col, at = c(0, 4),
             direction = 'horizontal', border = 'black', legend_width = unit(3, "cm"))
draw(lgd, x = unit(0.15, "npc"), y = unit(0.95, "npc"))

p = recordPlot()
plot.new()
png('18S/Figures/Heatmap_metazoa.png', height = 600, width = 800)
print(p)
dev.off()