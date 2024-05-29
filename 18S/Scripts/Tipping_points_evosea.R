library(gsveasyr)
library(ComplexHeatmap)
library(changepoint)
library(changepoint.geo)
library(gt)

## Read environmental data
env <- read.csv('18S/Data/env_16S.csv', row.names = 1)
env <- env[order(env$pH),] ## Sort the way you want the samples to be ordered on heatmap

## Load community data
otu = read.csv2('18S/Data/evosea_100.csv', row.names = 1)
tax = otu[,98:106]
otu = otu[1:97] |> t() |> as.data.frame()
env = env[rownames(env) %in% rownames(otu),]
otu <- otu[rownames(env),] ## sort otu table by env table

## PERMANOVA
perm = vegan::adonis2(otu ~ pH, env)
# gtsave(gt(perm), '18S/Results/PERMANOVA_evosea.docx')

## Calculate distances
hellinger_diversity = function(otu_table) {
  dist_matrix = dist(vegan::decostand(otu_table, method = 'hellinger'),
                     method = 'euc') %>% as.matrix() 
  return(dist_matrix)
}
beta_diversity = hellinger_diversity(otu)

############################### Change point analysis #########################
cpt = geomcp(beta_diversity)
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
plot(dist_cpt, ncpts = 3)
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

############################## Cluster analysis ###############################
set.seed(123)
clust = kmeans(as.dist(beta_diversity), centers=3, iter.max = 999)
clust = clust$cluster
clust = as.factor(clust)

################################ Build heatmap #################################
my_palette = colorRampPalette(c('#290133', '#672976', '#AE7DBB', 'white'))
col_fun = circlize::colorRamp2(c(3.7, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.5, 8),
                               c("#9e0142", "#d53e4f", "#f46d43", "#fdae61", "#fee08a",
                                 "#e6f598", "#aadda4", "#66a2a5", "#3288ad", "#5e4fa2"))
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

## Create HeatmapAnnotation
ha = HeatmapAnnotation(
  pH = env$pH,
  empty = anno_empty(border = FALSE, height = unit(3, "mm")),
  pH_labels = anno_text(names_for_pH, rot = 0),
  col = list(pH = col_fun),
  show_legend = FALSE,
  gp = gpar(col = "black")
)

ha1 = rowAnnotation(pH_labels = anno_text(names_for_pH, rot = 0), pH = env$pH,
                    col = list(pH = col_fun), show_legend = F, 
                    show_annotation_name = F, gp = gpar(col = "black"))

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
  show_annotation_name = F,
  cluster = clust,
  show_legend = F
)

Heatmap(
  beta_diversity,
  row_order = rownames(beta_diversity),
  column_order = env$Hoosfield.ID,
  col = my_palette(100),
  show_column_names = F,
  show_row_names = F,
  bottom_annotation = ha,
  left_annotation = ha1,
  top_annotation = ha2,
  show_heatmap_legend = FALSE,
  heatmap_width = unit(24, "cm"), 
  heatmap_height = unit(24, "cm"))

decorate_annotation("cluster", 
                    grid.text(
                      "K-means cluster",
                      x = unit(0, "npc"),
                      hjust = 1,
                      y = unit(0.5, "npc"),
                      gp = gpar(fontsize = 9)
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

# p = recordPlot()
# plot.new()
# png('18S/Figures/tipping_points_evosea.png', height = 700, width = 900)
# print(p)
# dev.off()
