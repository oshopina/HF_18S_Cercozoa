library(gsveasyr)
library(ggplot2)

otu = read.csv2('Cercozoa/Data/cerco_180.csv', row.names = 1)
env <- read.csv('18S/Data/env_16S.csv', row.names = 1)
env <- env[order(env$pH),] ## Sort the way you want the samples to be ordered on heatmap
env = env[rownames(env) %in% colnames(otu),]

collapsed_otu = aggregate(otu[,1:98], by = list(Morphology = otu$nutrition), FUN = 'sum')
rownames(collapsed_otu) = collapsed_otu$Morphology
collapsed_otu = collapsed_otu[,-1] |> t() |> as.data.frame()

collapsed_otu <- collapsed_otu[rownames(env),]
anova <- auto_aov_fixed(collapsed_otu, ~ pH, env_df = env)$Results
anova = subset(anova, str_detect(Parameter, 'pH'))[, c('Data', 'F_value', 'p_value', 'Signif')]
rownames(anova) = anova$Data
anova = anova[colnames(collapsed_otu),]
anova$F_value = round(anova$F_value, digits = 2)
anova$p_value = round(anova$p_value, digits = 4)
anova$Signif <- gsub("\\*+", "*", anova$Signif)

df = data.frame(collapsed_otu, pH = env$pH)

mypal.pH <-
  colorRampPalette(
    c(
      "#9e0142",
      "#d53e4f",
      "#f46d43",
      "#fdae61",
      "#fee08b",
      "#e6f598",
      "#abdda4",
      "#66c2a5",
      "#3288bd",
      "#5e4fa2"
    )
  )

ggplot(df, aes(x = pH, y = bacterivore, color = pH)) +
  geom_point(size = 2.5) +
  scale_color_gradientn(colours = mypal.pH(256)) +
  theme_classic(base_size = 15) +
  theme(legend.position = 'none') +
  xlab('pH') +
  ylab('Bacterivore') 
