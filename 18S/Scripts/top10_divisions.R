library(ggplot2)

otu = read.csv2('18S/Data/euk_6160.csv', row.names = 1)

collapsed_otu = aggregate(otu[,1:98], by = list(Subdivision = otu$Subdivision), FUN = 'sum')
rownames(collapsed_otu) = collapsed_otu$Subdivision
collapsed_otu = collapsed_otu[,-1]
collapsed_otu = collapsed_otu/6160 * 100
collapsed_otu = collapsed_otu[order(collapsed_otu$HF100, decreasing = T),]
collapsed_otu = collapsed_otu[1:10,]

env <- read.csv('18S/Data/env_16S.csv', row.names = 1)
env <- env[order(env$pH),] ## Sort the way you want the samples to be ordered on heatmap
env = env[rownames(env) %in% colnames(otu),]

sample45 = rownames(env[env$pH <= 4.5,])
sample6 = rownames(env[env$pH > 4.5 & env$pH <= 6,])
sample7 = rownames(env[env$pH > 6,])

otu_45 = collapsed_otu[,sample45]
otu_45$Mean = rowMeans(otu_45)

otu_6 = collapsed_otu[,sample6]
otu_6$Mean = rowMeans(otu_6)

otu_7 = collapsed_otu[,sample7]
otu_7$Mean = rowMeans(otu_7)

otu_45 = otu_45[order(otu_45$Mean, decreasing = T),]
otu_45$Division = rownames(otu_45)
otu_45$Division = factor(otu_45$Division, levels = rev(otu_45$Division))
p1 = ggplot(otu_45, aes(x = Division, y = Mean)) +
  geom_bar(stat="identity", width=.4, show.legend = F) +
  coord_flip() +
  xlab("") +
  theme_bw(base_size = 10) +
  ylim(0, 80) +
  ylab('Relative abundance, %')

otu_6 = otu_6[order(otu_6$Mean, decreasing = T),]
otu_6$Division = rownames(otu_6)
otu_6$Division = factor(otu_6$Division, levels = rev(otu_6$Division))
p2 = ggplot(otu_6, aes(x = Division, y = Mean)) +
  geom_bar(stat="identity", width=.4, show.legend = F) +
  coord_flip() +
  xlab("") +
  theme_bw(base_size = 10) +
  ylim(0, 80) +
  ylab('Relative abundance, %')

otu_7 = otu_7[order(otu_7$Mean, decreasing = T),]
otu_7$Division = rownames(otu_7)
otu_7$Division = factor(otu_7$Division, levels = rev(otu_7$Division))
p3 = ggplot(otu_7, aes(x = Division, y = Mean)) +
  geom_bar(stat="identity", width=.4, show.legend = F) +
  coord_flip() +
  xlab("") +
  theme_bw(base_size = 10) +
  ylim(0, 80) +
  ylab('Relative abundance, %')

library(patchwork)

p = p1 + p2 +p3

ggsave('18S/Figures/top10_divisions.png', p, width = 11, height = 3)
