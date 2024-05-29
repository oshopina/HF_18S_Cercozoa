library(tidyr)
library(vegan)

otu = read.csv('18S/Data/non_normalised_otu_18S_table_with_tax_TRAITS.csv', row.names = 1)
tax = separate_wider_delim(otu, taxonomy, delim = ';', 
                           names = c("Domain","Supergroup", "Division", 
                                      "Subdivision","Class","Order", 
                                      "Family", "Genus", "Species"),
                           too_few = "align_start",
                           too_many = 'debug') |> as.data.frame()
tax = tax[,c(99:107)]
rownames(tax) = rownames(otu)

unique(tax$Domain)
tax = tax[tax$Domain != "Unassigned",] 

euk = otu[rownames(tax),1:98]
r_level = colSums(euk) |> min() 
r_level = floor(r_level/10) * 10
euk = rrarefy(t(euk), r_level)
euk = t(euk) |> as.data.frame() |> cbind(tax)
traits = otu[rownames(tax), 100:103]
euk = cbind(euk, traits)
# write.csv2(euk, '18S/Data/euk_6160.csv')

protist_tax = tax[!(is.na(tax$Supergroup)),]
protist_tax = protist_tax[protist_tax$Class != 'Embryophyceae',]
protist_tax = protist_tax[protist_tax$Subdivision != 'Fungi',]
protist_tax = protist_tax[protist_tax$Subdivision != 'Metazoa',]
protist_tax = protist_tax[rowSums(is.na(protist_tax)) < ncol(protist_tax),]
protist = otu[rownames(protist_tax), 1:98]
r_level = colSums(protist) |> min() 
r_level = floor(r_level/10) * 10
protist = rrarefy(t(protist), r_level)
protist = t(protist) |> as.data.frame() |> cbind(protist_tax)
traits = otu[rownames(protist_tax), 100:103]
protist = cbind(protist, traits)
# write.csv2(protist, '18S/Data/protist_1060.csv')

embryo_tax = tax[!(is.na(tax$Class)),]
embryo_tax = embryo_tax[embryo_tax$Class == 'Embryophyceae',]
embryo = otu[rownames(embryo_tax), 1:98]
r_level = colSums(embryo) |> min()  ## really small number of reads

metazoa_tax = tax[!(is.na(tax$Subdivision)),]
metazoa_tax = metazoa_tax[metazoa_tax$Subdivision == 'Metazoa',]
metazoa = otu[rownames(metazoa_tax), 1:98]
r_level = 110
metazoa = metazoa[,colSums(metazoa) >= r_level]
metazoa = rrarefy(t(metazoa), r_level)
metazoa = t(metazoa) |> as.data.frame() |> cbind(metazoa_tax)
traits = otu[rownames(metazoa_tax), 100:103]
metazoa = cbind(metazoa, traits)
# write.csv2(metazoa, '18S/Data/metazoa_110.csv')

subdivisions = unique(tax$Subdivision)

gyrista_tax = tax[!(is.na(tax$Subdivision)),]
gyrista_tax = gyrista_tax[gyrista_tax$Subdivision == 'Gyrista',]
gyrista = otu[rownames(gyrista_tax), 1:98]
r_level = 110 
gyrista = gyrista[,colSums(gyrista) >= r_level]
gyrista = rrarefy(t(gyrista), r_level)
gyrista = t(gyrista) |> as.data.frame() |> cbind(gyrista_tax)
traits = otu[rownames(gyrista_tax), 100:103]
gyrista = cbind(gyrista, traits)
# write.csv2(gyrista, '18S/Data/gyrista_110.csv')

gyrista_tax = tax[!(is.na(tax$Subdivision)),]
gyrista_tax = gyrista_tax[gyrista_tax$Subdivision == 'Gyrista',]
gyrista = otu[rownames(gyrista_tax), 1:98]
r_level = 110 
gyrista = gyrista[,colSums(gyrista) >= r_level]
gyrista = rrarefy(t(gyrista), r_level)
gyrista = t(gyrista) |> as.data.frame() |> cbind(gyrista_tax)
traits = otu[rownames(gyrista_tax), 100:103]
gyrista = cbind(gyrista, traits)
# write.csv2(gyrista, '18S/Data/gyrista_110.csv')

evosea_tax = tax[!(is.na(tax$Subdivision)),]
evosea_tax = evosea_tax[evosea_tax$Subdivision == 'Evosea_X',]
evosea = otu[rownames(evosea_tax), 1:98]
sort(colSums(evosea))
r_level = 100 
evosea = evosea[,colSums(evosea) >= r_level]
evosea = rrarefy(t(evosea), r_level)
evosea = t(evosea) |> as.data.frame() |> cbind(evosea_tax)
traits = otu[rownames(evosea_tax), 100:103]
evosea = cbind(evosea, traits)
# write.csv2(evosea, '18S/Data/evosea_100.csv')

chlorophyta_tax = tax[!(is.na(tax$Subdivision)),]
chlorophyta_tax = chlorophyta_tax[chlorophyta_tax$Subdivision == 'Chlorophyta_X',]
chlorophyta = otu[rownames(chlorophyta_tax), 1:98]
sort(colSums(chlorophyta))
r_level = 120 
chlorophyta = chlorophyta[,colSums(chlorophyta) >= r_level]
chlorophyta = rrarefy(t(chlorophyta), r_level)
chlorophyta = t(chlorophyta) |> as.data.frame() |> cbind(chlorophyta_tax)
traits = otu[rownames(chlorophyta_tax), 100:103]
chlorophyta = cbind(chlorophyta, traits)
# write.csv2(chlorophyta, '18S/Data/chlorophyta_120.csv')


