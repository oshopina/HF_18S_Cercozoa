library(tidyr)
library(vegan)

otu = read.csv('18S/Data/non_normalised_otu_with_tax_18S.csv', row.names = 1)
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
# write.csv2(euk, '18S/Data/euk_6160.csv')

protist_tax = tax[!(is.na(tax$Supergroup)),]
protist_tax = protist_tax[tax$Class != 'Embryophyceae' &
                protist_tax$Subdivision != 'Fungi' &
                protist_tax$Subdivision != 'Metazoa',]
protist_tax = protist_tax[rowSums(is.na(protist_tax)) < ncol(protist_tax),]
protist = otu[rownames(protist_tax), 1:98]
r_level = colSums(protist) |> min() 
r_level = floor(r_level/10) * 10
protist = rrarefy(t(protist), r_level)
protist = t(protist) |> as.data.frame() |> cbind(protist_tax)
# write.csv2(protist, '18S/Data/protist_1230.csv')
