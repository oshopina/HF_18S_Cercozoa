library(tidyr)
library(vegan)

otu = read.csv('Cercozoa/Data/non_normalised_otu_CERCOZOA_table_with_tax_TRAITS.csv', row.names = 1)
tax = separate_wider_delim(otu, taxonomy, delim = ';', 
                           names = c("Domain","Supergroup", "Division", 
                                      "Subdivision","Class","Order", 
                                      "Family", "Genus", "Species"),
                           too_few = "align_start",
                           too_many = 'debug') |> as.data.frame()
tax = tax[,c(99:107)]
rownames(tax) = rownames(otu)

unique(tax$Domain)
tax = tax[tax$Subdivision == 'Cercozoa',] 
tax = tax[rowSums(is.na(tax)) < ncol(tax),]

euk = otu[rownames(tax),1:98]
r_level = colSums(euk) |> min() 
r_level = floor(r_level/10) * 10
euk = rrarefy(t(euk), r_level)
euk = t(euk) |> as.data.frame() |> cbind(tax)
traits = otu[rownames(tax), 100:103]
euk = cbind(euk, traits)
# write.csv2(euk, 'Cercozoa/Data/cerco_180.csv')

