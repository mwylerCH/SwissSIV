library(ggimage)
library(ggtree)
library(TDbook)
library(ape)
library(dplyr)
library(RColorBrewer)
library(data.table)
library(ggnewscale)

## classico splitted HA -----------------------------------

### leggi informazioni sample HA
tabella <- fread('//wsl.localhost/Ubuntu-22.04/home/mwyler/tempSIV/sampleDescription2.csv',
                 sep = ',', data.table = F) %>%
  filter(!is.na(Haid))
rownames(tabella) <- tabella$Haid
# read tree using ape
albero <- read.tree('//wsl.localhost/Ubuntu-22.04/home/mwyler/tempSIV/trimmed_HA.tree')
#root
albero <- root(albero, c('AB628080', 'JX138509'))

# add table data (to ggtree)
ggalbero <- ggtree(albero) %<+% tabella
ggalbero
alberoPulito <- ggalbero + 
  geom_tiplab(offset = 0.1, hjust = -0.2, size = 1.5, 
              align =TRUE, linesize = 0.1)
alberoPulito

# add reference or not
RefNames <- tabella %>% select(Country) %>% 
  mutate(Source = ifelse(Country != 'Switzerland', 'GenBank', 'Monitoring')) %>% 
  select(-Country)
p1 <- gheatmap(alberoPulito, RefNames, offset=.2, width=0.14, 
               colnames=T, legend_title="genotype", 
               font.size = 2, colnames_angle = 45, hjust = 1) +
  #geom_tiplab(offset = .2, hjust = .5)+
  scale_fill_manual(values=c('firebrick3', 'limegreen'), name="Source")
p1

# add country
paese <- tabella %>% select(Country) %>% mutate(Country = ifelse(Country == 'Reference', 'GenBank', Country))
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, paese, , offset=.3, width=0.14, 
               colnames=T, legend_title="genotype", 
               font.size = 2, colnames_angle = 45, hjust = 1) +
  scale_fill_brewer(palette = 'Pastel2', name="Country", na.value = "grey50")

p2

# add clade
clades <- tabella %>% select(Clade)
p3 <- p2 + new_scale_fill()
p3 <-gheatmap(p3, clades, , offset=.4, width=0.14, 
              colnames=T, legend_title="genotype", 
              font.size = 2, colnames_angle = 45, hjust = 1) +
  scale_fill_manual(values = c('olivedrab4', "#FC8D62", "#8DA0CB","#E78AC3", "#A6D854", "#FFD92F", "#E5C494"), name="Clade", na.value = "grey50")
p3

# add host
host <- tabella %>% select(Host)
p4 <- p3 + new_scale_fill()
p4 <-gheatmap(p4, host,, offset=.5, width=0.14, 
              colnames=T, legend_title="genotype", 
              font.size = 2, colnames_angle = 45, hjust = 1) +
  scale_fill_manual(values=c('dodgerblue2', 'palevioletred'), name="Host") +
  coord_cartesian(clip = "off")+
  theme(plot.margin = margin(2, 4, 2, 2, "cm"))
p4

p4 + guides(Country = guide_legend(order=1),
            Pathogenicity = guide_legend(order=2),
            Clade = guide_legend(order=3),
            Host = guide_legend(order=4))

ggsave('NewFigure7.pdf', width = 10, height = 10)

## classico splitted NA -----------------------------------

## leggi informazioni sample NA
tabella <- fread('//wsl.localhost/Ubuntu-22.04/home/mwyler/tempSIV/sampleDescription2.csv',
                 sep = ',', data.table = F) %>%
  filter(!is.na(Naid))
rownames(tabella) <- tabella$Naid
# read tree using ape
albero <- read.tree('//wsl.localhost/Ubuntu-22.04/home/mwyler/tempSIV/trimmed_NA.tree')
# root (error if also JX138511 is considered)
albero <- root(albero, c('AB628082'))


# add table data (to ggtree)
ggalbero <- ggtree(albero) %<+% tabella
ggalbero
alberoPulito <- ggalbero + 
  geom_tiplab(offset = 0.1, hjust = -0.2, size = 1.2, 
              align =TRUE, linesize = 0.1)
alberoPulito

# add reference or not
RefNames <- tabella %>% select(Country) %>% 
  mutate(Source = ifelse(Country != 'Switzerland', 'GenBank', 'Monitoring')) %>% 
  select(-Country)
p1 <- gheatmap(alberoPulito, RefNames, offset=.2, width=0.11, 
               colnames=T, legend_title="genotype", 
               font.size = 2, colnames_angle = 45, hjust = 1) +
  scale_fill_manual(values=c('firebrick3', 'limegreen'), name="Source")

p1

# add country
paese <- tabella %>% select(Country) %>% mutate(Country = ifelse(Country == 'Reference', 'GenBank', Country))
p2 <- p1 + new_scale_fill()
p2 <- gheatmap(p2, paese, , offset=.3, width=0.11, 
               colnames=T, legend_title="genotype", 
               font.size = 2, colnames_angle = 45, hjust = 1) +
  scale_fill_brewer(palette = 'Pastel2', name="Country", na.value = "grey50")
p2

# add clade
clades <- tabella %>% mutate(`NA Subclade` = as.factor(clusterNA)) %>% select(`NA Subclade`)
p3 <- p2 + new_scale_fill()
p3 <-gheatmap(p3, clades, , offset=.4, width=0.11, 
              colnames=T, legend_title="genotype", 
              font.size = 2, colnames_angle = 45, hjust = 1) +
  scale_fill_brewer(palette = 'Dark2', name="NA Subclade", na.value = "grey50")
p3

# add host
host <- tabella %>% select(Host)
p4 <- p3 + new_scale_fill()
p4 <-gheatmap(p4, host,, offset=.5, width=0.11, 
              colnames=T, legend_title="genotype", 
              font.size = 2, colnames_angle = 45, hjust = 1) +
  scale_fill_manual(values=c('dodgerblue2', 'palevioletred'), name="Host") +
  coord_cartesian(clip = "off")+
  theme(plot.margin = margin(2, 4, 2, 2, "cm"))
p4

p4 + guides(Country = guide_legend(order=1),
            Pathogenicity = guide_legend(order=2),
            Clade = guide_legend(order=3),
            Host = guide_legend(order=4))

ggsave('NewFigure8.pdf', width = 10, height = 10)
