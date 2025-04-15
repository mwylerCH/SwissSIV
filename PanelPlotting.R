library(data.table)
library(tidyverse)
library(viridis)
library(ggbiplot)
library(cowplot)

tabellaNA <- fread('~/tempSIV/trimmed_NOref_NA.DistMatrix')
matriceNA <- as.matrix(tabellaNA[,-1])
rownames(matriceNA) <- tabellaNA$V1


# PCA
pcaResNA <- cmdscale(matriceNA)
pcaResNA <- pcaResNA %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'nome')

# no color
pcaResNA %>% 
  ggplot(aes(x=V1, y=V2, label = nome))+
  geom_text()+
  theme_bw()+
  labs(title = 'NA')


# specify clades 
firstCladesNA <- rownames(pcaResNA[pcaResNA$V1 < 0, ])
secondCladesNA <- rownames(pcaResNA[pcaResNA$V1 > 0, ])
pcaResNA[rownames(pcaResNA) %in% firstCladesNA, 'clusterNA' ] <- 1
pcaResNA[rownames(pcaResNA) %in% secondCladesNA, 'clusterNA' ] <- 2
pcaResNA$clusterNA <- as.factor(pcaResNA$clusterNA)

pcaResNA %>% 
  rename(Subclade = clusterNA) %>% 
  ggplot(aes(x=V1, y=V2, label = nome, colour = Subclade))+
  geom_text(size = 1)+
  scale_color_manual(values = c('black', 'red3'))+
  theme_bw()+
  labs(title = 'NA Diversity',
       x='PC1', y='PC2')
ggsave('~/Downloads/SIV_CH/SuppFig1.pdf')

#### temporal pattern -----------------------------
SampleTable <- fread('~/Downloads/SIV_CH/sampleDescription2.csv')
doubleSeq <- SampleTable %>% 
  filter(!is.na(Clade) & !is.na(clusterNA))

# only HA clades
tempo <- SampleTable %>% 
  #filter(!grepl('EPI', Haid)) %>% 
  mutate(SmpYear = gsub('SIV_.+_(.+)_H1', '20\\1', Haid)) %>% 
  mutate(SmpYear = ifelse(Haid == 'EPI_ISL_85201_HA', 2010, SmpYear)) %>% 
  mutate(SmpYear = ifelse(Haid == 'EPI_ISL_98715_HA' | Haid == 'EPI_ISL_19031578_HA', 2011, SmpYear)) %>% 
  mutate(SmpYear = as.numeric(SmpYear)) %>% 
  filter(!is.na(SmpYear)) %>% 
  filter(!is.na(Clade)) %>% 
  group_by(Clade, SmpYear) %>% 
  summarise(Counts = n()) %>% 
  na.omit() %>% 
  ggplot(aes(x=SmpYear, y=Counts, fill = Clade)) +
  geom_bar(stat="identity", color="black")+
  theme_bw()+
  labs(title = 'Clade Occurence (HA)',
       x='')+
  theme(legend.position="bottom",
        legend.title = element_blank()) + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
tempo


# assign HA-NA groups
doubleSeq$Rearrangement <- ifelse(doubleSeq$Clade == 'C.2.1' & doubleSeq$clusterNA == 1, 'HA C.2.1 - NA Subclade 1', NA)
doubleSeq$Rearrangement <- ifelse(doubleSeq$Clade == 'C.2.2' & doubleSeq$clusterNA == 2, 'HA C.2.2 - NA Subclade 2', doubleSeq$Rearrangement)
doubleSeq$Rearrangement <- ifelse(doubleSeq$Clade == 'C.2.1' & doubleSeq$clusterNA == 2, 'HA C.2.1 - NA Subclade 2', doubleSeq$Rearrangement)
doubleSeq$Rearrangement <- ifelse(doubleSeq$Clade == 'C.2.2' & doubleSeq$clusterNA == 1, 'HA C.2.2 - NA Subclade 1', doubleSeq$Rearrangement)

# plot
tipi <- doubleSeq %>% 
  filter(grepl('SIV_', Haid)) %>% 
  mutate(SmpYear = gsub('SIV_.+_(.+)_H1', '20\\1', Haid)) %>% 
  mutate(SmpYear = as.numeric(SmpYear)) %>% 
  group_by(Rearrangement, SmpYear) %>% 
  summarise(Counts = n()) %>% 
  na.omit() %>% 
  ggplot(aes(x=SmpYear, y=Counts, fill=Rearrangement)) +
  geom_bar(stat="identity", color="black")+
  scale_fill_manual(values=c("blue4", "lightblue", "red2", "DarkRed"))+
  theme_bw()+
  labs(title = 'Reassortments (HA and NA)',
       x='')+
  theme(legend.position="bottom",
        legend.title = element_blank()) + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
tipi


# make panel
plot_grid(tempo, tipi, labels = c('A', 'B'), label_size = 12)
ggsave('~/Downloads/SIV_CH/panel.pdf', width = 9, height = 5)
