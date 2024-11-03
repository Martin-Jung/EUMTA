library(tidyr)
library(dplyr)
library(prioritzr)
library(ggplot2)


load("export/Species_Habitats_Range_protected.RData")
hist(spp_hab_rpz$prop_range_N2K)

##
# slope=90/(log10(250000)-log10(1000))
# range=runif(10000,50,350000)
# target=pmin(100,pmax(10,(100-(slope*(log10(range)-log10(1000))))))
# df=data.frame(range=range, target=target)
#
# p <- ggplot(df, aes(x = log10(range), y = target)) +
#   geom_point() +
#   geom_line()
#
# # Set the x-axis limits
# p + ylim(0,110)
# xlim



df <- separate(spp_hab_rpz, name, into = c("country", "bioregion", "genus", "species"), sep = "_")

spp_hab_rpz = spp_hab_rpz %>% filter(!grepl("Rattus|Tursiops|Alectoris_chukar|Psittacula|Mus_musculus",name))
test=spp_hab_rpz %>% arrange(name)
# species within the EU
spEU=spp_hab_rpz %>% filter(feat_type=='EU_species')

# soecies within the EU but only HD annexes
spEUhab=spEU %>% filter(!is.na(Art17_eu_status))


totalgap1=spEU %>% filter(prop_range_N2K<0.01)
totalgap2=spEU %>% filter(prop_range_allPA<0.01)

partialgap1=spEU %>% filter(prop_range_N2K<0.1)
partialgap2=spEU %>% filter(prop_range_allPA<0.1)

write.csv(partialgap1,'partialgapspeciesN2K.csv',row.names=FALSE)
# species within bioregion
spbioregion=spp_hab_rpz %>% filter(feat_type=='MS_BGR_species')

spbioregionsep <- separate(spbioregion, name, into = c("country", "bioregion", "genus", "species"), sep = "_") %>% arrange(genus,species)
totalSPB1=spbioregionsep %>% filter(prop_range_N2K<0.01) %>% arrange(genus)
totalSPB2=spbioregionsep %>% filter(prop_range_allPA<0.01)
partialSPB1=spbioregionsep %>% filter(prop_range_N2K<0.10) %>% arrange(genus)
partialSPB2=spbioregionsep %>% filter(prop_range_allPA<0.10)

# habitats within europe
hbeurope=spp_hab_rpz %>% filter(feat_type=='EU_habitats')

totalgapHB1=hbeurope %>% filter(prop_range_N2K<0.01)
totalgapHB2=hbeurope %>% filter(prop_range_allPA<0.01)
partialgapHB1=hbeurope %>% filter(prop_range_N2K<0.1)
partialgapHB2=hbeurope %>% filter(prop_range_allPA<0.1)


# habitats within bioregion
hbbioregion=spp_hab_rpz %>% filter(feat_type=='MS_BGR_habitats')
hbbioregionsep <- separate(hbbioregion, name, into = c("country", "bioregion", "habitat"), sep = "_") %>% arrange(habitat)
totalHABB1=hbbioregionsep %>% filter(prop_range_N2K<0.01)
totalHABB2=hbbioregionsep %>% filter(prop_range_allPA<0.01)

#rodriguez targets
slope=90/(log10(250000)-log10(1000))
spEU= spEU %>% mutate(rodrtargetperc= pmin(100,pmax(10,(100-(slope*(log10(total_range_size)-log10(1000)))))))
spEU=spEU %>% mutate(rodrtargetabs=rodrtargetperc*total_range_size/100)
spEU=spEU %>% mutate(shortfall=pmin(100,pmax(0,rodrtargetperc-prop_range_allPA*100)))
head(spEU)

hist(spEU$shortfall)




#plot spEU and spEUhab
ggplot(spEUhab, aes(x = 100-(100*shortfall/rodrtargetperc))) +
  geom_histogram(fill = "gray", color = "black", breaks = seq(0, 100, by = 10)) +
  ggtitle("Percentage of protection target met for 296 plants and animal species in Annex II of the HD in 2023") +
  xlab("Percentage of the protection target met") +
  ylab("Number of species")+
  scale_x_continuous(breaks = seq(0,100,10))

mean(spEUhab$shortfall/spEUhab$rodrtargetperc)
mean(spEU$shortfall/spEU$rodrtargetperc)
