#######use phyloseq######
library(phyloseq)
library(ggplot2)
library(dplyr)
library(scales)
library(grid)
library(reshape2)
library(vegan)
####import biom otu table with sample metadata and taxonomy######
setwd("/Users/Beto/Desktop/BIOSFERA")

bios_otu_table<-("./BIOSFERA_otu-table_n64_rare21712_barrel_collapsed_w_envsmd.biom")
bios_tree<-("./rep_set.tre")
bios_physeq_data<-import_biom(bios_otu_table,
                              bios_tree,
                              parseFunction = parse_taxonomy_greengenes)
bios_physeq_data
####make PCoA and ordination plot of unifrac distances####
theme_set(theme_bw())
bios_pcoa<-ordinate(bios_physeq_data, "PCoA", "unifrac", weighted=TRUE)
PCoA_unifrac_plot<-plot_ordination(bios_physeq_data,
                                   bios_pcoa,
                                   color = "Plant",
                                   shape = "Environment")
PCoA_unifrac_plot + geom_point(size=7, alpha=0.75) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "PCoA Axis 1 (38.8%)", y = "PCoA Axis 2 (16.5%)")


####NMDS of UniFrac distance and environmental variables as vectors####
bios_nmds<-ordinate(bios_physeq_data,
                    method = "NMDS",
                    distance = "wunifrac")
NMDS_bios_unifrac_plot<-plot_ordination(bios_physeq_data,
                                        bios_nmds,
                                        color = "Plant",
                                        shape = "Environment")
print(NMDS_bios_unifrac_plot)
env.data<-as.data.frame(bios_physeq_data@sam_data)
env.data<-env.data[c(1:18), c(1,3,6,7,10:11)]
vec.env<-envfit(bios_nmds, as.data.frame(t(env.data)), perm=100, na.rm = TRUE)

####NMDS of UniFrac distances pylum level####
bios_phylum<-tax_glom(bios_physeq_data, taxrank = "Phylum")
bios_phylum_nmds<-ordinate(bios_phylum,
                           method = "NMDS",
                           distance = "wunifrac")
NMDS_bios_phylum_unifrac_plot<-plot_ordination(bios_phylum,
                                               bios_phylum_nmds,
                                               color = "Plant",
                                               shape = "Environment") +
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  theme_bw()
print(NMDS_bios_phylum_unifrac_plot)
vec.sp<-envfit(bios_phylum_nmds, as.data.frame(t(bios_phylum@otu_table)), perm=1000)
vec.sp.df<-as.data.frame(vec.sp$vectors$arrows*sqrt(vec.sp$vectors$r))
vec.sp.df$species<-rownames(vec.sp.df)
taxa_phylum<-as.data.frame(bios_phylum@tax_table)
vec.sp.df$phylum<-taxa_phylum$Phylum
vec.sp.df$pval<-vec.sp$vectors$pvals
vec.sp.df<-filter(vec.sp.df, pval < 0.01)



NMDS_bios_phylum_unifrac_plot_df<-plot_ordination(bios_phylum, bios_phylum_nmds, justDF = TRUE)

ggplot(data = NMDS_bios_phylum_unifrac_plot_df, aes(NMDS1, NMDS2)) +
  geom_point(data = NMDS_bios_phylum_unifrac_plot_df, aes(colour = Plant, shape = Environment), size = 3, alpha = 0.75) +
  geom_segment(data = vec.sp.df, aes(x=0,xend=NMDS1*0.25,y=0,yend=NMDS2*0.25),
               arrow = arrow(length = unit(0.3, "cm")),
               colour="grey", inherit.aes = FALSE) +
  geom_text(data=vec.sp.df, aes(x=NMDS1*0.25,y=NMDS2*0.25,label=phylum),size=3) +
  theme_bw() +
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t")


###Stacked barplots###
#melt to long fomat table(ggplot format)
#prune phlya below 1% in each sample
bios_phylum<-bios_physeq_data %>%   
  tax_glom(taxrank = "Phylum") %>%  #aglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% #transform to relative abundance
  psmelt() %>% #melt to long format
  filter(Abundance > 0.01) %>% #filter low abundance taxa
  arrange(Phylum) #sort data frame alphabetically by phylum

phylum_colors <- c(
  "#CBD588", "#5F7FC7", "orange","#DA5724", "#508578", "#CD9BCD",
  "#AD6F3B", "#673770","#D14285", "#652926", "#C84248", 
  "#8569D5", "#5E738F","#D1A33D", "#8A7C64", "#599861"
)
ggplot(bios_phylum, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phylum_colors) +
  theme_bw() +
  facet_grid(~Plant, scale = "free", space = "free_x") +
  scale_y_continuous(labels = percent) + ylab("Abundace (%)")
