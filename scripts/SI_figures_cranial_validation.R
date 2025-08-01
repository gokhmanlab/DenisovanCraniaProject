################################################################################
# Statistical analysis of cranial matching data
################################################################################

# Load packages
library(dplyr)
library(ggplot2)
library(ggrepel)  # For non-overlapping labels
library(FSA)      # For Dunn's test


################################################################################
# Load and prepare data
################################################################################

calvarium <- read.csv('./data/foreahead_height_and_oxycephaly.csv', skip = 3, header = TRUE) %>%
  filter(name_in_code != 'bodo') %>%
  mutate(group = as.factor(group)) 

# Combine UPS and EHS into "AMH"
levels(calvarium$group) <- c(levels(calvarium$group), 'AMH')
calvarium[calvarium$group %in% c('UPS', 'EHS'), 'group'] <- 'AMH'

# Define color scheme
group_colors <- c("AMH" = '#3B94D1', "NE" = '#FBB040', "ERC" = "forestgreen", "MPH" = "red")


################################################################################
# Calvarial curvature and forehead height
################################################################################

calvarium_lat = calvarium %>%
  #filter(intact_yes_no == 'yes' & intact_yes_no.1 == 'yes')%>%
  filter(!Name_in_image == 'STW53')  #not homo erectus


calvarium_lat$flatness_error_SD = scale(calvarium_lat$SkFltnsCS_Diff)[,1]
calvarium_lat$flatness_outliers = abs(calvarium_lat$flatness_error_SD)>2

calvarium_lat$forehead_error_SD = scale(calvarium_lat$FhTYDiff_Diff)[,1]
calvarium_lat$forehead_outliers = abs(calvarium_lat$FhTYDiff_Diff)>2


# Calc corr and P val
flatness_data = calvarium_lat %>%
  filter(!flatness_outliers)

SkFltnsCS_cor = cor(flatness_data$SkFltnsCS, flatness_data$SkFltnsCS.1, method = 'pearson')
SkFltnsCS_cor_test <- cor.test(flatness_data$SkFltnsCS, flatness_data$SkFltnsCS.1, method = 'pearson')

# Plot the clean data:
ggplot(flatness_data, aes(x = SkFltnsCS, y = SkFltnsCS.1, color = group)) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linewidth = 1) + 
  geom_point() +
  geom_text_repel(aes(label = Name_in_graph),max.overlaps = Inf) +
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 15)) +
  labs(x = "Calvarial flatness in left view (unitless)", y = "Calvarial flatness in right view (unitless)")+
  scale_color_manual(values = group_colors)+
  guides(color = NULL) 

ggsave(filename = './results/supplementary_analyses/Calvarial Curevature validation.svg',device = 'svg',
       width = 10,height = 8)

################################################################################

# Calculate the final meanFhTYDiff values:


forehead_data = calvarium_lat %>%
  filter(!forehead_outliers)

# Calc corr and P val
FhTYDiff_cor = cor(forehead_data$FhTYDiff, forehead_data$FhTYDiff.1, method = 'pearson')
FhTYDiff_cor_test <- cor.test(forehead_data$FhTYDiff, forehead_data$FhTYDiff.1, method = 'pearson')

# Calculate the correlation between left and right

ggplot(forehead_data, aes(x = FhTYDiff, y = FhTYDiff.1, color = group)) +
  geom_abline(intercept = 0, slope = 1, color = "grey", linewidth = 1) + # Increased size of the abline
  geom_point() +
  geom_text_repel(aes(label = Name_in_graph),max.overlaps = Inf) +
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title = element_text(size = 15),
        legend.position = "none") +  # Remove background
  labs(x = "Forehead projected height in left view (unitless)",
       y = "Forehead projected height in right view (unitless)")+  # Change axis labels
  scale_color_manual(values = group_colors)+
  guides(color = NULL) 

ggsave(filename = './results/supplementary_analyses/Forehead height validation.svg',device = 'svg',
       width = 10,height = 8)
  

################################################################################
# Compare flatness to Ni et al.'s discrete classification
################################################################################


ggplot(flatness_data,
       aes(x = Forehead_elevation_Ni,
           y = SkFltnsCS_Final))+
  geom_boxplot()+
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 15))+  # Remove background
  guides(color = "none")+
  scale_color_manual(values = group_colors)+
  labs(x = "Glbaella concavity (taken from Ni et al. 2021)", y = "Glabellar curevature (unitless)")

ggsave(filename = './results/supplementary_analyses/Calvarial curvature cont vs disc.svg',device = 'svg',
       width = 10,height = 8)

# Perform Kruskal-Wallis test

model = kruskal.test(SkFltnsCS_Final ~ Forehead_elevation_Ni,
             data = flatness_data)

print(model)

dunnTest(SkFltnsCS_Final ~ Forehead_elevation_Ni,
         data = flatness_data)


################################################################################
# Glabellar curvature by group
################################################################################

calvarium_sup = calvarium

calvarium_sup = filter(calvarium_sup,!is.na(calvarium_sup$NewK))

ggplot(calvarium_sup,aes(x = group, y = NewK, color = group))+
  geom_boxplot(data = filter(calvarium_sup,!group == 'MPH'))+
  geom_point()+
  geom_text_repel(data = filter(calvarium_sup,group == 'MPH'), aes(label = Name_in_graph),
                  max.overlaps = Inf,hjust = "left",nudge_x = 0.05,direction = 'y',min.segment.length = 0)+
  scale_color_manual(values = group_colors)+
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))+  # Remove background
  scale_y_continuous(breaks = c(round(min(calvarium_sup$NewK),digits = 4), round(max(calvarium_sup$NewK),4)), 
                     labels = c(round(min(calvarium_sup$NewK),digits = 4), round(max(calvarium_sup$NewK),4))) +
  guides(color = FALSE)+
  scale_x_discrete(labels = c("H. erectus", "Neanderthals", "AMHs",'Middle Pleistocene Homo')) +
  labs(x = "Group", y = "Glabellar curevature (unitless)")

ggsave(filename = './results/supplementary_analyses/Glabellar curvature values.svg',device = 'svg',
       width = 10,height = 8)


################################################################################
# Glabellar curvature by Ni et al. categories
################################################################################

calvarium_sup = filter(calvarium_sup,!is.na(calvarium_sup$Glabella_concavity_Ni))


new_order = c('deep','shallow','absent')
calvarium_sup$Glabella_concavity_Ni <- factor(calvarium_sup$Glabella_concavity_Ni, levels = new_order)


ggplot(calvarium_sup,
       aes(x = Glabella_concavity_Ni,
           y = NewK))+
  geom_boxplot()+
  theme(panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 15))+  # Remove background
  guides(color = FALSE)+
  scale_y_continuous(breaks = c(round(min(calvarium_sup$NewK),digits = 3), round(max(calvarium_sup$NewK),3)), 
                     labels = c(round(min(calvarium_sup$NewK),digits = 3), round(max(calvarium_sup$NewK),3))) +
  scale_color_manual(values = group_colors)+
  labs(x = "Glbaella concavity (taken from Ni et al. 2021)", y = "Glabellar curevature (unitless)")


ggsave(filename = './results/supplementary_analyses/Glabellar curvature validation.svg',device = 'svg',
       width = 10,height = 8)

model = kruskal.test(NewK ~ Glabella_concavity_Ni,
                     data = calvarium_sup)

print(model)

dunnTest(NewK ~ Glabella_concavity_Ni,method = 'bh',
         data = calvarium_sup)

################################################################################
# Glenoid Fossa size
################################################################################

full_craniometric_data <- read.csv(
  './../processed_craniometric_data.csv')

source('./scripts/functions.R')


GF.mat = full_craniometric_data
GF.mat$specimen = fix.names(GF.mat$specimen)
GF.mat$group = as.character(GF.mat$group)
GF.mat[GF.mat$group %in% c('EHS','UPS'),'group'] = 'AMH'
GF.mat[GF.mat$group %in% c('ERC_EP','ERC_MP'),'group'] = 'ERC'
GF.mat[GF.mat$group %in% c('LMPEA','MPH'),'group'] = 'MPH'
GF.mat$group = as.factor(GF.mat$group)
GF.mat$group = factor(GF.mat$group, levels=c("MPH", "NE","ERC", "AMH","HOMO","ANT"))
GF.mat = filter(GF.mat,group%in%c("MPH", "NE","ERC", "AMH",NA))%>%
  slice(-c(1,2)) #Remove first two metadata rows

ggplot(GF.mat,
       aes(x = group,
           y = mandibular_fossa_area.neu,
           color = group))+
  geom_boxplot(data = filter(GF.mat,!group == 'MPH'))+
  geom_point()+
  geom_text_repel(data = filter(GF.mat,group == 'MPH'),aes(label = specimen),hjust = "left",nudge_x = 0.05,
                  direction = 'y',min.segment.length = 0)+
  scale_color_manual(values = group_colors)+
  theme(panel.background = element_blank(),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15))+  # Remove background
  scale_y_continuous(breaks = c(round(min(GF.mat$mandibular_fossa_area.neu),digits = 4), round(max(GF.mat$mandibular_fossa_area.neu),4)), 
                     labels = c(round(min(GF.mat$mandibular_fossa_area.neu),digits = 4), round(max(GF.mat$mandibular_fossa_area.neu),4))) +
  guides(color = FALSE)+
  scale_x_discrete(labels = c("Neanderthals", "H. erectus",  "AMHs",'Middle Pleistocene Homo')) +
  labs(x = "Group", y = expression("Glenoid fossa area (cm"^2*")"))

ggsave(filename = './results/supplementary_analyses/glenoid fossa size by group supplementary.svg',device = 'svg',
       width = 10,height = 8)








