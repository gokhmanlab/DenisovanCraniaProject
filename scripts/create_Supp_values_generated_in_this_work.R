
# Add Nasion angle and Prosthion angle
NAAPRA = read.csv(
  paste('./data/facial_prognathism_angles.csv',
        sep = '/'),colClasses = c('character',rep('NULL','6'),rep('numeric','1')))

# Add Forehead height and Calvarial sagital flatness----------------------------

lateral = read.csv('./data/foreahead_height_and_oxycephaly.csv',header = TRUE,skip = 3)

lateral = lateral[order(lateral$group),]%>%
  dplyr::select(!c(SkFltnsCS,FhTYDiff))%>% #These cols interfere with the followup code, so I remove them here
  dplyr::rename(specimen = name_in_code)%>%
  dplyr::rename(glab.prot = NewK)%>%
  dplyr::rename(SkFltnsCS = SkFltnsCS_Final)%>%
  dplyr::rename(FhTYDiff = FhTYDiff_Final)%>%
  dplyr::select(specimen,SkFltnsCS,FhTYDiff,glab.prot)

# Add malar flatenning results--------------------------------------------------

# import the data
malar = read_excel('./data/malar_flatenning_values.xlsx')

# process the data
malar = malar[order(malar$group),]%>%
  rename(trimmedR = CurvRPAngFit1,trimmedL = CurvLPAngFit1)%>%
  rename(RANSAC_R = CurvRPAngFit2,RANSAC_L = CurvLPAngFit2)%>%
  rename(specimen = name_in_code)%>%
  filter(include=='TRUE')%>%
  filter(Gadi_name!='M676573_Harbin ventral')%>% #remove the original Harbin version #
  #filter(Gadi_name!='M676573_Harbin_2_Edited')%>% #remove the additional Harbin version #
  filter(!Gadi_name %in% c('M683591_Cro-Magnon II ventral','M684240_Qafzeh IX EM2027 ventral'))%>% # remove specimens where both sides are problematic
  select(specimen,trimmedR,trimmedL)

# fix specimens where one side was not acquired correctly 
malar[malar['specimen']=='broken_hill','trimmedR'] = malar[malar['specimen']=='broken_hill','trimmedL'] # R->L
malar[malar['specimen']=='petralona_1','trimmedR'] = malar[malar['specimen']=='petralona_1','trimmedL'] # R->L
malar[malar['specimen']=='turkana','trimmedL'] = malar[malar['specimen']=='turkana','trimmedR'] # L->R
malar[malar['specimen']=='irhoud_1','trimmedL'] = malar[malar['specimen']=='irhoud_1','trimmedR'] # L->R
malar[malar['specimen']=='dali','trimmedL'] = malar[malar['specimen']=='dali','trimmedR'] # L->R

# for specimens with only one sufficiently preserved side replace the missing side to the correct side
malar[malar['specimen']=='gibraltar1','trimmedL'] = malar[malar['specimen']=='gibraltar1','trimmedR'] # L->R
malar[malar['specimen']=='steinheim_s11','trimmedL'] = malar[malar['specimen']=='steinheim_s11','trimmedR'] # L->R

# calculate PC1 of both sides
fdata = malar%>%
  select(trimmedR,trimmedL)

mdata = malar%>%
  select(all_of('specimen'))

pca_result <- prcomp(fdata, center = TRUE, scale. = TRUE)
PCA_data = cbind(mdata,pca_result$x)%>%
  rename(malar_flatenning_PC1 = PC1)%>%
  select(specimen,malar_flatenning_PC1)

# Add calculated glenoid fossa size and cranial base area-----------------------

mandibular_fossa_and_cranial_base_area <- read.csv(
  './data/craniometric_data.csv',
  colClasses = c('character', 'factor', rep('numeric', 324)))%>%
  slice(-1, -2)%>%
  mutate(cranial.base.elipse = pi*AUB*BNL)%>%
  mutate(fossa_temp_s = 0.5*(Ectoglenoid.entoglenoid.lengt + 
                               Postglenoid.ectoglenoid.lengt + 
                               Postglenoid.entoglenoid.lengt))%>%
  mutate(mandibular_fossa_area = sqrt(fossa_temp_s*
                                            (fossa_temp_s-Ectoglenoid.entoglenoid.lengt)*
                                            (fossa_temp_s-Postglenoid.ectoglenoid.lengt)*
                                            (fossa_temp_s-Postglenoid.entoglenoid.lengt)))%>%
  select(specimen,mandibular_fossa_area,cranial.base.elipse)


# Combine all datasets----------------------------------------------------------

datasets <- list(NAAPRA, lateral, PCA_data,mandibular_fossa_and_cranial_base_area)
merged_datastes <- reduce(datasets, full_join, by = "specimen")

# fix names---------------------------------------------------------------------

source('./scripts/functions.R')

merged_datastes <- merged_datastes %>%
  mutate(specimen_updated = fix.names(specimen)) %>%
  select(specimen_updated, everything())  


# Add group cols----------------------------------------------------------------

mapping <- c(
  "EHS" = "AMHs",
  "UPS" = "AMHs",
  "ERC" = "H. erectus",
  "MPH" = "Middle Pleistocene Homo",
  "NE " = "Neanderthals"
)

full_craniometric_data <- read.csv(
  './data/craniometric_data.csv',
  colClasses = c('character', 'factor', rep('numeric', 324))) %>%
  slice(-1, -2)%>%
  select(specimen,group)%>%
  filter(group %in% c("EHS", "UPS", "ERC", "NE", "MPH"))%>%
  mutate(group = recode(group, !!!mapping))

merged_datastes <- left_join(merged_datastes,full_craniometric_data,by = "specimen")


# Tidy up-----------------------------------------------------------------------

merged_datastes <- merged_datastes %>%
  select(-specimen)%>%
  rename(
    `specimen` = specimen_updated,
    `cranial base area` = cranial.base.elipse,
    `forehead height` = FhTYDiff,
    `glavellar protrusion` = glab.prot,
    `malar flattening` = malar_flatenning_PC1,
    `mandibular fossa size` = mandibular_fossa_area,
    `facial protrusion` = PRA.NAA.Prin1,
    `calvarial curevature` = SkFltnsCS
  )%>%
  select(specimen, group, `cranial base area`, `forehead height`,
         `glavellar protrusion`, `malar flattening`,`mandibular fossa size`,
         `facial protrusion`,`calvarial curevature`)%>%
  filter(!is.na(group))


write.csv(merged_datastes, "./results/Supp_values_generated_in_this_work.csv", row.names = FALSE)








