
# Add Nasion angle and Prosthion angle
if(NAAPRA){
  NAAPRA = read.csv(
    paste('./data/facial_prognathism_angles.csv',
          sep = '/'),colClasses = c('character',rep('NULL','6'),rep('numeric','1')))
  # Add to the matrix
  full_craniometric_data = left_join(full_craniometric_data, NAAPRA, by = 'specimen')
  full_craniometric_data[1,'PRA.NAA.Prin1'] = 2
  full_craniometric_data[2,'PRA.NAA.Prin1'] = 0
}

# Add Forehead height and Calvarial sagital flatness----------------------------

if(FH_height_Oxy){
  lateral = read.csv('./data/foreahead_height_and_oxycephaly.csv',header = TRUE,skip = 3)
  
  lateral = lateral[order(lateral$group),]%>%
    dplyr::select(!c(SkFltnsCS,FhTYDiff))%>% #These cols interfere with the followup code, so I remove them here
    dplyr::rename(specimen = name_in_code)%>%
    dplyr::rename(glab.prot = NewK)%>%
    dplyr::rename(SkFltnsCS = SkFltnsCS_Final)%>%
    dplyr::rename(FhTYDiff = FhTYDiff_Final)%>%
    dplyr::select(specimen,SkFltnsCS,FhTYDiff,glab.prot)
  
  
  # Add to the matrix  
  full_craniometric_data = left_join(full_craniometric_data, lateral, by = 'specimen')
  full_craniometric_data[1,c('SkFltnsCS','FhTYDiff','glab.prot')] = as.integer(1)
  full_craniometric_data[2,c('SkFltnsCS','FhTYDiff','glab.prot')] = as.integer(0)
  
}

# Add malar flatenning results--------------------------------------------------
if(malar_flatenning){
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

  # Add to the matrix  
  full_craniometric_data = left_join(full_craniometric_data, PCA_data, by = 'specimen')
  full_craniometric_data[1,c('malar_flatenning_PC1')] = as.integer(2)
  full_craniometric_data[2,c('malar_flatenning_PC1')] = as.integer(0)
}
