# 15. Load Datasets - Status

ps = readRDS("../Reference Files/phyloseq_object_leafyadj_absolute_counts.rds")
ps@sam_data$bristol = as.numeric(ps@sam_data$bristol)
ps@sam_data$Sex = as.factor(ps@sam_data$Sex)
ps@sam_data$ent_boolean = ps@sam_data$entacapone>0
# Sequencing depth
ps@sam_data$depth = sample_sums(ps)

# CLR Transform Data
ps = microbiome::transform(ps,'clr')

# Select only taxa of interest~~~~~~~~~~~
other_levels = colnames(ps@tax_table)[colnames(ps@tax_table) != 'Species']
psm = ps %>% psmelt()
# Univar
temp = readxl::read_xlsx('../Reference Files/DA Outputs/DA_Microbiome_Status_Univar_Summary.xlsx',sheet='Species') %>%
  filter(sig_in_at_least_2==1)
# Multivar
temp2 = readxl::read_xlsx('../Reference Files/DA Outputs/DA_Microbiome_Status_Multivar_Summary.xlsx',sheet='Species') %>%
  filter(sig_in_at_least_2==1, Variable=='StatusPD')
# Filter multivar to only include taxa significant in univar model
temp2 = temp2[temp2[['Species']] %in% temp[['Species']],]

# Also add F. prausnitzii since it drives a lot of pathway changes
sig = c(temp2[['Species']],'s__Faecalibacterium_prausnitzii')
temp = psm %>% filter(Species %in% sig) %>%
  select(Species,everything()) %>% select(-any_of(other_levels),-OTU)
names(temp)[1] = 'Species'
temp = temp %>% pivot_wider(names_from = Species, values_from = Abundance)

# Transform everything to 0-1 range
temp = temp %>% mutate_at(c(sig,all_of(c('bristol','depth'))), function(x) (x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T)))

df = list(data = temp,sig = sig)
rm(other_levels,psm,temp2)