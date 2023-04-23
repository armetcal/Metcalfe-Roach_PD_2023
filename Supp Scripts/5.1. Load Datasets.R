# Supplementary script for 5. Species vs Metabolites and Disease Progression

source('../id_functions.R')

# Load CRP and buutyrate
crp = read.csv('C:/Users/armetcal/OneDrive - UBC/Grad School/Data/PD data from Mihai/Reference Files/Master/CRP_cytokine_PD.csv') %>% rename(Sample = ID) %>% 
  mutate(Sample = change_ids_to_longitudinal(Sample))
but = read.csv('C:/Users/armetcal/OneDrive - UBC/Grad School/Data/PD data from Mihai/Reference Files/Master/butyrate_production.csv') %>% rename(Sample = Record.ID) %>% 
  mutate(Sample = change_ids_to_longitudinal(Sample))


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ps = readRDS("../Reference Files/phyloseq_object_leafyadj_absolute_counts.rds")

# Add CRP and butyrate
x = ps@sam_data %>% as.matrix %>% as.data.frame()
ps@sam_data = ps@sam_data %>% as.matrix %>% as.data.frame() %>% 
  rownames_to_column('Sample') %>% left_join(crp) %>% left_join(but) %>% 
  column_to_rownames('Sample') %>% sample_data()

ps@sam_data$bristol = as.numeric(ps@sam_data$bristol)
ps@sam_data$Sex = as.factor(ps@sam_data$Sex)
ps@sam_data$ent_boolean = ps@sam_data$entacapone>0
ps@sam_data$butyrate = as.numeric(ps@sam_data$butyrate)
ps@sam_data$Log2CRP = as.numeric(ps@sam_data$Log2CRP)
# Sequencing depth
ps@sam_data$depth = sample_sums(ps)

# CLR Transform Data
ps = microbiome::transform(ps,'clr')

# Select only taxa of interest~~~~~~~~~~~
other_levels = colnames(ps@tax_table)[colnames(ps@tax_table) != 'Species']
psm = ps %>% psmelt()
# Univar
getwd()
temp = readxl::read_xlsx('../Final Outputs/DA Results/DA_Microbiome_Status_All_Univar_Summary.xlsx',sheet='Species') %>%
  filter(sig_in_at_least_n==1)
# Multivar
temp2 = readxl::read_xlsx('../Final Outputs/DA Results/DA_Microbiome_Status_All_Multivar_Summary.xlsx',sheet='Species') %>%
  filter(sig_in_at_least_n==1, Variable=='StatusPD')
# Filter multivar to only include taxa significant in univar model
temp2 = temp2[temp2[['Species']] %in% temp[['Species']],]

# Also add F. prausnitzii since it drives a lot of pathway changes
sig = c(temp2[['Species']],'s__Faecalibacterium_prausnitzii')

temp = psm %>% filter(Species %in% sig) %>%
  select(Species,everything()) %>% select(-any_of(other_levels),-OTU)
names(temp)[1] = 'Species'
temp = temp %>% pivot_wider(names_from = Species, values_from = Abundance)

crpbut = temp

df = list(data = temp,sig = sig)