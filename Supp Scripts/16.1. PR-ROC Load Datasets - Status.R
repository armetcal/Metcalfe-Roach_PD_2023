source('../id_functions.R')

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prep_tsv = function(df){
  names(df)[1] = 'temp'
  return(df %>% filter(!(temp %in% c('UNMAPPED','UNINTEGRATED','UNGROUPED'))) %>% column_to_rownames('temp') %>%
           select(-contains('Atcc'),-contains('Zymo'),-contains('blnk')) %>%
           `colnames<-`(reformat_ids(colnames(.)) %>% change_ids_to_longitudinal()))
}
# Load RPK Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df = data.table::fread('C:/Users/armetcal/OneDrive - UBC/Grad School/Data/Metagenomics/RPK counts/merged_pathabundance_unstratified.tsv') %>% prep_tsv()

# Load total metagenomic reads~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cnt = read_delim('C:/Users/armetcal/OneDrive - UBC/Grad School/Data/Metagenomics/ORIGINAL DATA/mqc_fastqc_sequence_counts_plot_1.txt') %>%
  filter(!str_detect(Sample,'Atcc')&!str_detect(Sample,'blnk')&!str_detect(Sample,'Zymo')) %>% 
  filter(str_detect(Sample,'bmtagged_1')) %>% 
  mutate(Sample = reformat_ids(Sample) %>% change_ids_to_longitudinal()) %>%
  mutate(depth=`Unique Reads`+`Duplicate Reads`) %>%
  select(Sample,depth)

# Load metadata and convert to phyoseq object~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ps = readRDS("../Reference Files/phyloseq_object_leafyadj_absolute_counts.rds")
s = psmelt(ps) %>% select(-OTU,-Abundance,-any_of(c('Kingdom','Phylum','Class','Order',
                                                    'Family','Genus','Species'))) %>% unique() %>% 
  mutate(Status = as.factor(Status),
         bristol = as.numeric(bristol),
         Sex = as.factor(Sex),
         ent_boolean = entacapone>0) %>% 
  left_join(cnt) %>% `rownames<-`(.$Sample)

# Taxonomy tables
t = matrix(data = rownames(df), ncol = 1, dimnames = list(rownames(df),'Species'))

# Filter OTU tables to remove the two failed samples (no bacterial OTUs)
df = df %>% select(all_of(rownames(s))) %>% as.matrix

# Create phyloseq objects
ps = phyloseq(sample_data(s),otu_table(df,taxa_are_rows = T),tax_table(t))

# CLR Transform Data
ps = microbiome::transform(ps,'clr')

# Select only pathways of interest~~~~~~~~~~~
psm = ps %>% psmelt()
# Load significant pwys
sig_multi = readxl::read_xlsx('../Final Outputs/DA Results/DA_PWY_Status_All_Summary_3tools.xlsx',sheet='multivar') %>% 
  filter(Variable=='StatusPD',sig_in_at_least_n == T)
sig_uni = readxl::read_xlsx('../Final Outputs/DA Results/DA_PWY_Status_All_Summary_3tools.xlsx',sheet='univar') %>% 
  filter(sig_in_at_least_n == T) %>% filter(Species %in% sig_multi$Species)

sig = sig_uni$Species

temp = psm %>% filter(OTU %in% sig) %>%
  select(Species,everything()) %>% select(-OTU)
names(temp)[1] = 'Species'
temp = temp %>% pivot_wider(names_from = Species, values_from = Abundance)

# Transform everything to 0-1 range
temp = temp %>% mutate_at(c(sig,all_of(c('bristol','depth'))), function(x) (x-min(x,na.rm = T))/(max(x,na.rm = T)-min(x,na.rm = T)))

df = list(data = temp,sig = sig)
rm(psm,t,sig_multi,sig_uni,cnt,prep_tsv,s)
