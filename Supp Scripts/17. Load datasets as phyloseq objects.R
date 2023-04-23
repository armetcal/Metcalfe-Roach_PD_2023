source('../id_functions.R')

prep_tsv = function(df){
  names(df)[1] = 'temp'
  return(df %>% filter(!(temp %in% c('UNMAPPED','UNINTEGRATED','UNGROUPED'))) %>% column_to_rownames('temp') %>%
           select(-contains('Atcc'),-contains('Zymo'),-contains('blnk')) %>%
           `colnames<-`(reformat_ids(colnames(.)) %>% change_ids_to_longitudinal()))
}

# MetaCyc ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Load RPK Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df = fread('../../../Metagenomics/RPK counts/merged_pathabundance_unstratified.tsv') %>% prep_tsv()

# Load total metagenomic reads~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cnt = read_delim('../../../Metagenomics/ORIGINAL DATA/mqc_fastqc_sequence_counts_plot_1.txt') %>%
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
ps.meta = ps

###################################################################

# KO terms

# Load RPK Data~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
df = fread('../../../Metagenomics/RPK counts/merged_genefamilies_ko_renamed_unstratified.tsv') %>% prep_tsv()

# Load total metagenomic reads~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cnt = read_delim('../../../Metagenomics/ORIGINAL DATA/mqc_fastqc_sequence_counts_plot_1.txt') %>%
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
ps.ko = ps
