# Load dataset - Disease progression

ps = readRDS("C:/Users/armetcal/OneDrive - UBC/Grad School/Data/HUMAN COHORT/PD Functional Analysis/Reference Files/phyloseq_object_leafyadj_absolute_counts.rds")
ps@sam_data$bristol = as.numeric(ps@sam_data$bristol)
ps@sam_data$Sex = as.factor(ps@sam_data$Sex)
ps@sam_data$ent_boolean = ps@sam_data$entacapone>0
# Sequencing depth
ps@sam_data$depth = sample_sums(ps)

# CLR Transform Data
ps = microbiome::transform(ps,'clr')

long0 = read.csv('C:/Users/armetcal/OneDrive - UBC/Grad School/Data/HUMAN COHORT/Reference Files/Metadata/shayanandkiradata_nov22.csv') %>%
  select(Record.ID,redcap_event_name,Disease.duration,bdi_total,contains('mds'),fss_total,months_last_visit) %>% 
  rename(Sample = Record.ID) %>% 
  # Add an X to sample names to match the other datasets
  mutate(Sample = paste0('X_',Sample)) %>% 
  arrange(months_last_visit) %>% 
  group_by(Sample) %>% 
  # Set the baseline as 0 months
  mutate(months_last_visit = ifelse(redcap_event_name=='visit_1_arm_1',0,months_last_visit)) %>%
  # Select only those with longitudinal data available
  add_count() %>% filter(n>1) %>% select(-n) %>%
  ungroup

# This gives an estimate of progression for each person, calculated by taking the slope of mds.total~months_last_visit.
long_est = long0 %>% group_by(Sample) %>% 
  filter(!is.na(mds.total), !is.na(months_last_visit)) %>% 
  group_modify(~summary(glm(.$mds.total~.$months_last_visit))$coefficients %>% 
                 as.data.frame %>% rownames_to_column('Variable') %>% filter(Variable != '(Intercept)')) %>% 
  ungroup() %>% mutate(Variable = str_remove(Variable,'[.][$]'))

long_est_pd = long_est %>% left_join(ps@sam_data %>% as.matrix() %>% as.data.frame() %>% rownames_to_column('Sample') %>% select(Sample,Status)) %>% 
  filter(Status=='PD') %>% 
  mutate(Progression = ifelse(Estimate<=quantile(.$Estimate,1/3),'Slow',
                              ifelse(Estimate>quantile(.$Estimate,2/3),'Fast','Med'))) %>% 
  mutate(Progression = factor(Progression, levels = c('Slow','Med','Fast')))

hist(long_est_pd$Estimate)
table(long_est_pd$Progression)

long_all = long_est_pd %>% 
  left_join(temp %>% select(Sample,bristol,Sex,ent_boolean,depth,starts_with('s__'))) 
long = long_all %>% 
  filter(Progression != 'Med') %>% droplevels