# CRP correlations for Silke

a = ps %>% psmelt
a = a %>% select(-OTU,-Abundance) %>% select(-c(Kingdom:Species)) %>% unique

a %>% ggplot(aes(Status,Log2CRP,fill=Status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height=0,width=0.2) +
  theme_classic(base_size = 16) +
  facet_wrap("Sex") +
  ggpubr::stat_compare_means(comparisons = list(c('PD','Ctrl')),size=5)

a %>% mutate(Agegroup = ifelse(as.numeric(Age)<quantile(as.numeric(Age),0.6),'Young','Old')) %>% 
  ggplot(aes(Status,Log2CRP,fill=Status)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height=0,width=0.2) +
  theme_classic(base_size = 16) +
  facet_wrap(c("Sex",'Agegroup')) +
  ggpubr::stat_compare_means(comparisons = list(c('PD','Ctrl')),size=5)

a %>% filter(!is.na(Age.of.PD.onset), Status=='PD') %>% 
  mutate(Age.of.PD.onset = as.numeric(Age.of.PD.onset)) %>% 
  ggplot(aes(Log2CRP,Age.of.PD.onset)) +
  geom_point() + geom_smooth(method='lm',col='black') +
  theme_classic(base_size = 16) +
  facet_wrap("Sex") +
  ggpubr::stat_cor(method='spearman')

a %>%  
  mutate(Age = as.numeric(Age)) %>% 
  ggplot(aes(Age,Log2CRP,fill=Status,col=Status)) +
  geom_point() + geom_smooth(method='lm') +
  theme_classic(base_size = 16) +
  facet_wrap("Sex") +
  ggpubr::stat_cor(method='spearman')

summary(glm(Log2CRP~Age*Status,data=a %>% mutate(Age=as.numeric(Age)) %>% filter(Sex=='Female')))
summary(glm(Log2CRP~Age*Status,data=a %>% mutate(Age=as.numeric(Age)) %>% filter(Sex=='Male')))

a %>%  
  mutate(mindscore = as.numeric(mindscore)) %>% 
  ggplot(aes(mindscore,Log2CRP,fill=Status,col=Status)) +
  geom_point() + geom_smooth(method='lm') +
  theme_classic(base_size = 16) +
  facet_wrap("Sex") +
  ggpubr::stat_cor(method='spearman')
a %>%  
  mutate(mindscore = as.numeric(mindscore)) %>% 
  ggplot(aes(mindscore,Log2CRP)) +
  geom_point() + geom_smooth(method='lm') +
  theme_classic(base_size = 16) +
  facet_wrap("Sex") +
  ggpubr::stat_cor(method='spearman')