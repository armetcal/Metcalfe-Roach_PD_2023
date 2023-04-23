# LMER Functions

run_lmer = function(df,species,variable,log=F,pc = 1, meanscale=T){
  # df=pd %>% filter(clinstate==1);species = sig[1];variable='mds3.total';pc=1;meanscale=T;log=F
  names(df)[which(names(df)==variable)] = 'test'
  
  if(log){df$test = log10(df$test+pc)}
  
  # Set range to 0-1
  df$test = rescale(df$test)
  
  if(meanscale){df$test = df$test-mean(df$test,na.rm=T)}
  
  a = lmer(test ~ Abundance*months_last_visit + (months_last_visit|Sample), 
           data = df %>% filter(Species==species), na.action = 'na.omit',REML=F, 
           control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))) %>% 
    summary %>% .$coefficients %>% as.data.frame() %>% rownames_to_column('Variable') %>% mutate(Taxon = species,Test=variable)
  a$Variable = sapply(a$Variable,function(x)str_replace_all(x,'test',variable))
  allres <<- allres %>% rbind(a)
}

# Multivar for interaction
run_lmer_multi = function(df,species,variable,log=F,pc = 1, meanscale=T){
  names(df)[which(names(df)==variable)] = 'test'
  if(log){df$test = log10(df$test+pc)}
  # Set range to 0-1
  df$test = rescale(df$test)
  if(meanscale){df$test = df$test-mean(df$test,na.rm=T)}
  a = lmer(test ~ Sex + ent_boolean + bristol + depth + Abundance*months_last_visit + (months_last_visit|Sample), 
           data = df %>% filter(Species==species), na.action = 'na.omit',REML=F, 
           control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))) %>% 
    summary %>% .$coefficients %>% as.data.frame() %>% rownames_to_column('Variable') %>% mutate(Taxon = species,Test=variable)
  a$Variable = sapply(a$Variable,function(x)str_replace_all(x,'test',variable))
  allresmulti <<- allresmulti %>% rbind(a)
}

# For abundance, contains disease duration
run_lmer_multi2 = function(df,species,variable,log=F,pc = 1, meanscale=T){
  names(df)[which(names(df)==variable)] = 'test'
  if(log){df$test = log10(df$test+pc)}
  # Set range to 0-1
  df$test = rescale(df$test)
  if(meanscale){df$test = df$test-mean(df$test,na.rm=T)}
  a = lmer(test ~ Sex + ent_boolean + bristol + depth + Abundance*months_last_visit + (months_last_visit|Sample), 
           data = df %>% filter(Species==species), na.action = 'na.omit',REML=F, 
           control = lmerControl(optimizer="bobyqa", optCtrl=list(maxfun=100000))) %>% 
    summary %>% .$coefficients %>% as.data.frame() %>% rownames_to_column('Variable') %>% mutate(Taxon = species,Test=variable)
  a$Variable = sapply(a$Variable,function(x)str_replace_all(x,'test',variable))
  allresmulti <<- allresmulti %>% rbind(a)
}

edit_lmer = function(df,v){
  
  df = df %>% filter(Test %in% v)
  
  # Edit levels for plotting
  df$Test[df$Test=='mds.total'] = 'Total UPDRS'
  df$Test[df$Test=='mds1.total'] = 'UPDRS 1'
  df$Test[df$Test=='mds2.total'] = 'UPDRS 2'
  df$Test[df$Test=='mds3.total'] = 'UPDRS 3'
  df$Test[df$Test=='mds4.total'] = 'UPDRS 4'
  df$Test[df$Test=='bdi_total'] = 'Depression'
  df$Test[df$Test=='fss_total'] = 'Fatigue'
  df$Taxon = str_remove_all(df$Taxon,'s__') %>% str_replace_all('_',' ')
  
  df$Test = factor(df$Test,levels = c('Total UPDRS','UPDRS 1','UPDRS 2','UPDRS 3','UPDRS 4','Depression','Fatigue'))
  df = df %>% droplevels() %>% as.data.frame()
  
  return(df)
}

plot_lmer = function(df,v,bsize=20){
  # df = allres.p; v=c('mds.total','mds1.total','mds2.total','mds3.total','mds4.total')
  
  df = edit_lmer(df,v)
  
  # Arrange data
  temp = df %>% filter(Test==(levels(.$Test)[1])) %>% arrange(Estimate)
  df$Taxon = factor(df$Taxon, levels = temp$Taxon)
  
  # Plot
  p = df %>% 
    mutate(Qval = ifelse(qval>0.05,'',
                         ifelse(qval>0.01,'*',
                                ifelse(qval>0.001,'**',
                                       ifelse(qval<=0.001,'***',''))))) %>% 
    mutate(Qval2 = ifelse(qval>0.1,'',
                          ifelse(qval>0.05,'+',''))) %>% 
    ggplot(aes(Test,Taxon,col = Estimate)) +
    geom_point(size =14) +
    theme_classic(base_size = bsize) +
    scale_color_gradient2(low="darkblue", high="darkred", guide="colorbar") +
    geom_text(aes(label=Qval),size = 12, col = 'white',nudge_y = -0.25) +
    geom_text(aes(label=Qval2),size = 12, col = 'white',nudge_y = 0) +
    ylab(NULL) + xlab(NULL) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  return(p)
}

# Below is for td pigd plots
plot_lmer2 = function(df,v,bsize=16){
  # df = allres.p; v=c('mds.total','mds1.total','mds2.total','mds3.total','mds4.total')
  
  df = edit_lmer(df,v)
  
  # Arrange data
  temp = df %>% filter(td.pigd == 'PIGD',Test==(levels(.$Test)[1])) %>% arrange(Estimate)
  df$Taxon = factor(df$Taxon, levels = temp$Taxon)
  
  # Plot
  p = df %>% #filter(Test == 'Total UPDRS') %>% 
    mutate(Qval = ifelse(qval>0.05,'',
                         ifelse(qval>0.01,'*',
                                ifelse(qval>0.001,'**',
                                       ifelse(qval<=0.001,'***',''))))) %>% 
    mutate(Qval2 = ifelse(qval>0.1,'',
                          ifelse(qval>0.05,'+',''))) %>% 
    rename(Category = td.pigd) %>% 
    ggplot(aes(Category,Taxon,col = Estimate)) +
    geom_point(size =14) +
    theme_classic(base_size = bsize) +
    scale_color_gradient2(low="darkblue", high="darkred", guide="colorbar") +
    geom_text(aes(label=Qval),size = 12, col = 'white',nudge_y = -0.25) +
    geom_text(aes(label=Qval2),size = 12, col = 'white',nudge_y = 0) +
    ylab(NULL) + xlab(NULL) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  if(length(v)>5) c = 4 else c = 3
    p = p + facet_wrap('Test',ncol=c)
  return(p)
}