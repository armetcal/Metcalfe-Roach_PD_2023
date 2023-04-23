# PR-ROC Functions

convert_prroc = function(x){
  # x = pr_train
  # print(str(x))
  x2 = x$curve %>% as.data.frame %>% `names<-`(c('Recall','Precision','Cutoff'))
  # x2 = x2 %>% ggplot(aes(Recall,Precision,col=Cutoff)) +
  #   geom_line(aes(col=Cutoff),size=1) + 
  #   scale_color_viridis_b() + 
  #   theme_classic(base_size=18) + 
  #   ylim(0,1) + xlim(0,1)
  x2 = x2 %>% ggplot(aes(Recall,Precision,col='black')) +
    geom_line(size=1) + 
    theme_classic(base_size=18) + 
    ylim(0,1) + xlim(0,1)
  return(x2)
}

run_ROC = function(dataset,y,x, name = NA,testsize=0.8,testfamily = 'binomial',weights=F){
  # dataset=temp;y='Status';x=indiv_univar$Name[1:4];name=4;testsize=0.8;testfamily = 'binomial'
  require(ROCR)
  require(PRROC)
  set.seed(421)
  # Divide data
  index <- sample(nrow(dataset),nrow(dataset)*testsize)
  train = dataset[index,]
  test = dataset[-index,]
  # Run logistic regression on training dataset
  f = paste0(y,'~',paste(x,collapse = '+'))
  
  if(weights){
    model = glm(as.formula(f),family = testfamily,data = train, weights = train$weight)
  } else {
    model = glm(as.formula(f),family = testfamily,data = train)
  }
  
  # Within-dataset ROC (less important)
  pred_model_train<- predict(model, type="response")
  pred_model_test<- predict(model, newdata = test, type="response")
  # Precision and recall curve and its AUC is more appropriate for imbalanced data.
  # Extract y variable levels. 1st is the level being investigated, 2nd is the control/baseline.
  ylevels = rev(levels(dataset[[y]]))
  # Training dataset
  score1= pred_model_train[train[[y]]==ylevels[1]]
  score0= pred_model_train[train[[y]]==ylevels[2]]
  pr_train= pr.curve(score1, score0, curve = T)
  p_train_pr = convert_prroc(pr_train)
  # Test prediction: 
  score1.test= pred_model_test[test[[y]]==ylevels[1]]
  score0.test= pred_model_test[test[[y]]==ylevels[2]]
  pr_test= pr.curve(score1.test, score0.test, curve = T)
  p_test_pr = convert_prroc(pr_test)
  # Format Results
  ROC = tibble(Name = name,
               Test = c('PR-ROC','PR-ROC'),
               Data = c('Train','Test'),
               ROC = c(pr_train$auc.integral,pr_test$auc.integral),
               `n (Total)` = c(nrow(train),nrow(test)),
               `n1` = c(nrow(train[train[[y]]==ylevels[1],]),nrow(test[test[[y]]==ylevels[1],])),
               `n0` = c(nrow(train[train[[y]]==ylevels[2],]),nrow(test[test[[y]]==ylevels[2],])))
  names(ROC)[names(ROC)=='n1'] = paste0('n (',ylevels[1],')')
  names(ROC)[names(ROC)=='n0'] = paste0('n (',ylevels[2],')')
  Plots = list(`PR-ROC (Train)` = p_train_pr,`PR-ROC (Test)` = p_test_pr)
  return(list(Model = model,Formula = f, ROC = ROC, Plots = Plots))
}

prep_plots = function(data,plotvar,plotheights,comparisons=NULL){
  # data=temp;plotvar='Status';plotheights=c(5,3,5);comparisons=NULL
  # For nice plotting
  correct_names = tibble(Term = c('SexMale','bristol','depth'),correct=c('Sex (Male)','Bristol','Seq Depth'))
  
  # Univar ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ntaxa = comb_univar$Name[1] # the top combination of univar taxa
  taxa = indiv_univar$Name[1:ntaxa] %>% str_remove_all('`') # Taxa names
  
  # Estimates for each taxon
  m = list_comb_univar[[ntaxa]]$Model$coefficients %>% as.data.frame() %>% rownames_to_column('Species') %>%
    mutate(Species=str_remove_all(Species,'`'))
  names(m)[2] = 'Estimate'

  # Attach estimates, calculate predicted y value
  modeldf = data %>% select(Sample, !!sym(plotvar), all_of(taxa)) %>%
    pivot_longer(cols = -c(Sample,!!sym(plotvar)), names_to = 'Species',values_to = 'value') %>%
    left_join(m) %>%
    mutate(value = value*Estimate) %>% select(-Estimate) %>%
    group_by(Sample,!!sym(plotvar)) %>%
    summarize(value = sum(value),.groups='keep') %>% ungroup
  
  lfc1 = modeldf %>% group_by(!!sym(plotvar)) %>% summarize(value = mean(value)) %>% ungroup %>% 
    pivot_wider(names_from = !!sym(plotvar), values_from = value) %>% unname() %>% unlist()
  lfc1 = (lfc1[2]/lfc1[1]) %>% log2 %>% round(2)
  
  # Plot y value vs predicted value
  p1 = modeldf %>% ggplot(aes(!!sym(plotvar),value,fill = !!sym(plotvar))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.2) +
    theme_classic(base_size = 18) + 
    xlab(NULL) + ylab('Model Fit') + ggtitle('PWY') + ggeasy::easy_center_title() +
    ggpubr::stat_compare_means(label = 'p.format', size = 5.5,comparisons=comparisons) +
    theme(legend.position='none') +
    scale_y_continuous(expand = expansion(mult = 0.15))
  
  # PR-AUC for the selected model
  auc = list_comb_univar[[ntaxa]]$ROC %>% filter(Test == 'PR-ROC', Data == 'Test') %>% .$ROC %>% round(2)
  p111 = list_comb_univar[[ntaxa]]$Plots$`PR-ROC (Test)` + 
    annotate("text", label = paste0('PR AUC = ',auc), x = 0.5, y = 0.6, size = 6, colour = "black")
  
  # Estimates of each taxon
  p1111 = m %>% filter(Species!='(Intercept)') %>% arrange(-Estimate) %>% 
    mutate(Species=str_remove(Species,'s__') %>% str_wrap(30)) %>% 
    mutate(Species = factor(Species,levels = .$Species)) %>% 
    ggplot(aes(y = Estimate, x = Species, fill = Estimate<0)) + geom_col() +
    theme_classic(base_size=14) +
    theme(legend.position = 'none') +xlab(NULL) + ylab('Model Estimate') +
    theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1))
  
  # Multivar~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ntaxa = comb_multi$Name[1]
  taxa = indiv_multi$Name[1:ntaxa] %>% str_remove_all('`')
  
  m = list_comb_multi[[ntaxa]]$Model$coefficients %>% as.data.frame() %>% rownames_to_column('Term') %>%
    mutate(Term=str_remove_all(Term,'`'))
  names(m)[2] = 'Estimate'

  modeldf = data %>% select(Sample, !!sym(plotvar),
                            all_of(c('bristol','Sex','depth')),all_of(str_remove_all(taxa,'`'))) %>%
    mutate(SexMale = ifelse(Sex=='Female',0,1)) %>%
    select(-Sex) %>%
    pivot_longer(cols = -c(Sample,!!sym(plotvar)), names_to = 'Term',values_to = 'value') %>%
    left_join(m) %>%
    mutate(value = value*Estimate) %>% select(-Estimate) %>%
    group_by(Sample,!!sym(plotvar)) %>%
    summarize(value = sum(value),.groups='keep') %>% ungroup
  
  lfc2 = modeldf %>% filter(!is.na(value)) %>% group_by(!!sym(plotvar)) %>% 
    summarize(value = mean(value)) %>% ungroup %>% 
    pivot_wider(names_from = !!sym(plotvar), values_from = value) %>% unname() %>% unlist()
  lfc2 = (lfc2[2]/lfc2[1]) %>% log2 %>% round(2)
  
  p2 = modeldf %>% ggplot(aes(!!sym(plotvar),value,fill = !!sym(plotvar))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.2) +
    theme_classic(base_size = 18) + 
    xlab(NULL) + ylab('Model Fit') + ggtitle('PWY and Covar') + ggeasy::easy_center_title() +
    ggpubr::stat_compare_means(label = 'p.format', size = 5.5,comparisons=comparisons) +
    theme(legend.position='none') +
    scale_y_continuous(expand = expansion(mult = 0.15))
  
  auc = list_comb_multi[[ntaxa]]$ROC %>% filter(Test == 'PR-ROC', Data == 'Test') %>% .$ROC %>% round(2)
  p222 = list_comb_multi[[ntaxa]]$Plots$`PR-ROC (Test)` + 
    annotate("text", label = paste0('PR AUC = ',auc), x = 0.5, y = 0.6, size = 6, colour = "black")
  
  p2222 = m %>% filter(Term!='(Intercept)') %>% arrange(-Estimate) %>% 
    left_join(correct_names) %>% mutate(Term = ifelse(!is.na(correct),correct,Term)) %>% select(-correct) %>% 
    mutate(Term=str_remove(Term,'s__') %>% str_wrap(30)) %>% 
    mutate(Term = factor(Term,levels = .$Term)) %>% 
    ggplot(aes(y = Estimate, x = Term, fill = Estimate<0)) + geom_col() +
    theme_classic(base_size=14) +
    theme(legend.position = 'none') +xlab(NULL) + ylab('Model Estimate') +
    theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1))
  
  # Covars Only~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # The first model in the list is covariables only.
  m = list_indiv_multi[[1]]$Model$coefficients %>% as.data.frame() %>% rownames_to_column('Term')
  names(m)[2] = 'Estimate'
  
  modeldf = data %>% select(Sample, !!sym(plotvar),
                            all_of(c('bristol','Sex','depth'))) %>%
    mutate(SexMale = ifelse(Sex=='Female',0,1)) %>%
    # mutate(ent_booleanTRUE = ifelse(ent_boolean == T,1,0)) %>%
    select(-Sex) %>%
    pivot_longer(cols = -c(Sample,!!sym(plotvar)), names_to = 'Term',values_to = 'value') %>%
    left_join(m) %>%
    mutate(value = value*Estimate) %>% select(-Estimate) %>%
    group_by(Sample,!!sym(plotvar)) %>%
    summarize(value = sum(value),.groups='keep') %>% ungroup
  
  lfc3 = modeldf %>% filter(!is.na(value)) %>% group_by(!!sym(plotvar)) %>% 
    summarize(value = mean(value)) %>% ungroup %>% 
    pivot_wider(names_from = !!sym(plotvar), values_from = value) %>% unname() %>% unlist()
  lfc3 = (lfc3[2]/lfc3[1]) %>% log2 %>% round(2)
  
  p3 = modeldf %>% ggplot(aes(!!sym(plotvar),value,fill = !!sym(plotvar))) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.2) +
    theme_classic(base_size = 18) + 
    xlab(NULL) + ylab('Model Fit') + ggtitle('Covar') + ggeasy::easy_center_title() +
    ggpubr::stat_compare_means(label = 'p.format', size = 5.5,comparisons=comparisons) +
    theme(legend.position='none') +
    scale_y_continuous(expand = expansion(mult = 0.15))
  
  # Covars uniquely uses the list_INDIV_multi, as that is where the covar-only data is.
  auc = list_indiv_multi[[1]]$ROC %>% filter(Test == 'PR-ROC', Data == 'Test') %>% .$ROC %>% round(2)
  p333 = list_comb_univar[[ntaxa]]$Plots$`PR-ROC (Test)` + 
    annotate("text", label = paste0('PR AUC = ',auc), x = 0.5, y = 0.6, size = 6, colour = "black")
  
  p3333 = m %>% filter(Term!='(Intercept)') %>% arrange(-Estimate) %>% 
    left_join(correct_names) %>% mutate(Term = ifelse(!is.na(correct),correct,Term)) %>% select(-correct) %>%
    mutate(Term=str_remove(Term,'s__')) %>%
    mutate(Term = factor(Term,levels = .$Term)) %>% 
    ggplot(aes(y = Estimate, x = Term, fill = Estimate<0)) + geom_col() +
    theme_classic(base_size=14) +
    theme(legend.position = 'none') +xlab(NULL) + ylab('Model Estimate') +
    theme(axis.text.x = element_text(angle = 70, vjust = 1, hjust=1))
  
  # Boxplots
  g1 = ggarrange(plotlist = list(p3,p1,p2),ncol=3,nrow=1)
  # PR-ROC curves
  g3 = ggarrange(plotlist = list(p333,p111,p222),ncol=3,nrow=1, legend = 'none')
  # Estimates
  g4 = ggarrange(plotlist = list(p3333,p1111,p2222),ncol=3,nrow=1, legend = 'none')
  # Combined
  g = ggarrange(plotlist = list(g1,g3, g4), ncol = 1, nrow = 3, heights = plotheights)
  return(list(box=g1,roc=g3,est=g4,comb=g,lfc=c(pwy=lfc1,cov_pwy=lfc2,covar=lfc3)))
}