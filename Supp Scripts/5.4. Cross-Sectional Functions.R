# Cross-sectional functions

correlate_met = function(df,var,meth='spearman',grp){
  # df=crpbut;var='Log2CRP';meth='spearman';grp='All'
  test = apply(df[,sig],2,function(x){
    # x = df[,sig[1]] %>% unname() %>% unlist()
    ct = cor.test(x,df[[var]],method=meth)
    c = ct$p.value
    rho = ct$estimate
    return(c(c,rho))
  }) %>% as.data.frame(row.names = c('pval','rho')) %>% 
    t() %>% as.data.frame() %>% 
    mutate(method = meth,
           variable = var,
           group = grp,
           qval = p.adjust(pval,method='fdr')) %>% 
    arrange(pval) %>% 
    rownames_to_column('taxon') %>% 
    select(group,method,variable,taxon,rho,everything())
  if(meth=='pearson'){names(test)[names(test)=='rho']=='r'}
  return(test)
}

glm_met = function(df,var,grp){
  # df=met %>% filter(Status=='Ctrl');var='pcresol';grp='Ctrl'
  test = data.frame()
  for(i in 1:length(sig)){
    if(grp=='All'){
      g = summary(glm(rescale(df[[var]])~rescale(df[[sig[i]]])+df$Sex+df$Status+df$bristol+df$ent_boolean+df$depth))$coefficients
    } else if(grp=='PD'){
      g = summary(glm(rescale(df[[var]])~rescale(df[[sig[i]]])+df$Sex+df$bristol+df$ent_boolean+df$depth))$coefficients
    } else {
      g = summary(glm(rescale(df[[var]])~rescale(df[[sig[i]]])+df$Sex+df$bristol+df$depth))$coefficients
    }
    g = g %>% as.data.frame %>% rownames_to_column('Explanatory')
    g$Explanatory[g$Explanatory=='rescale(df[[sig[i]]])'] = sig[i]
    test = test %>% rbind(g %>% mutate(Taxon = sig[i]))
  }
  test = test %>% 
    rename(pval = 'Pr(>|t|)') %>% 
    mutate(Explanatory = str_remove(.$Explanatory,'df$'),
           Response = var,group = grp) %>% 
    group_by(Response) %>% 
    mutate(qval = p.adjust(pval,method = 'fdr')) %>% ungroup() %>% 
    arrange(pval) %>%
    select(group,Taxon,Response,Explanatory,everything())
  return(test)
}