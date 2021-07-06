library(nnet)

########## multinomial model
t=quantile(data_cross %>% filter(dur_positive>0) %>% pull(dur_positive),0.95) %>% as.numeric()

data_mlr<-data_sum %>% 
  mutate(class3=as.factor(class3)) %>% 
  mutate(sympt=ifelse(symptom=="No","No symptoms",ifelse(symptom_classic=="No","Other symptoms","Classic symptoms"))%>% factor(levels=c("No symptoms","Other symptoms","Classic symptoms"))) %>% 
  mutate(dur_21=ifelse(dur_positive>t,t,dur_positive)) 

m<-nnet::multinom(class3~ns(age,knots = c(50),Boundary.knots = c(20,80))+Sex+ethnicity+
                    lthc+hcw+ct+sympt+multi_positive+dur_21,
                  data=data_mlr)

summary(m)

coef=t(coef(m)) %>% as.data.frame() %>% rename(class2=`2`,class3=`3`)
confint=confint(m) %>% as.data.frame() %>% 
  rename(class2l=`2.5 %.2`,class3l=`2.5 %.3`,
         class2u=`97.5 %.2`,class3u=`97.5 %.3`)

z<-summary(m)$coefficients/summary(m)$standard.errors
p<- t((1-pnorm(abs(z),0,1))*2) %>% as.data.frame() %>% rename(p2=`2`,p3=`3`)

OR=exp(cbind(coef, confint)) %>% 
  as.data.frame() %>% 
  round(2) %>% 
  rownames_to_column() %>% 
  mutate(class2CI=paste0(class2l,"-",class2u),
         class3CI=paste0(class3l,"-",class3u)) %>% 
  cbind(p) %>% 
  select(rowname,class2,class2CI,p2,class3,class3CI,p3)


################### logistic regression
data_lr=data_cross %>% 
  filter(class3!=2) %>% 
  mutate(class=ifelse(class3==1,0,1) %>% as.factor()) %>% 
  mutate(evidence=ifelse(ct<=32&ctp!=1,1,0))

m<-glm(class~ns(age,knots = c(50),Boundary.knots = c(20,80))+Sex+ethnicity+
         lthc+hcw+fever+muscle_ache_myalgia+fatigue_weakness+
         sore_throat+cough+shortness_of_breath+
         headache+nausea_vomiting+abdominal_pain+
         diarrhoea+loss_of_smell+loss_of_taste,
       data=data_lr, family=binomial)
summary(m)

pvalue=summary(m)$coefficients[,4]

OR=exp(cbind(coef(m), confint(m))) %>% 
  as.data.frame() %>% 
  round(2) %>% 
  rownames_to_column() %>% 
  mutate(CI=paste0(`2.5 %`,"-",`97.5 %`)) %>% 
  cbind(p=pvalue) 


################### logistic regression (restricted)
m<-glm(class~ns(age,knots = 50,Boundary.knots = c(20,80))+Sex+ethnicity+
         lthc+hcw+fever+muscle_ache_myalgia+fatigue_weakness+
         sore_throat+cough+shortness_of_breath+
         headache+nausea_vomiting+abdominal_pain+
         diarrhoea+loss_of_smell+loss_of_taste,
       data=data_lr %>% filter(!(evidence==0&class3==3)), family=binomial)
summary(m)

pvalue=summary(m)$coefficients[,4]

OR=exp(cbind(coef(m), confint(m))) %>% 
  as.data.frame() %>% 
  round(2) %>% 
  rownames_to_column() %>% 
  mutate(CI=paste0(`2.5 %`,"-",`97.5 %`)) %>% 
  cbind(p=pvalue) 

