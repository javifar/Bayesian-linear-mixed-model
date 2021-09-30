library(brms)

############################
start=56
end=180
dur=end-start

df<-lcmm1 %>% 
  filter(class3==1) %>% 
  mutate(time=as.numeric(time)) %>% 
  group_by(participant_id) %>%
  mutate(Before=ifelse(any(time<0&assay>=23),1,0)) %>% 
  filter(Before==0) %>% 
  ungroup() %>% 
  select(participant_id, class3, assay, assay_mabs, assay_log, time, age, sex, ethnicity, symptom,
         ct, lthc) %>% 
  filter(time>=start & time<=end) %>% 
  mutate(time=time-start) %>% 
  filter(assay>=23) %>% 
  group_by(participant_id) %>% 
  mutate(neg=ifelse(time<=14&assay<42,1,0)) %>% 
  filter(neg==0) %>% 
  mutate(ID=cur_group_id()) %>% 
  mutate(ID=as.factor(ID)) %>% 
  ungroup() %>% 
  mutate(cen1=ifelse(assay_mabs>800,1,0)) %>% 
  mutate(age2=age-43,
         ct2=ct-22,
         ethnicity=as.factor(ethnicity)) %>% 
  mutate(age10=age2/10)


###########fit a censored baseline model
mean_assay=list(Intercept=8)
inits<-list(mean_assay,mean_assay,mean_assay,mean_assay)

priors=c(
  set_prior("normal(7.4,2.8)",class="Intercept"),
  set_prior("normal(0,0.1)",coef="time",class="b"),
  set_prior("normal(0,0.01)",coef="time",class="sd",group="ID"),
  set_prior("normal(0,1)",coef="Intercept",class="sd",group="ID"),
  set_prior("lkj_corr_cholesky(2)",class="L")
)

fit<- brm(formula=log2(assay)|cens(cen1)~1+time+(1+time|ID),
                 data=df,cores=4,family=gaussian(),
                 prior =priors,
                 chains = 4, iter=4000, warmup = 2000, seed=42,inits=inits,
                 init_r=10, control = list(adapt_delta=0.99,max_treedepth=20))

a=ggpredict(fit,c("time [all]")) %>% mutate(pred=2^predicted,lc=2^`conf.low`,uc=2^`conf.high`)
plot<-ggplot(a)+
  geom_line(data=df, aes(time, assay,group=ID),color="gray",alpha=0.3)+
  geom_point(data=df, aes(time, assay),alpha=0.2)+
  geom_ribbon(aes(x=x,ymin=lc,ymax=uc),alpha=1,fill="red")+
  geom_line(aes(x,pred),color="black",lwd=0.6)+
  geom_hline(yintercept = 42, color="black", linetype=2)+
  geom_hline(yintercept = 28, color="black", linetype=3)+
  scale_y_log10(limits = c(1.8,800), breaks=c(2,6,15,28,42,100,300,800))+
  scale_x_continuous(limits = c(0,dur), breaks = seq(0,dur,30))+
  labs(title="", x="Time from peak level (days)",y="Predicted anti-spike IgG levels (ng/ml)")+
  theme_light()


###################################age model
priors=c(
  set_prior("normal(7.4,2.8)",class="Intercept"),
  set_prior("normal(0,0.1)",coef="time",class="b"),
  set_prior("normal(0,2)",coef="age10",class="b"),
  set_prior("normal(0,0.1)",coef="time:age10",class="b"),
  set_prior("normal(0,0.01)",coef="time",class="sd",group="ID"),
  set_prior("normal(0,1)",coef="Intercept",class="sd",group="ID"),
  set_prior("lkj_corr_cholesky(2)",class="L",group="ID")
)

fit_age <- brm(formula=log2(assay)|cens(cen1)~1+time*age10+(1+time|ID),
               data=df,cores=4,family=gaussian(),
               prior =priors,
               chains = 4, iter=4000, warmup = 2000, seed=42,inits=inits,
               init_r=10, control = list(adapt_delta=0.99,max_treedepth=20))


####################################sex model
priors=c(
  set_prior("normal(7.4,2.8)",class="Intercept"),
  set_prior("normal(0,0.1)",coef="time",class="b"),
  set_prior("normal(0,6)",coef="sex",class="b"),
  set_prior("normal(0,0.1)",coef="time:sex",class="b"),
  set_prior("normal(0,0.01)",coef="time",class="sd",group="ID"),
  set_prior("normal(0,1)",coef="Intercept",class="sd",group="ID"),
  set_prior("lkj_corr_cholesky(2)",class="L",group="ID")
  
)

fit_sex <- brm(formula=log2(assay)|cens(cen1)~1+time*sex+(1+time|ID),
               data=df,cores=4,family=gaussian(),
               prior = priors,
               chains = 4, iter=4000, warmup = 2000, seed=42,inits=inits,
               init_r=10, control = list(adapt_delta=0.99,max_treedepth=20))


###################################ethnicity model
priors=c(
  set_prior("normal(7.4,2.8)",class="Intercept"),
  set_prior("normal(0,0.1)",coef="time",class="b"),
  set_prior("normal(0,10)",coef="ethnicity",class="b"),
  set_prior("normal(0,0.2)",coef="time:ethnicity",class="b"),
  set_prior("normal(0,0.01)",coef="time",class="sd",group="ID"),
  set_prior("normal(0,1)",coef="Intercept",class="sd",group="ID"),
  set_prior("lkj_corr_cholesky(2)",class="L",group="ID")
  
)

fit_ethnicity <- brm(formula=log2(assay)|cens(cen1)~1+time*ethnicity+(1+time|ID),
                     data=df,cores=4,family=gaussian(),
                     prior = priors,
                     chains = 4, iter=4000, warmup = 2000, seed=42,inits=inits,
                     init_r=10, control = list(adapt_delta=0.99,max_treedepth=20))

###################################lthc model
priors=c(
  set_prior("normal(7.4,2.8)",class="Intercept"),
  set_prior("normal(0,0.1)",coef="time",class="b"),
  set_prior("normal(0,7)",coef="lthc1",class="b"),
  set_prior("normal(0,0.2)",coef="time:lthc1",class="b"),
  set_prior("normal(0,0.01)",coef="time",class="sd",group="ID"),
  set_prior("normal(0,1)",coef="Intercept",class="sd",group="ID"),
  set_prior("lkj_corr_cholesky(2)",class="L",group="ID")
)

fit_lthc <- brm(formula=log2(assay)|cens(cen1)~1+time*lthc+(1+time|ID),
                data=df,cores=4,family=gaussian(),
                prior = priors,
                chains = 4, iter=4000, warmup = 2000, seed=42,inits=inits,
                init_r=10, control = list(adapt_delta=0.99,max_treedepth=20))


################################ct model
priors=c(
  set_prior("normal(7.4,2.8)",class="Intercept"),
  set_prior("normal(0,0.1)",coef="time",class="b"),
  set_prior("normal(0,1)",coef="ct2",class="b"),
  set_prior("normal(0,0.1)",coef="time:ct2",class="b"),
  set_prior("normal(0,0.01)",coef="time",class="sd",group="ID"),
  set_prior("normal(0,1)",coef="Intercept",class="sd",group="ID"),
  set_prior("lkj_corr_cholesky(2)",class="L",group="ID")
  
)

fit_ct <- brm(formula=log2(assay)|cens(cen1)~1+time*ct2+(1+time|ID),
              data=df,cores=4,family=gaussian(),
              prior = priors,
              chains = 4, iter=4000, warmup = 2000, seed=42,inits=inits,
              init_r=10, control = list(adapt_delta=0.99,max_treedepth=20))


####################################symptom model
priors=c(
  set_prior("normal(7.4,2.8)",class="Intercept"),
  set_prior("normal(0,0.1)",coef="time",class="b"),
  set_prior("normal(0,6)",coef="symptom1",class="b"),
  set_prior("normal(0,0.1)",coef="time:symptom1",class="b"),
  set_prior("normal(0,0.01)",coef="time",class="sd",group="ID"),
  set_prior("normal(0,1)",coef="Intercept",class="sd",group="ID"),
  set_prior("lkj_corr_cholesky(2)",class="L",group="ID")
  
)
fit_sympt <- brm(formula=log2(assay)|cens(cen1)~1+time*symptom+(1+time|ID),
                 data=df,cores=4,family=gaussian(),
                 prior = priors,
                 chains = 4, iter=4000, warmup = 2000, seed=42,inits=inits,
                 init_r=10, control = list(adapt_delta=0.99,max_treedepth=20))


#######################################multivariable model
priors=c(
  set_prior("normal(7.4,2.8)",class="Intercept"),
  set_prior("normal(0,0.1)",coef="time",class="b"),
  set_prior("normal(0,2)",coef="age10",class="b"),
  set_prior("normal(0,6)",coef="sex1",class="b"),
  set_prior("normal(0,10)",coef="ethnicity1",class="b"),
  set_prior("normal(0,7)",coef="lthc1",class="b"),
  set_prior("normal(0,1)",coef="ct2",class="b"),
  set_prior("normal(0,6)",coef="symptom1",class="b"),
  
  set_prior("normal(0,0.1)",coef="time:age10",class="b"),
  set_prior("normal(0,0.1)",coef="time:sex1",class="b"),
  set_prior("normal(0,0.2)",coef="time:ethnicity1",class="b"),
  set_prior("normal(0,0.2)",coef="time:lthc1",class="b"),
  set_prior("normal(0,0.1)",coef="time:ct2",class="b"),
  set_prior("normal(0,0.1)",coef="time:symptom1",class="b"),
  
  set_prior("normal(0,0.01)",coef="time",class="sd",group="ID"),
  set_prior("normal(0,1)",coef="Intercept",class="sd",group="ID"),
  set_prior("lkj_corr_cholesky(2)",class="L",group="ID")
)

fit_all <- brm(formula=log2(assay)|cens(cen1)~1+time*age10+time*sex+time*ethnicity+time*lthc+
                 time*ct2+time*symptom+(1+time|ID),
               data=df,cores=4,family=gaussian(),
               prior = priors,
               chains = 4, iter=4000, warmup = 2000, seed=42,inits=inits,
               init_r=10, control = list(adapt_delta=0.99,max_treedepth=20))


############bi-exponential model

df=df %>% 
  mutate(time2=ifelse(time<=28,0,time-28)) %>% 
  mutate(t1=(time<28)*time+(time>=28)*28,
         t2=(time>=28)*(time-28))

mean_assay=list(Intercept=8)
inits<-list(mean_assay,mean_assay,mean_assay,mean_assay)

priors=c(
  set_prior("normal(7.5,1)",class="Intercept"),
  set_prior("normal(0,0.1)",coef="t1",class="b"),
  set_prior("normal(0,0.1)",coef="t2",class="b"),
  set_prior("cauchy(0,0.01)",coef="t1",class="sd",group="ID"),
  set_prior("cauchy(0,0.01)",coef="t2",class="sd",group="ID"),
  set_prior("cauchy(0,1)",coef="Intercept",class="sd",group="ID"),
  set_prior("cauchy(0,0.5)",class="sigma"),
  set_prior("lkj_corr_cholesky(2)",class="L")
)

fit_bi<- brm(formula=log2(assay)|cens(cen1)~1+t1+t2+(1+t1+t2|ID),
                 data=df,cores=4,family=gaussian(),
                 prior =priors,
                 chains = 4, iter=4000, warmup = 2000, seed=42,inits=inits,
                 init_r=10, control = list(adapt_delta=0.99))


#################spline model
mean_assay=list(Intercept=8)
inits<-list(mean_assay,mean_assay,mean_assay,mean_assay)

priors=c(
  set_prior("normal(7.4,2.8)",class="Intercept"),
  set_prior("normal(0,0.1)",coef="nstimeknotsEQc3070Boundary.knotsEQc51101",class="b"),
  set_prior("normal(0,0.1)",coef="nstimeknotsEQc3070Boundary.knotsEQc51102",class="b"),
  set_prior("normal(0,0.1)",coef="nstimeknotsEQc3070Boundary.knotsEQc51103",class="b"),

  set_prior("cauchy(0,0.01)",coef="nstimeknotsEQc3070Boundary.knotsEQc51101",class="sd",group="ID"),
  set_prior("cauchy(0,0.5)",coef="nstimeknotsEQc3070Boundary.knotsEQc51102",class="sd",group="ID"),
  set_prior("cauchy(0,0.5)",coef="nstimeknotsEQc3070Boundary.knotsEQc51103",class="sd",group="ID"),
  set_prior("cauchy(0,1)",coef="Intercept",class="sd",group="ID"),
  set_prior("cauchy(0,0.5)",class="sigma"),
  set_prior("lkj_corr_cholesky(2)",class="L",group="ID")
)

fit_spline <- brm(formula=log2(assay)|cens(cen1)~
                 1+ns(time,knots =c(30,70),Boundary.knots = c(5,110))+
                 (1+ns(time,knots =c(30,70),Boundary.knots = c(5,110))|ID),
               data=df,cores=4,family=gaussian(),prior = priors,
               chains = 4, iter=4000, warmup = 2000, seed=42,inits=inits,
               init_r=10, control = list(adapt_delta=0.95,max_treedepth=15))


###########compare model fit
loo=loo(fit)
loo_spline=loo(fit_spline)
loo_bi=loo(fit_bi)

loo_compare(loo,loo_spline,loo_bi)