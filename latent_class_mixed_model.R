library(tidyverse)
library(lcmm)

####################
lcmm1<-data %>% 
  select(participant_id,time,cross,source,tnt_before,swab,visit_date,result_mk,assay_mabs,assay_log,
         assay,age,sex,ethnicity,hcw,lthc,index_date,first_cis_pos_date,ct,ctp,ct_pattern,symptom,
         symptom_classic,ever_negative_before,dur_negative_before,dur_positive,dur_positive_21,dur_tnt) %>% 
  group_by(participant_id) %>%
  mutate(ID=cur_group_id()) %>% 
  mutate(ID=as.factor(ID)) %>% 
  ungroup()  %>% 
  as.data.frame()

# add in time:
time <- data.frame(splines::ns(lcmm1$time,knots = c(-10,30,60),Boundary.knots = c(-60,140)))
names(time) <- c("TIME1","TIME2","TIME3","TIME4")
lcmm1 <- cbind(lcmm1,time)

# add in age: 
age <- data.frame(splines::ns(lcmm1$age,knots = c(50),Boundary.knots = c(20,80)))
names(age) <- c("AGE1","AGE2")
lcmm1 <- cbind(lcmm1,age)

lcmm1$ID<-as.numeric(lcmm1$ID)

model1 <- lcmm(assay_log ~ TIME1 + TIME2 + TIME3 + TIME4,
               random = ~ TIME1 + TIME2 + TIME3 + TIME4, 
               subject='ID',maxiter = 1000,
               data = lcmm1)

model2 <- lcmm(assay_log ~ TIME1 + TIME2 + TIME3 + TIME4, 
               mixture = ~ TIME1 + TIME2 + TIME3 + TIME4,
               classmb = ~ AGE1 + AGE2 +lthc+ct+symptom,
               random = ~ TIME1 + TIME2 + TIME3 + TIME4, maxiter=10000,
               subject='ID', ng=2, B=model1,
               data = lcmm1)

model3 <- lcmm(assay_log ~ TIME1 + TIME2 + TIME3 + TIME4, 
               mixture = ~ TIME1 + TIME2 + TIME3 + TIME4,
               classmb = ~ AGE1 + AGE2 +lthc+ct+symptom, 
               random = ~ TIME1 + TIME2 + TIME3 + TIME4, maxiter=10000,
               subject='ID', ng=3, B=model1,
               data = lcmm1)

t <- unique(lcmm1$time) %>% sort()
a <- unique(lcmm1$age) %>% sort()

plot_data <- expand.grid(time=t, age=a)
time_spline <- lcmm1 %>% select(time, TIME1:TIME4) %>% unique() 
time_spline$time <- as.numeric(time_spline$time)
age_spline <- lcmm1 %>% select(age, AGE1:AGE2) %>% unique()

plot_data <- plot_data %>% 
  left_join(time_spline) %>% 
  left_join(age_spline) %>% 
  unique()
plot_data$lthc <- 0
plot_data$ct<-26
plot_data$symptom<-1

pred <- predictY(model3, plot_data[plot_data$age==40,], draws = T)
pred3 <- cbind(pred$pred,pred$times) %>% 
  pivot_longer(cols=1:9) %>%
  separate(col="name",into=c("a","m","class"),sep="_") %>% 
  mutate(type=ifelse(m==50,"pred",ifelse(m==2.5,"ll","ul"))) %>% 
  select(-a,-m) %>% 
  pivot_wider(names_from = type,values_from=value) %>% 
  mutate(Class=ifelse(class=="class1","1",ifelse(class=="class2","2","3"))) %>% 
  mutate(assay=10^pred,
         lower=10^ll,
         upper=10^ul)

cbpalette=c("#B24745FF","#DF8F44FF", "#00A1D5FF","#374E55FF","#F0E442", "#6A6599FF")

plot<-ggplot(pred3,aes(time,assay))+
  geom_line(aes(color=Class))+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=Class),alpha=0.2)+
  geom_hline(yintercept = 42, color="black", linetype=2)+
  geom_hline(yintercept = 28, color="black", linetype=3)+
  scale_y_log10(limits = c(1.8,800), breaks=c(2,6,15,28,42,100,300,800))+
  scale_x_continuous(limits = c(-90,180), breaks = seq(-90,180,30))+
  labs(title="", x="Time from first swab positive (days)",y="Predicted anti-spike IgG levels (ng/ml)")+
  theme_light()+
  scale_fill_manual(values=cbpalette, labels=c("1: seroconverted","2: possible late/reinfection","3: seronegative non-responders"))+
  scale_color_manual(values=cbpalette, labels=c("1: seroconverted","2: possible late/reinfection","3: seronegative non-responders"))+
  theme(legend.position = "bottom")


people <- as.data.frame(model3$pprob) %>% 
  mutate(prob=ifelse(class==1,prob1,ifelse(class==2, prob2, prob3))*100)

lcmm1 <- lcmm1 %>% left_join(people, by="ID") %>% rename(class3=class)
lcmm1$class3 <- as.factor(lcmm1$class3)


data_sum <- lcmm1 %>% 
 mutate(Sex=ifelse(sex==1,"Male","Female") %>% as.factor(),
         hcw=ifelse(hcw==1,"Yes","No") %>% as.factor(),
         lthc=ifelse(lthc==1,"Yes","No") %>% as.factor(),
         ethnicity=ifelse(ethnicity==0,"White","Non-white") %>% factor(levels=c("White","Non-white")),
         cross=factor(cross,levels = c("No symptom&Ct>=30","No symptom&Ct<30","Symptom&Ct>=30","Symptom&Ct<30")),
         multi_positive=ifelse(dur_positive==0,"No","Yes") %>% as.factor(),
         symptom=ifelse(symptom==1,"Yes","No") %>% as.factor(),
         symptom_classic=ifelse(symptom_classic==1,"Yes","No") %>% as.factor(),
         dur_tnt=ifelse(dur_tnt==0,NA,dur_tnt),
         dur_positive=ifelse(dur_positive==0,NA,dur_positive),
         ctp=ifelse(ctp==1,"1/2/3",ifelse(ctp==2,"4","5/6/7"))%>% factor(levels = c("5/6/7","4","1/2/3")),
         ever_negative_before=ifelse(ever_negative_before==1,"Yes","No") %>% factor()%>% relevel(ref="Yes")) %>% 
  select(participant_id,ID,age,source,tnt_before, class3, Sex, ethnicity,lthc, hcw, ct,ctp,symptom,symptom_classic,
         cross,multi_positive,dur_positive,ever_negative_before,dur_negative_before,cross,dur_tnt,prob) %>% 
  unique()



