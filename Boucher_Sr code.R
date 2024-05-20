#updated fit for models - cleaned up 
library(ggplot2)
library(readxl)
library(ggpubr)
library(ggbeeswarm)
library(MuMIn)
library(glmmTMB)
library(TMB)
library(fitdistrplus)
library(metRology)
library(DHARMa)
library(emmeans)
library(effects)
library(carData)
library(tidyverse)
library(dplyr)
library(tidyr)
library(tibble)
library(devtools)
library(broom)
library(mvabund)
library(tweedie)
library(RcppGSL)

Sr_Chimps_3 <- read_excel("/Users/Renee/Desktop/Sr Paper/Tai_Chimps_Sr-3.xlsx")
Tai_SrEco <- read_excel("/Users/Renee_1/Desktop/Sr Paper/Tai_SrEco.xlsx")

#updated groups
Sr_Chimps_4 <- read_excel("/Users/Renee/Desktop/Sr Paper/Sr-chimps-4.xlsx") #BIJ, GOM post-1980 #does best!##

Sr_Chimps_7 <- read_excel("/Users/Renee_1/Desktop/Sr Paper/Sr-chimps-7.xlsx")

#differences in variance between groups
bartlett.test(Sr87_Sr86 ~ Time_Slice2, data = Tai_SrEco)


#we can accept the null that they are the same
Bartlett test of homogeneity of variances

data:  Sr87_Sr86 by Time_Slice2
Bartletts K-squared = 11.152, df = 6, p-value = 0.08378

#linear mixed models for Sr and teeth

#check data
str(Sr_Chimps_3)
hist(Sr_Chimps_3$Sr87_Sr86, breaks  = 50) # not normal distribution
plot(density(Sr_Chimps_3$Sr87_Sr86)) # not normal distribution

Sr_Chimps_4$Origin=as.factor(Sr_Chimps_4$Origin)
Sr_Chimps_4$Sex=as.factor(Sr_Chimps_4$Sex)
Sr_Chimps_4$Time_Slice2=as.factor(Sr_Chimps_4$Time_Slice2)
Sr_Chimps_4$Age_year=as.numeric(Sr_Chimps_4$Age_year)
Sr_Chimps_4$Teeth_Sampled=as.factor(Sr_Chimps_4$Teeth_Sampled)

Sr_Chimps_7$Origin=as.factor(Sr_Chimps_7$Origin)
Sr_Chimps_7$Age_year=as.numeric(Sr_Chimps_7$Age_year)
Sr_Chimps_7$Teeth_Sampled=as.factor(Sr_Chimps_7$Teeth_Sampled)

str(Sr_residents)
Sr_residents

Sr_residents$Origin=as.factor(Sr_residents$Origin)
Sr_residents$Sex=as.factor(Sr_residents$Sex)
Sr_residents$Time_Slice2=as.factor(Sr_residents$Time_Slice2)
Sr_residents$Age_year=as.numeric(Sr_residents$Age_year)
Sr_residents$Teeth_Sampled=as.factor(Sr_residents$Teeth_Sampled)

ggplot(Sr_Chimps_3, aes(Teeth_Sampled, Sr87_Sr86)) +
  geom_beeswarm() +
  theme_classic()

ggplot(Sr_Chimps_3, aes(Origin_Home_Range , Sex)) +
  geom_point() +
  theme_classic()

ggplot(Sr_Chimps_3, aes(Time_Slice2 , Sex )) +
  geom_point() +
  theme_classic()

ggplot(Sr_Chimps_3, aes(Origin_Home_Range , Time_Slice2 )) +
  geom_point() +
  theme_classic()

##lots of redundancy between Sex, Origin_Home_Range, and Time_Slice2. Models won't work with all three together. Consider leaving out Sex and/or one of the other two

exampledist1 <- fitdist(as.numeric(na.omit(Sr_Chimps_3$Sr87_Sr86)), "norm", method = "mle")
exampledist2 <- fitdist(as.numeric(na.omit(Sr_Chimps_3$Sr87_Sr86)), "gamma", method = "mme")
#does not work# exampledist3 <- fitdist(as.numeric(na.omit(Sr_Chimps_3$Sr87_Sr86)), "t", start=list(mean=mean(as.numeric(na.omit(Sr_Chimps_3$Sr87_Sr86))),sd=sd(as.numeric(na.omit(Sr_Chimps_3$Sr87_Sr86))), df=3), lower=c(0, 0,0))
exampledist4 <- fitdist(as.numeric(na.omit(Sr_Chimps_4$Sr87_Sr86)), "beta", method = "mme")
exampledist5 <- fitdist(as.numeric(na.omit(Sr_Chimps_5$Sr87_Sr86)), "beta", method = "mme")


qqcomp(list(exampledist1)) #Gaussian fits good
qqcomp(list(exampledist2))
qqcomp(list(exampledist3))
qqcomp(list(exampledist4)) #beta fits good
qqcomp(list(exampledist5))


# run a model on all groups to see if we can speak of differences across chimpanzee group territories if we exclude all adult females
# (while controlling for tooth type!), make subset without adult females (n=22)

#females <- filter(Sr_Chimps_3, age_sex == "female")

#red.Tai_Sr <-filter(Tai_Sr,!age_sex%in%c("female"))
red.Tai_Sr=Sr_Chimps_3[-which(as.character(Sr_Chimps_3$Sex)=="f"),]
str(red.Tai_Sr)

#groups with new delineations of pre- and post- 1980
red.Tai_Sr2=Sr_Chimps_7[-which(as.character(Sr_Chimps_7$age_sex) %in% c("female", "subadult female")),]

#original model - don't use
#full_Tai_groups = lmer (Sr.sc ~   Territory  + (1|tooth) , data= red.Tai_Sr, REML=F)
#null_Tai_groups = lmer (Sr.sc ~    (1|tooth), data= red.Tai_Sr, REML=F) 
full_Tai_groups = glmmTMB (Sr87_Sr86  ~   Time_Slice2  + (1|Teeth_Sampled) , data= red.Tai_Sr, family=beta_family(link="logit"))
null_Tai_groups = glmmTMB (Sr87_Sr86  ~ (1|Teeth_Sampled), data= red.Tai_Sr, family=beta_family(link="logit")) 

#groups with pre and post- 1980 
full_Tai_groups2 = glmmTMB (Sr87_Sr86  ~   Origin  + (1|Teeth_Sampled) , data= red.Tai_Sr2, family=beta_family(link="logit"))
null_Tai_groups2 = glmmTMB (Sr87_Sr86  ~ (1|Teeth_Sampled), data= red.Tai_Sr2, family=beta_family(link="logit")) 

################including resident subadult females in chimp baseline model "Sr_residents" n =23 ###################### 

full_Tai_groups3 = glmmTMB (Sr87_Sr86  ~   Origin  + (1|Teeth_Sampled) , data= Sr_residents, family=beta_family(link="logit"))
null_Tai_groups3 = glmmTMB (Sr87_Sr86  ~ (1|Teeth_Sampled), data= Sr_residents, family=beta_family(link="logit")) 

MuMIn::AICc(full_Tai_groups, null_Tai_groups) #full tai groups wins
MuMIn::AICc(full_Tai_groups2, null_Tai_groups2) #full tai groups2 wins, lower AICc than null
MuMIn::AICc(full_Tai_groups3, null_Tai_groups3) #null does better, lower AICc

summary(full_Tai_groups)
#summary 
Family: beta  ( logit )
Formula:          Sr87_Sr86 ~ Time_Slice2 + (1 | Teeth_Sampled)
Data: red.Tai_Sr

AIC      BIC   logLik deviance df.resid 
-177.7   -171.8     95.8   -191.7       10 

Random effects:
  
  Conditional model:
  Groups        Name        Variance Std.Dev. 
Teeth_Sampled (Intercept) 3.49e-13 5.908e-07
Number of obs: 17, groups:  Teeth_Sampled, 4

Dispersion parameter for beta family (): 2.71e+05 

Conditional model:
  Estimate Std. Error z value Pr(>|z|)    
(Intercept)        0.928347   0.003014  308.04  < 2e-16 ***
  Time_Slice2Middle  0.025412   0.005239    4.85 1.23e-06 ***
  Time_Slice2North 1 0.016751   0.003420    4.90 9.68e-07 ***
  Time_Slice2North 2 0.013560   0.003895    3.48 0.000499 ***
  Time_Slice2South   0.019582   0.003696    5.30 1.17e-07 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(full_Tai_groups2)

Family: beta  ( logit )
Formula:          Sr87_Sr86 ~ Origin + (1 | Teeth_Sampled)
Data: red.Tai_Sr2

AIC      BIC   logLik deviance df.resid 
-177.3   -171.5     95.7   -191.3       10 

Random effects:
  
  Conditional model:
  Groups        Name        Variance  Std.Dev.
Teeth_Sampled (Intercept) 1.571e-06 0.001253
Number of obs: 17, groups:  Teeth_Sampled, 4

Dispersion parameter for beta family (): 2.86e+05 

Conditional model:
  Estimate Std. Error z value Pr(>|z|)    
(Intercept)                 0.928214   0.003063  303.08  < 2e-16 ***
  OriginMiddle reference      0.026261   0.005547    4.73 2.20e-06 ***
  OriginNorth natal post-1980 0.015491   0.003418    4.53 5.83e-06 ***
  OriginNorth natal pre-1980  0.018063   0.004155    4.35 1.38e-05 ***
  OriginSouth reference       0.018948   0.004000    4.74 2.17e-06 ***
  ---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(full_Tai_groups3)

simOut3 <-  simulateResiduals(fittedModel = full_Tai_groups3 , n = 500)
plot(simOut3) #no problems detected
testDispersion(simOut2) 
plotResiduals(simOut2, full_Tai_groups2$Origin, xlab = "Origin", main=NULL) #looks good

#making dataframe to write into table
res_male_baseline = data.frame(summary(full_Tai_groups))
res_male_baseline <- data.frame(unclass(summary(full_Tai_groups)))

write.csv(summary(full_Tai_gorups), file=“file.csv”)

write.csv(summary, file=“file.csv”)

library(emmeans) #treated with bonferroni b/c we have known groups we want to see true diff between
pairwisecomparisons1 <- emmeans(full_Tai_groups, pairwise ~ Time_Slice2, type = "response", adjust = "bonferroni")
summary(pairwisecomparisons1)

#with new groups
pairwisecomparisons1 <- emmeans(full_Tai_groups2, pairwise ~ Origin, type = "response", adjust = "bonferroni")
summary(pairwisecomparisons1)


#comparing pairwise p values
pwpp(pairwisecomparisons1, adjust = "bonferroni") #p-values adjusted by bonferroni 


library(effects)

groups=allEffects(full_Tai_groups2)
plot(groups)

#plotting model outputs

res_males <- as.data.frame(groups$Origin)

res_males$Origin <- factor(res_males$Origin, levels = c("North natal pre-1980",
                                                        "North natal post-1980",
                                                        "South reference",
                                                        "East reference",
                                                        "Middle reference"))


results_full_Tai_groups=data.frame(confint(full_Tai_groups, full = TRUE))
colnames(results_full_Tai_groups)=c("lower","upper", "Estimate")
results_full_Tai_groups$coef=rownames(results_full_Tai_groups)
results_full_Tai_groups_plot=ggplot(results_full_Tai_groups,aes(x=Estimate,y=coef))+geom_point()+geom_linerange(aes(xmin=upper,xmax=lower))+
  geom_vline(xintercept = 0)+theme_classic()
results_full_Tai_groups_plot

##Plot without unnecessary ones
results_full_Tai_groups2=results_full_Tai_groups[-which(results_full_Tai_groups$coef=="sigma"),]
results_full_Tai_groups2=results_full_Tai_groups2[-which(results_full_Tai_groups2$coef=="cond.Std.Dev.(Intercept)|Teeth_Sampled"),]

results_full_Tai_groups_plot2=ggplot(results_full_Tai_groups2,aes(x=Estimate,y=coef))+geom_point()+geom_linerange(aes(xmin=upper,xmax=lower))+
  geom_vline(xintercept = 0)+theme_classic()
results_full_Tai_groups_plot2 ##Territory significant for resident males


# anova(null_Tai_groups, full_Tai_groups, test="Chisq")
# # suggests significant differences in Sr between territories
# #                npar    AIC    BIC  logLik deviance Chisq Df Pr(>Chisq)  
# #null_Tai_groups    3 51.908 55.181 -22.954   45.908                      
# #full_Tai_groups    6 48.941 55.487 -18.470   36.941 8.967  3    0.02973 *

simOut2 <-  simulateResiduals(fittedModel = full_Tai_groups , n = 500)
plot(simOut2) #quantile deviation detected
testDispersion(simOut2) 
plotResiduals(simOut2, full_Tai_groups$Time_Slice2, xlab = "Time Slice", main=NULL)
#qqplots look good, no under/over distribution -- model fit is OK, no cause for concern from stack overflow

#qqplot, histogramm, residuals  look okay, a little skewed but not too bad
## check vifs #######################
vif_full_Tai_groups = glmmTMB (Sr87_Sr86 ~    Teeth_Sampled , data= red.Tai_Sr, family=beta_family(link="logit"))

#results (results can be plotted in boxplots as Renee already did, but now removing adult females here)

pairwisecomparisons_vif_full_Tai_groups <- emmeans(vif_full_Tai_groups , specs = pairwise ~ Teeth_Sampled, type = "response")
pairwisecomparisons_vif_full_Tai_groups # not significant for tooth type
plot(pairwisecomparisons_vif_full_Tai_groups) #no difference between tooth type

teeth=allEffects(vif_full_Tai_groups)
plot(teeth)

results_vif_full_Tai_groups=data.frame(confint(vif_full_Tai_groups, full = TRUE))
colnames(results_vif_full_Tai_groups)=c("lower","upper", "Estimate")
results_vif_full_Tai_groups$coef=rownames(results_vif_full_Tai_groups)
results_results_vif_full_Tai_groups_plot=ggplot(results_vif_full_Tai_groups,aes(x=Estimate,y=coef))+geom_point()+geom_linerange(aes(xmin=upper,xmax=lower))+
  geom_vline(xintercept = 0)+theme_classic()
results_results_vif_full_Tai_groups_plot

##Plot without unnecessary ones
results_vif_full_Tai_groups2=results_vif_full_Tai_groups[-which(results_vif_full_Tai_groups$coef=="sigma"),]

results_results_vif_full_Tai_groups_plot2=ggplot(results_vif_full_Tai_groups2,aes(x=Estimate,y=coef))+geom_point()+geom_linerange(aes(xmin=upper,xmax=lower))+
  geom_vline(xintercept = 0)+theme_classic()
results_results_vif_full_Tai_groups_plot2 ##Nothing significant 


###########################################################################################################
### Run a model ONLY on NORTH GROUP now, to see how sex interacting with age, as well as territory shifts 
#("time slice2" influence Sr ratios across individuals (all while controlling for tooth type!) , n = 27

N.Tai_Sr <-Sr_Chimps_3[which(as.character(Sr_Chimps_3$Territory)=="North"),]

N.Tai_Sr4 <-Sr_Chimps_4[which(as.character(Sr_Chimps_4$Territory)=="North"),]
N.Tai_Sr5 <-Sr_Chimps_5[which(as.character(Sr_Chimps_5$Territory)=="North"),]
N.Tai_Sr6 <-Sr_Chimps_6[which(as.character(Sr_Chimps_6$Territory)=="North"),]

#exclude resident females 
N.Tai_Sr_resmales=N.Tai_Sr[-which(as.character(N.Tai_Sr$Origin_group)=="resident female"),]

#model without res females and sex as fixed effect

#hypothesis - **** model 3 in Sr manuscript ****
red_N.Tai = glmmTMB (Sr87_Sr86 ~   Origin_Home_Range + Age_year + Time_Slice2 + (1|Teeth_Sampled) , data= N.Tai_Sr, family=beta_family(link="logit"))
red_N.Taidrop = glmmTMB (Sr87_Sr86 ~   Origin_Home_Range + (1|Teeth_Sampled) , data= N.Tai_Sr, family=beta_family(link="logit"))
null_N.Tai = glmmTMB (Sr87_Sr86 ~ (1|Teeth_Sampled), data= N.Tai_Sr, family=beta_family(link="logit"))  

#updated models with new group temporal cuts  **** model 3 in Sr manuscript ****
red_N.Taidrop4 = glmmTMB (Sr87_Sr86 ~   Origin + (1|Teeth_Sampled) , data= N.Tai_Sr4, family=beta_family(link="logit"))
red_N.Taidrop5 = glmmTMB (Sr87_Sr86 ~   Origin + (1|Teeth_Sampled) , data= N.Tai_Sr5, family=beta_family(link="logit"))
red_N.Taidrop6 = glmmTMB (Sr87_Sr86 ~   Origin + (1|Teeth_Sampled) , data= N.Tai_Sr6, family=beta_family(link="logit"))
null_N.Tai4 = glmmTMB (Sr87_Sr86 ~ (1|Teeth_Sampled), data= N.Tai_Sr4, family=beta_family(link="logit"))  
null_N.Tai5 = glmmTMB (Sr87_Sr86 ~ (1|Teeth_Sampled), data= N.Tai_Sr5, family=beta_family(link="logit"))  
null_N.Tai6 = glmmTMB (Sr87_Sr86 ~ (1|Teeth_Sampled), data= N.Tai_Sr6, family=beta_family(link="logit"))  


#enacted temporal fusion for North 1 and 2
resf_N.Tai = glmmTMB (Sr87_Sr86 ~   Origin_Home_Range + Age_year + Time_Slice2 + (1|Teeth_Sampled) , data= N.Tai_Sr_resmales, family=beta_family(link="logit"))
resf_N.Taidrop = glmmTMB (Sr87_Sr86 ~  Origin_Home_Range + Time_Slice2 + (1|Teeth_Sampled) , data= N.Tai_Sr_resmales, family=beta_family(link="logit"))
nullresf_N.Tai = glmmTMB (Sr87_Sr86 ~ (1|Teeth_Sampled), data= N.Tai_Sr_resmales, family=beta_family(link="logit"))  

MuMIn::AICc(red_N.Taidrop4, null_N.Tai4,red_N.Taidrop5, null_N.Tai5,red_N.Taidrop6, null_N.Tai6) #red_N Tai drop4 wins, 
#lowest AICc,GOM, BIJ post-1980 
summary(red_N.Taidrop4)

Family: beta  ( logit )
Formula:          Sr87_Sr86 ~ Origin + (1 | Teeth_Sampled)
Data: N.Tai_Sr4

AIC      BIC   logLik deviance df.resid 
-244.2   -236.4    128.1   -256.2       21 

Random effects:
  
  Conditional model:
  Groups        Name        Variance  Std.Dev.
Teeth_Sampled (Intercept) 1.047e-05 0.003236
Number of obs: 27, groups:  Teeth_Sampled, 4

Dispersion parameter for beta family (): 4.89e+04 

Conditional model:
  Estimate Std. Error z value Pr(>|z|)    
(Intercept)                     0.941296   0.003215  292.78  < 2e-16 ***
  Originnatal pre-1980            0.005584   0.007744    0.72    0.471    
Originunknown female post-1980  0.021056   0.005067    4.16 3.25e-05 ***
  Originunknown female pre-1980  -0.007677   0.005869   -1.31    0.191    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

simOut5 <-  simulateResiduals(fittedModel = red_N.Taidrop4, n = 500)
plot(simOut5) #some deviation, overall not significant 
testDispersion(simOut5)
plotResiduals(simOut5, red_N.Taidrop4$Origin, xlab = "Origin", main=NULL) #Venus as outlier

pcomp_red_N.Tai4_Origin <- emmeans(red_N.Taidrop4, specs = pairwise ~ Origin, type = "response", adjust = "bonferroni")
pcomp_red_N.Tai4_Origin
plot(pcomp_red_N.Tai4_Origin)

$emmeans
Origin                   response        SE  df asymp.LCL asymp.UCL
natal post-1980            0.7194 0.0006490 Inf    0.7181    0.7206
natal pre-1980             0.7205 0.0015008 Inf    0.7175    0.7234
unknown female post-1980   0.7236 0.0009248 Inf    0.7218    0.7254
unknown female pre-1980    0.7178 0.0011097 Inf    0.7156    0.7200

Confidence level used: 0.95 
Intervals are back-transformed from the logit scale 

$contrasts
contrast                                               odds.ratio      SE  df null
(natal post-1980) / (natal pre-1980)                        0.994 0.00770 Inf    1
(natal post-1980) / (unknown female post-1980)              0.979 0.00496 Inf    1
(natal post-1980) / (unknown female pre-1980)               1.008 0.00591 Inf    1
(natal pre-1980) / (unknown female post-1980)               0.985 0.00804 Inf    1
(natal pre-1980) / (unknown female pre-1980)                1.013 0.00883 Inf    1
(unknown female post-1980) / (unknown female pre-1980)      1.029 0.00660 Inf    1
z.ratio p.value
-0.721  1.0000
-4.155  0.0002 #significant, natal and unknown post-1980
1.308  1.0000 
-1.895  0.3482
1.521  0.7691
4.482  <.0001 #significant, unknown female post-1980 and pre-1980

P value adjustment: bonferroni method for 6 tests 
Tests are performed on the log odds ratio scale 

MuMIn::AICc(red_N.Taidrop, null_N.Tai) #red_N Tai drop wins
summary(red_N.Taidrop)
#comparing to North baseline w/ res females wins
df      AICc
resf_N.Taidrop  6 -179.7152
nullresf_N.Tai  3 -178.7197
red_N.Taidrop   6 -235.9190 #wins 
null_N.Tai      3 -178.7197

simOut3 <-  simulateResiduals(fittedModel = red_N.Taidrop, n = 500)
plot(simOut3) #model fit looks ok

#summary 
Family: beta  ( logit )
Formula:          Sr87_Sr86 ~ Origin_Home_Range + (1 | Teeth_Sampled)
Data: N.Tai_Sr

AIC      BIC   logLik deviance df.resid 
-241.9   -235.4    125.9   -251.9       22 

Random effects:
  
  Conditional model:
  Groups        Name        Variance  Std.Dev. 
Teeth_Sampled (Intercept) 3.137e-13 5.601e-07
Number of obs: 27, groups:  Teeth_Sampled, 4

Dispersion parameter for beta family (): 3.88e+04 

Conditional model:
  Estimate Std. Error z value Pr(>|z|)    
(Intercept)                              0.931281   0.006512  143.01  < 2e-16 ***
  Origin_Home_Rangebirth community unknown 0.026791   0.007648    3.50  0.00046 ***
  Origin_Home_Rangenatal group             0.010594   0.007099    1.49  0.13561    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
> 
  
  pcomp_red_N.Tai_Origin <- emmeans(red_N.Taidrop, specs = pairwise ~ Origin_Home_Range, type = "response", adjust = "bonferroni")
pcomp_red_N.Tai_Origin
plot(pcomp_red_N.Tai_Origin)

pwpp(pcomp_red_N.Tai_Origin, adjust = "bonferroni")

pcomp_red_N.Tai_Time_Slice2 <- emmeans(red_N.Taidrop, specs = pairwise ~ Time_Slice2, type = "response", adjust = "bonferroni")
pcomp_red_N.Tai_Time_Slice2 #not significant between time slice
plot(pcomp_red_N.Tai_Time_Slice2)

North=allEffects(red_N.Taidrop)
plot(North)

effect_N <- emmeans(red_N.Tai, pairwise ~ Origin_Home_Range, type = "response", adjust = "bonferroni")
effect_N
plot(effect_N)

plot(effect_N, comparisons = TRUE)

pwpp(pcomp_red_N.Tai_Time_Slice2, adjust = "bonferroni")

#model with res females
str(N.Tai_Sr)
N.Tai_Sr$Territory


#full_N.Tai = lmer (Sr.sc ~   Origin + Sex * age + Time_Slice2   + (1|tooth) , data= N.Tai_Sr, REML=F)
red_N.Tai = glmmTMB (Sr87_Sr86 ~   Origin_Home_Range + Sex + Age_year + Time_Slice2   + (1|Teeth_Sampled) , data= N.Tai_Sr, family=beta_family(link="logit"))
#anova(red_N.Tai, full_N.Tai, test="Chisq") # interaction of age and sex insignificant, can be dropped
red_N.Taidrop = glmmTMB (Sr87_Sr86 ~   Origin_Home_Range + Sex + Time_Slice2 + (1|Teeth_Sampled) , data= N.Tai_Sr, family=beta_family(link="logit"))
null_N.Tai = glmmTMB (Sr87_Sr86 ~ (1|Teeth_Sampled), data= N.Tai_Sr, family=beta_family(link="logit"))  

MuMIn::AICc(red_N.Taidrop, null_N.Tai) #red_N Tai drop wins
summary(red_N.Taidrop)

#summary of model
Family: beta  ( logit )
Formula:          Sr87_Sr86 ~ Origin_Home_Range + Time_Slice2 + (1 | Teeth_Sampled)
Data: N.Tai_Sr

AIC      BIC   logLik deviance df.resid 
-240.1   -232.3    126.1   -252.1       21 

Random effects:
  
  Conditional model:
  Groups        Name        Variance  Std.Dev. 
Teeth_Sampled (Intercept) 2.467e-13 4.967e-07
Number of obs: 27, groups:  Teeth_Sampled, 4

Dispersion parameter for beta family (): 3.91e+04 

Conditional model:
  Estimate Std. Error z value Pr(>|z|)    
(Intercept)                          0.931281   0.006485  143.61  < 2e-16 ***
  Origin_Home_Rangebirthplace unknown  0.026791   0.007616    3.52 0.000436 ***
  Origin_Home_Rangenatal group         0.011365   0.007254    1.57 0.117159    
Time_Slice2North 2                  -0.003087   0.006498   -0.48 0.634689  


pcomp_red_N.Tai_Origin <- emmeans(red_N.Taidrop, specs = pairwise ~ Origin_Home_Range, type = "response", adjust = "bonferroni")
pcomp_red_N.Tai_Origin
plot(pcomp_red_N.Tai_Origin)

pwpp(pcomp_red_N.Tai_Origin, adjust = "bonferroni")


#results for North group only

$emmeans
Origin_Home_Range                   response        SE  df asymp.LCL asymp.UCL
birth community & territory unknown   0.7173 0.0013204 Inf    0.7147    0.7199
birth community unknown               0.7227 0.0008038 Inf    0.7212    0.7243
natal group                           0.7195 0.0005704 Inf    0.7184    0.7206

Confidence level used: 0.95 
Intervals are back-transformed from the logit scale 

$contrasts
contrast                                                      odds.ratio      SE  df null z.ratio p.value
birth community & territory unknown / birth community unknown      0.974 0.00745 Inf    1  -3.503  0.0014 #significant
birth community & territory unknown / natal group                  0.989 0.00702 Inf    1  -1.492  0.4068
birth community unknown / natal group                              1.016 0.00499 Inf    1   3.301  0.0029 #significant

P value adjustment: bonferroni method for 3 tests 
Tests are performed on the log odds ratio scale 


Results are averaged over the levels of: Origin_Home_Range 
Tests are performed on the log odds ratio scale 

North=allEffects(red_N.Taidrop4)
plot(North)

#plotting model

N_ts <- as.data.frame(North$Time_Slice2)

N_origin <- as.data.frame(North$Origin)

N_ts$Time_Slice2 <- factor(N_ts$Time_Slice2, levels = c("North 1", "North 2"))

N_origin$Origin <- factor(N_origin$Origin, levels = c("natal pre-1980","unknown female pre-1980",
                                                      "natal post-1980",
                                                      "unknown female post-1980"))




ggsave("/Users/Renee/Desktop/Sr Paper/Figures/Fig-S3_origin.pdf",width=14, height = 9, units="in", encoding = "MacRoman")


effect_N <- emmeans(red_N.Taidrop, pairwise ~ Origin_Home_Range, type = "response", adjust = "bonferroni")
effect_N
plot(effect_N)

plot(effect_N, comparisons = TRUE)

pwpp(pcomp_red_N.Tai_Time_Slice2, adjust = "bonferroni")


results_null_N.Tai=data.frame(confint(red_N.Tai, full = TRUE))
colnames(results_null_N.Tai)=c("lower","upper", "Estimate")
results_null_N.Tai$coef=rownames(results_null_N.Tai)
results_null_N.Tai_plot=ggplot(results_null_N.Tai,aes(x=Estimate,y=coef))+geom_point()+geom_linerange(aes(xmin=upper,xmax=lower))+
  geom_vline(xintercept = 0)+theme_classic()
results_null_N.Tai_plot

##Plot without unnecessary ones
results_null_N.Tai2=results_null_N.Tai[-which(results_null_N.Tai$coef=="sigma"),]
results_null_N.Tai2=results_null_N.Tai2[-which(results_null_N.Tai2$coef=="cond.Std.Dev.(Intercept)|Teeth_Sampled"),]

results_null_N.Tai_plot2=ggplot(results_null_N.Tai2,aes(x=Estimate,y=coef))+geom_point()+geom_linerange(aes(xmin=upper,xmax=lower))+
  geom_vline(xintercept = 0)+theme_classic()
results_null_N.Tai_plot2 ##Nothing significant 

#are there differences between Sr baselines for territories? (plants)

str(Tai_SrEco)
hist(Tai_SrEco$Sr87_Sr86, breaks  = 50) # not normal distribution
plot(density(Tai_SrEco$Sr87_Sr86)) # normal density


Tai_SrEco$project_site=as.factor(Tai_SrEco$project_site)
Tai_SrEco$Time_Slice2=as.factor(Tai_SrEco$Time_Slice2)
Tai_SrEco$Group=as.factor(Tai_SrEco$Group)
Tai_SrEco$PEMA_No=as.factor(Tai_SrEco$No)
Tai_SrEco$sample_type=as.factor(Tai_SrEco$sample_type)
Tai_SrEco$habitat=as.factor(Tai_SrEco$habitat)
Tai_SrEco$stratigraphy=as.factor(Tai_SrEco$stratigraphy)
Tai_SrEco$age=as.factor(Tai_SrEco$age)

ggplot(Tai_SrEco, aes(Time_Slice2, Sr87_Sr86)) +
  geom_beeswarm() +
  theme_classic()

exampledist5 <- fitdist(as.numeric(na.omit(Tai_SrEco$Sr87_Sr86)), "norm", method = "mle")
exampledist6 <- fitdist(as.numeric(na.omit(Tai_SrEco$Sr87_Sr86)), "gamma", method = "mme")
exampledist7 <- fitdist(as.numeric(na.omit(Tai_SrEco$Sr87_Sr86)), "t", start=list(mean=mean(as.numeric(na.omit(Tai_SrEco$Sr87_Sr86))),sd=sd(as.numeric(na.omit(Tai_SrEco$Sr87_Sr86))), df=3), lower=c(0, 0,0))
exampledist8 <- fitdist(as.numeric(na.omit(Tai_SrEco$Sr87_Sr86)), "beta", method = "mme")

qqcomp(list(exampledist5)) #Gaussian fits good
qqcomp(list(exampledist6))
qqcomp(list(exampledist7))
qqcomp(list(exampledist8)) #beta fits good

model_s <- glmmTMB(Sr87_Sr86 ~ Time_Slice2 + sample_type, data = Tai_SrEco, family=beta_family(link="logit"))
model_r <- glmmTMB(Sr87_Sr86 ~ Time_Slice2 + (1|sample_type), data = Tai_SrEco, family=beta_family(link="logit")) #best fit with AICc
model_t <- glmmTMB(Sr87_Sr86 ~ (1|sample_type), data = Tai_SrEco, family=beta_family(link="logit"))

MuMIn::AICc(model_s, model_r, model_t) #model_r wins, df 10

model_r_output <- summary(model_r)

simOut4 <-  simulateResiduals(fittedModel = model_r, n = 500)
plot(simOut4) #model fit looks ok

#summary model_r
Family: beta  ( logit )
Formula:          Sr87_Sr86 ~ Time_Slice2 + (1 | sample_type)
Data: Tai_SrEco

AIC      BIC   logLik deviance df.resid 
-282.0   -266.5    151.0   -302.0       25 

Random effects:
  
  Conditional model:
  Groups      Name        Variance  Std.Dev. 
sample_type (Intercept) 3.688e-13 6.073e-07
Number of obs: 35, groups:  sample_type, 3

Dispersion parameter for beta family (): 1.92e+04 

Conditional model:
  Estimate Std. Error z value Pr(>|z|)    
(Intercept)          0.925748   0.007998  115.75  < 2e-16 ***
  Time_Slice2Middle    0.006074   0.010736    0.57  0.57155    
Time_Slice2North 1   0.005968   0.013864    0.43  0.66684    
Time_Slice2North 2   0.004467   0.013861    0.32  0.74723    
Time_Slice2North 3   0.066341   0.012321    5.38 7.27e-08 ***
  Time_Slice2Northeast 0.027747   0.010759    2.58  0.00991 ** 
  Time_Slice2Northwest 0.040349   0.009638    4.19 2.84e-05 ***
  Time_Slice2South     0.007634   0.010738    0.71  0.47714    
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

simOutr <-  simulateResiduals(fittedModel = model_r , n = 500)
plot(simOutr) #model fit looks ok

##Plot CI of coefficients, significant if they exclude 0
results_r=data.frame(confint(model_r, full = TRUE))
colnames(results_r)=c("lower","upper", "Estimate")
results_r$coef=rownames(results_r)
results_r_plot=ggplot(results_r,aes(x=Estimate,y=coef))+geom_point()+geom_linerange(aes(xmin=upper,xmax=lower))+
  geom_vline(xintercept = 0)+theme_classic()
results_r_plot

##Plot without unnecessary ones
results_r2=results_r[-which(results_r$coef=="sigma"),]
results_r2=results_r2[-which(results_r2$coef=="cond.Std.Dev.(Intercept)|sample_type"),]

results_r_plot2=ggplot(results_r2,aes(x=Estimate,y=coef))+geom_point()+geom_linerange(aes(xmin=upper,xmax=lower))+
  geom_vline(xintercept = 0)+theme_classic()
results_r_plot2 #significant  for Northwest, North 3, Northeast, South

pairwisecomparisons_model_r <- emmeans(model_r, specs = pairwise ~ Time_Slice2, type = "response", adjust = "bonferroni")
pairwisecomparisons_model_r
plot(pairwisecomparisons_model_r) #significant
plot(pairwisecomparisons_model_r, comparisons = TRUE)

$emmeans
Time_Slice2 response       SE  df asymp.LCL asymp.UCL
East          0.7162 0.001626 Inf    0.7130    0.7194
Middle        0.7174 0.001452 Inf    0.7146    0.7203
North 1       0.7174 0.002296 Inf    0.7129    0.7219
North 2       0.7171 0.002297 Inf    0.7126    0.7216
North 3       0.7295 0.001849 Inf    0.7259    0.7331
Northeast     0.7218 0.001445 Inf    0.7190    0.7246
Northwest     0.7243 0.001074 Inf    0.7222    0.7264
South         0.7178 0.001452 Inf    0.7149    0.7206

Confidence level used: 0.95 
Intervals are back-transformed from the logit scale 

$contrasts
contrast              odds.ratio      SE  df null z.ratio p.value
East / Middle              0.994 0.01067 Inf    1  -0.566  1.0000
East / North 1             0.994 0.01378 Inf    1  -0.430  1.0000
East / North 2             0.996 0.01380 Inf    1  -0.322  1.0000
East / North 3             0.936 0.01153 Inf    1  -5.384  <.0001
East / Northeast           0.973 0.01046 Inf    1  -2.579  0.2775
East / Northwest           0.960 0.00926 Inf    1  -4.186  0.0008
East / South               0.992 0.01066 Inf    1  -0.711  1.0000
Middle / North 1           1.000 0.01340 Inf    1   0.008  1.0000
Middle / North 2           1.002 0.01342 Inf    1   0.120  1.0000
Middle / North 3           0.942 0.01111 Inf    1  -5.109  <.0001
Middle / Northeast         0.979 0.00994 Inf    1  -2.134  0.9185
Middle / Northwest         0.966 0.00866 Inf    1  -3.826  0.0036
Middle / South             0.998 0.01012 Inf    1  -0.154  1.0000
North 1 / North 2          1.002 0.01604 Inf    1   0.094  1.0000
North 1 / North 3          0.941 0.01384 Inf    1  -4.107  0.0011
North 1 / Northeast        0.978 0.01313 Inf    1  -1.623  1.0000
North 1 / Northwest        0.966 0.01211 Inf    1  -2.742  0.1709
North 1 / South            0.998 0.01338 Inf    1  -0.124  1.0000
North 2 / North 3          0.940 0.01382 Inf    1  -4.210  0.0007
North 2 / Northeast        0.977 0.01311 Inf    1  -1.735  1.0000
North 2 / Northwest        0.965 0.01209 Inf    1  -2.863  0.1176
North 2 / South            0.997 0.01336 Inf    1  -0.236  1.0000
North 3 / Northeast        1.039 0.01228 Inf    1   3.266  0.0305
North 3 / Northwest        1.026 0.01109 Inf    1   2.405  0.4526
North 3 / South            1.060 0.01251 Inf    1   4.976  <.0001
Northeast / Northwest      0.987 0.00887 Inf    1  -1.403  1.0000
Northeast / South          1.020 0.01036 Inf    1   1.981  1.0000
Northwest / South          1.033 0.00926 Inf    1   3.651  0.0073

P value adjustment: bonferroni method for 28 tests 
Tests are performed on the log odds ratio scale 


SrEco=allEffects(model_r)
plot(SrEco)

effect_SrEco <- emmeans(model_r, specs = pairwise ~ Time_Slice2, type = "response", adjust = "bonferroni")
effect_SrEco
plot(effect_SrEco)

#plotting model outputs
territories <- as.data.frame(SrEco$Time_Slice2)

territories$Time_Slice2 <- factor(territories$Time_Slice2, levels = c("North 1", "North 2", 
                                                                      "Northwest", "Northeast", "South", "East", "Middle"))


#regression with res males and eco baseline

write_csv(red.Tai_Sr, "/Users/Renee/Desktop/Sr Paper/Res_Males.csv")

baseline <- read_excel("/Users/Renee/Desktop/Sr Paper/baseline_reg.xlsx")

baseline_avg <- read_excel("/Users/Renee/Desktop/Sr Paper/all_baseline-avg.xlsx")

ggplot(baseline_avg, aes(res, eco)) + 
  geom_smooth(method="lm") +
  geom_point() + theme_classic()

lm()

my.formula <- y ~ x

ylab(expression(paste(""^{87},"Sr/"^86,"Sr")))

#ordering dataframe

baseline_up <- read_excel("/Users/Renee/Desktop/Sr Paper/all_baseline.xlsx")
baseline_up$territory <- factor(baseline_up$territory, levels = c("North", "South", "East", "Middle"))

baselineavg_lm <- lm(res ~ eco, data = baseline_avg)

print(summary(baselineavg_lm))

ggplot(baseline)_up, aes(res, eco)) + 
  geom_smooth(method="lm") +
  geom_point() + theme_classic()

#ordering dataframe
baseline_avg$territory <- factor(baseline_avg$territory, levels = c("North", "South", "East", "Middle"))

#cant figure out why the order of the names is not matching the colors
# baseline %>% #reorder names#
#   #mutate(territory = fct_relevel(territory, 
#                             "north 1", "north 2", 
#                             "south", "east", 
#                             "middle")) %>%
#checking that resident males and females accurately represent baseline for North

n_baseline <- read_excel("/Users/Renee/Desktop/Sr Paper/baseline_reg.xlsx")

hist(n_baseline$res)
hist(n_baseline$eco)

plot(eco ~ res, data = n_baseline)


ggplot(n_baseline, aes(eco, res)) + 
  geom_smooth(method="lm") +
  geom_point() + theme_classic()

baseline_lm <- lm(res ~ eco, data = n_baseline)

print(summary(baseline_lm))

Call:
  lm(formula = res ~ eco, data = n_baseline)

Residuals:
  1          2          3          4 
-0.0005423 -0.0007143 -0.0001034  0.0013600 

Coefficients:
  Estimate Std. Error t value Pr(>|t|)  
(Intercept)  -1.3669     0.6995  -1.954   0.1899  
eco           2.9095     0.9754   2.983   0.0964 .
---
  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.001154 on 2 degrees of freedom
Multiple R-squared:  0.8165,	Adjusted R-squared:  0.7247 
F-statistic: 8.899 on 1 and 2 DF,  p-value: 0.0964

my.formula <- y ~ x

citation("stats")

#example 
library(fitdistrplus)

data = read.csv("Sr-chimps_resident baseline copy.csv")

data2 = data[which(data$Territory=="North"),]
data3 = data[-which(data$Territory=="North"),]

fit = fitdist(data2$Sr87_Sr86, "beta",method = "mme")
summary(fit)
plot(fit)

plot(density(data2$Sr87_Sr86))

newdataNorth = rbeta(10000,shape1 = 44804.36, shape2 = 17468.86)
newdataNorth2 = dbeta(newdataNorth,shape1 = 44804.36, shape2 = 17468.86)

newdataother = dbeta(data3$Sr87_Sr86,shape1 = 44804.36, shape2 = 17468.86)
newdataNorth = data.frame(cbind(newdataNorth, newdataNorth2))

library(ggplot2)
first  = ggplot(newdataNorth, aes(x = newdataNorth, y = newdataNorth2)) +
  geom_line() +
  geom_vline(data = data3, aes(xintercept = Sr87_Sr86, color = factor(Territory))) +
  ylab("density") + xlab("Sr87_Sr86")+ theme_classic()

second =  ggplot(newdataNorth, aes(x = newdataNorth)) +
  geom_density() + ggtitle("Empirical density") + theme_classic()

ggpubr::ggarrange(first,second, legend = "bottom")

#figuring out small sample sizes 
library(fitdistrplus)

#north natal, split by sex
North_natal <- read_excel("/Users/Renee_1/Desktop/Sr Paper/Sr-chimps_North natal.xlsx")

#territory and north group 
data2 = north, data3 = non-north

North = Sr_residents[which(Sr_residents$Territory=="North"),]
North_m_post = North_natal[which(North_natal$Origin== "North natal male post-1980"),]
natal_pre1980 = Sr_Chimps_7[which(Sr_Chimps_7$Origin=="North natal pre-1980"),]
North_f_post = North_natal[which(North_natal$Origin== "North natal female post-1980"),]
non_North = Sr_residents[-which(Sr_residents$Territory=="North"),]

#just North males 
North_m = fitdist(North_m_post$Sr87_Sr86, "beta",method = "mme")
summary(North_m)
plot(North_m)

North_m = rbeta(10000,shape1 = 95803.09, shape2 = 37290.00)
North2_m = dbeta(North_m,shape1 = 95803.09, shape2 = 37290.00)

newNorth_m = data.frame(cbind(North_m, North2_m))

quantile(newNorth_m$North_m, probs = .025)
#0.7174175

quantile(newNorth_m$North_m, probs = .975)
#0.722241 

library(ggplot2)
first  = ggplot(newNorth_m, aes(x = North_m, y = North2_m)) +
  geom_line() +
  geom_vline(data = non_North, aes(xintercept = Sr87_Sr86, color = factor(Territory))) +
  geom_vline(xintercept = c(.7158617,0.7229387 ), linetype  = "dotted")+
  ylab("density") + xlab("Sr87_Sr86")+ theme_classic()

#north natal males and females
fit = fitdist(North$Sr87_Sr86, "beta",method = "mme")
summary(fit)
plot(fit)

plot(density(North$Sr87_Sr86))

newdataNorth = rbeta(10000,shape1 = 44804.36, shape2 = 17468.86)
newdataNorth2 = dbeta(newdataNorth,shape1 = 44804.36, shape2 = 17468.86)

newdataNorth = data.frame(cbind(newdataNorth, newdataNorth2))

quantile(newdataNorth$newdataNorth, probs = .025)

quantile(newdataNorth$newdataNorth, probs = .975)

library(ggplot2)
first  = ggplot(newdataNorth, aes(x = newdataNorth, y = newdataNorth2)) +
  geom_line() +
  geom_vline(data = non_North, aes(xintercept = Sr87_Sr86, color = factor(Territory))) +
  geom_vline(xintercept = c(.7158617,0.7229387 ), linetype  = "dotted")+
  ylab("density") + xlab("Sr87_Sr86")+ theme_classic()

#origin model 

natal_post1980 = Sr_Chimps_7[which(Sr_Chimps_7$Origin=="North natal post-1980"),]
natal_pre1980 = Sr_Chimps_7[which(Sr_Chimps_7$Origin=="North natal pre-1980"),]
non_North = Sr_residents[-which(Sr_residents$Territory=="North"),]

shapiro.test(natal_post1980$Sr87_Sr86) #normal

plot(density(natal_post1980$Sr87_Sr86))

fit2 = fitdist(natal_post1980$Sr87_Sr86, "beta",method = "mme")
summary(fit2)
plot(fit2)

plot(density(North$Sr87_Sr86))

natal_post1980 = rbeta(10000,shape1 = 40719.09, shape2 = 15886.31)
natal_post1980_2 = dbeta(natal_post1980,shape1 = 40719.09, shape2 = 15886.31)

North_post1980 = data.frame(cbind(natal_post1980, natal_post1980_2))

quantile(North_post1980$natal_post1980, probs = .025)

quantile(North_post1980$natal_post1980, probs = .975)

library(ggplot2)
first2  = ggplot(North_post1980, aes(x = natal_post1980, y = natal_post1980_2)) +
  geom_line() +
  geom_vline(data = natal_pre1980, aes(xintercept = Sr87_Sr86, color = Origin)) +
  scale_colour_manual(name = "Origin",
                      values = c("North natal pre-1980" = "orchid")) +
  geom_vline(xintercept = c(0.715666,0.7230852), linetype  = "dotted")+
  ylab("density") + xlab("Sr87_Sr86")+ theme_minimal()

second =  ggplot(newdataNorth, aes(x = newdataNorth)) +
  geom_density() + ggtitle("Empirical density") + theme_minimal()

ggpubr::ggarrange(first2, legend = "bottom") #figure with densities

ggsave("/Users/Renee_1/Desktop/Sr revision docs/Revised figures/pre and post 1980 density.pdf",width=8, height = 8, units="in", encoding = "MacRoman")

#reference males vs. north natal with combine pre and post 1980
North_combined <- read_excel("/Users/Renee_1/Desktop/Sr Paper/North natal males.xlsx")

North = Sr_residents[which(Sr_residents$Territory=="North"),]
Ncomb_m_post = North_combined[which(North_combined$Origin== "North natal male post-1980"),]
natal_pre1980 = Sr_Chimps_7[which(Sr_Chimps_7$Origin=="North natal pre-1980"),]
North_f_post = North_natal[which(North_natal$Origin== "North natal female post-1980"),]
non_North = Sr_residents[-which(Sr_residents$Territory=="North"),]

#testing normality 
shapiro.test(Ncomb_m_post$Sr87_Sr86)

#just North males 
Ncomb = fitdist(Ncomb_m_post$Sr87_Sr86, "beta",method = "mme")
summary(Ncomb)
plot(Ncomb)

Ncomb_m = rbeta(10000,shape1 = 114593.32, shape2 = 44578.48)
Ncomb2_m = dbeta(Ncomb_m,shape1 = 114593.32, shape2 = 44578.48)

newNcomb_m = data.frame(cbind(Ncomb_m, Ncomb2_m))

quantile(newNcomb_m$Ncomb_m, probs = .025)
#0.7177086

quantile(newNcomb_m$Ncomb_m, probs = .975)
#0.7221183 

library(ggplot2)
#changing line colors
first3 = ggplot(newNcomb_m, aes(x = Ncomb_m, y = Ncomb2_m)) +
  geom_line() + geom_vline(data = non_North, aes(xintercept = Sr87_Sr86, color = Territory)) +
  scale_colour_manual(name = "Territory",
                      values = c("South" = "cyan2",
                                 "East" = "chartreuse2", "Middle" = "black")) +
  geom_vline(xintercept = c(0.7177086,0.7221183), linetype  = "dotted") +
  ylab("density") + xlab("Sr87_Sr86")+ theme_minimal()

ggpubr::ggarrange(first3, legend = "bottom")

ggsave("/Users/Renee_1/Desktop/Sr revision docs/Revised figures/resident baseline only males.pdf",width=8, height = 8, units="in", encoding = "MacRoman")


#north males pre-1980 on north males post-1980 
first2  = ggplot(newNorth_m, aes(x = North_m, y = North2_m)) +
  geom_line() +
  geom_vline(data = natal_pre1980, aes(xintercept = Sr87_Sr86, color = Origin)) +
  scale_colour_manual(name = "Origin",
                      values = c("North natal pre-1980" = "orchid")) +
  geom_vline(xintercept = c(0.7174175,0.722241), linetype  = "dotted")+
  ylab("density") + xlab("Sr87_Sr86")+ theme_minimal()

ggpubr::ggarrange(first2, legend = "bottom")

ggsave("/Users/Renee_1/Desktop/Sr revision docs/Revised figures/pre 1980 vs post 1980 only males.pdf",width=8, height = 8, units="in", encoding = "MacRoman")


ggpubr::ggarrange(first2,first3, legend = "bottom") #combined figure with north pre and post 1980 and then resident male baselines

ggsave("/Users/Renee_1/Desktop/Sr revision docs/Revised figures/only resident male baseline_revised.pdf",width=12, height = 8, units="in", encoding = "MacRoman")


