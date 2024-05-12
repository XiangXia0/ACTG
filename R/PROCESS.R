library(BART)  # The package that contains ACTG175 data set
library(dplyr)
library(ggplot2)
library(tidyverse)
library(naniar)
library(reshape2)
library(causalweight)
library(tmle)
library(dbarts)
library(chngpt)

# Read in data
data(ACTG175)
mydata<-ACTG175

# Data clean
cols_to_factor<-c("hemo","homo","drugs","race","gender","offtrt","arms")
mydata[cols_to_factor]<-lapply(mydata[cols_to_factor],factor)

# Rows and Columns
nrow(mydata)
ncol(mydata)

# Missing variables
vis_miss(mydata)

# Missing patterns
gg_miss_var(as.data.frame(mydata[c("cd496","age")]),age,show_pct=T)
gg_miss_var(as.data.frame(mydata[c("cd496","karnof")]),karnof,show_pct=T)
gg_miss_var(as.data.frame(mydata[c("cd496","arms")]),arms,show_pct=T)


hemo_pro<-mydata %>% group_by(hemo) %>% summarise(missingrate=sum(r)/length(r)) %>% as.matrix()
homo_pro<-mydata %>% group_by(homo) %>% summarise(missingrate=sum(r)/length(r)) %>% as.matrix()
drugs_pro<-mydata %>% group_by(drugs) %>% summarise(missingrate=sum(r)/length(r)) %>% as.matrix()
race_pro<-mydata %>% group_by(race) %>% summarise(missingrate=sum(r)/length(r)) %>% as.matrix()
gender_pro<-mydata %>% group_by(gender) %>% summarise(missingrate=sum(r)/length(r)) %>% as.matrix()
offtrt_pro<-mydata %>% group_by(offtrt) %>% summarise(missingrate=sum(r)/length(r)) %>% as.matrix()
combined_df<-rbind(hemo_pro,homo_pro,drugs_pro,race_pro,gender_pro,offtrt_pro) %>%
           cbind(c("Hemo","Hemo","Homo","Homo","Drugs","Drugs","Race","Race","Gender","Gender","Offtrt","Offtrt")) %>%
           data.frame()
colnames(combined_df)<-c("variable","proportion","divide")
combined_df$proportion<-combined_df$proportion %>% as.numeric()
ggplot(combined_df,aes(x=interaction(variable,divide),y=proportion,fill=divide)) +
  geom_bar(stat="identity",position="dodge") + theme_minimal() +
  labs(x="Variable value",y="Observed proportion",fill="Groups") +
  scale_x_discrete(labels=c("0","1","0","1","0","1","0","1","0","1","0","1")) +
  scale_y_continuous(breaks=seq(0,1,by=0.1)) +
  scale_fill_manual(values=c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854","#ffd92f"),
                    labels=c("Hemo","Homo","Drugs","Race","Gender","Offtrt")) +
  guides(fill=guide_legend(title="Groups"))

# Main variable
summary(mydata$cd496)
summary(mydata$cd40)

# Treatment arms
# Data clean
mydata2<-mydata %>% mutate(Gender=fct_recode(gender,'Female'='0','Male'='1')) %>%
                   mutate(Arms=fct_recode(arms,'ZDV only'='0','ZDV+ddI'='1','ZDV+ddC'='2','ddI only'='3')) %>%
                   mutate(Race=fct_recode(race,'White'='0','Non-white'='1'))

ggplot(mydata2,aes(Gender,fill=Arms))+geom_bar()

ggplot(mydata2,aes(Race,fill=Arms))+geom_bar()

ggplot(mydata2 %>% filter(r==1),aes(x=Arms,y=cd40,fill="CD4_0")) +
  geom_boxplot(position=position_nudge(x=-0.2),width=0.3) +
  geom_boxplot(aes(x=Arms,y=cd496,fill="CD4_96"),
               position=position_nudge(x=0.2),width=0.3) +
  labs(x="Arms",y="CD4 counts",fill="") +
  scale_fill_manual(values=c("CD4_0"="blue","CD4_96"="red")) +
  theme_minimal()

ggplot(mydata2 %>% filter(r==1) %>% mutate(cd4dif=cd496-cd40),aes(x=Arms,y=cd4dif,fill="CD4_gap"))+ geom_boxplot() +
  labs(x="Arms",y="CD4 counts gap",fill="") + theme_minimal()

ggplot(mydata2 %>% filter(r==1),aes(x=cd40,y=cd496)) + geom_point(alpha=0.5,size=1,shape=20) +
  geom_smooth(method="loess",se=T) +
  labs(x="CD4_0",y="CD4_96") +
  theme_minimal()

ggplot(mydata2 %>% filter(r==1) %>% filter(arms==0),aes(x=cd40,y=cd496)) + geom_point(alpha=0.5,size=1,shape=20) +
  geom_smooth(method="loess",se=T) +
  labs(x="CD4_0",y="CD4_96") +
  theme_minimal()

ggplot(mydata2 %>% filter(r==1) %>% filter(arms==1),aes(x=cd40,y=cd496)) + geom_point(alpha=0.5,size=1,shape=20) +
  geom_smooth(method="loess",se=T) +
  labs(x="CD4_0",y="CD4_96") +
  theme_minimal()

ggplot(mydata2 %>% filter(r==1) %>% filter(arms==2),aes(x=cd40,y=cd496)) + geom_point(alpha=0.5,size=1,shape=20) +
  geom_smooth(method="loess",se=T) +
  labs(x="CD4_0",y="CD4_96") +
  theme_minimal()

ggplot(mydata2 %>% filter(r==1) %>% filter(arms==3),aes(x=cd40,y=cd496)) + geom_point(alpha=0.5,size=1,shape=20) +
  geom_smooth(method="loess",se=T) +
  labs(x="CD4_0",y="CD4_96") +
  theme_minimal()

corre<-cor((mydata %>% filter(r==1))[c("cd40","cd420","cd496","cd80","cd820","preanti","age","wtkg","karnof","hemo","homo","drugs","race","offtrt")])
heatmap<-melt(corre)
ggplot(heatmap,aes(Var1,Var2,fill=value)) + geom_tile() + scale_fill_gradient2(low="blue",mid="white",high="red",
                       midpoint=0,limits=c(-1,1),breaks=c(-1,-0.5,0,0.5,1),labels=c("-1","-0.5","0","0.5","1"),guide="colorbar") +
  theme_minimal() + theme(text=element_text(size=12),plot.title=element_text(face="bold",size=14),axis.text.x=element_blank())



# Estimate missing probability
data.used<-mydata
Y_full<-data.used$cd496
C_full<-1-is.na(Y_full)
Cov_full<-data.used[c("age","wtkg","race","offtrt","hemo","drugs")]
weight_fit<-glm(C_full~.,data=cbind(C_full,Cov_full),family='binomial')
summary(weight_fit)


# Pairwise treatment effect estimation
data.obs<-data.used[which(C_full==1),]
pro_est<-predict(weight_fit,type="response")[which(C_full==1)]
weight_est<-1/pro_est
cov<-c("age","wtkg","race","gender","offtrt","drugs")
#Pair: 0 vs other 
result<-tmle(data.obs$cd496,as.numeric(data.obs$arms!=0),data.obs[cov],obsWeights=weight_est)
summary(result)
#Pair: 0 vs 1
weight_est_01<-weight_est[which(data.obs$arms==0|data.obs$arms==1)]
data.obs_01<-data.obs[which(data.obs$arms==0|data.obs$arms==1),]
result_01<-tmle(data.obs_01$cd496,as.numeric(data.obs_01$arms==1),data.obs_01[cov],obsWeights=weight_est_01)
summary(result_01)
#Pair: 0 vs 2
weight_est_02<-weight_est[which(data.obs$arms==0|data.obs$arms==2)]
data.obs_02<-data.obs[which(data.obs$arms==0|data.obs$arms==2),]
result_02<-tmle(data.obs_02$cd496,as.numeric(data.obs_02$arms==2),data.obs_02[cov],obsWeights=weight_est_02)
summary(result_02)
#Pair: 0 vs 3
weight_est_03<-weight_est[which(data.obs$arms==0|data.obs$arms==3)]
data.obs_03<-data.obs[which(data.obs$arms==0|data.obs$arms==3),]
result_03<-tmle(data.obs_03$cd496,as.numeric(data.obs_03$arms==3),data.obs_03[cov],obsWeights=weight_est_03)
summary(result_03)
# Changepoint estimation
# threshold linear regression

par(mfrow=c(1,3))
types=c("segmented","M12","M21","M22")
for (type in types) {
  fit=chngptm(formula.1=cd496~age+wtkg+race+gender+offtrt+drugs,formula.2=~cd40,data.obs,type=type,family="gaussian")
  print(summary(fit))
  print(AIC(fit))
  for (i in 1:3) plot(fit,which=i,main=type)
}
data.obs_0<-data.obs[which(data.obs$arms==0),]
par(mfrow=c(1,3))
types=c("segmented","M12","M21","M22")
for (type in types) {
  fit=chngptm(formula.1=cd496~age+wtkg+race+gender+offtrt+drugs,formula.2=~cd40,data.obs_0,type=type,family="gaussian")
  print(summary(fit))
  print(AIC(fit))
  for (i in 1:3) plot(fit,which=i,main=type)
}
data.obs_1<-data.obs[which(data.obs$arms==1),]
par(mfrow=c(1,3))
types=c("segmented","M12","M21","M22")
for (type in types) {
  fit=chngptm(formula.1=cd496~age+wtkg+race+gender+offtrt+drugs,formula.2=~cd40,data.obs_1,type=type,family="gaussian")
  print(summary(fit))
  print(AIC(fit))
  for (i in 1:3) plot(fit,which=i,main=type)
}
data.obs_2<-data.obs[which(data.obs$arms==2),]
par(mfrow=c(1,3))
types=c("segmented","M12","M21","M22")
for (type in types) {
  fit=chngptm(formula.1=cd496~age+wtkg+race+gender+offtrt+drugs,formula.2=~cd40,data.obs_2,type=type,family="gaussian")
  print(summary(fit))
  print(AIC(fit))
  for (i in 1:3) plot(fit,which=i,main=type)
}
data.obs_3<-data.obs[which(data.obs$arms==3),]
par(mfrow=c(1,3))
types=c("segmented","M12","M21","M22")
for (type in types) {
  fit=chngptm(formula.1=cd496~age+wtkg+race+gender+offtrt+drugs,formula.2=~cd40,data.obs_3,type=type,family="gaussian")
  print(summary(fit))
  print(AIC(fit))
  for (i in 1:3) plot(fit,which=i,main=type)
}
