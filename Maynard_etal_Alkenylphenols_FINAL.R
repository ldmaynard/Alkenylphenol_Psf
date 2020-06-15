rm(list=ls())


# Load Libraries ----------------------------------------------------------

library(dplyr)
library(tidyr)
library(Rmisc)
library(reshape)
library(multcomp)
library(MASS)
library(RColorBrewer)
library(ggplot2)
library(multcompView)
library(emmeans)
library(lme4)
library(AICcmodavg)
library(betareg)
library(car)
library(viridis)
library(knitr)
library(kableExtra)
library(tinytex)
library(ggpubr)
library(pals)


# Objective 2: Quantitative chemistry -------------------------------------

df_all <- read.csv(file="Maynard_etal_AlkenylphenolQuantChem.csv",head=TRUE)
#removing columns unnec. for analysis
df_all <- subset(df_all, select = -c(1, 6:32))
df_all[df_all == "#VALUE!"] <- NA
df_all[df_all == "NA"] <- 0

##DATA CLEANING
df_all$A_pdw<-as.numeric(as.character(df_all$A_pdw))
df_all$B_pdw<-as.numeric(as.character(df_all$B_pdw))
df_all$C_pdw<-as.numeric(as.character(df_all$C_pdw))
df_all$D_pdw<-as.numeric(as.character(df_all$D_pdw))
df_all$E_pdw<-as.numeric(as.character(df_all$E_pdw))
df_all$F_pdw<-as.numeric(as.character(df_all$F_pdw))
df_all$G_pdw<-as.numeric(as.character(df_all$G_pdw))
df_all$H_pdw<-as.numeric(as.character(df_all$H_pdw))
df_all$I_pdw<-as.numeric(as.character(df_all$I_pdw))
df_all$J_pdw<-as.numeric(as.character(df_all$J_pdw))

df_all$tissue<-as.factor(df_all$tissue)
df_all$plant<-as.factor(df_all$plant)
df_all$stage<-as.factor(df_all$stage)

#gather rows into columns, more data cleaning
datall <- gather(df_all, "compound", "per_dry_wt", 5:14)
datall$per_dry_wt<-as.numeric(datall$per_dry_wt)
datall$compound<-as.character(datall$compound)
datall$stage<-as.character(datall$stage)

#dry weight % -> proportions
datall$props<-(datall$per_dry_wt/100)

##splitting between tissue types
##only leaves
dfl <- filter(datall, tissue == "L")

##only pulp
dat <- filter(datall, tissue == "P") 


#labeling
datall <- datall[order(datall$tissue),]
datall <- datall[order(datall$stage),]
datall$tissue<-as.character(datall$tissue)
datall$plant<-as.character(datall$plant)
  

#sep stages within tissues 
datall$tissue[1:200]<-"R"
datall$tissue[261:470]<-"U2"
datall$tissue[471:660]<-"U3"
datall$tissue[661:870]<-"F4"
datall$tissue[871:1020]<-"F5"
datall$tissue[1021:1130]<-"F6"

datall$tissue[datall$tissue=="R"]="Ripe pulp (1)"
datall$tissue[datall$tissue=="U2"]="Unripe pulp (2)"
datall$tissue[datall$tissue=="U3"]="Unripe pulp (3)"
datall$tissue[datall$tissue=="F4"]="Flowers (4)"
datall$tissue[datall$tissue=="F5"]="Flowers (5)"
datall$tissue[datall$tissue=="F6"]="Dev. flowers (6)"
datall$tissue[datall$tissue=="S"]="Seeds"
datall$tissue[datall$stage=="M"]="Mature leaves"

datall$compound[datall$compound=="A_pdw"]="A"
datall$compound[datall$compound=="B_pdw"]="B"
datall$compound[datall$compound=="C_pdw"]="C"
datall$compound[datall$compound=="D_pdw"]="D"
datall$compound[datall$compound=="E_pdw"]="E"
datall$compound[datall$compound=="F_pdw"]="F"
datall$compound[datall$compound=="G_pdw"]="G"
datall$compound[datall$compound=="H_pdw"]="H"
datall$compound[datall$compound=="I_pdw"]="I"

#ANALYSIS
a5<-aggregate(props~tissue+plant, data=datall, FUN=sum) 

beta1<-betareg(props~tissue, data=a5)
summary(beta1)
Anova(beta1)

d1<-emmeans(beta1,pairwise~tissue, type="response")
CLD(d1$emmeans,  Letters ='ABCDEFGHIJKLMNOPQRS')


##plot
tissueplot2<-ggplot(a5, aes(x=tissue, y=props)) +
	geom_boxplot() + geom_point(aes(color=as.character(plant)))+
	labs(x=" ", y="Total alkenylphenols (% dry wt)")+
	theme_minimal()+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds","Exp. leaves", "Mature leaves"))+
	scale_color_manual(values=as.vector(kelly(21)), name = "Plant ID")
tissueplot2

##too busy, selecting representative stages for each tissue

##selecting representative stages for each tissue
##"Ripe pulp (1)","Unripe pulp (2)", "Flowers (4)","Seeds","Mature leaves"
datall <- datall[order(datall$stage),]
dat_select<-datall[c(1:470,661:870,1131:1170),]

#ANALYSIS
a10<-aggregate(props~tissue+plant, data=dat_select, FUN=sum)

beta10<-betareg(props ~ tissue, data=a10)
summary(beta10)

d10<-emmeans(beta10,pairwise~tissue, type="response")
d10
CLD(d10$emmeans,  Letters ='ABCDEFGHIJKLMNOPQRS')


##plot
a10$tissue[a10$tissue=="Unripe pulp (2)"]="Late unripe pulp"
a10$tissue[a10$tissue=="Flowers (4)"]="Late flowers"
a10$tissue[a10$tissue=="Ripe pulp (1)"]="Ripe pulp"


tissueplot<-ggplot(a10, aes(x=tissue, y=props)) +
	geom_boxplot(outlier.shape = NA) + geom_jitter(position=position_jitter(width = 0.04), alpha=0.25)+
	labs(x=" ", y="Total alkenylphenols (prop. dw)")+
	theme_classic()+
	scale_x_discrete(limits=c("Mature leaves","Late flowers", "Late unripe pulp","Ripe pulp",
							  "Seeds"))+
	stat_summary(geom = 'text', label = c("B","A","C","A","C"),
				 fun.y = max, vjust = -0.8)+
	scale_y_continuous(limits = c(0, 0.120))+
	theme(text = element_text(size=15),
		  axis.text.x = element_text(angle=45, hjust=1))+
	stat_summary(fun.y=mean, colour="black", geom="point", 
				 shape=18, size=4,show.legend = FALSE)
tissueplot

#EXPORT PLOT
tiff('Maynard_etal_Fig3.tiff', units="in", width=6, height=5, res=400)
tissueplot
dev.off()

#barplot (means and SEs)
a10_se <- summarySE(a10, measurevar="props", groupvars=c("tissue"))
tissue_bar<-ggplot(a10_se, aes(x=tissue, y=props, fill=tissue)) +
	geom_bar(stat = "summary")+
	labs(x=" ", y="Total alkenylphenols (prop. dw)")+
	theme_classic()+
	scale_x_discrete(limits=c("Mature leaves","Late flowers", "Late unripe pulp","Ripe pulp",
							  "Seeds"))+
	scale_y_continuous(limits = c(0, 0.06))+
	theme(text = element_text(size=15),
		  axis.text.x = element_text(angle=45, hjust=1), legend.position = "none")+
	geom_errorbar(aes(ymin=props-se, ymax=props+se), width=.1)+
	scale_fill_viridis(discrete = T, option = "D")
tissue_bar

#SUMMARY STATS
library(plyr)
tissue.tab <- ddply(a10, c("tissue"), summarise,
               N    = length(props),
               mean = mean(props),
               sd   = sd(props),
               se   = sd / sqrt(N))
tissue.tab

#EXPORT TABLE
#write.table(tissue.tab, file = "TableS4.csv", sep = ",", quote = FALSE, row.names = F)

 ##Alkenylphenols over fruit ripening------------------------------------

#creating columns for stage as continuous variable and stage^2 (quadratic term)
dat[is.na(dat)] <- 0
dat$stage2 = (as.numeric(dat$stage))^2
dat$stage.num = as.numeric(dat$stage)

#aggregating dataset so each point=one plant and isn't sep by compound
ag_dat<-aggregate(props~plant+stage.num+stage2,data=dat,FUN=sum)
ag_dat$plant<-as.character(ag_dat$plant)

#betaregression with quadratic term
beta12<-betareg(props~stage.num+stage2 , data=ag_dat)
summary(beta12)#stage and quad. term are significant
#since quad. term is significant, suggests the data are not linear
#stage p=0.0234, quad term p<0.001

#prediction plot
ag_dat$yhat_12<-predict(beta12)
predplot12<-ggplot(ag_dat)+
	geom_point(aes(x=stage.num, y=props))+
	geom_point(aes(x=stage.num, y=yhat_12), color="red", size=2)+
	geom_line(aes(x=stage.num, y=yhat_12) ,color="red", size=1)
predplot12

#regular betaregression
beta13<-betareg(props~stage.num, data=ag_dat)
summary(beta13)#stage is significant

#prediction plot
ag_dat$yhat_13<-predict(beta13)
predplot13<-ggplot(ag_dat)+
	geom_point(aes(x=stage.num, y=props))+
	geom_point(aes(x=stage.num, y=yhat_13), color="red", size=2)+
	geom_line(aes(x=stage.num, y=yhat_13) ,color="red", size=1)
predplot13

#prediction plot with all three models
predplot_all<-ggplot(ag_dat)+
	geom_point(aes(x=stage.num,y=props))+
	geom_line(aes(x=stage.num, y=yhat_12),color="red",size=1)+
	geom_line(aes(x=stage.num, y=yhat_13),color="blue",size=1)
predplot_all

null.mod<-betareg(props~1, data=ag_dat)

modcomp<-aictab(cand.set=list(beta12,beta13,null.mod),
	   modnames=c("nonlinear","linear","null"))#AIC table

modcomp #top model is nonlinear (mod with quad term)


stage_line_plot<-ggplot(ag_dat)+
	geom_line(aes(x=stage.num,y=props,color=plant),alpha=0.5,size=1,show.legend = F)+
	stat_smooth(aes(x=stage.num,y=props),method = "lm", formula = y ~ x + I(x^2), size = 1.5,
				linetype="solid", color="black")+
	labs(x=" ", y="Total alkenylphenols (prop. dw)")+
	theme_classic()+
	scale_color_viridis(discrete = T, option = "D")+
	scale_x_reverse(breaks=c(6,5,4,3,2,1), expand=c(0,0))+
	theme(text = element_text(size = 15), axis.text.x.bottom = element_blank())+
	scale_y_continuous(expand=c(0.0,0.0))+
	coord_cartesian(xlim=c(1.0,6.0), ylim=c(0.0,.1))
stage_line_plot


##EXPORT PLOT
tiff('Maynard_etal_Fig4.tiff', units="in", width=8, height=5, res=400)
stage_line_plot
dev.off()

#SUMMARY STATS
library(plyr)
dev.tab <- ddply(ag_dat, c("stage.num"), summarise,
                    N    = length(props),
                    mean = mean(props),
                    sd   = sd(props),
                    se   = sd / sqrt(N))
dev.tab

#EXPORT TABLE
#write.table(dev.tab, file = "TableS5.csv", sep = ",", quote = FALSE, row.names = F)

##summary stats
library(plyr)
cdata <- ddply(datall, c("tissue", "compound"), summarise,
			   N    = length(props),
			   mean = mean(props),
			   sd   = sd(props),
			   se   = sd / sqrt(N))

cdata <- cdata[order(cdata$tissue),]
cdata[is.na(cdata)] <- 'nd'
head(cdata)

#EXPORT SUMMARY STATS TABLE
#write.table(cdata, file = "stats.csv", sep = ",", quote = FALSE, row.names = F)

##SUMMARY PLOT
datall <- datall[order(datall$compound),]
datall[is.na(datall)] <- 0
a.dat<-slice(datall, 1:117)
b.dat<-slice(datall, 118:234)
c.dat<-slice(datall, 235:351)
d.dat<-slice(datall, 352:468)
e.dat<-slice(datall, 469:585)
f.dat<-slice(datall, 586:702)
g.dat<-slice(datall, 703:819)
h.dat<-slice(datall, 820:936)
i.dat<-slice(datall, 937:1053)
j.dat<-slice(datall, 1054:1170)

a.plot<-ggplot(a.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (prop. dw)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))+
  coord_flip()

a.plot

b.plot<-ggplot(b.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (prop. dw)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))+
  coord_flip()

c.plot<-ggplot(c.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (prop. dw)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))+
  coord_flip()

d.plot<-ggplot(d.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (prop. dw)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))+
  coord_flip()

e.plot<-ggplot(e.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (prop. dw)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15))+
  coord_flip()

f.plot<-ggplot(f.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (prop. dw)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15), axis.text.y=element_blank())+
  coord_flip()

g.plot<-ggplot(g.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (prop. dw)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15), axis.text.y=element_blank())+
  coord_flip()

h.plot<-ggplot(h.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (prop. dw)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15), axis.text.y=element_blank())+
  coord_flip()

i.plot<-ggplot(i.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (prop. dw)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15), axis.text.y=element_blank())+
  coord_flip()

j.plot<-ggplot(j.dat, aes(x=tissue,y=props))+
	geom_boxplot(aes(fill=tissue)) + geom_point()+
	labs(x=" ", y="Total alkenylphenols (prop. dw)")+
	scale_x_discrete(limits=c("Dev. flowers (6)",
							  "Flowers (5)","Flowers (4)","Unripe pulp (3)", "Unripe pulp (2)",
							  "Ripe pulp (1)","Seeds", "Mature leaves"))+
	theme_minimal()+
	scale_fill_viridis(discrete = T, option = "D")+ 
	theme(legend.position="none",text = element_text(size=15), axis.text.y=element_blank())+
  coord_flip()

sum.plot1<-ggarrange(a.plot, f.plot,b.plot,g.plot,c.plot,h.plot,d.plot,i.plot,e.plot,j.plot, 
                     labels = c("A", "F","B","G","C","H","D","I","E","J"),heights = c(20, 20),
                     ncol = 2, nrow = 5, align = "v",
                     font.label = list(size = 22))
sum.plot1


# Objective 3: Fungal bioassays --------------------------------------------------------

datf<- read.csv(file="Maynard_etal_FungalBioassays.csv",head=TRUE,fill=T)
datf <- datf[1:81,] 

datf$Conc<-as.numeric(datf$Conc)
datf$fungi<-as.character(datf$fungi)
datf$nfungi[datf$fungi=="R3"]="Microdochium lycopodinum"
datf$nfungi[datf$fungi=="R23"]="Fusarium A"
datf$nfungi[datf$fungi=="R26"]="Fusarium B"

#average triplicates
ag.fun<-aggregate(abs_corr~Conc+nfungi+well, data=datf, FUN=mean) 
mod.fun<-lm(abs_corr~Conc*nfungi, data=ag.fun)
summary(mod.fun)

mod.add<-lm(abs_corr~Conc+nfungi, data=ag.fun)
mod.con<-lm(abs_corr~Conc, data=ag.fun)
mod.fung<-lm(abs_corr~nfungi, data=ag.fun)
mod.null<-lm(abs_corr~1, data=ag.fun)

modcomp.fung<-aictab(cand.set=list(mod.fun,mod.add,mod.con,mod.fung,mod.null),
				modnames=c("Interaction","Add","Concentration", "Fung. sp.", "Null"))#AIC table

modcomp.fung #top model is interactive


Anova(mod.fun)#two-way ANOVA

library(effects)
summary(allEffects(mod.fun))
plot(allEffects(mod.fun))

datf <- datf[order(datf$fungi),]
ag.fun<-ag.fun[order(ag.fun$nfungi),]

R23<-slice(ag.fun, 1:9) #Fusarium A
R26<-slice(ag.fun, 10:18) #Fusarium B
R3<-slice(ag.fun, 19:27) #M. lychopodinum

mod23<-lm(data = R23, abs_corr~Conc)
mod26<-lm(data = R26, abs_corr~Conc)
mod3<-lm(data = R3, abs_corr~Conc)

summary(mod23)#t=-7.03,p=0.0002,r2=0.88
summary(mod26)#t=-1.0, p=0.351, r2=0.13
summary(mod3)#t=-5.213, p=0.00124, r2=0.80

mod23#m= -0.4339 
#-0.43  Fusarium A dec every 1 mg/mL increase

mod3 #m= -0.4985
#0.50 Ml decrease every 1 mg/mL increase


R23$yhat<-predict(mod23)
predplot23<-ggplot(R23)+
	geom_smooth(aes(x=Conc, y=yhat), color="red")+
	geom_point(aes(x=Conc, y=abs_corr))
predplot23

R26$yhat<-predict(mod26)
predplot26<-ggplot(R26)+
	geom_smooth(aes(x=Conc, y=yhat), color="red")+
	geom_point(aes(x=Conc, y=abs_corr))
predplot26

R3$yhat<-predict(mod3)
predplot3<-ggplot(R3)+
	geom_smooth(aes(x=Conc, y=yhat), color="red")+
	geom_point(aes(x=Conc, y=abs_corr))
predplot3

##back-calculating alkenylphenol concentrations in wells
#ripe fruit mean (10mg) = 0.2925986
#average concentration of alkenylphenols in 10mg ripe fruit
0.2925986/10

#alkenylphenols in 1mg ripe fruit=0.02925986
#alkenylphenols in 5g (5g starting plant material)
0.02925986*5000

#in the vial, there was 146.2993 mg of alkenylphenols

#146.2993mg resuspended in 2mL=73.14965mg/mL
146.2993/2

#Heather's dilutions
73.14965/80

#highest concentration in well is 0.9143706 mg/ml
0.9143706/3000
#~0.0003047902 prop fresh weight of fruit


##DATA PLOT
###load italics for legend
leg_fung <- c(expression(paste(italic("Fusarium "), "A")),
			 expression(paste(italic("Fusarium "), "B")), 
			 expression(paste(italic("Microdochium lycopodinum"))))

#plot
plota<-ggplot(ag.fun, aes(x=Conc, y=abs_corr, group=nfungi))+
	geom_smooth(aes(color=nfungi), method = "lm", se=T)+
	geom_point(aes(color=nfungi))+
	theme_classic()+
	scale_color_viridis(discrete = T, option = "D", labels=leg_fung)+
	labs(x="Alkenylphenol concentration (mg/mL)", y="Average absorbance (OD)", color=" ")+
	theme(legend.text.align = 0, text = element_text(size=18),legend.position="top",
		  plot.margin = margin(10, 30, 10, 10))+ 
	scale_x_continuous(expand = c(0, 0), limits = c(0.0,0.92), 
					   breaks = c(0.0,0.23,0.46,0.69, 0.92))+
	annotate("text", x = 0.8, y = 0.68,
			 label = "paste(italic(R) ^ 2, \" = 0.88\")", parse = TRUE, size =5)+ 
	annotate("text", x = 0.8, y = 0.3,
			 label = "paste(italic(R) ^ 2, \" = 0.80\")", parse = TRUE, size =5)
#to include R^2 for Fusarium B,
#annotate("text", x = 0.8, y = 0.465, label = "paste(italic(R) ^ 2, \" = 0.13\")", parse = TRUE, size =5)

plota
#Sometimes error message "polygon edge not found", just run it again
plota

#EXPORT PLOT
tiff('Maynard_etal_Fig5.tiff', units="in", width=8, height=5, res=400)
plota
dev.off()


# Objective 4b: Animal preference trials ------------------------------------------------

animal <- read.csv(file="Maynard_etal_AlkenylphenolAnimalTrials.csv",head=TRUE)

#only animals that participated
animal<-animal%>%
	filter(participate==1)

#separating bats and birds 
bat<-slice(animal, 1:31)
bird<-slice(animal, 32:58)


birdg <- gather(bird, "treatment", "amount_eaten", 17:18)
birdg$amount_eaten<-as.numeric(birdg$amount_eaten)
birdg$amount_eaten[birdg$amount_eaten>3] <- 3
birdg$treatment<-as.character(birdg$treatment)
bird_ag<-aggregate(amount_eaten~ID+treatment, data=birdg, FUN=mean) 
bird_1<-spread(data=bird_ag, treatment, amount_eaten)

shapiro.test(bird_1$c.eaten)#norm dist
shapiro.test(bird_1$t.eaten)#norm dist

batg <- gather(bat, "treatment", "amount_eaten", 17:18)
batg$amount_eaten<-as.numeric(batg$amount_eaten)
batg$treatment<-as.character(batg$treatment)
bat_ag<-aggregate(amount_eaten~ID+treatment, data=batg, FUN=mean) 
bat_1<-spread(data=bat_ag, treatment, amount_eaten)

shapiro.test(bat_1$t.eaten)#not norm dist
shapiro.test(bat_1$c.eaten)#norm dist

t.test(bird_1$c.eaten, bird_1$t.eaten, paired = T)#t=0.24,df=9,p=0.8104
t.test(bat_1$c.eaten, bat_1$t.eaten, paired = T)#t=3.8959,df=15,p=0.00143

#change names for graph
bat_ag$treatment[bat_ag$treatment=="t.eaten"]="Treatment"
bat_ag$treatment[bat_ag$treatment=="c.eaten"]="Control"
bird_ag$treatment[bird_ag$treatment=="t.eaten"]="Treatment"
bird_ag$treatment[bird_ag$treatment=="c.eaten"]="Control"

batpref<-ggplot(bat_ag,aes(x=treatment,y=amount_eaten))+geom_boxplot()+
	geom_jitter(position=position_jitter(width = 0.03), alpha=0.25, size=1.7)+
	theme_classic()+
	labs(x=" ",y="Average amount eaten (g)")+
	scale_y_continuous(limits =  c(0,3))+
	theme(text = element_text(size = 18))+
	stat_summary(fun.y=mean, colour="black", geom="point", 
				 shape=18, size=5,show.legend = FALSE)
batpref

birdpref<-ggplot(bird_ag,aes(x=treatment,y=amount_eaten))+geom_boxplot()+
	geom_jitter(position=position_jitter(width = 0.03), alpha=0.25,size=1.7)+
	theme_classic()+
	labs(x=" ",y="")+
	scale_y_continuous(limits =  c(0,3.01))+
	theme(text = element_text(size = 18))+
	stat_summary(fun.y=mean, colour="black", geom="point", 
				 shape=18, size=5,show.legend = FALSE)
birdpref

animal.plot<-ggarrange(batpref, birdpref, 
                       labels = c("a", "b"),heights = c(2, 2),
                       ncol = 2, nrow = 1, align = "v")
animal.plot

#EXPORT PLOT
tiff('Maynard_etal_Fig6.tiff', units="in", width=8, height=4, res=400)
animal.plot
dev.off()

##barplot
bat_agse <- summarySE(bat_ag, measurevar="amount_eaten", groupvars=c("treatment"))

batpref_bar<-ggplot(bat_agse,aes(x=treatment,y=amount_eaten, fill=treatment))+geom_bar(stat = "summary")+
	theme_classic()+
	labs(x=" ",y="Average amount eaten (g)")+
	scale_y_continuous(limits =  c(0,3))+
	theme(text = element_text(size = 18), legend.position = "none")+
	geom_errorbar(aes(ymin=amount_eaten-se, ymax=amount_eaten+se), width=.1)+
	scale_fill_viridis(discrete = T, option = "D")
batpref_bar

bird_agse <- summarySE(bird_ag, measurevar="amount_eaten", groupvars=c("treatment"))
birdpref_bar<-ggplot(bird_agse,aes(x=treatment,y=amount_eaten, fill=treatment))+
	geom_bar(stat = "summary")+
	theme_classic()+
	labs(x=" ",y="")+
	scale_y_continuous(limits =  c(0,3))+
	theme(text = element_text(size = 18), legend.position = "none")+
	geom_errorbar(aes(ymin=amount_eaten-se, ymax=amount_eaten+se), width=.1)+
	scale_fill_viridis(discrete = T, option = "D")
birdpref_bar

animal.plot.bar<-ggarrange(batpref_bar, birdpref_bar, 
						   labels = c("a", "b"),heights = c(2, 2),
						   ncol = 2, nrow = 1, align = "v")
animal.plot.bar

#EXPORT PLOT
tiff('animal_pref_bar.tiff', units="in", width=8, height=4, res=500)
animal.plot.bar
dev.off()

bat.tab <- ddply(bat_ag, c("treatment"), summarise,
                 N    = length(amount_eaten),
                 mean = mean(amount_eaten),
                 sd   = sd(amount_eaten),
                 se   = sd / sqrt(N))
bat.tab

bird.tab <- ddply(bird_ag, c("treatment"), summarise,
                 N    = length(amount_eaten),
                 mean = mean(amount_eaten),
                 sd   = sd(amount_eaten),
                 se   = sd / sqrt(N))
bird.tab


# Objective 4a: Removal study####
remo <- read.csv(file="Maynard_etal_PiperRemovalStudy.csv",head=TRUE,fill=T)

remo<-remo[-c(71),]

remo$Ripeness <- NA
for(i in 1:length(remo$Ripeness)){
	if(remo$ripeness[i]=="1"){remo$Ripeness[i]="Ripe"}
	if(remo$ripeness[i]=="2"){remo$Ripeness[i]="Unripe"}
}

remo$Time <- NA
for(i in 1:length(remo$Time)){
	if(remo$time.of.day[i]=="1"){remo$Time[i]="Nocturnal"}
	if(remo$time.of.day[i]=="2"){remo$Time[i]="Diurnal"}
}

tbl = table(remo$Ripeness, remo$Time) 
tbl
chisq.test(tbl)

tbl.1 = table(remo$Time, remo$Ripeness)
tbl.1
chisq.test(tbl.1)

tbl2<-as.data.frame.matrix(tbl)
tbl2

tbl3<-as.data.frame(tbl)
tbl3

bp <- ggplot(tbl3, aes(x=Var2, y=Freq, fill=Var1))+
	geom_bar(aes(fill=Var1),stat = "identity", position = "stack")+
	theme_classic()+
	labs(x="", y="No. fruits removed",fill="")+
	theme(text = element_text(size = 18))+
	scale_fill_manual(values=c('#006837','#9EBE7D'))+
	coord_flip()+
	theme(legend.position="top")
	
bp

# Natural history information (Objective 4)####
birb <- read.csv(file="Maynard_etal_Birds_Psf.csv",head=TRUE,fill=T)
birb[is.na(birb)] <- 0
birb$gleaning<-as.numeric(birb$gleaning)
birb$frugivory<-as.numeric(birb$frugivory)
birb$cover.perch<-as.numeric(birb$cover.perch)
birb$defense<-as.numeric(birb$defense)
birb$calling<-as.numeric(birb$calling)
birb$socializing<-as.numeric(birb$socializing)
birb$parental.care<-as.numeric(birb$parental.care)


#summary table of bird activities
#Warning, sum not meaningful for factors. Not sure why it isn't recognizing the commands above...
birb$no.visits<-aggregate(gleaning+frugivory+cover.perch+defense+
						calling+socializing+parental.care~sp, data=birb, FUN=sum) 

head(birb.sum)

birb.melt<-melt(birb, id.vars="sp", measure.vars=c('gleaning', 'frugivory', 'cover.perch',
												   'defense','calling','socializing','parental.care'))

colnames(birb.melt$gleaning) <- birb.melt$Gleaning

birb.plot<-ggplot(birb.melt) +
	geom_bar(aes(x=sp, y=value, fill=as.factor(variable)), stat = "identity")+
	labs(y="No. visits", fill = "Activity", x="")+
	scale_fill_manual(labels=c("Gleaning","Frugivory", "Cover/perch", "Defense", "Calling",
							   "Socializing", "Parental care"), 
					  values = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f", "#e5c494"))+
	theme_classic()+
	theme(text = element_text(size=13),
											 axis.text.x = element_text(angle=45, hjust=1))
birb.plot

#EXPORT PLOT
tiff('Maynard_etal_FigS2.tiff', units="in", width=8, height=4, res=500)
birb.plot
dev.off()
