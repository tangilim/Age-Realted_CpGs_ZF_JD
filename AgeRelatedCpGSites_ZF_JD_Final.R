#Age-related CpG sites in zebra finch AND jackdaw 
#Marianthi Tangili
#last updated 19/03/2024

setwd("C:/Users/MWP-/Dropbox/ZF_JD_AgeRelatedSites/Draft/Data_Submission")
setwd("C:/Users/MWP-/Dropbox/ZF_JD_AgeRelatedSites/Figures")

#load necessary libraries
library(ggrepel)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(readr)
library(janitor)
library(readxl)
library(gridExtra)
library(grid)
library(ggchicklet)
library(tidyverse)

#load graph theme
graph_theme<- theme_update(
  panel.grid.major=element_line(colour=NA),
  panel.grid.minor=element_line(colour=NA),
  panel.background = element_rect(colour =NA ,fill=NA,size=0.5),
  axis.title.x=element_text(size=15,hjust=0.5,vjust=0.5,angle=0),
  axis.title.y=element_text(size=15,hjust=0.5,vjust=0.5,angle=90),
  axis.text.x=element_text(colour="black",angle=0,size=10),
  axis.text.y=element_text(colour="black",angle=0,size=10),
  axis.ticks=element_line(colour="black",size=0.5),
  axis.line.x=element_line(size=0.5),
  axis.line.y=element_line(size=0.5))

###### ZEBRA FINCH #######


# ZF_Age-related CpG sites ------------------------------------------------

#data frame from pipeline run in Habrok (Cor_ZF_noMT_Mar24_W.R)
s_zf <- read_csv("ZF_AgeRelated_CpG_noMT_Final.csv")

length(unique(s_zf$Pos))

#separate Pos to get Chromosome
s_zf <- s_zf %>%
  separate(Pos, into = c("Chromosome", "Position"), sep = "_", remove = FALSE)

#subset to only keep sites in W
w<- subset(s_zf, s_zf$Chromosome=="W")

#how many samples each position is shared by
counts <- table(s_zf$Pos)

#how many samples each position is shared by (only W)
counts_w<- table(w$Pos)

#subset the cor data frame to only keep sites shared by >=15 samples
cor_subset <- counts[counts >= 15]

#subset the cor data frame to only keep sites shared by >=6 samples for W
cor_subset_w<- counts_w[counts_w >= 6]

# View the subsetted data frame
cor_subset<- as.data.frame(cor_subset)

cor_subset_w<- as.data.frame(cor_subset_w)

#merge the subset for w and rest of chromosomes
final<- rbind(cor_subset, cor_subset_w)

length(unique(final$Var1))

#change colnames
colnames(final)<- c("Pos", "Freq")

#separate Pos to get chromosome
final<- final %>%
  separate(Pos, into = c("Chromosome", "Position"), sep = "_", remove = FALSE)

#write.csv(final[,2:3], "ZF_Positions_AR_Final.csv")

#how many sitees per Chromosome
ft<- tabyl(final$Chromosome)
#change colnames
colnames(ft)<- c("Chromosome", "n", "percent")

# ZF-AR-CpGs-TotalCpGs+Length---------------------------------------------
#corelate with total CpGs in each chromosome
p_zf <- read_csv("ZF_PositionsPerChr_All_Final.csv")

ftp<- merge(ft, p_zf, by="Chromosome")

ftp<- ftp[,c(1:3,5)]

names(ftp)[names(ftp) == 'Unique_Pos_Count'] <- 'TotalCpGs'

#add info on chromosome length 
length_zf <- read_csv("ZF_Chromosome_Length.csv")

ftpl<- merge(ftp, length_zf, by="Chromosome")

ftpl$Length_kbp<- ftpl$Length/1000

#find what percentage of all CpGs/chromosome are age-related
ftpl$percentage<- (ftpl$n*100)/ftpl$TotalCpGs

#write.csv(ftpl, "FinalPosLength_ZF.csv")

ftpl$p<- (0.006*ftpl$TotalCpGs)/100

(sum(ftpl$n)*100)/sum(ftpl$TotalCpGs)

s_all<- ggplot(ftpl) +
  geom_point(aes(x = TotalCpGs, y = n, color=Chromosome, size=Length_kbp)) +
  geom_text(data = subset(ftp, Chromosome %in% c(10,32, 35, "Z", "W")),
            aes(x = TotalCpGs, y = n,label = Chromosome),
            hjust = -0.4, vjust = 0, size = 6) +
  scale_color_manual(values = c("W" = "orange", "Z" = "turquoise1", "10"="black", "32"="black","35"="black"))+
  ylab(NULL) +
  xlab(NULL) +
  labs(size = "Chromosome Length (kbp)")+
  ylim(0,220)+
  guides(color = FALSE)+theme(legend.title = element_text(size=16), legend.text = element_text(size=16), aspect.ratio = 0.5)+
  theme(
    legend.position = c(.95, .95),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    panel.border = element_rect(colour = "black", fill=NA, size=1)
  )+
  scale_size_continuous(name = "Chromosome Length (kbp)", labels = scales::comma)+
  scale_x_continuous(labels = scales::comma, limits = c(0,2100000)) +
  geom_line(aes(x=TotalCpGs, y=p),linetype = "dotted", size=0.6) +
  geom_point(aes(x = TotalCpGs, y = n, color=Chromosome, size=Length_kbp), data=subset(ftpl, ftpl$Chromosome=="10"))
s_all


# ZF-Increasing/Decreasing DNAm -------------------------------------------
#correlation information for CpG sites
info_zf <- read_csv("ZF_CpG_Info_Final.csv")

zf_cor<- merge(final, info_zf, by="Pos")

zf_cor$Chromosome <- factor(zf_cor$Chromosome, levels=c('1', '1A', '2', '3', '4','4A', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '34', '35', '37','W', 'Z'))

zf_cor_count <- zf_cor %>%
  group_by(Chromosome) %>%
  summarise(Increasing = sum(mCor > 0),
            Decreasing = sum(mCor < 0))

result_long <- tidyr::gather(zf_cor_count, key = "Change", value = "Count", -Chromosome)

result_long$Chromosome <- factor(result_long$Chromosome, levels=c('1', '1A', '2', '3', '4','4A', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '34', '35', '37','W', 'Z'))

#count of age-related CpG sites increasing/decreasing in DNAm with age
result_long %>%
  group_by(Change) %>%
  summarise(Sum_Count = sum(Count))

#stacked bar graph increasing/decreasing per chromosome
zf_incdec_stack<- ggplot(result_long, aes(x = Chromosome, y = Count, fill = Change)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = NULL,
       y = NULL,
       fill = "Change") +
  scale_fill_manual(values = c("Increasing" = "azure4", "Decreasing" = "black"))+
  scale_y_continuous(expand = c(0,0))
zf_incdec_stack

# ZF-location AR-CpGs -----------------------------------------------------
zf_cor_length<- merge(zf_cor, length_zf, by="Chromosome")
zf_cor_length$Position<- as.numeric(zf_cor_length$Position)
zf_cor_length$Chromosome <- factor(zf_cor_length$Chromosome, levels=c('1', '1A', '2', '3', '4','4A', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '34', '35', '37','W', 'Z'))

u<- zf_cor_length %>% distinct(Chromosome, Length)
u$Length_kbp<- u$Length/1000
zf_cor_length$Length_kbp<- zf_cor_length$Length/1000
zf_cor_length$Position_kbp<- zf_cor_length$Position/1000

#plot with lines on chromosomes for each age-related CpG site
u_zf<- ggplot(u) +
  geom_chicklet(aes(x = Chromosome, y = Length_kbp), col="black", fill="white")  +
  scale_y_continuous(expand = c(0, 0),limits=c(0,160000),labels = scales::comma)+
  geom_errorbar(data=zf_cor_length, aes(x = as.numeric(Chromosome), ymin = Position_kbp, ymax = Position_kbp), color = "navyblue", size = 0.2, width = 0.85) +
  labs(x = NULL, y = "Length (kbp)") 

u_zf


# ZF_Lenient --------------------------------------------------------------
#data frame from pipeline run in Habrok (Cor_ZF_noMT_Oct23_W.R)
s_zf_len<- read_csv("ZF_AgeRelated_CpG_noMT_LENIENT_Final.csv")
View(s_zf_len)

#separate Pos to get chromosome
s_zf_len <- s_zf_len %>%
  separate(Pos, into = c("Chromosome", "Position"), sep = "_", remove = FALSE)

#subset to only keep sites in W
w_len<- subset(s_zf_len, s_zf_len$Chromosome=="W")

#how many samples each position is shared by
counts_len <- table(s_zf_len$Pos)
#how many samples each position is shared by (only W)
counts_w_len<- table(w_len$Pos)

#subset the cor data frame to only keep sites shared by >=15 samples
cor_subset_len <- counts_len[counts_len >= 15]
#subset the cor data frame to only keep sites shared by >=6 samples for W
cor_subset_w_len<- counts_w_len[counts_w_len >= 6]
# View the subsetted data frame
cor_subset_len<- as.data.frame(cor_subset_len)
cor_subset_w_len<- as.data.frame(cor_subset_w_len)
#merge the subset for w and rest of chromosomes
final_len<- rbind(cor_subset_len, cor_subset_w_len)
#change colnames
colnames(final_len)<- c("Pos", "Freq")
#separate Pos to get chromosome
final_len<- final_len %>%
  separate(Pos, into = c("Chromosome", "Position"), sep = "_", remove = FALSE)
#write.csv(final_len[,2:3], "ZF_Positions_AR_Final_Len.csv")
#how many sitees per Chromosome
ft_len<- tabyl(final_len$Chromosome)
#change colnames
colnames(ft_len)<- c("Chromosome", "n", "percent")

#corelate with total CpGs in each chromosome

ftp_len<- merge(ft_len, p_zf, by="Chromosome")

ftp_len<- ftp_len[,c(1:3,5)]

names(ftp_len)[names(ftp_len) == 'StartPosition'] <- 'TotalCpGs'

ftp_zf<- merge(ftp, ftp_len, by="Chromosome")

sum(ftp_len$n)

ggplot(data=ftp_zf)+geom_point(aes(x=n.y, y=n.x))

#correlation between default and lenient pipelines
cor.test(ftp_zf$n.x,ftp_zf$n.y)


# ZF_Genomic_Annotation ---------------------------------------------------

zf_anot<- readRDS("ZF_Annotation_All.rds")

zf_loc<- merge(zf_anot, final, by="Pos")

colnames(zf_loc)[colnames(zf_loc) == 'category'] <- 'Location'

zf_loc[zf_loc=="prom"]<- "Promoter"
zf_loc[zf_loc=="exon"]<- "Exon"
zf_loc[zf_loc=="intron"]<- "Intron"
zf_loc[zf_loc=="intergenic"]<- "Intergenic"

zf_loc$Location<- factor(zf_loc$Location, levels=c("Promoter", "Exon", "Intron", "Intergenic"))

prop_zf_anot <- prop.table(table(zf_loc$Location))
prop_zf_anot<- as.data.frame(prop_zf_anot)
prop_zf_anot$ObsExp<- "Observed"
colnames(prop_zf_anot)<- c("Location", "Proportion", "ObsExp")

zf_anot[zf_anot=="prom"]<- "Promoter"
zf_anot[zf_anot=="exon"]<- "Exon"
zf_anot[zf_anot=="intron"]<- "Intron"
zf_anot[zf_anot=="intergenic"]<- "Intergenic"


prop_zf_anot_exp<-  prop.table(table(zf_anot$category))
prop_zf_anot_exp<- as.data.frame(prop_zf_anot_exp)
prop_zf_anot_exp$ObsExp<- "Expected"
colnames(prop_zf_anot_exp)<- c("Location", "Proportion", "ObsExp")

zf_loc_all<- rbind(prop_zf_anot, prop_zf_anot_exp)

zf_loc_all$ObsExp<- factor(zf_loc_all$ObsExp, levels=c("Observed", "Expected"))


loc_zf<- ggplot(zf_loc_all, aes(x=ObsExp, y=Proportion,fill=Location))+
  geom_bar( position = "fill", stat="identity", width=0.5)+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_brewer(palette="RdBu")+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.6))+
  xlab(NULL)+  ylab(NULL)
loc_zf

loc_zf1<- ggplot(data=subset(zf_loc_all, zf_loc_all$ObsExp=="Observed"))+
  geom_bar(aes(x=Location, fill=Location, y=Proportion), stat="identity",  width=0.75)+
  scale_y_continuous(expand = c(0, 0), limits=c(0,0.6))+
  scale_fill_brewer(palette="RdBu")+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.6))+
  xlab(NULL)+ ylab("Proportion of AR-CpG sites")
loc_zf1

#write.csv(zf_loc_all, "ZF_CpG_Annotation.csv")

#Per chromosome type

zf_loc$type <- ifelse(zf_loc$chromossome == "W", "W", 
                      ifelse(zf_loc$chromossome == "Z", "Z", "Autosomes"))

prop_zf_type<- prop.table(table(zf_loc$Location,zf_loc$type))
prop_zf_type<- as.data.frame(prop_zf_type)
colnames(prop_zf_type)<- c("Location", "Chr", "Proportion")
prop_zf_type$ObsExp<- "Observed"

zf_anot$type <- ifelse(zf_anot$chromossome == "W", "W", 
                       ifelse(zf_anot$chromossome == "Z", "Z", "Autosomes"))

prop_zf_type_exp<- prop.table(table(zf_anot$category,zf_anot$type))
prop_zf_type_exp<- as.data.frame(prop_zf_type_exp)
colnames(prop_zf_type_exp)<- c("Location", "Chr", "Proportion")
prop_zf_type_exp$ObsExp<- "Expected"

zf_aut<- rbind(prop_zf_type, prop_zf_type_exp)

zf_aut$Location<- factor(zf_aut$Location, levels=c("Promoter", "Exon", "Intron", "Intergenic"))

aut_zf<- ggplot(prop_zf_type, aes(x=Chr, y=Proportion,fill=Location))+
  geom_bar( position = "fill", stat="identity", width=0.65)+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_brewer(palette="RdBu")+
  scale_color_brewer(palette="RdBu")+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.6),
         legend.title = element_text(size = 16),
         legend.text =  element_text(size = 16))+
  xlab(NULL)+  ylab(NULL)

aut_zf



# ZF_Promoters ------------------------------------------------------------

#subset for promoters
zf_anot_prom<- subset(zf_anot, zf_anot$prom=="1")


#make variable Pos
zf_anot_prom <- zf_anot_prom %>%
  mutate(Pos = paste(chromossome, site, sep = "_"))

#merge final CpGs of analysis with information on increasing/decreasing 
zf_mrged<- merge(final, info_zf, by="Pos")
#merge that with CpGs in promoters
zf_merged<- merge(zf_mrged, zf_anot_prom, by="Pos")

zf_merged_count<- zf_merged %>%
  group_by(Chromosome) %>%
  summarise(Increasing = sum(mCor > 0),
            Decreasing = sum(mCor < 0))

change_count_zf <- tidyr::gather(zf_merged_count, key = "Change", value = "Count", -Chromosome)

change_count_zf$Chromosome <- factor(change_count_zf$Chromosome, levels=c('1', '1A', '2', '3', '4','4A', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28', '29', '30', '31', '32', '34', '35', '37','W', 'Z'))

change_count_zf<- merge(change_count_zf, length_zf, by="Chromosome")

prom_zf<- ggplot(change_count_zf, aes(x = factor(Chromosome, 
                                                 levels = rev(unique(Chromosome[order(Length)]))), y = Count, fill = Change)) +
  geom_bar( stat = "identity") +
  labs(x = "Chromosome",
       y = NULL,
       fill = "Change") +
  scale_fill_manual(values = c("Increasing" = "azure4", "Decreasing" = "black"))+
  ylab("Count of AR-sites on promoters")+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.6))+
  scale_y_continuous(expand = c(0, 0))
prom_zf

sum(subset(change_count_zf,change_count_zf$Change=="Increasing")$Count)
sum(subset(change_count_zf,change_count_zf$Change=="Decreasing")$Count)


# ZF_All_IncDec -----------------------------------------------------------
zf_all_incdec<- merge(zf_loc,info_zf, by="Pos" )

zf_all_incdec_m <- zf_all_incdec%>%
  group_by(Location,type) %>%
  summarise(Increasing = sum(mCor > 0),
            Decreasing = sum(mCor < 0))

zf_all_incdec_m$Location<- factor(zf_all_incdec_m$Location, levels=c("Promoter", "Exon", "Intron", "Intergenic"))

zf_all_incdec_count <- tidyr::gather(zf_all_incdec_m, key = "Change", value = "Count", -Location, -type)

zf_all_incdec_plot<- ggplot(zf_all_incdec_count, aes(x=Location,y=Count,fill = Change)) +
  geom_bar(stat="identity", position = "fill") +
  labs(x = NULL,
       y = NULL,
       fill = "Change") +
  scale_fill_manual(values = c("Increasing" = "azure4", "Decreasing" = "black"))+
  ylab("Proportion of AR-CpG sites")+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.6))+
  scale_y_continuous(expand = c(0, 0))+facet_wrap(~type)+
  theme(strip.text.x = element_text(
    size = 14), 
    strip.background = element_rect(
      color="black", fill="white", size=0.75, linetype="solid"
    ),legend.title=element_text(size=14),legend.text=element_text(size=14))
zf_all_incdec_plot

###### JACKDAW #######

# JD_Age-related CpG sites ------------------------------------------------

#data frame from pipeline run in Habrok (Cor_JD_noMT_Oct23_W.R)
jd25 <- read_csv("JD_AgeRelated_CpG_noMT_Final.csv")

#add info on chromosome length 
length_jd<- read_csv("JD_Chromosome_Length.csv")

length(unique(jd25$Pos))

#separate Pos to get chromosome
jd25 <- jd25 %>%
  separate(Pos, into = c("Chromosome", "Position"), sep = "_", remove = FALSE)

#subset to only keep sites in W
w_jd25<- subset(jd25, jd25$Chromosome=="W")

#how many samples each position is shared by
counts_jd25 <- table(jd25$Pos)
#how many samples each position is shared by (only W)
counts_w_jd25<- table(w_jd25$Pos)

#subset the cor data frame to only keep sites shared by >=15 samples
cor_subset_jd25 <- counts_jd25[counts_jd25 >= 15]

#subset the cor data frame to only keep sites shared by >=6 samples for W
cor_subset_w_jd25<- counts_w_jd25[counts_w_jd25 >= 8]

# View the subsetted data frame
cor_subset_jd25<- as.data.frame(cor_subset_jd25)
cor_subset_w_jd25<- as.data.frame(cor_subset_w_jd25)

#merge the subset for w and rest of chromosomes
final_jd25<- rbind(cor_subset_jd25, cor_subset_w_jd25)

#change colnames
colnames(final_jd25)<- c("Pos", "Freq")

#separate Pos to get chromosome
final_jd25<- final_jd25 %>%
  separate(Pos, into = c("Chromosome", "Position"), sep = "_", remove = FALSE)

#write.csv(final_jd25[,2:3], "JD_Positions_AR_Final.csv")
#how many sitees per Chromosome
ft_jd25<- tabyl(final_jd25$Chromosome)

#change colnames
colnames(ft_jd25)<- c("Chromosome", "n", "percent")

sum(ft_jd25$n)

# JD-AR-CpGs-TotalCpGs+Length---------------------------------------------
#data frame of age-related CpGs per chromosome
pos_jd<- read_csv("JD_PositionsPerChr_All_Final.csv")

ftp_jd25<- merge(ft_jd25, pos_jd, by="Chromosome")

ftpl_jd25<- merge(ftp_jd25, length_jd, by="Chromosome")

ftpl_jd25<- ftpl_jd25[,c(1:2,5:6)]

colnames(ftpl_jd25)<- c("Chromosome", "ARsites", "TotalCpGs", "Length")

ftpl_jd25$Length_kbp<- ftpl_jd25$Length/1000

#write.csv(ftpl_jd25, "FinalPosLength_JD_Cov25.csv")

ftpl_jd25$p<- (0.006*ftpl_jd25$TotalCpGs)/100

jd_all_cov25<- ggplot(ftpl_jd25) +
  geom_point(aes(x = TotalCpGs, y = ARsites, col=Chromosome, size=Length_kbp), fill="black")  +
  geom_text_repel(data = subset(ftpl_jd25, Chromosome %in% c(16,18,35,"Z", "W")),
                  aes(x = TotalCpGs, y = ARsites,label = Chromosome),
                  hjust = -0.5, vjust = 0, size = 6) +
  ylab(NULL) +
  xlab("Total number of CpG sites") +
  scale_color_manual(values = c("W" = "orange", "Z" = "turquoise1", "16"="black", "18"="black","35"="black"))+ 
  labs(size = "Chromosome Length (kbp)")+
  ylim(0,220)+
  guides(color = FALSE)+theme(legend.title = element_text(size=16), legend.text = element_text(size=16), aspect.ratio = 0.5)+
  theme(
    legend.position = c(0.95, 1.2),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.margin = margin(6, 6, 6, 6),
    panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_size_continuous(name = "Chromosome Length (kbp)", labels = scales::comma)+
  scale_x_continuous(labels = scales::comma, limits = c(0,2100000)) +
  geom_line(aes(x=TotalCpGs, y=p),linetype = "dotted", size=0.6)
jd_all_cov25

#find what percentage of all CpGs/chromosome are age-related
ftpl_jd25$percentage<- (ftpl_jd25$ARsites*100)/ftpl_jd25$TotalCpGs

(sum(ftpl_jd25$ARsites)*100)/sum(ftpl_jd25$TotalCpGs)

info_jd <- read_csv("JD_CpG_Info_Final.csv")

jd_cor<- merge(final_jd25, info_jd, by="Pos")

jd_cor$Chromosome <- factor(jd_cor$Chromosome, levels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15','16','17', '18', '19', '20', '21', '22', '24', '25', '26', '27', '28', '29', '30', '31', '33',  '35', '37','39', 'W', 'Z'))


# JD-Increasing/Decreasing DNAm ------------------------------------------
jd_cor_count <- jd_cor %>%
  group_by(Chromosome) %>%
  summarise(Increasing = sum(mCor > 0),
            Decreasing = sum(mCor < 0))

result_long_jd <- tidyr::gather(jd_cor_count, key = "Change", value = "Count", -Chromosome)

result_long_jd$Chromosome <- factor(result_long_jd$Chromosome, levels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15','16','17', '18', '19', '20', '21', '22', '24', '25', '26', '27', '28', '29', '30', '31', '33',  '35', '37','39', 'W', 'Z'))

jd_incdec_stack<- ggplot(result_long_jd, aes(x = Chromosome, y = Count, fill = Change)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Chromosome",
       y = NULL,
       fill = "Change") +
  scale_fill_manual(values = c("Increasing" = "azure4", "Decreasing" = "black"))+
  scale_y_continuous(expand = c(0,0))

result_long_jd %>%
  group_by(Change) %>%
  summarise(Sum_Count = sum(Count))

# JD-Location AR-CpGs -----------------------------------------------------
jd_cor_length<- merge(jd_cor, length_jd, by="Chromosome")
jd_cor_length$Position<- as.numeric(jd_cor_length$Position)
jd_cor_length$Chromosome <- factor(jd_cor_length$Chromosome, levels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15','16','17', '18', '19', '20', '21', '22', '24', '25', '26', '27', '28', '29', '30', '31', '33',  '35', '37','39', 'W', 'Z'))

u_jd<-jd_cor_length %>% distinct(Chromosome, Length)

u_jd$Length_kbp<- u_jd$Length/1000
jd_cor_length$Position_kbp<- jd_cor_length$Position/1000


u_jd_p<-ggplot(u_jd) +
  geom_chicklet(aes(x = Chromosome, y = Length_kbp), col="black", fill="white")  +
  scale_y_continuous(expand = c(0, 0), limits=c(0,127000),labels = scales::comma)+
  geom_errorbar(data=jd_cor_length, aes(x = as.numeric(Chromosome), ymin = Position_kbp, ymax = Position_kbp), color = "navyblue", size = 0.2, width = 0.85) +
  labs(x = "Chromosome", y = "Length (kbp)") 

u_jd_p


# JD_Lenient --------------------------------------------------------------
jd_len <- read_csv("JD_AgeRelated_CpG_noMT_LENIENT_Final.csv")

#separate Pos to get chromosome
jd_len <- jd_len %>%
  separate(Pos, into = c("Chromosome", "Position"), sep = "_", remove = FALSE)

#subset to only keep sites in W
w_jd_len<- subset(jd_len, jd_len$Chromosome=="W")

#how many samples each position is shared by
counts_jd_len <- table(jd_len$Pos)
#how many samples each position is shared by (only W)
counts_w_jd_len<- table(w_jd_len$Pos)

#subset the cor data frame to only keep sites shared by >=15 samples
cor_subset_jd_len <- counts_jd_len[counts_jd_len >= 15]
#subset the cor data frame to only keep sites shared by >=6 samples for W
cor_subset_w_jd_len<- counts_w_jd_len[counts_w_jd_len >= 8]

# View the subsetted data frame
cor_subset_jd_len<- as.data.frame(cor_subset_jd_len)
cor_subset_w_jd_len<- as.data.frame(cor_subset_w_jd_len)
#merge the subset for w and rest of chromosomes
final_jd_len<- rbind(cor_subset_jd_len, cor_subset_w_jd_len)
#change colnames
colnames(final_jd_len)<- c("Pos", "Freq")
#separate Pos to get chromosome
final_jd_len<- final_jd_len %>%
  separate(Pos, into = c("Chromosome", "Position"), sep = "_", remove = FALSE)
#write.csv(final_jd_len[,2:3], "JD_Positions_AR_Final_Len.csv")
#how many sitees per Chromosome
ft_jd_len<- tabyl(final_jd_len$Chromosome)
#change colnames
colnames(ft_jd_len)<- c("Chromosome", "n", "percent")

sum(ft_jd_len$n)

ftp_jd_len<- merge(ft_jd_len, pos_jd, by="Chromosome")


ftp_jd25$Chromosome <- factor(ftp_jd25$Chromosome, levels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15','16','17', '18', '19', '20', '21', '22','23', '24', '25', '26', '27', '28', '29', '30', '31', '33','32', '34', '35','36', '37','39', 'W', 'Z'))

ppc_jd_25<-ggplot(data=ftp_jd25)+
  geom_bar(aes(x=Chromosome, y=n), stat="identity") +
  ylab("Age-related CpG sites")+
  xlab("Chromosome" )
ppc_jd_25

ftp_all<- rbind(ftp_jd_len, ftp_jd25)


ftp_all$Chromosome <- factor(ftp_all$Chromosome, levels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15','16','17', '18', '19', '20', '21', '22','23', '24', '25', '26', '27', '28', '29', '30', '31', '33','32', '34', '35','36', '37','39', 'W', 'Z'))


ftp_jd<- merge(ftp_jd25, ftp_jd_len, by="Chromosome")

ggplot(ftp_jd)+geom_point(aes(x=n.y, y=n.x))

cor.test(ftp_jd$n.x,ftp_jd$n.y)

# JD_Genomic_Annotation ------------------------------------------------------
jd_anot<- readRDS("JD_Annotation_All.rds")

jd_loc <- read_excel("C:/Users/MWP-/Dropbox/ZF_JD_AgeRelatedSites/Stats/JD_CpG_Annotation.xlsx")

jd_loc$Location<- factor(jd_loc$Location, levels=c("Promoter", "Exon", "Intron", "Intergenic"))

jd_loc$ObsExp<- factor(jd_loc$ObsExp, levels=c("Observed", "Expected"))

loc_jd<- ggplot(jd_loc, aes(x=ObsExp, y=Percentage,fill=Location))+
  geom_bar( position = "fill", stat="identity", width=0.5)+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_brewer(palette="RdBu")+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.6))+
  xlab(NULL)+  ylab(NULL)

loc_jd

loc_jd1<- ggplot(data=subset(jd_loc, jd_loc$ObsExp=="Observed"))+
  geom_bar(aes(y=(Percentage/100),x=Location, fill=Location), stat="identity", width=0.75)+
  scale_y_continuous(expand = c(0, 0), limits=c(0,0.6))+
  scale_fill_brewer(palette="RdBu")+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.6))+
  xlab(NULL)+  ylab("Proportion of AR-CpG sites")

loc_jd1

jd_aut <- read_excel("C:/Users/MWP-/Dropbox/ZF_JD_AgeRelatedSites/Stats/JD_CpG_Annotation.xlsx",  sheet = "Sheet2")

jd_aut$Location<- factor(jd_aut$Location, levels=c("Promoter", "Exon", "Intron", "Intergenic"))

aut_jd<- ggplot(jd_aut, aes(x=Chr, y=Percentage,fill=Location))+
  geom_bar( position = "fill", stat="identity", width=0.65)+
  scale_y_continuous(expand = c(0, 0))+
  scale_fill_brewer(palette="RdBu")+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.6),
         legend.title = element_text(size = 16),
         legend.text =  element_text(size = 16))+
  xlab(NULL)+  ylab(NULL)
aut_jd


# JD_Promoters ------------------------------------------------------------
#subset sites in promoters
jd_anot_prom<- subset(jd_anot, jd_anot$prom=="1")

#make variable Pos
jd_anot_prom <- jd_anot_prom %>%
  mutate(Pos = paste(chromossome, site, sep = "_"))

#merge final CpGs with information on increasing/decreasing 
jd_mrged<- merge(final_jd25, info_jd, by="Pos")
#merge that with CpGs in promoters
jd_merged<- merge(jd_mrged, jd_anot_prom, by="Pos")

#write.csv(jd_merged, "JD_ARCpGs_Promoters.csv")

jd_merged<- read.csv("C:/Users/MWP-/Dropbox/ZF_JD_AgeRelatedSites/Scripts/JD_ARCpGs_Promoters.csv")

jd_merged_count<- jd_merged %>%
  group_by(Chromosome) %>%
  summarise(Increasing = sum(mCor > 0),
            Decreasing = sum(mCor < 0))

change_count_jd <- tidyr::gather(jd_merged_count, key = "Change", value = "Count", -Chromosome)

change_count_jd$Chromosome<- factor(change_count_jd$Chromosome, levels=c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15','16','17', '18', '19', '20', '21', '22','23', '24', '25', '26', '27', '28', '29', '30', '31', '33','32', '34', '35','36', '37','39', 'W', 'Z'))

change_count_jd<- merge(change_count_jd, length_jd, by="Chromosome")

ggplot(change_count_jd, aes(x = Chromosome, y = Count, fill = Change)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Chromosome",
       y = NULL,
       fill = "Change") +
  scale_fill_manual(values = c("Increasing" = "azure4", "Decreasing" = "black"))+
  ylab("Proportion of AR-sites on promoters")

prom_jd<- ggplot(change_count_jd, aes(x = factor(Chromosome, 
                                                 levels = rev(unique(Chromosome[order(Length)]))), y = Count, fill = Change)) +
  geom_bar( stat = "identity") +
  labs(x = "Chromosome",
       y = NULL,
       fill = "Change") +
  scale_fill_manual(values = c("Increasing" = "azure4", "Decreasing" = "black"))+
  ylab("Count of AR-sites on promoters")+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.6))+
  scale_y_continuous(expand = c(0, 0))
prom_jd

sum(change_count_jd$Count)
sum(subset(change_count_jd,change_count_jd$Change=="Increasing")$Count)
sum(subset(change_count_jd,change_count_jd$Change=="Decreasing")$Count)

# JD_All_IncDec -----------------------------------------------------------
jd_loc_all<- merge(jd_anot, final_jd25, by="Pos")

jd_all_incdec<- merge(jd_loc_all,info_jd, by="Pos" )

jd_all_incdec$type<- NA
jd_all_incdec$type[jd_all_incdec$Chromosome=="W"]<- "W"
jd_all_incdec$type[jd_all_incdec$Chromosome=="Z"]<- "Z"
jd_all_incdec$type[is.na(jd_all_incdec$type)] = "Autosomes"

jd_all_incdec[jd_all_incdec=="prom"]<- "Promoter"
jd_all_incdec[jd_all_incdec=="exon"]<- "Exon"
jd_all_incdec[jd_all_incdec=="intron"]<- "Intron"
jd_all_incdec[jd_all_incdec=="intergenic"]<- "Intergenic"

jd_all_incdec_m <- jd_all_incdec%>%
  group_by(category,type) %>%
  summarise(Increasing = sum(mCor > 0),
            Decreasing = sum(mCor < 0))

jd_all_incdec_m$category<- factor(jd_all_incdec_m$category, levels=c("Promoter", "Exon", "Intron", "Intergenic"))

jd_all_incdec_count <- tidyr::gather(jd_all_incdec_m, key = "Change", value = "Count", -category, -type)

jd_all_incdec_plot<-ggplot(jd_all_incdec_count, aes(x=category,y=Count,fill = Change)) +
  geom_bar(stat="identity", position = "fill") +
  labs(x = NULL,
       y = NULL,
       fill = "Change") +
  scale_fill_manual(values = c("Increasing" = "azure4", "Decreasing" = "black"))+
  ylab("Proportion of AR-CpG sites")+
  theme( panel.border = element_rect(colour = "black", fill=NA, size=0.6))+
  scale_y_continuous(expand = c(0, 0))+facet_wrap(~type)+
  theme(strip.text.x = element_text(
    size = 14), 
    strip.background = element_rect(
      color="black", fill="white", size=0.75, linetype="solid"
    ),legend.title=element_text(size=14),legend.text=element_text(size=14))
jd_all_incdec_plot


# STATS -------------------------------------------------------------------
# ZF_PER_CHROMOSOME_STATS -------------------------------------------------
zf_mod<- ftpl[,c(1:2,4:5)]

colnames(zf_mod)<- c("Chromosome", "ARsites", "TotalCpGs","Length")
# Calculate the expected number of age-related CpGs
# If 0.006% of all CpGs in a chromosome were age-related
zf_mod<- zf_mod %>%
  mutate(
    prop_exp=(TotalCpGs*0.006)/sum(TotalCpGs*0.006),
    prop_obs=ARsites/sum(ARsites)
  )

# Make sure it adds up to 1
sum(zf_mod$prop_exp)
sum(zf_mod$prop_obs)

with(zf_mod, chisq.test(x = ARsites, p = prop_exp))# Permutation version
with(zf_mod, chisq.test(x = ARsites, p = prop_exp, simulate.p.value = TRUE))

N <- sum(zf_mod$ARsites)

# For each category...
with(zf_mod, map2(ARsites, prop_exp, function(ARsites, prop_exp) {
  
  # Perform a binomial test: is the no. of CpG in that category different than expected?
  binom.test(ARsites, N, p = prop_exp, alternative = "two.sided")
  
  # Note: the test is two-sided to identify both under and over-representation
  
}))

# We can do the same and attach the results to our original data
zf_mod <- zf_mod %>%
  mutate(
    p_value = map2_dbl(
      ARsites, prop_exp,
      ~ binom.test(.x, N, .y, alternative = "two.sided")$p.value
    )
  )

# Correction for multiple testing
zf_mod <- zf_mod %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

zf_mod$sig <- ifelse(zf_mod$p_adj <= 0.05, "Yes", "No")

#write.csv(zf_mod, "ZF_ARsites_PerChr_Stats.csv")

tabyl(zf_mod$sig)

#check correlation prop_obs-prop_exp vs. chromosome length
zf_mod<- zf_mod %>%
  mutate(
    abs_dif=abs(prop_obs-prop_exp),
    dif=prop_obs-prop_exp
  )

#count how many chromosomes with p<0.05 were over and under-represented
filtered_zf_mod<- zf_mod %>% 
  filter(p_adj < 0.05)

# Count positive and negative dif values
positive_count <- sum(filtered_zf_mod$dif > 0)
negative_count <- sum(filtered_zf_mod$dif < 0)


# ZF_Observed vs. Expected STATS ------------------------------------------
#all CpG sites detected is 18906252
#subset expected
zf_loc_exp<- subset(zf_loc_all, zf_loc_all$ObsExp=="Expected")
#find counts expected
zf_loc_exp$Expected_Count <- round(zf_loc_exp$Proportion *18906252)
#make sure they sum to 18906252
sum(zf_loc_exp$Expected_Count)

#subset observed
zf_loc_obs<- subset( zf_loc_all,  zf_loc_all$ObsExp=="Observed")

# Calculate the count of observed and expected based on the total count of 1120
zf_loc_obs$Observed_Count <- round(zf_loc_obs$Proportion * 1085)

sum(zf_loc_obs$Observed_Count)

# Calculate the count of observed and expected based on the total count of 1120
#zf_loc_exp$Expected_Count <- round(zf_loc_exp$Proportion * 1121)

sum(zf_loc_exp$Expected_Count)

zf_loc_obs<- zf_loc_obs[, c(1:2,4)]
zf_loc_exp<- zf_loc_exp[,c(2,4)]

zf_loc_count<- cbind(zf_loc_obs, zf_loc_exp)
colnames(zf_loc_count)<- c("Location", "prop_obs", "count_obs", "prop_exp", "count_exp")

sum(zf_loc_count$prop_obs)
sum(zf_loc_count$prop_exp)

chisq.test(x=zf_loc_count$count_obs, p=zf_loc_count$prop_exp)

with(zf_loc_count, map2(count_obs, prop_exp, function(count_obs, prop_exp) {
  
  # Perform a binomial test
  binom.test(count_obs, sum(zf_loc_count$count_obs), p = prop_exp, alternative = "two.sided")
  
  # Note: the test is two-sided to identify both under and over-representation
  
}))



# JD_PER_CHROMOSOME_STATS -------------------------------------------------
jd_mod<- ftpl_jd25[,c(1:4)]

jd_mod<- jd_mod %>%
  mutate(
    prop_exp=(TotalCpGs*0.006)/sum(TotalCpGs*0.006),
    prop_obs=ARsites/sum(ARsites)
  )

# Make sure it adds up to 1
sum(jd_mod$prop_exp)
sum(jd_mod$prop_obs)

with(jd_mod, chisq.test(x = ARsites, p = prop_exp))

# Permutation version
with(jd_mod, chisq.test(x = ARsites, p = prop_exp, simulate.p.value = TRUE))

N <- sum(jd_mod$ARsites)

# For each category...
with(jd_mod, map2(ARsites, prop_exp, function(ARsites, prop_exp) {
  
  # Perform a binomial test: is the no. of CpG in that category different than expected?
  binom.test(ARsites, N, p = prop_exp, alternative = "two.sided")
  
  # Note: the test is two-sided to identify both under and over-representation
  
}))

# We can do the same and attach the results to our original data
jd_mod <- jd_mod %>%
  mutate(
    p_value = map2_dbl(
      ARsites, prop_exp,
      ~ binom.test(.x, N, .y, alternative = "two.sided")$p.value
    )
  )

# Correction for multiple testing
jd_mod <- jd_mod %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

jd_mod$sig <- ifelse(jd_mod$p_adj <= 0.05, "Yes", "No")

#write.csv(jd_mod, "JD_ARsites_PerChr_Stats.csv")

tabyl(jd_mod$sig)
#check correlation prop_obs-prop_exp vs. chromosome length
jd_mod<- jd_mod %>%
  mutate(
    abs_dif=abs(prop_obs-prop_exp),
    dif=prop_obs-prop_exp
  )

#count how many chromosomes with p<0.05 were over and under-represented
filtered_jd_mod<- jd_mod %>% 
  filter(p_adj < 0.05)

# Count positive and negative dif values
positive_count <- sum(filtered_jd_mod$dif > 0)
negative_count <- sum(filtered_jd_mod$dif < 0)


# JD_Observed vs. Expected STATS ------------------------------------------
###EXPECTED###
#subset expected
jd_loc_exp<- subset(jd_loc, jd_loc$ObsExp=="Expected")
#find counts expected
jd_loc_exp$Expected_Count <- round((jd_loc_exp$Percentage)/100 *21835349)
#make sure they sum to 18901869
sum(jd_loc_exp$Expected_Count)
#find proportion exp
jd_loc_exp$prop_exp <- jd_loc_exp$Percentage / 100

###OBSERVED###
#subset observed
jd_loc_obs<- subset(jd_loc, jd_loc$ObsExp=="Observed")
#proportion observed
jd_loc_obs$Proportion <- jd_loc_obs$Percentage / 100

# Calculate the count of observed and expected based on the total count of 1120
jd_loc_obs$count_obs<- round(jd_loc_obs$Proportion * 1299)

# Calculate the count of observed and expected based on the total count of 1120
#subset expected


sum(jd_loc_exp$prop_exp)

jd_loc_obs<- jd_loc_obs[, c(1,4:5)]
jd_loc_exp<- jd_loc_exp[,4:5]

jd_loc_count<- cbind(jd_loc_obs, jd_loc_exp)
colnames(jd_loc_count)<- c("Location", "prop_obs", "count_obs", "count_exp", "prop_exp")

sum(jd_loc_count$prop_obs)
sum(jd_loc_count$prop_exp)

chisq.test(x=jd_loc_count$count_obs, p=jd_loc_count$prop_exp)

with(jd_loc_count, map2(count_obs, prop_exp, function(count_obs, prop_exp) {
  
  # Perform a binomial test
  binom.test(count_obs, sum(jd_loc_count$count_obs), p = prop_exp, alternative = "two.sided")
  
  # Note: the test is two-sided to identify both under and over-representation
  
}))

# FIGURES ------------------------------------------------
#merge graphs together

setwd("C:/Users/MWP-/Dropbox/ZF_JD_AgeRelatedSites/Figures")

###FIGURE 2###
#Increasing/decreasing age-relatedCpGs
tiff("Fig2_Prop_AR-CpGsites_PerChromosome.tiff", units="in", width=10, height=7, res=300)
incdec<- ggarrange(zf_incdec_stack,ggplot() + theme_void(),jd_incdec_stack, nrow=3, common.legend = T, labels = c("(a)","","(b)"), label.x = c(0.03,0.03, 0.03), label.y = c(1.1,0.92, 1.1), heights = c(1.3,0.4,1.4))
incdec <- incdec + theme(plot.margin = margin(l = 20, t=30))
annotate_figure(incdec,left = textGrob("Proportion of age-related CpG sites", rot = 90, vjust = 1, gp = gpar(cex = 1.6)))  
dev.off()


###FIGURE 3###
#Age-related CpG location
tiff("Fig3_AR-CpGsites_PerChromosome.jpeg.tiff", units="in", width=10, height=8, res=300)
ggarrange(ggplot() + theme_void(),u_zf,  ggplot() + theme_void(), u_jd_p, nrow=4,labels = c("", "(a)","", "(b)"), common.legend = T, label.y = c(1,1.1,1,1.15), heights=c(0.1,1,0.1, 1))
dev.off()

###FIGURE 4###
#Age-related CpGs/Total CpGs
svg("Fig3_AR-CpGsites_AllCpGSites.svg", width=10, height=8)
art<- ggarrange(s_all, jd_all_cov25, 
                nrow=2, 
                common.legend = T, 
                align="v", 
                labels = c("A", "B"),
                widths=c(1,1),
                label.x = c(0.1,0.1), 
                label.y = c(1,1), 
                heights = c(1,1.07), 
                legend = "top")
annotate_figure(art,left = textGrob("Number of age-related CpG sites", rot = 90, vjust = 6, gp = gpar(cex = 1.3)))
dev.off()

###FIGURE 5###
p1<- ggarrange(loc_zf1, aut_zf, loc_zf, nrow = 1, ncol = 3, align = "h", common.legend = T,  legend = "right",labels=c("A", "B", "C"), label.y = c(1.1,1.1,1.1),label.x = c(0.1,0.1,0.1))+
  theme(plot.margin = margin(1,0.1,0.1,0.1, "cm"))
p1
p2<- ggarrange(prom_zf, nrow = 1, align = "h", labels="D", label.x=0.025, label.y=1.1)+
  theme(plot.margin = margin(0.6,0.1,0.1,0.1, "cm"))
p2
p3<- ggarrange(loc_jd1, aut_jd, loc_jd, nrow = 1, ncol = 3, align = "h", common.legend = T, legend = "right",labels=c("E", "F", "G"), label.y = c(1.1,1.1,1.1),label.x = c(0.1,0.1,0.1))+
  theme(plot.margin = margin(0.6,0.1,0.1,0.1, "cm"))
p3
p4<-ggarrange(prom_jd, nrow = 1, align = "h", labels="H", label.x=0.025, label.y=1.1)+
  theme(plot.margin = margin(0.6,0.1,0.1,0.1, "cm"))
p4
svg("Fig4.svg", width=10, height=14)
ggarrange(p1,p2,p3,p4, nrow=4, align="v")
dev.off()


##SUPPLEMETARY INFORMATION##
#FIGURE S2#
plot_len_zf <- ggplot(ftp_len, aes(x = Chromosome, y = n)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  ylab("Number of age-related CpG sites") +
  xlab(NULL) +  
  ylim(0,900)+
  theme_minimal()+
  geom_text(label="Step 2: Top 8% correlation, Coverage>=17 \n Step 3:Average DNAm >=30% & <=70%", x=16, y=600, size=6)
plot_len_zf

plot_zf <- ggplot(ftp, aes(x = Chromosome, y = n)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  ylab("Number of age-related CpG sites") +
  xlab("Chromosome")  +  
  ylim(0,230)+
  theme_minimal() +  
  geom_text(label="Step 2: Top 5% correlation, Coverage>=20 \n Step 3:Average DNAm >=40% & <=60%", x=16, y=150, size=6)
plot_zf

tiff("FigS2_ZF_Len.tiff", units="in", width=10, height=8, res=300)
ggarrange( plot_zf, plot_len_zf,nrow=2, labels=c("(a)","(b)"), label.x = c(0.02,0.02))
dev.off()

# Plotting for ftp_jd_len
plot_len <- ggplot(ftp_jd_len, aes(x = Chromosome, y = n)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  ylab("Number of age-related CpG sites") +
  xlab("Chromosome") +  
  ylim(0,900)+
  geom_text(label="Step 2: Top 8% correlation, Coverage>=20 \n Step 3:Average DNAm >=30% & <=70%", x=16, y=600, size=6)+
  theme_minimal()
plot_len

# Plotting for ftp_jd25
plot_jd25 <- ggplot(ftp_jd25, aes(x = Chromosome, y = n)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +
  ylab("Number of age-related CpG sites") +
  xlab(NULL)  +  
  geom_text(label="Step 2: Top 5% correlation, Coverage>=25 \n Step 3:Average DNAm >=40% & <=60%", x=16, y=150, size=6)+
  theme_minimal()
plot_jd25

tiff("FigS3_JD_Len.tiff", units="in", width=10, height=8, res=300)
ggarrange(plot_jd25,plot_len,  nrow=2, labels=c("(a)","(b)"), label.x = c(0.02,0.02))
dev.off()

svg("FigS4_IncDec_ZF_JD_Len.svg", , width=10, height=8)
ggarrange(zf_all_incdec_plot, jd_all_incdec_plot, nrow=2, labels=c("A","B"), common.legend = T)
dev.off()

# Figures_Merged_(Science) ------------------------------------------------

u_zf<- ggplot(u) +
  geom_chicklet(aes(x = Chromosome, y = Length_kbp), col="black", fill="white")  +
  scale_y_continuous(expand = c(0, 0),limits=c(0,160000),labels = scales::comma)+
  geom_errorbar(data=zf_cor_length, aes(x = as.numeric(Chromosome), ymin = Position_kbp, ymax = Position_kbp), color = "navyblue", size = 0.2, width = 0.85) +
  labs(x = NULL, y ="Length (kbp)") 

u_zf

zf_incdec_stack<- ggplot(result_long, aes(x = Chromosome, y = Count, fill = Change)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Chromosome",
       y = "Proportion of AR-CpG sites",
       fill = "Change") +
  scale_fill_manual(values = c("Increasing" = "azure4", "Decreasing" = "black"))+
  scale_y_continuous(expand = c(0,0))
zf_incdec_stack



u_jd_p<-ggplot(u_jd) +
  geom_chicklet(aes(x = Chromosome, y = Length_kbp), col="black", fill="white")  +
  scale_y_continuous(expand = c(0, 0), limits=c(0,127000),labels = scales::comma)+
  geom_errorbar(data=jd_cor_length, aes(x = as.numeric(Chromosome), ymin = Position_kbp, ymax = Position_kbp), color = "navyblue", size = 0.2, width = 0.85) +
  labs(x = NULL, y = NULL) 

u_jd_p

jd_incdec_stack<- ggplot(result_long_jd, aes(x = Chromosome, y = Count, fill = Change)) +
  geom_bar(position = "fill", stat = "identity") +
  labs(x = "Chromosome",
       y = NULL,
       fill = "Change") +
  scale_fill_manual(values = c("Increasing" = "azure4", "Decreasing" = "black"))+
  scale_y_continuous(expand = c(0,0))

svg(file = "Fig2_Merged.svg", width = 18, height = 10)
ggarrange(u_zf, u_jd_p, zf_incdec_stack, jd_incdec_stack, nrow=2, ncol=2, common.legend = T, align="v", legend="bottom", labels=c("A", "B", "C", "D"), label.x = c(0.001,0.001,0.01,0.01))
dev.off()



