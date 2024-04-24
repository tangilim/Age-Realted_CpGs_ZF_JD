#title: "Identification of age-related CpG sites in two avian species"
#author: "Marianthi Tangili & Joanna Sudyka"
#date: "24/04/2024"

#this code was run using the University of Groningen hpc 'Habrok'

library(readr)
library(data.table)
library(dplyr)
library(janitor)
library(tidyr)
library(ggplot2)


# merge .cov files including DNAm of all samples --------------------------

#ZEBRA FINCH
cov_data_list <- list()

# Loop through files, read data, add Sample column, and set column names
for (i in 1:20) {
  file_path <- paste0(directory_path, "S", i, ".deduplicated.bismark.cov.gz")
  cov_data <- read.delim(file_path)
  cov_data$Sample <- as.character(i)
  
  # Set column names
  colnames(cov_data) <- c("Chromosome", "StartPosition", "EndPosition", "MethylationPercentage", "CountMethylated", "CountNonMethylated", "Sample")
  
  cov_data_list[[paste0("Cov", i)]] <- cov_data
}

# Combine the data frames into a single data frame
CovN <- do.call(rbind, cov_data_list)

#Replace Chromosome names
CovN[CovN == "CM012081.2"] <- "1"
CovN[CovN == "CM012082.2"] <- "1A"
CovN[CovN == "CM012083.2"] <- "2"
CovN[CovN == "CM012084.2"] <- "3"
CovN[CovN == "CM012085.2"] <- "4"
CovN[CovN == "CM012086.2"] <- "4A"
CovN[CovN == "CM012087.2"] <- "5"
CovN[CovN == "CM012088.2"] <- "6"
CovN[CovN == "CM012089.2"] <- "7"
CovN[CovN == "CM012090.2"] <- "8"
CovN[CovN == "CM012091.2"] <- "9"
CovN[CovN == "CM012092.2"] <- "10"
CovN[CovN == "CM012093.2"] <- "11"
CovN[CovN == "CM012094.2"] <- "12"
CovN[CovN == "CM012095.2"] <- "13"
CovN[CovN == "CM012096.2"] <- "14"
CovN[CovN == "CM012097.2"] <- "15"
CovN[CovN == "CM012098.2"] <- "16"
CovN[CovN == "CM012099.2"] <- "17"
CovN[CovN == "CM012100.2"] <- "18"
CovN[CovN == "CM012101.2"] <- "19"
CovN[CovN == "CM012102.2"] <- "20"
CovN[CovN == "CM012103.2"] <- "21"
CovN[CovN == "CM012104.2"] <- "22"
CovN[CovN == "CM012105.2"] <- "23"
CovN[CovN == "CM012106.2"] <- "24"
CovN[CovN == "CM012107.2"] <- "25"
CovN[CovN == "CM012108.2"] <- "26"
CovN[CovN == "CM012110.2"] <- "27"
CovN[CovN == "CM012112.2"] <- "28"
CovN[CovN == "CM012109.2"] <- "29"
CovN[CovN == "CM012111.2"] <- "30"
CovN[CovN == "CM025860.2"] <- "31"
CovN[CovN == "CM025861.2"] <- "32"
CovN[CovN == "CM025862.2"] <- "33"
CovN[CovN == "CM025863.2"] <- "34"
CovN[CovN == "CM025864.2"] <- "35"
CovN[CovN == "CM025865.2"] <- "36"
CovN[CovN == "CM030979.1"] <- "37"
CovN[CovN == "CM012113.2"] <- "Z"
CovN[CovN == "CM018260.2"] <- "W"
CovN[CovN == "NC_007897.1"] <- "MT"

CovN<- subset(CovN, CovN$Chromosome!="MT")

CovN$Coverage<- CovN$CountMethylated+CovN$CountNonMethylated

saveRDS(CovN, file="CovN_ZF_noMT.rds")

names<- unique(CovN$Chromosome)

write.csv(names, "ChrNamesZF_noMT.csv")

CovN$SampleL<- NA
CovN$SampleL[CovN$Sample=="1" | CovN$Sample=="18"]<- 1
CovN$SampleL[CovN$Sample=="2" | CovN$Sample=="11"]<- 2
CovN$SampleL[CovN$Sample=="3" | CovN$Sample=="20"]<- 3
CovN$SampleL[CovN$Sample=="4" | CovN$Sample=="13"]<- 4
CovN$SampleL[CovN$Sample=="5" | CovN$Sample=="14"]<- 5
CovN$SampleL[CovN$Sample=="6" | CovN$Sample=="10"]<- 6
CovN$SampleL[CovN$Sample=="15" | CovN$Sample=="7"]<- 7
CovN$SampleL[CovN$Sample=="16" | CovN$Sample=="8"]<- 8
CovN$SampleL[CovN$Sample=="17" | CovN$Sample=="9"]<- 9
CovN$SampleL[CovN$Sample=="19" | CovN$Sample=="12"]<- 10

Dage <- read_csv("/scratch/p301533/Pipeline/ZF_MT/Dage_ZF.csv")

df<- merge(CovN, Dage, by=c("Sample", "SampleL"))

df$Pos<-paste(df$Chromosome,df$StartPosition,sep="_")

saveRDS(df, file = "MethDage_ZF_noMT.rds")

#JACKDAW
# Loop through files, read data, add Sample column, and set column names
for (i in 1:22) {
  file_path <- paste0(directory_path, "S", i, ".deduplicated.bismark.cov.gz")
  cov_data <- read.delim(file_path)
  cov_data$Sample <- as.character(i)
  
  # Set column names
  colnames(cov_data) <- c("Chromosome", "StartPosition", "EndPosition", "MethylationPercentage", "CountMethylated", "CountNonMethylated", "Sample")
  
  cov_data_list[[paste0("Cov", i)]] <- cov_data
}

# Combine the data frames into a single data frame
CovN <- do.call(rbind, cov_data_list)

#replace chromosome names
CovN[CovN == "NC_063213.1"] <- "1"
CovN[CovN == "NC_063214.1"] <- "2"
CovN[CovN == "NC_063215.1"] <- "3"
CovN[CovN == "NC_063216.1"] <- "4"
CovN[CovN == "NC_063217.1"] <- "5"
CovN[CovN == "NC_063218.1"] <- "6"
CovN[CovN == "NC_063219.1"] <- "7"
CovN[CovN == "NC_063220.1"] <- "8"
CovN[CovN == "NC_063221.1"] <- "9"
CovN[CovN == "NC_063222.1"] <- "10"
CovN[CovN == "NC_063223.1"] <- "11"
CovN[CovN == "NC_063224.1"] <- "12"
CovN[CovN == "NC_063225.1"] <- "13"
CovN[CovN == "NC_063226.1"] <- "14"
CovN[CovN == "NC_063227.1"] <- "15"
CovN[CovN == "NC_063228.1"] <- "16"
CovN[CovN == "NC_063229.1"] <- "17"
CovN[CovN == "NC_063230.1"] <- "18"
CovN[CovN == "NC_063231.1"] <- "19"
CovN[CovN == "NC_063232.1"] <- "20"
CovN[CovN == "NC_063233.1"] <- "21"
CovN[CovN == "NC_063234.1"] <- "22"
CovN[CovN == "NC_063235.1"] <- "23"
CovN[CovN == "NC_063236.1"] <- "24"
CovN[CovN == "NC_063237.1"] <- "25"
CovN[CovN == "NC_063238.1"] <- "26"
CovN[CovN == "NC_063239.1"] <- "27"
CovN[CovN == "NC_063240.1"] <- "28"
CovN[CovN == "NC_063241.1"] <- "29"
CovN[CovN == "NC_063242.1"] <- "30"
CovN[CovN == "NC_063243.1"] <- "31"
CovN[CovN == "NC_063244.1"] <- "32"
CovN[CovN == "NC_063245.1"] <- "33"
CovN[CovN == "NC_063246.1"] <- "34"
CovN[CovN == "NC_063247.1"] <- "35"
CovN[CovN == "NC_063248.1"] <- "36"
CovN[CovN == "NC_063249.1"] <- "37"
CovN[CovN == "NC_063250.1"] <- "38"
CovN[CovN == "NC_063251.1"] <- "39"
CovN[CovN == "NC_063252.1"] <- "40"
CovN[CovN == "NC_063253.1"] <- "41"
CovN[CovN == "NC_063255.1"] <- "Z"
CovN[CovN == "NC_063254.1"] <- "W"
CovN[CovN == "NC_026783.1"] <- "MT"

#save data frame
saveRDS(CovN, file="CovJD.rds")

#Save data frame
saveRDS(CovN, file="CovJD_HC__Oct23_onlyMT.rds")


CovN$SampleL<- NA
CovN$SampleL[CovN$Sample=="1" | CovN$Sample=="10"]<- 1
CovN$SampleL[CovN$Sample=="2" | CovN$Sample=="11"]<- 2
CovN$SampleL[CovN$Sample=="3" | CovN$Sample=="12"]<- 3
CovN$SampleL[CovN$Sample=="4" | CovN$Sample=="13"]<- 4
CovN$SampleL[CovN$Sample=="5" | CovN$Sample=="14"]<- 5
CovN$SampleL[CovN$Sample=="6" | CovN$Sample=="15"]<- 6
CovN$SampleL[CovN$Sample=="7" | CovN$Sample=="16"]<- 7
CovN$SampleL[CovN$Sample=="8" | CovN$Sample=="17"]<- 8
CovN$SampleL[CovN$Sample=="9" | CovN$Sample=="18"]<- 9
CovN$SampleL[CovN$Sample=="19" | CovN$Sample=="21"]<- 10
CovN$SampleL[CovN$Sample=="20" | CovN$Sample=="22"]<- 10

Dage <- read_csv("/scratch/p301533/Pipeline/JD_MT/DAge_JD.csv")

df<- merge(CovN, Dage, by=c("Sample", "SampleL"))

df$Pos<-paste(df$Chromosome,df$StartPosition,sep="_")

saveRDS(df, file = "MethDageJD__Oct23_onlyMT.rds")


# pipeline steps ----------------------------------------------------------
# ZEBRA FINCH   
#insert data frame with Dage information
d<- readRDS("/scratch/p301533/Pipeline/ZF_MT/MethDage_ZF_noMT.rds") 

#find avg coverage per Pos
avgcov<- d%>%
  group_by(Pos) %>%
  summarise(MeanC=mean(Coverage),SDC=sd(Coverage))

#make data table
dt <- data.table(d)

#find correlations of % DNAm and Dage per position
dtCor <- dt[, .(mCor = cor(MethylationPercentage,Dage)), by=Pos]

cor<- as.data.frame(dtCor)

#merge average coverage and correlation information per position
ab1<- merge(avgcov,cor, by="Pos")

#make correlation absolute
ab1$abs<- abs(ab1$mCor)

#find absolute correlation per 20 coverage
b1<- subset(ab1, ab1$MeanC>=0 & ab1$MeanC<=20)
b1s<-top_frac(b1, 0.08)

b2<- subset(ab1, ab1$MeanC>21 & ab1$MeanC<=40)
b2s<-top_frac(b2, 0.08)

b3<- subset(ab1, ab1$MeanC>41 & ab1$MeanC<=60)
b3s<-top_frac(b3, 0.08)

b4<- subset(ab1, ab1$MeanC>61 & ab1$MeanC<=80)
b4s<-top_frac(b4, 0.08)

b5<- subset(ab1, ab1$MeanC>81 & ab1$MeanC<=100)
b5s<-top_frac(b5, 0.08)

b6<- subset(ab1, ab1$MeanC>101 & ab1$MeanC<=120)
b6s<-top_frac(b6, 0.08)

b7<- subset(ab1, ab1$MeanC>121 & ab1$MeanC<=140)
b7s<-top_frac(b7, 0.08)

b8<- subset(ab1, ab1$MeanC>141 & ab1$MeanC<=160)
b8s<-top_frac(b8, 0.08)

b9<- subset(ab1, ab1$MeanC>161 & ab1$MeanC<=180)
b9s<-top_frac(b9, 0.08)

b10<- subset(ab1, ab1$MeanC>181 & ab1$MeanC<=200)
b10s<-top_frac(b10, 0.08)

b11<- subset(ab1, ab1$MeanC>201 & ab1$MeanC<=220)
b11s<-top_frac(b11, 0.08)

b12<- subset(ab1, ab1$MeanC>221 & ab1$MeanC<=240)
b12s<-top_frac(b12, 0.08)

b13<- subset(ab1, ab1$MeanC>241 & ab1$MeanC<=260)
b13s<-top_frac(b13, 0.08)

b14<- subset(ab1, ab1$MeanC>261 & ab1$MeanC<=280)
b14s<-top_frac(b14, 0.08)

b15<- subset(ab1, ab1$MeanC>281 & ab1$MeanC<=300)
b15s<-top_frac(b15, 0.08)

t<- rbind(b1s, b2s, b3s, b4s, b5s, b6s, b7s, b8s, b9s, b10s, b11s, b12s, b13s, b14s, b15s)

#subset for mean coverage >=20
t1<- subset(t, t$MeanC>=18)

#find mean %DNAm per position
mean<-aggregate(x=d$MethylationPercentage,
                by=list(d$Pos), FUN=mean)

colnames(mean)<- c("Pos", "Mean")

#only keep top 5% correlation positions
t1$MeanM <- mean$Mean[match(t1$Pos, mean$Pos)]

#subset for mean %DNAm

t2<- subset(t1, t1$MeanM>=35 & t1$MeanM<=65)

t2<-separate(t2, col=Pos, into=c('Chromosome', 'Position'), sep='_')

j3<- tabyl(t2, Chromosome)

#make position variable from chromosome and startposition
t2$Pos<-paste(t2$Chromosome,t2$Position,sep="_")

#merge information from initial dataframe and subsetted pos
merged_df <- merge(d, t2, by ="Pos")

# Keep only the rows where Pos is present in t2
subset_df <- subset(merged_df, Pos %in% t2$Pos)

# find by how many samples each position is shared
s <- aggregate(Sample ~ Pos, subset_df, function(x) length(unique(x)))

merged_df <- merge(s, d, by = "Pos")

#Keep only the rows where Pos is present in t2
c <- subset(merged_df, Pos %in% t2$Pos)

#JACKDAW 
d<- readRDS("/scratch/p301533/Pipeline/JD_MT/MethDageJD__Oct23_noMT.rds") 

#find avg coverage per Pos
avgcov<- d%>%
  group_by(Pos) %>%
  summarise(MeanC=mean(Coverage),SDC=sd(Coverage))

#make data table
dt <- data.table(d)

#find correlations of % DNAm and Dage per position
dtCor <- dt[, .(mCor = cor(MethylationPercentage,Dage)), by=Pos]

cor<- as.data.frame(dtCor)

#merge average coverage and correlation information per position
ab1<- merge(avgcov,cor, by="Pos")

ab1$abs<- abs(ab1$mCor)

b1<- subset(ab1, ab1$MeanC>=0 & ab1$MeanC<=20)
b1s<-top_frac(b1, 0.05)

b2<- subset(ab1, ab1$MeanC>21 & ab1$MeanC<=40)
b2s<-top_frac(b2, 0.05)

b3<- subset(ab1, ab1$MeanC>41 & ab1$MeanC<=60)
b3s<-top_frac(b3, 0.05)

b4<- subset(ab1, ab1$MeanC>61 & ab1$MeanC<=80)
b4s<-top_frac(b4, 0.05)

b5<- subset(ab1, ab1$MeanC>81 & ab1$MeanC<=100)
b5s<-top_frac(b5, 0.05)

b6<- subset(ab1, ab1$MeanC>101 & ab1$MeanC<=120)
b6s<-top_frac(b6, 0.05)

b7<- subset(ab1, ab1$MeanC>121 & ab1$MeanC<=140)
b7s<-top_frac(b7, 0.05)

b8<- subset(ab1, ab1$MeanC>141 & ab1$MeanC<=160)
b8s<-top_frac(b8, 0.05)

b9<- subset(ab1, ab1$MeanC>161 & ab1$MeanC<=180)
b9s<-top_frac(b9, 0.05)

b10<- subset(ab1, ab1$MeanC>181 & ab1$MeanC<=200)
b10s<-top_frac(b10, 0.05)

b11<- subset(ab1, ab1$MeanC>201 & ab1$MeanC<=220)
b11s<-top_frac(b11, 0.05)

b12<- subset(ab1, ab1$MeanC>221 & ab1$MeanC<=240)
b12s<-top_frac(b12, 0.05)

b13<- subset(ab1, ab1$MeanC>241 & ab1$MeanC<=260)
b13s<-top_frac(b13, 0.05)

b14<- subset(ab1, ab1$MeanC>261 & ab1$MeanC<=280)
b14s<-top_frac(b14, 0.05)

b15<- subset(ab1, ab1$MeanC>281 & ab1$MeanC<=300)
b15s<-top_frac(b15, 0.05)

t<- rbind(b1s, b2s, b3s, b4s, b5s, b6s, b7s, b8s, b9s, b10s, b11s, b12s, b13s, b14s, b15s)

t_distinct <- length(unique(t$Pos))

write.csv(t_distinct, "AllPosAllSamplesJD_Post5%_noMT_Cov30_W..csv")

#subset for mean coverage >=30
t1<- subset(t, t$MeanC>=25)

t1_distinct <- length(unique(t1$Pos))

write.csv(t1_distinct, "AllPosAllSamplesJD_PostCov_noMT_Cov25_W.csv")


#find mean %DNAm per position
mean<-aggregate(x=d$MethylationPercentage,
                by=list(d$Pos), FUN=mean)

colnames(mean)<- c("Pos", "Mean")

#only keep age-related positions
t1$MeanM <- mean$Mean[match(t1$Pos, mean$Pos)]

write.csv(t1, "Top5Cor_JD_AvgDNAm_noMT__Oct23_Cov25_W.csv")

#subset for mean %DNAm

t2<- subset(t1, t1$MeanM>=40 & t1$MeanM<=60)

t2_distinct <- length(unique(t2$Pos))

write.csv(t2_distinct, "Pos_JD_PostAvgDNAm_noMT__Oct23_Cov25_W.csv")

t2<-separate(t2, col=Pos, into=c('Chromosome', 'Position'), sep='_')

write.csv(t2, "AgeRelatedPos_JD_noMT_Oct23_Cov25_W.csv")


j3<- tabyl(t2, Chromosome)

#find correlation information about selected positions
t2$Pos<-paste(t2$Chromosome,t2$Position,sep="_")

#find by how many samples CpG sites are shared by
merged_df <- merge(d, t2, by ="Pos")

# Keep only the rows where Pos is present in f
subset_df <- subset(merged_df, Pos %in% t2$Pos)

# View the resulting data frame
s <- aggregate(Sample ~ Pos, subset_df, function(x) length(unique(x)))

merged_df <- merge(s, d, by = "Pos")

#Keep only the rows where Pos is present in t2
c <- subset(merged_df, Pos %in% t2$Pos)




