setwd("C:\\Users\\Javan\\Desktop\\UFS_collaboration") #set working directory
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggVennDiagram)
library(BioVenn)
library(cowplot)
#Load the first wave data
firstwave = read.csv("firstwave_March2020-October2020.csv",header = T, sep = ',')
firstwave

firstwave = select(firstwave,Lineage) #grab the Lineage column

df = table(firstwave$Lineage) #create the table with frequency of observations

df = as.data.frame(df) #convert the table into a dataframe

df = df[-c(11), ]

names(df)[1] = "Variant_type"
names(df)[2] = "Infected_individuals"

#write.csv(df,"first_wave.csv") #Write the summarized variant types

# Barplot
ggplot(df, aes(x=Variant_type, y=Infected_individuals,fill=Variant_type)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  theme_cowplot()+
  ggtitle("First wave variant types in Free state")+
  xlab("Variant types") +
  ylab("Number of infected individuals") 

#Second wave dataframe
secondwave = read.csv("secondwave_November2020-March2021.csv",header = T, sep = ',')
secondwave

secondwave = select(secondwave,Lineage) #grab the Lineage column

df2 = table(secondwave$Lineage) #create the table with frequency of observations

df2 = as.data.frame(df2) #convert the table into a dataframe

df2 = df2[-c(5), ]

names(df2)[1] = "Variant_type"
names(df2)[2] = "Infected_individuals"

#write.csv(df2, "second_wave.csv")
# Barplot
ggplot(df2, aes(x=Variant_type, y=Infected_individuals,fill=Variant_type)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  theme_cowplot()+
  ggtitle("Second wave variant types in Free state")+
  xlab("Variant types") +
  ylab("Number of infected individuals") 

#Third wave
thirdwave = read.csv("thirdwave_April 2021_September2021.csv",header = T, sep = ',')
thirdwave

thirdwave = select(thirdwave,Lineage) #grab the Lineage column

df3 = table(thirdwave$Lineage) #create the table with frequency of observations

df3 = as.data.frame(df3) #convert the table into a dataframe

df3 = df3[-c(12), ]

names(df3)[1] = "Variant_type"
names(df3)[2] = "Infected_individuals"
#write.csv(df3,"third_wave.csv")

# Barplot
ggplot(df3, aes(x=Variant_type, y=Infected_individuals,fill=Variant_type)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  theme_cowplot()+
  ggtitle("Third wave variant types in Free state")+
  xlab("Variant types") +
  ylab("Number of infected individuals") 

#Getting variant trend
vaNtdata = read.csv("first_wave.csv",header = T,sep = ',')

# minimal horizontal grid theme
ggplot(vaNtdata, aes(Infected_individuals, fill = Variant_type)) + 
  geom_density(alpha = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_cowplot()+
  ggtitle("Variant type classification by waves")+
  xlab("Number of infected individuals") +
  ylab("Density") 


#variants classification by waves
ggplot(vaNtdata, aes(Infected_individuals, Variant_type, color = Wave,size = Infected_individuals)) + 
  geom_point() +
  theme_cowplot() +
  ggtitle("Variant type classification by waves")+
  xlab("Number of infected individuals") +
  ylab("Variant types") 

#====================================
p1 <- ggplot(vaNtdata, aes(Infected_individuals,Variant_type,fill=Variant_type)) + 
  geom_point() + theme_cowplot()
p2 <- ggplot(vaNtdata, aes(Wave,Variant_type,fill=Variant_type)) +
  geom_point() + theme_cowplot()

plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)

#What variants are common and different
# List of items
x <- list(First_wave = as.integer(df$Variant_type), 
          Second_wave = as.integer(df2$Variant_type), 
          Third_wave = as.integer(df3$Variant_type))
# 3D Venn diagram
ggVennDiagram(x)

First_wave = df$Variant_type
Second_wave = df2$Variant_type
Third_wave = df3$Variant_type

biovenn <- draw.venn(First_wave, Second_wave, Third_wave, 
                     subtitle=" ", nrtype="abs",
                     title = "Free state variants of concern",
                     xtitle = "First wave",
                     ytitle = "Second wave",
                     ztitle = "Third wave")

#Variant name and variant class analysis
setwd("C:\\Users\\Javan\\Desktop\\UFS_collaboration\\variantYpes")
f1varname = read.csv("thirdwavevariantTypes.csv",header = T, sep = ',')
f1varname

f1varname = select(f1varname,varname) #grab the varname column

df = table(f1varname$varname) #create the table with frequency of observations

df = as.data.frame(df) #convert the table into a dataframe

df = subset(df,Freq >= 200 ) ;dim(df)

names(df)[1] = "Variant_name"
names(df)[2] = "Frequency"

#write.csv(df,"first_wave.csv") #Write the summarized variant types

# Barplot
ggplot(df, aes(x=Variant_name, y=Frequency,fill=Variant_name)) + 
  geom_bar(stat = "identity") + coord_flip() + 
  theme_cowplot()+
  ggtitle("Third wave viral mutation types")+
  xlab("Viral protein mutations") +
  ylab("Mutations frequency")+
  theme(axis.text = element_text(size = 8)) 



