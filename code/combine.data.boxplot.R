library(ggplot2)
library(reshape2)

#1st raw data
input=read.csv("/Users/user/Desktop/splin403/Parkinson's Disease/2020.4.6.miRNA/combine.table/combine.miRNA.table.num.csv",row.names=1)

input.names=read.csv("/Users/user/Desktop/splin403/Parkinson's Disease/2020.4.6.miRNA/combine.table/combine.miRNA.table.csv",row.names=1)

input.types=t(input.names[1,])

pseudoCount = log2(input + 1)
pseudoCount = as.data.frame(pseudoCount)
pseudoCount = as.data.frame(input)

temp_names=sub("^X[0-9]{2}_","1st.",sub("^D[0-9]{2}","2nd",colnames(pseudoCount)))

for(i in 1:length(temp_names)){
	temp_names[i]=paste(temp_names[i],input.types[i],sep="_", collapse = NULL)
}
colnames(pseudoCount)=temp_names

#sub(".*_","",temp_names)

#try first column
hist(pseudoCount[,1])
ggplot(pseudoCount, aes(x = X18_APD163.raw.)) + geom_histogram(fill = "#525252", binwidth = 0.3)

df = melt(pseudoCount, variable.name = "Samples", value.name = "count")
df = data.frame(df, Condition = sub(".*_","",df$Samples))

ggplot(df, aes(x = Samples, y = count, fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count + 1))) + scale_fill_manual(values = c("#619CFF", "#F564E3","forestgreen","darkseagreen1","dodgerblue4","gold"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1,size=5,vjust=0.5))

#text = element_text(size=5),

#---TMM---

input_1st_TMM=read.csv("/Users/user/Desktop/splin403/Parkinson's Disease/2020.3.17.raw.data/YR made/1st_miRNA_raw&TMM_normalized_csv/1st_TMM_normalized_num.csv",row.names=1)

input_1st_TMM_names=read.csv("/Users/user/Desktop/splin403/Parkinson's Disease/2020.3.17.raw.data/YR made/1st_miRNA_raw&TMM_normalized_csv/1st_TMM_normalized.csv",row.names=1)

input_1st_TMM.names=t(input_1st_TMM_names[1,][,-1])

pseudoCount_TMM = log2(input_1st_TMM + 1)
pseudoCount_TMM = as.data.frame(pseudoCount_TMM)

temp_names_TMM=sub(".*_","",sub("\\..*","",colnames(pseudoCount_TMM)))

for(i in 1:length(temp_names_TMM)){
	temp_names_TMM[i]=paste(temp_names_TMM[i],input_1st_TMM.names[i],sep="_", collapse = NULL)
}
colnames(pseudoCount_TMM)=temp_names_TMM

df_1st_TMM = melt(pseudoCount_TMM, variable.name = "Samples", value.name = "count")
df_1st_TMM = data.frame(df_1st_TMM, Condition = sub(".*_","",df_1st_TMM$Samples))

ggplot(df_1st_TMM, aes(x = Samples, y = count, fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count + 1))) + scale_fill_manual(values = c("#619CFF", "#F564E3","forestgreen","darkseagreen1","dodgerblue4","gold"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

