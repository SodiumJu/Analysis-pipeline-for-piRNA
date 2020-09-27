library(edgeR)
library(ggplot2)
library(reshape2)
library(sva)
library(ggrepel)

#function for exact t test
get.exactTest <- function(Type,Control,dataf) {
	test.result <- exactTest(dataf, pair=c(Control,Type)) #Control first
	test.table=test.result$table[,c(1,3)]

	temp=colnames(test.table)
	for(i in 1:length(temp)){
		temp[i]=paste(paste(sub("-","_",Type),Control,sep="."),temp[i],sep=".", collapse = NULL)
	}
	colnames(test.table)=temp	
	return(test.table)
}

#function for finding row names of low expressed values for each columns
#Usage: Low_cut(dataframe,the number of low expressed miRNA we want to cut)
#The function return: A list with each element stored the name list of low expressed miRNA of one patient
Low_cut <- function(df,cut_num){
	low.cut.list=list()
	for(i in 1:ncol(df)){
	
		temp=as.factor(rownames(head(df[order(df[,i]),],cut_num)))
		temp2=c()
		for(j in 1:nrow(df)){
			if(df[j,i]==0){
				temp2=c(temp2,rownames(df[j,]))
			}
		}
		temp2=as.factor(temp2)
		temp.union=union(temp,temp2)
		low.cut.list[[colnames(df)[i]]]<-as.factor(temp.union)		
	}
	return(low.cut.list)
}

#input datasets
TMM.merge.num=read.csv("/Users/user/Desktop/splin403/Parkinson's Disease/2020.4.6.miRNA/combine.table/combine.miRNA.table.num.csv",header=T,sep=',',row.names=1)

#input group info
TMM.merge.group=as.vector(t(read.csv("/Users/user/Desktop/splin403/Parkinson's Disease/2020.4.6.miRNA/combine.table/combine.miRNA.table.csv",header=T,sep=',',row.names=1)[1,]))

#Change Names with batch information
colnames(TMM.merge.num)=sub("_",".",sub("^D","G2.",sub("^X","G1.",colnames(TMM.merge.num))))


##--- changing name with type---
temp_names=colnames(TMM.merge.num)
for(i in 1:length(temp_names)){
	#temp.group=sub("-","",TMM.merge.group[i])
	temp_names[i]=paste(temp_names[i],TMM.merge.group[i],sep="_", collapse = NULL)
}
colnames(TMM.merge.num)=temp_names
##--- changing name with type---

#filiter only miRNA
TMM.merge.num.miRNA= TMM.merge.num[-c(grep('hsa_piR.*',rownames(TMM.merge.num))),]
#get 2546 miRNAs

#first batch has 75(column 1~75) patients
TMM.merge.num.miRNA.1st=TMM.merge.num.miRNA[,c(1:75)]
#second batch has 99(column 76~174) patients
TMM.merge.num.miRNA.2nd=TMM.merge.num.miRNA[,c(76:174)]

#Using cut function to cut low expression miRNA
#Cut around 10% low expressed miRNA
low.cut.list.1st<-Low_cut(TMM.merge.num.miRNA.1st,255)
low.cut.list.2nd<-Low_cut(TMM.merge.num.miRNA.2nd,255)

#Get interesection of low expressed miRNA among different samples
intersect.1st<-intersect(as.factor(as.vector(unlist(low.cut.list.1st[1]))),as.factor(as.vector(unlist(low.cut.list.1st[2]))))
for(i in 3:length(low.cut.list.1st)){
	intersect.1st<-intersect(as.factor(as.vector(unlist(low.cut.list.1st[i]))),intersect.1st)
}

intersect.2nd<-intersect(as.factor(as.vector(unlist(low.cut.list.2nd[1]))),as.factor(as.vector(unlist(low.cut.list.2nd[2]))))
for(i in 3:length(low.cut.list.2nd)){
	intersect.2nd<-intersect(as.factor(as.vector(unlist(low.cut.list.2nd[i]))),intersect.2nd)
}

#Output dropped miRNA
cut.union=union(intersect.1st,intersect.2nd)
write.table(cut.union,file="/Users/user/Desktop/splin403/Parkinson's Disease/2020.4.6.miRNA/output/dropped miRNA/drop_miRNA.txt",sep=" ")

#Trim the union of low expressed miRNA
TMM.merge.num.miRNA.cut=TMM.merge.num.miRNA[!(row.names(TMM.merge.num.miRNA) %in% cut.union),]
TMM.merge.num.miRNA.cut.0=TMM.merge.num.miRNA.cut
#----------------------imputation of missing values--------------------------------

for(i in 1:ncol(TMM.merge.num.miRNA.cut)){
	
	temp=TMM.merge.num.miRNA.cut[,i]
	temp.min=min(temp[temp>0])/2
	for(j in 1:nrow(TMM.merge.num.miRNA.cut)){
		if(TMM.merge.num.miRNA.cut[j,i]==0){
			TMM.merge.num.miRNA.cut[j,i]=temp.min
		}
	}
}

#substitute missing values with 1/2 minimum value of each column (patient)

#----------------------imputation of missing values--------------------------------

#create batches information for sva noralization
input.batch=sub("G1.*",1,sub("G2.*",2,colnames(TMM.merge.num)))

input.combat = ComBat(as.matrix(TMM.merge.num.miRNA.cut), batch=input.batch, mod=NULL, par.prior=TRUE, prior.plots=FALSE)

input.combat.trans=input.combat-(round(min(input.combat),1)-0.1)

#----------------------imputation of missing values again--------------------------------

for(i in 1:ncol(input.combat.trans)){
	
	temp=input.combat.trans[,i]
	temp.min=min(temp[temp>0])/2
	for(j in 1:nrow(input.combat.trans)){
		if(TMM.merge.num.miRNA.cut.0[j,i]==0){
			input.combat.trans[j,i]=temp.min
			#input.combat.trans[j,i]=0 #imputation 0 for missing values 
		}
	}
}

#substitute missing values with 1/2 minimum value of each column (patient)

#----------------------imputation of missing values again--------------------------------



#transfer log2 based counts~
input.combat.trans.base=round(2^input.combat.trans)

#Create dataframe to draw boxplot
df = melt(input.combat.trans, variable.name = "Samples", value.name = "count")
df = data.frame(df, Condition = sub(".*_","",df$Var2))

#ggplot draw boxplots
ggplot(df, aes(x = Var2, y = count, fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count + 1))) + scale_fill_manual(values = c("#619CFF", "#F564E3","forestgreen","darkseagreen1","dodgerblue4","gold"))+ theme(axis.text.x = element_text(angle = 90, hjust = 1,size=8,vjust=0.5))
#---draw boxplot---check

#output sva merged table
write.table(input.combat.trans,file="/Users/user/Desktop/splin403/Parkinson's Disease/2020.4.6.miRNA/output/final.5.13/PD.TMM.sva.miRNA.csv",sep=",")

#input dataframe to EdgeR
y=DGEList(input.combat.trans.base, group=TMM.merge.group)
y <- estimateCommonDisp(y, verbose=TRUE)
y <- estimateTagwiseDisp(y)

#do exact t-test for every type of patients vs HC
for(i in c("PDND","PD-MCI","PDD","MSA-C","MSA-P")){
	assign(paste0(sub("-","_",i),".HC.table"),get.exactTest(i,"HC",y))
}

#merge all exact t-test results
merge.ttest.table=cbind(row.names=rownames(PDND.HC.table),PDND.HC.table,PD_MCI.HC.table,PDD.HC.table,MSA_C.HC.table,MSA_P.HC.table)

#Output feature selection table
write.table(merge.ttest.table,file="/Users/user/Desktop/splin403/Parkinson's Disease/2020.4.6.miRNA/output/final.5.21/merge.ttest.table.half.csv",sep=",")

#-------End-----------

#select miRNA
logFC.target=".HC.logFC"
PValue.target=".HC.PValue"

for(i in c("PDND","PD-MCI","PDD","MSA-C","MSA-P")){
	type.target=sub("-","_",i)
	assign(paste0(type.target,".HC.select"),subset(merge.ttest.table,(get(paste0(type.target,PValue.target))<=0.01) & (abs(get(paste0(type.target,logFC.target)))>=1)))
}

#test

t(t(as.vector(rownames(PDND.HC.select))))
rownames(PD_MCI.HC.select)
rownames(PDD.HC.select)
rownames(MSA_C.HC.select)
rownames(MSA_P.HC.select)
