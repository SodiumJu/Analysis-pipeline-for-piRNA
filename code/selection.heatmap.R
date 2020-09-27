library(ComplexHeatmap) #Import the package
library(circlize)
#Change Colors
library(circlize)
col_fun=col_fun = colorRamp2(c(0, 4, 6), c("blue", "white", "red"))

#---input all the tables---
input=read.csv("/Users/user/Desktop/splin403/Parkinson's Disease/2020.4.6.miRNA/output/combine.miRNA.TMM.cpm.csv",header=T,sep=',',row.names=1)

input.pvalue=read.csv("/Users/user/Desktop/splin403/Parkinson's Disease/2020.4.6.miRNA/output/Exact.t.test/merge.ettest.table.csv",header=T,sep=',',row.names=1)

input.group=read.csv("/Users/user/Desktop/splin403/Parkinson's Disease/2020.4.6.miRNA/output/diagnosis.combine.miRNA.TMM.cpm.csv",header=T,sep=',',row.names=1)[1,]

#---input all the tables---

input.log=(log2(input+1))

input.merge=cbind(row.names=rownames(input.log),input.log,input.pvalue)
colnames(input.merge)=sub("X","",colnames(input.merge))

#select data

HC.vector=c(grep("HC",as.vector(t(input.group))))
PDND.vector=c(grep("PDND",as.vector(t(input.group))))
annot.vector=c(HC.vector,PDND.vector)

#(PDND.HC.logFC>=1.5)|(PDND.HC.logFC<=-1.5)) &

PDND.HC.select=subset(input.merge,select=annot.vector,PDND.HC.PValue<=0.05)

PDND.HC.select.miRNA=PDND.HC.select[-c(grep('hsa_piR.*',rownames(PDND.HC.select))),]

input.matrix=data.matrix(PDND.HC.select.miRNA)

#Annotations
annot_df <- data.frame(Type=t(input.group[annot.vector]))
#col_annot <- list(Type=c("HC"="darkolivegreen1","PDND"="cyan1","PD-MCI"="deepskyblue1","PDD"="dodgerblue4",,"MSA-C"="indianred1","MSA-P"="firebrick3"))

col_annot <- list(Type=c("HC"="darkolivegreen1","PDND"="cyan1"))
#,"PD-MCI"="#3399FF"
annot <- HeatmapAnnotation(df=annot_df, col = col_annot,show_annotation_name = FALSE)





H=Heatmap(input.matrix, name = "",column_names_side = "bottom",show_column_names = TRUE,column_names_rot = 90,na_col = "blue",top_annotation = annot,column_title_gp = gpar(fontsize = 12, fontface = "bold"),column_title = "HC vs PDND",column_names_gp = gpar(fontsize = 10),row_names_gp = gpar(fontsize = 10)) #Draw the heatmap
#Default setting is using complete clustering
draw(H,heatmap_legend_side = "left")