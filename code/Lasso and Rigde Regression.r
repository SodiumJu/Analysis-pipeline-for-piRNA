#For Health Control vs PDND
df=read.csv('/Users/user/Desktop/splin403/MSA/PDND_HC_T.csv',header=T,sep=',')

rownames(df)=df$Name
df=df[,-1] #Cut out the column names

df_num=as.matrix(df[,-1])

###Using 2/3 data###

#number of the rows in dataframe
num_row=nrow(df)
train_rows=sample(1:num_row,0.66*num_row) #Make 2/3 data to be training data

x=df[train_rows,] #training data 
x.train=model.matrix(object=Type~.,data=x)[,-1]
x.test=model.matrix(object=Type~.,data=df[-train_rows,])[,-1]

y=as.matrix(df$Type)
y.train=as.factor(y[train_rows,])
y.test=y[-train_rows,]

###Using 2/3 data###

###Using all data to build model
x=df[,-1]
x.test=as.matrix(df[,-1])
y=df$Type 

###Using all data to build model

fit<-cv.glmnet(x=x.train, y=y.train, alpha=0,type.measure="deviance",family="binomial")
#type.measure=deviance is default setting which is for logistic regression
#If we are going to fit a linear model, then thetype.measure="mse".
#min(fit$cvm)
#plot(fit)

#Define a new list to store elastic net models

library(glmnet)
list.fits=list()

for(i in 0:10){
	fit.name = paste0("alpha",i/10)
	
	list.fits[[fit.name]]=cv.glmnet(x=x.train, y=y.train, alpha=i/10,type.measure="class",family="binomial",grouped=FALSE)
	
}

#Ridge Regression alpha=0
#Lasso Regression alpha=1

library(pROC) #To draw AUC

list.model=list()
list.auc=list()
y.binary=as.numeric(as.factor(y.test))-1 #For character answers data

for(i in 0:10){
	fit.name = paste0("alpha",i/10)
	print(fit.name)
	list.model[[fit.name]]=predict(list.fits[[fit.name]],x.test,s=list.fits[[fit.name]]$lambda.1se,type="response")
	#print(list.fits[[fit.name]]$lambda.1se)
	temp=roc(y.binary,as.numeric(list.model[[fit.name]]),plot=FALSE)
	print(temp$auc)
	list.auc[[fit.name]]=as.numeric(gsub("^Area under the curve: ","",temp$auc))	
}



p=predict(fit,x.test,s=fit$lambda.1se,type="response")
y.binary=as.numeric(as.factor(y.test))-1 #For character answers data
roc(y.binary,as.numeric(list.model$alpha0),plot=TRUE) #If U want to see the plot result

#train_x=model.matrix(object=Type~.,data=df)[,-1]
#train_y=as.factor(df$Type)

#glmmod <- cv.glmnet(x=train_x, y=train_y, family="binomial")
#lasso_tuning<-cv.glmnet(x=train_x, y=train_y, alpha=1)

#logit1=glm(Type~.,family=binomial,data=df)

#logit1=glm(Type~hsa\-\miR\-\3937,family=binomial,data=df)

#Output the coefficient of 
#Ridge Regression alpha=0

df2=as.matrix(coef(list.fits$alpha0))

write.table(df2, file = "/Users/user/Desktop/coef.txt",sep = " ", quote = FALSE, na = "NA")

#Output the coefficient of 
#Lasso Regression alpha=1

df3=as.matrix(coef(list.fits$alpha1))

#Select non zero variables
vec=vector()
for(i in 2:length(df3)){
	if(df3[i]!=0){
		vec=c(vec,i)
		
	}
}

df4=as.matrix(df3[vec,])
write.table(df4, file = "/Users/user/Desktop/coef.txt",sep = " ", quote = FALSE, na = "NA")









