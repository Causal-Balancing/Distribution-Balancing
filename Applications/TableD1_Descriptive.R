#read data
data<-read.csv("Applications/cate-birthdata-cleaned-1stkid-white-detailed.csv")
dataset<-data[data$wt_gain!=0&data$year==2002,]

mean(dataset$bweight[dataset$smoke==0])
mean(dataset$bweight[dataset$smoke==1])

sum(dataset$smoke==1);

#skewness
library('moments')
skewness(dataset$bweight[dataset$smoke==0])
skewness(dataset$bweight[dataset$smoke==1])

#kurtosis
kurtosis(dataset$bweight[dataset$smoke==0])
kurtosis(dataset$bweight[dataset$smoke==1])

#symmetry.test
library('lawstat')
symmetry.test(dataset$bweight[dataset$smoke==0])
symmetry.test(dataset$bweight[dataset$smoke==1])

D<-dataset$smoke
y<-dataset$bweight

X<-data.frame(dataset$mage,dataset$bweight,dataset$medu,dataset$X1st_prenatal,dataset$num_prenatal,
         dataset$male,dataset$married,dataset$fagemiss,dataset$diabetes,dataset$hyperpr,dataset$amnio,
         dataset$ultra,dataset$terms,dataset$drink)

datanew<-cbind(y,D,X)

results<-rbind(colMeans(datanew),apply(datanew,2,sd),apply(datanew,2,median),apply(datanew,2,quantile))

row_names<-c('Mean','sd','Median','Min','25%','50%','75%','Max')
rownames(results) <- row_names

col_names<-c('Birth_weight','Smoking','Mage','Weight gain','Medu','1st Prenatal',
             'num_prenatal','Male','Married','fagemiss','diabetes','hyperpr','amnio',
             'ultra','terms','drink')
colnames(results) <- col_names


