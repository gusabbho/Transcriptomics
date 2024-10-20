3*5
p53=8
p53
p53+1
p53=16
log2(p53)
a=log2(p53)
p53=log2(p53)
p54=p53^2


gg.ggg5=log2(-16)

x=100:1
x
x[c(26,27,38,42,51,69,74,96)]
x[1:3]
x[43]
x[-25]




x[-4:-8]
x[c(-17,-36,-50,-93)]
x[39]=40.5
x=c(x,14.3)
length(x)

y=sort(x)
y=sort(x,decreasing = T)
y=sort(x,decreasing = F)

tumorp53=9,8.3,7.5,10.1,9.6,8.4,8.3
tumorp53=c(9,8.3,7.5,10.1,9.6,8.4,8.3)

mean(tumorp53)
average(tumorp53)

gene.name=c("p53","oct4","sox2","sry","p53","oct4","sox2","sry")

class(gene.name)

class(x)

x=gene.name
fx=factor(x=gene.name)
fx
gene.name

setwd("D:/Transcriptomics")
getwd()
dir()
data1=read.delim(file ="data1.txt",header = T,
                 sep = "\t",quote = "\"'",dec = ".",)
View(data1)


data2=read.csv(file = "Data1.csv",header = T,sep = ",",row.names = 1 )


class(data2)
class(data1)


head(data2)
head(data2,2)
top10<-head(data2,10)
tail(data2,10)

dim(data2)
ncol(data2)
nrow(data2)


data2[1,1]
data2[1,1:5]
data2[1:10,1:5]
data2[,]

data4=data2[1:6136, c(2,4,6,7)]


colnames(data4)
data4[,2]
data4$tumor2
data4[ , "tumor2"]

row.names(data2)
###row.names(data2)=data2$id

data2[,-1]
data2[-1:-5,]

colnames(data2)
colnames(data2)=c("t1","t2","t3","t4","n1","n2","n3")



data3=data2[order(data2$t1),]       # it shows dataframe's t1 rows ordered based on ascending order of t1 col
data3=data2[order(-data2$t1),]      # like the above code but in descending order

data3=data2[ ,order(-data2[1,])]    # columns of the first row are ordered
data3=data2[data2$t1>1000, ]        # rows of the t1 column are greater than 1000

data3=data2[data2$t1==200.25,]
data3=data2[data2$t1!=200.25,]

data3=data2[c(data2$t1<1000 & data2$t2<1000), ]
data3=data2[c(data2$t1<1000 | data2$t2<1000), ]

###Analyze with R

data2=log2(data2)
boxplot(data2)

boxplot(data2,col = "green",xlab="grups",
        ylab="vaue", main="mydata",ylim=c(4,10))


meanc=colMeans(data2)
class(meanc)
meanc=data.frame(meanc)

meanr=rowMeans(data2)
class(meanr)
meanr=data.frame(meanr)

meanrt=rowMeans(data2[,1:4])
meanrn=rowMeans(data2[,5:7])

meanrn=data.frame(meanrn)
meanrt=data.frame(meanrt)

means=cbind(meanrt,meanrn)


data3=as.matrix(data2)

data4<-t(data2)
class(data2)
class(data4)


###apply...T-test
apply(data2,1,mean)
apply(data2,2,mean)
apply(data2[,1:4],1,mean)
apply(data2,1,sd)
apply(data2,1,sum)
apply(data2,1,median)


###function....T-test & P.VAUE
f(x)=x+1    f(1)

fx=function(x){x+1}



as.numeric(data2[x,y])
fm1=function(x){mean(as.numeric(data2[x,]))}
fm2=function(y){mean(as.numeric(data2[,y]))}


meanrt=sapply(1:6136,fm1)
meanrt=data.frame(meanrt)




###FC

setwd("D:/Transcriptomics")
getwd()
dir()
data=read.delim(file ="data1.txt",header = T,
                sep = "\t",quote = "\'",dec = ".",row.names = 1)
data2=log2(data)
head(data2)

boxplot(data2)

tumor=data2[,1:4]        # make tumor data
norm=data2[,5:7]         # make normal data


tmean=rowMeans(tumor)     # tumor
class(tmean)
tmean=data.frame(tmean)   # we need data frame to proceed with operations


nmean=rowMeans(norm)      # normal
class(nmean)
nmean=data.frame(nmean)


tsd=apply(tumor,1,sd)     # tumor std, 1 means rows, and if we had 2 mean columns
tsd=data.frame(tsd)

nsd=apply(norm,1,sd)      # normal std
nsd=data.frame(nsd)

log2fc=tmean-nmean
colnames(log2fc)="log2fc"

fc=2^log2fc
colnames(fc)="fold"                   

data3=cbind.data.frame(data2,tmean,nmean,tsd,nsd,log2fc,fc)


#####P-Value

fp1=function(x){t.test(data2[x,1:4],data2[x,5:7])$p.value}
p.val1=sapply(1:6136,fp1)
p.val1=data.frame(p.val1)




fp2=function(i){t.test(i[1:4],i[5:7])$ p.value}
fp(2)
p.val2=apply(data2,1,fp2)
p.val2=data.frame(p.val2)


#####adj-p.value

class(data2)
matr=as.matrix(data2)        # using matrix instead of dataframe is another way to work with that data
p.val3=apply(matr,1,fp1)
p.val3=data.frame(p.val3)

data3=cbind(data3,p.val2)
head(data3)

gup=data3[data3$log2fc>=0.5,]
gdown=data3[data3$log2fc<= -0.5,]

gup=subset(data3,(log2fc>=0.5& p.val2<=0.05))
gdown=subset(data3,(log2fc<(-0.5)& p.val2<=0.05))

genes=rbind(gup,gdown)
adj.p=p.adjust(p = genes$p.val2,method = "fdr")
adj.p=data.frame(adj.p)


adj.p=p.adjust(p = data3$p.val2,method = "fdr")
head(genes)
data3=cbind(data3,adj.p)

p.adjust.methods


#### Annotation
setwd("D:/Transcriptomics/")
dir()
annot=read.delim("annot.txt" ,sep = "\t")

genes$id = row.names(genes)


finald=merge(genes,annot,by="id",all = F)
#finald=final.d[,-15:-16]


write.table(x = finald,file = "fdata.txt",quote = F,sep = "\t")
getwd()
