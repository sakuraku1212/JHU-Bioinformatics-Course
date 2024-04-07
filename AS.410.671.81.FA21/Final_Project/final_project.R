#####  load the dataset from GEO, and add the feature to the column name
library(GEOquery)
dataset <- getGEO("GDS4345", GSEMatrix =TRUE, getGPL=FALSE)
dat.original<-dataset@dataTable@table
dat <- subset(dat.original, select = -c(ID_REF, IDENTIFIER))
rownames(dat)<-(dat.original[,1])
sample<-dataset@dataTable@columns[["sample"]]
specimen<-dataset@dataTable@columns[["specimen"]]
colnames(dat)<-paste(sample,specimen)



#####The student should first test for outlier samples and provide visual proof.  Remove these outliers.
##Correlation plot (heat map)
library(corrplot)
cor_dat<-cor(dat, use= "pairwise.complete.obs", method= "pearson")
corrplot(cor(cor_dat), method = "color", title =" correlation matrix of Quadriceps muscle \nfrom cancer patients before and after surgery", 
         mar=c(2,0,3,0))

##Hierarchical clustering dendrogram
dat.tree<- t(dat)
dat.dist <- dist(dat.tree,method="euclidean")
dat.clust <- hclust(dat.dist,method="single")
plot(dat.clust,labels=names(dat),cex=0.75, main="Hierarchical Clustering Dendrogram of Quadriceps muscle \nfrom cancer patients before and after surgery")

##CV vs. mean plot
dat.mean <- apply(log2(dat),2,mean) # calculate mean for each sample
dat.sd <- sqrt(apply(log2(dat),2,var)) # calculate st.deviation for each sample
dat.cv <- dat.sd/dat.mean #calculate cv
plot(dat.mean,dat.cv,main="Quadriceps muscle from cancer patients \nbefore and after surgery CV vs. Mean", xlab="Mean", ylab="CV", col='blue', cex=1.5, type="n")
points(dat.mean,dat.cv,bg="lightblue",col=1,pch=21)
text(dat.mean,dat.cv,label=dimnames(dat)[[2]],pos=1,cex=0.5)

##Average correlation plot
dat.avg <- apply(cor_dat,1,mean)
par(oma=c(3,0.1,0.1,0.1))
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avg r",main="Avg correlation of Quadriceps muscle \nfrom cancer patients before and after surgery",axes=F)
points(dat.avg,bg="red",col=1,pch=21,cex=1.25)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,62.5,1),col="grey")

##Remove the outliers: GSM842017 and GSM842022
dat_clean <- dat[,-c(6, 11)]



#####filter out genes that have lowest 10% expression values
m<-rowMeans(dat_clean)
dat_clean<-dat_clean[m>quantile(m, probs=.10),]



ex <- exprs(gset)
# log2 transform
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)
if (LogC) { ex[which(ex <= 0)] <- NaN
ex <- log2(ex) }

# box-and-whisker plot
dev.new(width=3+ncol(gset)/6, height=5)
par(mar=c(7,4,2,1))
title <- paste ("GSE16382", "/", annotation(gset), sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
dev.off()

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE16382", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

# mean-variance trend
ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE16382")

# UMAP plot (multi-dimensional scaling)
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)