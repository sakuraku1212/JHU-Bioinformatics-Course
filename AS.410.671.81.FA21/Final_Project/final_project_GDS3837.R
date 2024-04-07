#####  load the dataset from GEO, and add the feature to the column name
#BiocManager::install("GEOquery")
library(GEOquery)

dataset <- getGEO("GDS3837", GSEMatrix =TRUE, getGPL=FALSE)
dat.original<-dataset@dataTable@table
dat <- subset(dat.original, select = -c(ID_REF, IDENTIFIER))
rownames(dat)<-(dat.original[,1])
sample<-dataset@dataTable@columns[["sample"]]
factor<-dataset@dataTable@columns[["disease.state"]]
colnames(dat)<-paste(sample,factor)



#####The student should first test for outlier samples and provide visual proof.  Remove these outliers.
##Correlation plot (heat map)
#BiocManager::install("corrplot")
library(corrplot)

#par(mar=c(1,1,1,1))

layout(matrix(c(1,1,1,1,1,1,1,1,2,2), 4, 2, byrow = TRUE))
par(oma=c(5,5,5,5))

#mar=c(2,0,3,0),

cor_dat<-cor(dat, use= "pairwise.complete.obs", method= "pearson")
corrplot(cor(cor_dat), method = "color",  tl.cex=0.5,mar=c(2,0,3,0),
         title =" correlation matrix of lung cancer samples and control samples")

##Hierarchical clustering dendrogram
dat.tree<- t(dat)
dat.dist <- dist(dat.tree,method="euclidean")
dat.clust <- hclust(dat.dist,method="single")

plot(dat.clust,labels=names(dat),cex=0.75, 
     main="Hierarchical Clustering Dendrogram of lung cancer samples and control samples")

##CV vs. mean plot
dat.mean <- apply(na.omit(log2(dat)),2,mean) # calculate mean for each sample
dat.sd <- sqrt(apply(na.omit(log2(dat)),2,var)) # calculate st.deviation for each sample
dat.cv <- dat.sd/dat.mean #calculate cv
plot(dat.mean,dat.cv,main="lung cancer samples and control samples CV vs. Mean", 
     xlab="Mean", ylab="CV", col='blue', cex=1.5, type="n")
points(dat.mean,dat.cv,bg="lightblue",col=1,pch=21)
text(dat.mean,dat.cv,label=dimnames(dat)[[2]],pos=1,cex=0.5)

##Average correlation plot
dat.avg <- apply(cor_dat,1,mean)
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avg r",
     main="Avg correlation of lung cancer samples and control samples",axes=F)
points(dat.avg,bg="red",col=1,pch=21,cex=1.25)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,120,1),col="grey")


###density plot
plot(density(na.omit(log2(dat[,1]))), lwd=2, col=1, lty =1, main="Density plots of log-ratios M",xlab="Log-raios")
for (i in 2:119) {
  lines(density(na.omit(log2(dat[,i]))), col=i, lty =1, lwd=2)
}


##Remove the outliers: GSM494594, GSM494591, GSM494571, GSM494572, GSM494590, GSM494596, GSM494654, GSM494657, GSM494561
dat_clean <- dat[,-c(2, 5, 17, 21, 44, 49, 59, 62, 102)]
sample<-sample[-c(2, 5, 17, 21, 44, 49, 59, 62, 102)]
factor<-factor[-c(2, 5, 17, 21, 44, 49, 59, 62, 102)]
colnames(dat_clean)<-paste(sample,factor)

#####filter out genes that have lowest 10% expression values
m<-rowMeans(dat_clean, na.rm=TRUE)
dat_clean<-dat_clean[m>quantile(m, probs=.10),]


#####conduct some method of feature selection with a statistical test or other machine learning method.  
###Student's t-test
t.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}


cancer<-factor[1:53]
control<-factor[54:111]


pv<- apply(dat_clean,1,t.test.all.genes,s1=c(1:53),s2=c(54:111))

hist(pv,col="lightblue",xlab="p-values",cex.main=0.9,
     main="P-value dist'n between\nlung cancer samples and control samples")
abline(v=.05,col=2,lwd=2)
hist(-log10(pv),col="lightblue",xlab="log10(p-values)",cex.main=0.9,
     main="-log10(pv) dist'n between\nlung cancer samples and control samples")
abline(v= -log10(.05),col=2,lwd=2)




# calculate multiple adjusted p-values with Holm methods
library(multtest)
adj.p.holm<-p.adjust(pv, method = 'holm', n=length(pv))
sum(adj.p.holm<.001)
sum(adj.p.holm<.01)
sum(adj.p.holm<.05)


plot(c(1, length(pv)), range(adj.p.holm), type='n',main='Raw p-values vs Holm P-values', xlab='Genes', ylab='p-values')
points(sort(adj.p.holm), type='p', pch=8,col='red', cex=0.25)
lines(sort(adj.p.holm), col='red', lwd=2)
points(sort(pv), type='p', pch=8,col='blue', cex=0.25)
lines(sort(pv), col='blue', lwd=2)
legend('bottomright', inset=0.02, legend=c('Raw', 'Holm'), lty =c(1, 1),lwd=c(2,2),col=c('blue','red'),cex=0.8)



# calculate multiple adjusted p-values with Bonferroni methods
adj.p.bonferroni<-p.adjust(pv, method = 'bonferroni', n=length(pv))
sum(adj.p.bonferroni<.001)
sum(adj.p.bonferroni<.01)
sum(adj.p.bonferroni<.05)
sum(adj.p.bonferroni<.1)



plot(c(1, length(pv)), range(adj.p.holm), type='n',main='Raw p-values vs Bonferroni P-values', xlab='Genes', ylab='p-values')
points(sort(adj.p.holm), type='p', pch=8,col='red', cex=0.25)
lines(sort(adj.p.holm), col='red', lwd=2)
points(sort(pv), type='p', pch=8,col='blue', cex=0.25)
lines(sort(pv), col='blue', lwd=2)
legend('bottomright', inset=0.02, legend=c('Raw', 'Bonferroni'), lty =c(1, 1),lwd=c(2,2),col=c('blue','red'),cex=0.8)


####Top genes (p-value=0.05)
top_gene<-subset(adj.p.bonferroni,adj.p.bonferroni<=0.05)  #7195
length((top_gene))
length(names(top_gene))


hist(top_gene,col="lightblue",xlab="p-values",cex.main=0.9,
     main="P-value dist'n between\nGastric cancer samples and normal samples")
abline(v=0.05,col=2,lwd=2)
hist(-log10(top_gene),col="lightblue",xlab="log10(p-values)",cex.main=0.9,
     main="-log10(pv) dist'n between\nGastric cancer samples and normal samples")
abline(v= -log10(0.05),col=2,lwd=2)




####Top genes (p-value=0.001 and |fold change|>2)

cancer.m <- apply(log2(dat_clean[,c(1:53)]),1,mean,na.rm=T)
control.m <- apply(log2(dat_clean[,c(54:111)]),1,mean,na.rm=T)
fold <- control.m - cancer.m
summary(fold)


probes <- dimnames(dat_clean[(adj.p.bonferroni<=0.1 & abs(fold)>log2(2)),])[[1]]
probes1 <- dimnames(dat_clean[(pv<(.05/length(pv)) & abs(fold)>log2(2)),])[[1]]





fold_abs<-abs(2^fold)
probesets<-subset(dat_clean, adj.p.bonferroni<=0.05 & abs(fold)>log2(2))
length(probesets) #29
row.names(probesets)


###volcano plot
p.trans <- -1 * log(pv)
par(mar=c(1,1,1,1))
plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',
     main='Volcano Plot\ncancer and control group differences')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(.05)&fold>log2(2))],fold[(p.trans> -log10(.05)&fold>log2(2))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(.05)&fold< -log2(2))],fold[(p.trans> -log10(.05)&fold< -log2(2))],col=1,bg=3,pch=21)
abline(v= -log10(.05))
abline(h= -log2(2))
abline(h=log2(2))



p.trans <- -1 * log10(pv)
x.line <- -log10(.05)	#p-value=0.05
y.line <- log2(2)	#fold change=2
plot(range(p.trans),range(fold),type='n',xlab=expression(paste("-",log[10]," (p-value)")),
     ylab=expression(paste(log[2]," fold change")),main='Volcano Plot\nRat diet study')
points(p.trans,fold,col='black',pch=16)
points(p.trans[(p.trans>x.line&fold>y.line)],fold[(p.trans>x.line&fold>y.line)],col=1,pch=21,bg='red')
points(p.trans[(p.trans>x.line&fold<(-1*y.line))],fold[(p.trans>x.line&fold< (-1*y.line))],col=1,pch=21,bg='green')
abline(v=x.line)
abline(h=y.line)
abline(h=(-1* y.line))







gset <- getGEO("GSE19804", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

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
title <- paste ("GSE19804", "/", annotation(gset), sep ="")
boxplot(ex, boxwex=0.7, notch=T, main=title, outline=FALSE, las=2)
dev.off()

# expression value distribution plot
par(mar=c(4,4,2,1))
title <- paste ("GSE19804", "/", annotation(gset), " value distribution", sep ="")
plotDensities(ex, main=title, legend=F)

# mean-variance trend
ex <- na.omit(ex) # eliminate rows with NAs
plotSA(lmFit(ex), main="Mean variance trend, GSE19804")

# UMAP plot (multi-dimensional scaling)
ex <- ex[!duplicated(ex), ]  # remove duplicates
ump <- umap(t(ex), n_neighbors = 15, random_state = 123)
plot(ump$layout, main="UMAP plot, nbrs=15", xlab="", ylab="", pch=20, cex=1.5)
library("maptools")  # point labels without overlaps
pointLabel(ump$layout, labels = rownames(ump$layout), method="SANN", cex=0.6)



gset <- getGEO("GSE19804", GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

ex <- exprs(gset)
ex <- na.omit(ex) # eliminate rows with NAs
ex <- ex[!duplicated(ex), ]  # remove duplicates
ex<-as.data.frame(ex)
ex_clean <- ex[,-c(2, 5, 17, 21, 44, 49, 59, 62, 102)]


#####filter out genes that have lowest 10% expression values
ex_m<-rowMeans(ex_clean, na.rm=TRUE)
ex_clean<-ex_clean[ex_m>quantile(ex_m, probs=.10),]


#####conduct some method of feature selection with a statistical test or other machine learning method.  
###Student's t-test
t.test.all.genes <- function(x,s1,s2) {
  x1 <- x[s1]
  x2 <- x[s2]
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  t.out <- t.test(x1,x2, alternative="two.sided",var.equal=T)
  out <- as.numeric(t.out$p.value)
  return(out)
}


name<-gset@phenoData@data[["tissue:ch1"]]
name <- name[-c(2, 5, 17, 21, 44, 49, 59, 62, 102)]


ex.pv<- apply(ex_clean,1,t.test.all.genes,s1=name=='paired normal adjacent',s2=name=='lung cancer')

hist(pv,col="lightblue",xlab="p-values",cex.main=0.9,
     main="P-value dist'n between\nlung cancer samples and control samples")
abline(v=.05,col=2,lwd=2)
hist(-log10(pv),col="lightblue",xlab="log10(p-values)",cex.main=0.9,
     main="-log10(pv) dist'n between\nlung cancer samples and control samples")
abline(v= -log10(.05),col=2,lwd=2)




# calculate multiple adjusted p-values with Holm methods
library(multtest)
adj.p.holm<-p.adjust(ex.pv, method = 'holm', n=length(ex.pv))
sum(adj.p.holm<.001)
sum(adj.p.holm<.01)
sum(adj.p.holm<.05)


plot(c(1, length(pv)), range(adj.p.holm), type='n',main='Raw p-values vs Holm P-values', xlab='Genes', ylab='p-values')
points(sort(adj.p.holm), type='p', pch=8,col='red', cex=0.25)
lines(sort(adj.p.holm), col='red', lwd=2)
points(sort(pv), type='p', pch=8,col='blue', cex=0.25)
lines(sort(pv), col='blue', lwd=2)
legend('bottomright', inset=0.02, legend=c('Raw', 'Holm'), lty =c(1, 1),lwd=c(2,2),col=c('blue','red'),cex=0.8)



# calculate multiple adjusted p-values with Bonferroni methods
adj.p.bonferroni<-p.adjust(ex.pv, method = 'bonferroni', n=length(ex.pv))
sum(adj.p.bonferroni<.001)
sum(adj.p.bonferroni<.01)
sum(adj.p.bonferroni<.05)
sum(adj.p.bonferroni<.1)



plot(c(1, length(pv)), range(adj.p.holm), type='n',main='Raw p-values vs Bonferroni P-values', xlab='Genes', ylab='p-values')
points(sort(adj.p.holm), type='p', pch=8,col='red', cex=0.25)
lines(sort(adj.p.holm), col='red', lwd=2)
points(sort(pv), type='p', pch=8,col='blue', cex=0.25)
lines(sort(pv), col='blue', lwd=2)
legend('bottomright', inset=0.02, legend=c('Raw', 'Bonferroni'), lty =c(1, 1),lwd=c(2,2),col=c('blue','red'),cex=0.8)


####Top genes (p-value=0.05)
top_gene<-subset(adj.p.bonferroni,adj.p.bonferroni<=0.05)  #7195
length((top_gene))
length(names(top_gene))


hist(top_gene,col="lightblue",xlab="p-values",cex.main=0.9,
     main="P-value dist'n between\nGastric cancer samples and normal samples")
abline(v=0.05,col=2,lwd=2)
hist(-log10(top_gene),col="lightblue",xlab="log10(p-values)",cex.main=0.9,
     main="-log10(pv) dist'n between\nGastric cancer samples and normal samples")
abline(v= -log10(0.05),col=2,lwd=2)




####Top genes (p-value=0.001 and |fold change|>2)
is.data.frame(ex_clean)
dat.log2<-log2(ex_clean)
dat.log2<-as.numeric(dat.log2)
cancer<-colnames(ex_clean)[1:53]
control<-colnames(ex_clean)[54:111]

control.m <- apply(dat.log2[,control],1,mean,na.rm=T)
cancer.m <- apply(dat.log2[,cancer],1,mean,na.rm=T)
fold <- control.m - cancer.m
summary(fold)

fold_abs<-abs(2^fold)
probesets<-subset(ex_clean, adj.p.bonferroni<=0.05 & fold_abs>2)
length(probesets) #2
row.names(probesets)


###volcano plot
p.trans <- -1 * log(pv)
par(mar=c(1,1,1,1))
plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',
     main='Volcano Plot\ncancer and control group differences')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(.05)&fold>log2(2))],fold[(p.trans> -log10(.05)&fold>log2(2))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(.05)&fold< -log2(2))],fold[(p.trans> -log10(.05)&fold< -log2(2))],col=1,bg=3,pch=21)
abline(v= -log10(.05))
abline(h= -log2(2))
abline(h=log2(2))



p.trans <- -1 * log10(pv)
x.line <- -log10(.05)	#p-value=0.05
y.line <- log2(2)	#fold change=2
plot(range(p.trans),range(fold),type='n',xlab=expression(paste("-",log[10]," (p-value)")),
     ylab=expression(paste(log[2]," fold change")),main='Volcano Plot\nRat diet study')
points(p.trans,fold,col='black',pch=16)
points(p.trans[(p.trans>x.line&fold>y.line)],fold[(p.trans>x.line&fold>y.line)],col=1,pch=21,bg='red')
points(p.trans[(p.trans>x.line&fold<(-1*y.line))],fold[(p.trans>x.line&fold< (-1*y.line))],col=1,pch=21,bg='green')
abline(v=x.line)
abline(h=y.line)
abline(h=(-1* y.line))