#####  load the dataset from GEO, and add the feature to the column name
#BiocManager::install("GEOquery")
library(GEOquery)
dataset <- getGEO("GDS4222", GSEMatrix =TRUE, getGPL=FALSE)
dat.original<-dataset@dataTable@table
dat <- subset(dat.original, select = -c(ID_REF, IDENTIFIER))
rownames(dat)<-(dat.original[,1])

dat.1 <- subset(dat.original, select = -c(ID_REF))
rownames(dat.1)<-(dat.original[,1])



#####The student should first test for outlier samples and provide visual proof.  Remove these outliers.
##Correlation plot (heat map)
#BiocManager::install("corrplot")
library(corrplot)
cor_dat<-cor(dat, use= "pairwise.complete.obs", method= "pearson")
corrplot(cor(cor_dat), method = "color",  tl.cex=0.5,mar=c(2,0,3,0),
         title =" correlation matrix of Hodgkins lymphoma samples\nfemale vs male")

##Hierarchical clustering dendrogram
dat.tree<- t(dat)
dat.dist <- dist(dat.tree,method="euclidean")
dat.clust <- hclust(dat.dist,method="single")

plot(dat.clust,labels=names(dat),cex=0.75, 
     main="Hierarchical Clustering Dendrogram of Hodgkins lymphoma samples\nfemale vs male")

##CV vs. mean plot
dat.mean <- apply(na.omit(log2(dat)),2,mean) # calculate mean for each sample
dat.sd <- sqrt(apply(na.omit(log2(dat)),2,var)) # calculate st.deviation for each sample
dat.cv <- dat.sd/dat.mean #calculate cv
plot(dat.mean,dat.cv,main="Hodgkins lymphoma samples female and male\nCV vs. Mean", 
     xlab="Mean", ylab="CV", col='blue', cex=1.5, type="n")
points(dat.mean,dat.cv,bg="lightblue",col=1,pch=21)
text(dat.mean,dat.cv,label=dimnames(dat)[[2]],pos=1,cex=0.5)

##Average correlation plot
par(mfrow=c(1,1))
dat.avg <- apply(cor_dat,1,mean)
plot(c(1,length(dat.avg)),range(dat.avg),type="n",xlab="",ylab="Avg r",
     main="Avg correlation of Hodgkins lymphoma samples female and male",axes=F)
points(dat.avg,bg="red",col=1,pch=21,cex=1.25)
axis(1,at=c(1:length(dat.avg)),labels=dimnames(dat)[[2]],las=2,cex.lab=0.4,cex.axis=0.6)
axis(2)
abline(v=seq(0.5,131,1),col="grey")



##Remove the outliers: GSM447654, GSM447628, GSM447646, GSM447731
sample<-dataset@dataTable@columns[["sample"]]
gender<-dataset@dataTable@columns[["gender"]]
dat_clean <- dat[,-c(11,84, 109, 126)]
sample<-sample[-c(11, 84, 109, 126)]
gender<-gender[-c(11, 84, 109, 126)]
colnames(dat_clean)<-paste(sample,gender)


#check if need Normalization
boxplot(log2(dat_clean), main = "Boxplot of raw data (log2 transformed)",
        ylab = "log2(Intensity)")


#####filter out genes that have lowest 10% expression values
par(mfrow=c(1,2))
m<-rowMeans(dat_clean, na.rm=TRUE)
dat_clean<-dat_clean[m>quantile(m, probs=.10),]
hist(m,breaks = 100, main ='Mean of expression, all genes',xlab = "Row mean expression")
abline(v=quantile(m, probs=.10), col=2, lwd=2)

m_clean<-rowMeans(dat_clean, na.rm=TRUE)
hist(m_clean,breaks = 100, main ='Mean of expression, filtered genes',xlab = "Row mean expression")
abline(v=quantile(m, probs=.10), col=2, lwd=2)


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


#colnames(dat_clean) #female 1:56, male 57:124

pv<- apply(dat_clean,1,t.test.all.genes,s1=c(1:56),s2=c(57:126))

par(mfrow=c(1,2))
par(mar=c(5,5,5,5))
hist(pv,col="lightblue",xlab="p-values",cex.main=0.9, 
     main="P-value dist'n between\nHodgkins lymphoma samples female vs male")
abline(v=.05,col=2,lwd=2)
hist(-log10(pv),col="lightblue",xlab="log10(p-values)",cex.main=0.9, 
     main="-log10(pv) dist'n between\nHodgkins lymphoma samples female vs male")
abline(v= -log10(.05),col=2,lwd=2)



# calculate multiple adjusted p-values with Holm methods
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


plot(c(1, length(pv)), range(adj.p.holm), type='n',main='Raw p-values vs Bonferroni P-values', xlab='Genes', ylab='p-values')
points(sort(adj.p.holm), type='p', pch=8,col='red', cex=0.25)
lines(sort(adj.p.holm), col='red', lwd=2)
points(sort(pv), type='p', pch=8,col='blue', cex=0.25)
lines(sort(pv), col='blue', lwd=2)
legend('bottomright', inset=0.02, legend=c('Raw', 'Bonferroni'), lty =c(1, 1),lwd=c(2,2),col=c('blue','red'),cex=0.8)



####Top genes (p-value with Bonferroni<0.05 and |fold change|>log2(2))
#female 1:55, male 56:124
female.m <- apply(log2(dat_clean[,c(1:55)]),1,mean,na.rm=T)
male.m <- apply(log2(dat_clean[,c(56:124)]),1,mean,na.rm=T)
fold <- female.m - male.m
summary(fold)

top_gene<-subset(adj.p.bonferroni, (adj.p.bonferroni<0.05 & abs(fold)>log2(2)))
names(probesets)
par(mfrow=c(1,1))
par(mar=c(5,5,5,5))
hist(top_gene,col="lightblue",xlab="p-value",cex.main=0.9, breaks=20,
     main="Distribution of p-values for differentially expressed genes")



###volcano plot
par(mfrow=c(1,1))
p.trans <- -1 * log10(pv)
par(mar=c(5,5,5,5))
plot(range(p.trans),range(fold),type='n',xlab='-1*log10(p-value)',ylab='fold change',
     main='Volcano Plot between\nHodgkins lymphoma samples female vs male')
points(p.trans,fold,col='black',pch=21,bg=1)
points(p.trans[(p.trans> -log10(.05)&fold>log2(2))],fold[(p.trans> -log10(.05)&fold>log2(2))],col=1,bg=2,pch=21)
points(p.trans[(p.trans> -log10(.05)&fold< -log2(2))],fold[(p.trans> -log10(.05)&fold< -log2(2))],col=1,bg=3,pch=21)
abline(v= -log10(.05))
abline(h= -log2(2))
abline(h=log2(2))



#clustering
#female 1:56, male 57:126

library(limma)
design <-cbind(Grp1=1,Grp2vs1=c(rep(0,56),
                                rep(1,70)))
fit <- lmFit(dat_clean,design)
fit <- eBayes(fit)
#p-values
p.value<-fit$p.value[,2]


#Cluster Analysis
dat_new<-dat_clean[names(top_gene),] #subset the gene that selected as a new dataset
dis.cols <- dist(t(dat_new), method = 'manhattan')
hclust.cols <- hclust(dis.cols, method = 'median')
par(mar=rep(2,4))
plot(hclust.cols, main='hierarchical clustering of Hodgkins lymphoma samples')


hm.rg <- c("#FF0000","#CC0000","#990000","#660000",
           "#330000","#000000","#000000","#0A3300",
           "#146600","#1F9900","#29CC00","#33FF00")
heatmap(as.matrix(dat_new),col=hm.rg,
        main='heatmap of Hodgkins lymphoma samples')


dat.pca <- prcomp(t(dat_new))
dat.loadings <- dat.pca$x[,1:2]
cl <- kmeans(dat.loadings, centers=2, iter.max=20)

par(mar=c(5,5,5,5))
plot(dat.loadings, type="n",
     main = "PCA Cluster Plot of Hodgkins lymphoma samples with Classification Labels")
points(kclust$centers, col = 1:2, pch = "*",cex=2.5)
text(dat.pca.2, col = kclust$cluster, labels=names(kclust$cluster), cex=0.85, pos=3)


#Classification
library(MASS)
clas <- names(dat_new)
clas[grep("female",clas)] <- rep("F",length(clas[grep("female",clas)]))
clas[grep("male",clas)] <- rep("M",length(clas[grep("male",clas)]))


datx <- as.data.frame(t(dat_new))
datx <- data.frame(clas,datx)

# create train/test sets on 60%/40% of data
#female 1:56, male 57:126
train<-rbind(datx[1:34,],datx[57:98,])
test<-rbind(datx[35:56,],datx[99:126,])

dat.lda <- lda(clas~.,train)
dat.pred <- predict(dat.lda,test)
table(dat.pred$class,as.vector(test[['clas']]))

par(mfrow=c(1,1))
par(mar=c(5.5, 5.5, 4, 5))
plot(dat.pred$x,col=as.numeric(factor(test[['clas']])),pch=19,
     main="Discriminant function for Hodgkins lymphoma samples\nfemale vs male")
legend("bottomright",  legend=c("Female", "Male"),
       pch=19, col=c(1,2),cex=0.8, bty="n", xpd = TRUE)


#####top 5 discriminant genes (positive and negative direction) 
neg_5<-sort(fold[names(top_gene)])[1:5]
pos_5<-sort(fold[names(top_gene)])[11:15]
dat_clean[names(top_gene),] 
dat.1[names(sort(fold,decreasing = TRUE)[16:20]),1]
names(sort(fold,decreasing = TRUE)[16:20])

