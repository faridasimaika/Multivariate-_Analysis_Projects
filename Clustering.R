#clustering - Project 3 Farida Simaika 900201753

df=read.csv("wine.csv",header = TRUE)
head(df)
dim(df)
sum(is.na(df))
df=scale(df)
#install.packages('robustX')
require(robustX); library(robustbase);
library(MASS)
x=df[,-13] #removing label 
class=df[,13]

# Hierarchical Clustering with euclidean and complete linkage
d=dist(x,method="euclidean") 
hc=hclust(d, method = "complete")
plot(hc)
clusters=cutree(hc, k=4) 
rect.hclust(hc, k=4, border="red")

cm=table(clusters,class); print(cm); #  plot(x,pch=19,col=hd.clust)
cat(" Wrong error.rate =",1-sum(diag(cm))/sum(cm),"\n")


#Goodness of fit
# Computes R2: 0.673
R2=function(x,clusters,k){
  n=nrow(x); tss=var(x); tss=(n-1)*sum(diag(tss));
  wss=0
  for(j in 1:k){
    cj=x[clusters==j,]; nj=nrow(cj);
    vj=var(cj); wssj=0
    if(is.matrix(cj)) wssj=(nj-1)*sum(diag(vj));
    wss=wss+wssj
  }
  r2=1-wss/tss; cat("R2 = ",r2,"\n")
  return(r2)
}
r2=R2(as.matrix(x),clusters,4)


# Hierarchical Clustering with manhattan and average linkage 246 outlier
d=dist(x,method="manhattan") 
hc=hclust(d, method = "average")
plot(hc)
clusters=cutree(hc, k=4) 
rect.hclust(hc, k=4, border="red")

cm=table(clusters,class); print(cm); #  plot(x,pch=19,col=hd.clust)
cat(" Wrong error.rate =",1-sum(diag(cm))/sum(cm),"\n")


R2=function(x,clusters,k){
  n=nrow(x); tss=var(x); tss=(n-1)*sum(diag(tss));
  wss=0
  for(j in 1:k){
    cj=x[clusters==j,]; nj=nrow(cj);
    vj=var(cj); wssj=0
    if(is.matrix(cj)) wssj=(nj-1)*sum(diag(vj));
    wss=wss+wssj
  }
  r2=1-wss/tss; cat("R2 = ",r2,"\n")
  return(r2)
}
r2=R2(as.matrix(x),clusters,4)


# Hierarchical Clustering with euclidean and average linkage
d=dist(x,method="euclidean") 
hc=hclust(d, method = "average")
plot(hc)
clusters=cutree(hc, k=4) 
rect.hclust(hc, k=4, border="red")
cm=table(clusters,class); print(cm); #  plot(x,pch=19,col=hd.clust)
cat(" Wrong error.rate =",1-sum(diag(cm))/sum(cm),"\n")


R2=function(x,clusters,k){
  n=nrow(x); tss=var(x); tss=(n-1)*sum(diag(tss));
  wss=0
  for(j in 1:k){
    cj=x[clusters==j,]; nj=nrow(cj);
    vj=var(cj); wssj=0
    if(is.matrix(cj)) wssj=(nj-1)*sum(diag(vj));
    wss=wss+wssj
  }
  r2=1-wss/tss; cat("R2 = ",r2,"\n")
  return(r2)
}
r2=R2(as.matrix(x),clusters,4)


#Clustering complete and manhattan

d=dist(x,method="manhattan") 
hc=hclust(d, method = "complete")
plot(hc)
clusters=cutree(hc, k=4) 
rect.hclust(hc, k=4, border="red")
cm=table(clusters,class); print(cm); #  plot(x,pch=19,col=hd.clust)
cat(" Wrong error.rate =",1-sum(diag(cm))/sum(cm),"\n")

#Goodness of fit
# Computes R2
R2=function(x,clusters,k){
  n=nrow(x); tss=var(x); tss=(n-1)*sum(diag(tss));
  wss=0
  for(j in 1:k){
    cj=x[clusters==j,]; nj=nrow(cj);
    vj=var(cj); wssj=0
    if(is.matrix(cj)) wssj=(nj-1)*sum(diag(vj));
      wss=wss+wssj
  }
  r2=1-wss/tss; cat("R2 = ",r2,"\n")
  return(r2)
}
r2=R2(as.matrix(x),clusters,4)

# kMeans Clustering K=2

set.seed(1949)
kmc = kmeans(x, centers=2);    clusters=kmc$cluster

cm=table(clusters,class); print(cm); #  plot(x,pch=19,col=hd.clust)
cat(" Wrong error.rate =",1-sum(diag(cm))/sum(cm),"\n")

table(clusters)


plot(df, pch=19, col = kmc$cluster)

R2(as.matrix(x),clusters,2)


#K-means k=5
set.seed(1949)
kmc = kmeans(x, centers=2);    clusters=kmc$cluster

cm=table(clusters,class); print(cm); #  plot(x,pch=19,col=hd.clust)
cat(" Wrong error.rate =",1-sum(diag(cm))/sum(cm),"\n")

table(clusters)


plot(df, pch=19, col = kmc$cluster)

R2(as.matrix(x),clusters,4)


#Hungarian method:

match.label=minWeightBipartiteMatching(clusters,class)
id=c(match.label); cm=cm[,id];  print(cm)
cat("Error Rate:",1-sum(diag(cm))/sum(cm),"\n")


# Determining the number of clusters
x=scale(x, center=TRUE, scale=TRUE); # scale x
wss = (nrow(x)-1)*sum(apply(x,2,var))
for (i in 2:10) {
  wss[i] = sum(kmeans(x,centers=i)$withinss)}
plot(wss, type="b", pch=19, xlab="k",
     ylab="WSS", main="The L-Curve")

kmc = kmeans(x, centers=5);    clusters=kmc$cluster

cm=table(clusters,class); print(cm); #  plot(x,pch=19,col=hd.clust)
cat(" Wrong error.rate =",1-sum(diag(cm))/sum(cm),"\n")

table(clusters)


#Hungarian method:

match.label=minWeightBipartiteMatching(5,2)
id=c(match.label); cm=cm[,id];  print(cm)
cat("Error Rate:",1-sum(diag(cm))/sum(cm),"\n")
plot(df, pch=19, col = kmc$cluster)

R2=function(x,clusters,k){
  n=nrow(x); tss=var(x); tss=(n-1)*sum(diag(tss));
  wss=0
  for(j in 1:k){
    cj=x[clusters==j,]; nj=nrow(cj);
    vj=var(cj); wssj=0
    if(is.matrix(cj)) wssj=(nj-1)*sum(diag(vj));
    wss=wss+wssj
  }
  r2=1-wss/tss; cat("R2 = ",r2,"\n")
  return(r2)
}
r2=R2(as.matrix(x),clusters,5)



