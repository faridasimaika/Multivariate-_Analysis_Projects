## Project 2 - Discriminant Analysis 
df=read.csv("wine.csv",header = TRUE)
head(df)
dim(df)
sum(is.na(df))
#install.packages('robustX')
require(robustX); library(robustbase);

library(MASS)
df1=df[,-13] #removing label to plot numerical variables: will be put back again
#df1=scale(df1)

#dividing the data into G1 and G2

G1 <- df[df$style == 'white',]
G2<- df[df$style == 'red',]

op <- par(mfrow = c(2, 4),
          pty = "s")
for (i in 1:ncol(df1)){
  names=c()
  hist(G1[,i],main = colnames(df1)[i],
       xlab = colnames(df1)[i],col='pink')
} 
par(op)

G1$residual_sugar<-(1/G1$residual_sugar)
G1$chlorides<-log(G1$chlorides)
G1$density<-log(G1$density)
G1$alcohol<-sqrt(G1$alcohol)

op <- par(mfrow = c(2, 4),
          pty = "s")
for (i in 1:ncol(df1)){
  names=c()
  hist(G1[,i],main = colnames(df1)[i],
       xlab = colnames(df1)[i],col="blue")
} 
par(op)

op <- par(mfrow = c(2, 4),
          pty = "s")
for (i in 1:ncol(df1)){
  names=c()
  hist(G2[,i],main = colnames(df1)[i],
       xlab = colnames(df1)[i],col='darkmagenta')
} 
par(op)

G2$citric_acid<-sqrt(G2$citric_acid)
G2$residual_sugar<-log(G2$residual_sugar)
G2$chlorides<-log(G2$chlorides)
G2$alcohol<-log(G2$alcohol)
G2$total_sulfur_dioxide<-sqrt(G2$total_sulfur_dioxide)
G2$quality<-log(G2$quality)

op <- par(mfrow = c(2, 4),
          pty = "s")
for (i in 1:ncol(df1)){
  names=c()
  hist(G2[,i],main = colnames(df1)[i],
       xlab = colnames(df1)[i],col="blue")
} 
par(op)

#
df_t=rbind(G1,G2)
df1=df_t[,-13]
df1=scale(df1)

#Setting the label (categorical variable) 
class=df[,13] 
#df1=cbind(df1,label)

#splitting the data into feature variables and class label
features = df1


# Fisher Linear Discriminant Analysis;
pairs(features,col=factor(class),pch=19)
plot(features,pch=19,col=factor(class)); windows()

library(MASS)
flda=function(x,class) {
  cat("Fisher Linear Discriminant:\n")
  a = lda(x, class); d = predict(a)
  t=table(class, d$class); print(t)
  er=100*(sum(t)-sum(diag(t)))/nrow(x)
  cat("Error Rate =",er,"%\n")
  return(d)
}
rslt=flda(features,class)

# External Validation
loo=function(x,class){
  n=length(class)
  rslt={}
  for(i in 1:n){
    a = lda(x[-i,], class[-i])
    b = predict(a,x[i,])
    rslt[i]=b$class #[i]==class[i]
  }
   return(rslt)
}

a=lda(features,class)
b = predict(a,features)
cat("\nInternal validation:")
print(table(class,b$class))

loo=function(x,class){
  n=length(class)
  rslt={}
  for(i in 1:n){
    a = lda(x[-i,], class[-i])
    b = predict(a,x[i,])
    rslt[i]=b$class #[i]==class[i]
  }
  return(rslt)
}

result=loo(features,class)
cat("\nExternal validation:\n")
print(table(class,result))

#Implementing FLDA 2
flda2=function(x, class) { 
  if(ncol(x)!=2) {cat("Data should be 2-dimensional\n" ); return()}
  t=factor(class);   level=levels(t);   
  if(length(level)!=2) {cat("Data should have only two groups\n" ); return()}
  y=x[class==level[1],];     x=x[class==level[2],]
  n1=nrow(x);   n2=nrow(y); n=n1+n2;  xcenter=colMeans(x);  ycenter=colMeans(y)
  xcov=cov(x); ycov=cov(y);   sp=(n1-1)*xcov+(n2-1)*ycov; sp=sp/(n-2)
  d=xcenter-ycenter; m=(xcenter+ycenter)/2;  a=solve(sp)%*%d;  
  class=c(rep(1,n1),rep(2,n2));      p=1;      
  z=rbind(x,y);       pred=z-matrix(m,ncol=2,nrow=n, byrow=T); pred=as.matrix(pred)
  pred=(pred%*%a)<log(p);    C=(class!=pred+1);     ce=sum(C)
  cat("--------------------------------------------------\n")
  cat("           Correct         Incorrect\n Class  Classification   Classification    Total\n")  
  cat("--------------------------------------------------\n")
  cd1=n1-sum(C[1:n1]);         cat("  1          ",cd1,"             ", n1-cd1,"            ",n1,"\n")
  cd2=n2-sum(C[(n1+1):n]);  cat("  2          ",cd2,"             ", n2-cd2,"            ",n2,"\n")
  cat("--------------------------------------------------\n")
  cat(" Total:     ",cd1+cd2,"             ", n-(cd1+cd2),"           ",n,"\n")  
  cat("--------------------------------------------------\n")
  cat("Error Rate = ",100*(ce/n),"%\n");  const=(sum(a*m)+log(p))/a[2]; slope=-a[1]/a[2]
  z=rbind(x,y);   print(rbind(xcenter[1:2],ycenter[1:2]))
  plot(z[,c(1,2)],col=class,pch=19);  abline(const,slope, col = 'blue'); 
  points(rbind(xcenter[1:2],ycenter[1:2]),pch=19,col=3,cex=1.5); 
  segments(xcenter[1], xcenter[2], ycenter[1], ycenter[2],col=3)
  list(xcenter=xcenter[2:1],ycenter=ycenter[2:1],xcov=xcov,ycov=ycov,sp=sp,a=a,slope=slope,
       const=const,ce=ce,m=m,z=z) }

a=flda2(df1[,1:2],class)

#Projection Method Plot
plot(b$x, pch=19, col=as.factor(class), ylab="LD1")

#Multinomial Method Implementation
library(nnet)
mn = multinom(class ~ ., data = df[,-13])
results = predict(mn)
t=table(class, results)
er=100*(sum(t)-sum(diag(t)))/nrow(df[,-13])
cat("Error Rate =",er,"%\n")

looMN=function(x_subset,class){
  n=length(class)
  rslt={}
  for(i in 1:n){
    y=x_subset[-i,]
    a=multinom(class[-1]~.,x_subset[-i,])
    b=predict(a,x_subset[i,])
    rslt[i]=b
  }
  return(rslt)
}
rslt=looMN(df[,-13],class)
cat("Multinomial")
cat("\nExternal validation:\n")
print(table(class,rslt))
t2=table(class,rslt)
er2=100*(sum(t2)-sum(diag(t2)))/nrow(df[,-13])
cat("Error Rate =",er2,"%\n")