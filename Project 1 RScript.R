df=read.csv("wine.csv",header = TRUE)
head(df)
dim(df)
sum(is.na(df))
#install.packages('robustX')
require(robustX); library(robustbase);

library(MASS)
df1=df[,-13] #removing label to plot numerical variables: will be put again before bacon

#HISTOGRAMS to check Normality Assumptions
#x
op <- par(mfrow = c(2, 4),
          pty = "s")
for (i in 1:ncol(df1)){
  names=c()
  hist(x[,i],main = colnames(df1)[i],
       xlab = colnames(df1)[i],col='pink')
} 
par(op)


#transforming variables of  data
df1$volatile_acidity<-log(df1$volatile_acidity)
df1$residual_sugar<-log(df1$residual_sugar)
df1$chlorides<-log(df1$chlorides)
df1$density<-log(df1$density)
df1$alcohol<-log(df1$alcohol)
df1$sulphates<-log(df1$sulphates)
df1$free_sulfur_dioxide<-log(df1$free_sulfur_dioxide)
df1$total_sulfur_dioxide<-sqrt(df1$total_sulfur_dioxide)
df1$citric_acid<-sqrt(df1$citric_acid)

#HISTOGRAMS to check Normality Assumptions after transformation
#x
op <- par(mfrow = c(2, 4),
          pty = "s")
for (i in 1:ncol(df1)){
  names=c()
  hist(x[,i],main = colnames(df1)[i],
       xlab = colnames(df1)[i],col="blue")
} 
par(op)

#Putting the label (categorical variable) again in the data before BACON
label=df[,13] 
df1=cbind(df1,label)

#dividing data
x=df1[1:300,]
head(x)
dim(x)


###Detecting Outliers in 1st group X


bacon=function(x, v=2, bss=1, known.out=0, c=5, alpha=.05, pr.opt=0,fname="OutputFile.txt") 
{

# Checking user error
if(!is.numeric(x)){stop("x should be a numeric matrix.")}
if(!is.numeric(v)){stop("v should be 1 or 2.")}
if(!is.numeric(bss)){stop("bss should be 1 or 2.")}
if(!is.numeric(known.out)){stop("known.out should be numeric.")}
if(!is.numeric(c) ||  c < 1){stop("c should be >= 1.")}
if(!is.numeric(alpha) || alpha<=0 || alpha >= 1){stop("alpha should be between 0 and 1.")}
if(!is.numeric(pr.opt)){pr.opt=0}; colp=1 # colp: parameter to color the outliers

if(v==1){
    cat("Original BACON: Version 1,")
    cat("Original BACON: Version 1, ", file=fname,append=FALSE,sep="")
    }  
else{
    cat("Original BACON: Version 2,")
    cat("Original BACON: Version 2, ", file=fname,append=FALSE,sep="")
    }
if(bss!=1) {
    cat(" bss = c x p,")
    cat(" bss = c x p,", file=fname,append=TRUE,sep="")
    }     
else{
    cat(" bss = p + 1,")
    cat(" bss = p + 1,", file=fname,append=TRUE,sep="")
    }

cat("\nc = ",c,"    alpha = ",alpha)
cat("\nc = ",c,"    alpha = ",alpha, file=fname,append=TRUE,sep="")

x<- as.matrix(x); n<-nrow(x); p<-ncol(x); id<-c(1:n); h<- floor((n+p+1)/2)
cat("    n = ",n,"    p = ",p,"\n")
cat("    n = ",n,"    p = ",p,"\n", file=fname,append=TRUE,sep="")
#
chi2<-qchisq(alpha/n, df=p, lower.tail = FALSE, log.p = FALSE) # we could change alpha/n to alpha

md=try(mahalanobis(x, colMeans(x), var(x)), silent = TRUE)
if(!is.numeric(md)){stop("The covariance matrix is singular.")}

par(mfrow=c(2,2))   
qqplot(qchisq(ppoints(n), df=p), md, pch=20,
     main = expression("Q-Q plot of Squared MD vs. quantiles of " *  ~chi[p]^2))
     abline(0, 1, col = 'gray')

plot(id, sqrt(md), pch=20,cex=0.5, xlab = "Index", ylab = "MD",main="Index plot of Mahalanobis distance", 
   col=colp, abline(h=sqrt(qchisq(1-alpha,p)),col=4,lty=2))
if(length(known.out)==n){points(which(known.out==1),sqrt(md)[known.out==1],col=2,pch=20,cex=0.5)}

if(v==1){  # ---------- Version 1 ----------
   dd=md
   bso<-order(md) # is it possible to find the smallest c*p distances without ordering 
}
else{      # ---------- Version 2 ----------
   med<-apply(x,2,median)
   z<-scale(x, center = med, scale = FALSE)
   dd=d<-apply(z*z,1,sum)
   bso<-order(d) # is it possible to find the smallest c*p distances without ordering the entire vector d?
}
# -----------

singsub=0; m=round(c*p)
if(bss!=1) {bss = m} else {bss=p+1}
b=bso[1:bss]; xb<-x[b,]; xb<-as.matrix(xb) 
if(pr.opt==1){
  cat("\nThe mean vector of the basic subset is:\n"); print(colMeans(xb)); 
  cat("\nThe covariance matrix of the basic subset is:\n"); print(var(xb)); 
  cat("\nThe distances for the basic subset are:\n");print(matrix(cbind(b,dd[b]),ncol=2))
}
#
while(bss<=m){          #****************************************
  bss=bss+1;
  if(pr.opt==1){cat("------------------------------------------------\n")}
  d=try(mahalanobis(x, colMeans(xb), var(xb)), silent = TRUE)
  if(is.character(d)) {
     singsub=singsub+1  
     if(bss>m) m=m+1
     b=bso[1:bss]; xb<-x[b,]; xb<-as.matrix(xb) 
     if(pr.opt==1){ cat("Basic subset is singular. Size = ",length(b)," m = ",m,"\n")}
  }
  else{  
    dd=d
    new.order=order(d); bso=new.order;
    b=new.order[1:bss]; xb<-x[b,]; xb<-as.matrix(xb) 
    if(pr.opt==1){cat("Basic subset size = ",length(b)," m = ",m ,"\n")}
    }

  if(pr.opt==1){
    cat("\nThe mean vector of the basic subset is:\n"); print(colMeans(xb)); 
    cat("\nThe covariance matrix of the basic subset is:\n"); print(var(xb)); 
    cat("\nThe distances for the basic subset are:\n");print(matrix(cbind(b,dd[b]),ncol=2))
    }
}

if(pr.opt==1){cat("======  End of MVDBS. ======\n")}

bss=0; singsub=0; itr=0; flag= 1
while(bss!= length(b)  &&  flag){   #****************************************
  bss=length(b)
  Cnp=1.0+(2.0/(n-1.0-3.0*p))+1.0*(p+1)/(n-p)+1.0*(h-bss)/(h+bss);
  if (bss>=h) Cnp=1.0+(2.0/(n-1.0-3.0*p))+1.0*(p+1)/(n-p)
# Cnp=1;
  Cnp=Cnp*Cnp; cutoff=chi2*Cnp
 cat("\nChi2 is:\n"); print(chi2); 
cat("\nCnp is:\n"); print(Cnp); 
cat("\nCutoff is:\n"); print(cutoff); 


#  
  d=try(mahalanobis(x, colMeans(xb), var(xb)), silent = TRUE)
  if(is.character(d)){
    if(pr.opt==1){cat("Basic subset is singular. Size = ",length(b),"\n")}
    b=bso[1:(1+length(b))]; xb=x[b,]; xb<-as.matrix(xb);   singsub=singsub+1 
  }
  else{ b<-id[d< cutoff]; xb<-x[b,]; xb<-as.matrix(xb) ; 
    if (bss > length(b)) flag = 0
  }
#
  itr<-itr+1; 
  if(is.numeric(d)) {
    if(pr.opt==1){
      cat("------------------------------------------------ Iteration = ",itr,"\n")
      cat("The ",length(b)," obs. in the basic subset are:\n ",b,"\n")
      cat("\nThe mean vector of the basic subset is:\n"); print(colMeans(xb)); 
      cat("\nThe covariance matrix of the basic subset is:\n"); print(var(xb)); 
      cat("\nThe distances for the basic subset are:\n");print(matrix(cbind(b,d[b]),ncol=2))
    } 
  }
}
# ----------

cb=colMeans(xb); sb<-var(xb)
if (singsub!=0&& pr.opt==1)  cat("Number of singular subsets = ", singsub,"\n")

b.out<-d>= cutoff
bacon.out<- id[b.out]

qqplot(qchisq(ppoints(n), df=p), d,  pch=20,
   main = expression("Q-Q plot of Squared BD vs. quantiles of " *  ~chi[p]^2))
   abline(0, 1, col = 'gray')

d<-sqrt(d)
plot(id, d, pch=20,cex=0.5, xlab = "Index", ylab = "BD",main="Index plot of BACON distance", col=colp, abline(h=sqrt(cutoff),col=4,lty=2))
if(length(known.out)==n){points(which(known.out==1),d[known.out==1],col=2,pch=20,cex=0.5) }
dev.copy2eps(file="graphs.eps")
d<-matrix(d,ncol=1)
par(mfrow=c(1,1))
#-------------
fn<-0;fp<-0;n.out=sum(known.out); n.in=n-n.out
if((sum(known.out==0)+sum(known.out==1))==n){
  fn=sum(b.out[known.out==1]==0) 
  fp=sum(b.out[known.out==0]==1)
  tp=n.out-fn;   tn=n.in-fp;
#  if(pr.opt==1)
   {
     cat("\n tp =   ",tp,"   fp = ",fp,"\n");   cat("fn = ",fn," tn = ",tn,"\n")
     cat("\n tp =   ",tp,"   fp = ",fp,"    ", file=fname,append=TRUE,sep="") 
     cat("fn = ",fn,"    tn = ",tn,"\n", file=fname,append=TRUE,sep="")
   }
 }

if(length(bacon.out)>0){
#  if(pr.opt==1)
   { 
     cat("\nThe following are the ",length(bacon.out)," outliers in the data:\n")
     cat("\nThe following are the ",length(bacon.out)," outliers in the data:\n", file=fname,append=TRUE,sep="")
     aux=matrix(cbind(bacon.out,d[b.out]),ncol=2);    print(aux)  
     write(t(aux), ncolumns = 2, file=fname,append=TRUE)
   }
  }
else{ 
 #  if(pr.opt==1)
    {
      cat("There are no outliers in the data.\n")
      cat("There are no outliers in the data.\n", file=fname,append=TRUE,sep="")
    }  
  }

#  if(pr.opt==1){
  cat("\nBACON estimate of the mean vector is: \n")
  cat("\nBACON estimate of the mean vector is: \n", file=fname,append=TRUE,sep="")
  aux=colMeans(xb); print(aux)
  write(t(aux), ncolumns = p, file=fname,append=TRUE)

  cat("\nBACON estimate of the covariance matrix is: \n")
  cat("\nBACON estimate of the covariance matrix is: \n", file=fname,append=TRUE,sep="")
  aux=var(xb);  print(aux)
  write(t(aux), ncolumns = p, file=fname,append=TRUE)
 #}

if((sum(known.out==0)+sum(known.out==1))==n){
  list(dist=d,outliers=bacon.out,center=cb,cov=sb,iterations=itr,singsub=singsub,n.out=n.out,fn=fn,fp=fp)
  }
else{list(dist=d,outliers=bacon.out,center=cb,cov=sb,cutoff=sqrt(cutoff),iterations=itr,singsub=singsub,n.out=length(bacon.out))
  }
#----------------

}

df=x[,-13]
head(df)
d=as.matrix(df)
out=bacon(d, v=2, bss=1, known.out=0, c=5, alpha=.05, pr.opt=0,fname="OutputFile.txt")

output= mvBACON(df); 
y = cbind(1:nrow(df),output$dis)
colnames(y) <- c("Index","Distance");
plot(y, pch=19, main = "BACON Distances Red Wine Group")
abline(h=output$limit,  col= "red",  lty=2)
points(y[ ! output$subset, ], pch = 4, col = 2, cex = 1.5)
identify(y)

#Detecting outliers in second group (white wine)
x=df1[301:600,]
head(x)
dim(x)

bacon=function(x, v=2, bss=1, known.out=0, c=5, alpha=.05, pr.opt=0,fname="OutputFile.txt") 
{
  
  # Checking user error
  if(!is.numeric(x)){stop("x should be a numeric matrix.")}
  if(!is.numeric(v)){stop("v should be 1 or 2.")}
  if(!is.numeric(bss)){stop("bss should be 1 or 2.")}
  if(!is.numeric(known.out)){stop("known.out should be numeric.")}
  if(!is.numeric(c) ||  c < 1){stop("c should be >= 1.")}
  if(!is.numeric(alpha) || alpha<=0 || alpha >= 1){stop("alpha should be between 0 and 1.")}
  if(!is.numeric(pr.opt)){pr.opt=0}; colp=1 # colp: parameter to color the outliers
  
  if(v==1){
    cat("Original BACON: Version 1,")
    cat("Original BACON: Version 1, ", file=fname,append=FALSE,sep="")
  }  
  else{
    cat("Original BACON: Version 2,")
    cat("Original BACON: Version 2, ", file=fname,append=FALSE,sep="")
  }
  if(bss!=1) {
    cat(" bss = c x p,")
    cat(" bss = c x p,", file=fname,append=TRUE,sep="")
  }     
  else{
    cat(" bss = p + 1,")
    cat(" bss = p + 1,", file=fname,append=TRUE,sep="")
  }
  
  cat("\nc = ",c,"    alpha = ",alpha)
  cat("\nc = ",c,"    alpha = ",alpha, file=fname,append=TRUE,sep="")
  
  x<- as.matrix(x); n<-nrow(x); p<-ncol(x); id<-c(1:n); h<- floor((n+p+1)/2)
  cat("    n = ",n,"    p = ",p,"\n")
  cat("    n = ",n,"    p = ",p,"\n", file=fname,append=TRUE,sep="")
  #
  chi2<-qchisq(alpha/n, df=p, lower.tail = FALSE, log.p = FALSE) # we could change alpha/n to alpha
  
  md=try(mahalanobis(x, colMeans(x), var(x)), silent = TRUE)
  if(!is.numeric(md)){stop("The covariance matrix is singular.")}
  
  par(mfrow=c(2,2))   
  qqplot(qchisq(ppoints(n), df=p), md, pch=20,
         main = expression("Q-Q plot of Squared MD vs. quantiles of " *  ~chi[p]^2))
  abline(0, 1, col = 'gray')
  
  plot(id, sqrt(md), pch=20,cex=0.5, xlab = "Index", ylab = "MD",main="Index plot of Mahalanobis distance", 
       col=colp, abline(h=sqrt(qchisq(1-alpha,p)),col=4,lty=2))
  if(length(known.out)==n){points(which(known.out==1),sqrt(md)[known.out==1],col=2,pch=20,cex=0.5)}
  
  if(v==1){  # ---------- Version 1 ----------
    dd=md
    bso<-order(md) # is it possible to find the smallest c*p distances without ordering 
  }
  else{      # ---------- Version 2 ----------
    med<-apply(x,2,median)
    z<-scale(x, center = med, scale = FALSE)
    dd=d<-apply(z*z,1,sum)
    bso<-order(d) # is it possible to find the smallest c*p distances without ordering the entire vector d?
  }
  # -----------
  
  singsub=0; m=round(c*p)
  if(bss!=1) {bss = m} else {bss=p+1}
  b=bso[1:bss]; xb<-x[b,]; xb<-as.matrix(xb) 
  if(pr.opt==1){
    cat("\nThe mean vector of the basic subset is:\n"); print(colMeans(xb)); 
    cat("\nThe covariance matrix of the basic subset is:\n"); print(var(xb)); 
    cat("\nThe distances for the basic subset are:\n");print(matrix(cbind(b,dd[b]),ncol=2))
  }
  #
  while(bss<=m){          #****************************************
    bss=bss+1;
    if(pr.opt==1){cat("------------------------------------------------\n")}
    d=try(mahalanobis(x, colMeans(xb), var(xb)), silent = TRUE)
    if(is.character(d)) {
      singsub=singsub+1  
      if(bss>m) m=m+1
      b=bso[1:bss]; xb<-x[b,]; xb<-as.matrix(xb) 
      if(pr.opt==1){ cat("Basic subset is singular. Size = ",length(b)," m = ",m,"\n")}
    }
    else{  
      dd=d
      new.order=order(d); bso=new.order;
      b=new.order[1:bss]; xb<-x[b,]; xb<-as.matrix(xb) 
      if(pr.opt==1){cat("Basic subset size = ",length(b)," m = ",m ,"\n")}
    }
    
    if(pr.opt==1){
      cat("\nThe mean vector of the basic subset is:\n"); print(colMeans(xb)); 
      cat("\nThe covariance matrix of the basic subset is:\n"); print(var(xb)); 
      cat("\nThe distances for the basic subset are:\n");print(matrix(cbind(b,dd[b]),ncol=2))
    }
  }
  
  if(pr.opt==1){cat("======  End of MVDBS. ======\n")}
  
  bss=0; singsub=0; itr=0; flag= 1
  while(bss!= length(b)  &&  flag){   #****************************************
    bss=length(b)
    Cnp=1.0+(2.0/(n-1.0-3.0*p))+1.0*(p+1)/(n-p)+1.0*(h-bss)/(h+bss);
    if (bss>=h) Cnp=1.0+(2.0/(n-1.0-3.0*p))+1.0*(p+1)/(n-p)
    # Cnp=1;
    Cnp=Cnp*Cnp; cutoff=chi2*Cnp
    cat("\nChi2 is:\n"); print(chi2); 
    cat("\nCnp is:\n"); print(Cnp); 
    cat("\nCutoff is:\n"); print(cutoff); 
    
    
    #  
    d=try(mahalanobis(x, colMeans(xb), var(xb)), silent = TRUE)
    if(is.character(d)){
      if(pr.opt==1){cat("Basic subset is singular. Size = ",length(b),"\n")}
      b=bso[1:(1+length(b))]; xb=x[b,]; xb<-as.matrix(xb);   singsub=singsub+1 
    }
    else{ b<-id[d< cutoff]; xb<-x[b,]; xb<-as.matrix(xb) ; 
    if (bss > length(b)) flag = 0
    }
    #
    itr<-itr+1; 
    if(is.numeric(d)) {
      if(pr.opt==1){
        cat("------------------------------------------------ Iteration = ",itr,"\n")
        cat("The ",length(b)," obs. in the basic subset are:\n ",b,"\n")
        cat("\nThe mean vector of the basic subset is:\n"); print(colMeans(xb)); 
        cat("\nThe covariance matrix of the basic subset is:\n"); print(var(xb)); 
        cat("\nThe distances for the basic subset are:\n");print(matrix(cbind(b,d[b]),ncol=2))
      } 
    }
  }
  # ----------
  
  cb=colMeans(xb); sb<-var(xb)
  if (singsub!=0&& pr.opt==1)  cat("Number of singular subsets = ", singsub,"\n")
  
  b.out<-d>= cutoff
  bacon.out<- id[b.out]
  
  qqplot(qchisq(ppoints(n), df=p), d,  pch=20,
         main = expression("Q-Q plot of Squared BD vs. quantiles of " *  ~chi[p]^2))
  abline(0, 1, col = 'gray')
  
  d<-sqrt(d)
  plot(id, d, pch=20,cex=0.5, xlab = "Index", ylab = "BD",main="Index plot of BACON distance", col=colp, abline(h=sqrt(cutoff),col=4,lty=2))
  if(length(known.out)==n){points(which(known.out==1),d[known.out==1],col=2,pch=20,cex=0.5) }
  dev.copy2eps(file="graphs.eps")
  d<-matrix(d,ncol=1)
  par(mfrow=c(1,1))
  #-------------
  fn<-0;fp<-0;n.out=sum(known.out); n.in=n-n.out
  if((sum(known.out==0)+sum(known.out==1))==n){
    fn=sum(b.out[known.out==1]==0) 
    fp=sum(b.out[known.out==0]==1)
    tp=n.out-fn;   tn=n.in-fp;
    #  if(pr.opt==1)
    {
      cat("\n tp =   ",tp,"   fp = ",fp,"\n");   cat("fn = ",fn," tn = ",tn,"\n")
      cat("\n tp =   ",tp,"   fp = ",fp,"    ", file=fname,append=TRUE,sep="") 
      cat("fn = ",fn,"    tn = ",tn,"\n", file=fname,append=TRUE,sep="")
    }
  }
  
  if(length(bacon.out)>0){
    #  if(pr.opt==1)
    { 
      cat("\nThe following are the ",length(bacon.out)," outliers in the data:\n")
      cat("\nThe following are the ",length(bacon.out)," outliers in the data:\n", file=fname,append=TRUE,sep="")
      aux=matrix(cbind(bacon.out,d[b.out]),ncol=2);    print(aux)  
      write(t(aux), ncolumns = 2, file=fname,append=TRUE)
    }
  }
  else{ 
    #  if(pr.opt==1)
    {
      cat("There are no outliers in the data.\n")
      cat("There are no outliers in the data.\n", file=fname,append=TRUE,sep="")
    }  
  }
  
  #  if(pr.opt==1){
  cat("\nBACON estimate of the mean vector is: \n")
  cat("\nBACON estimate of the mean vector is: \n", file=fname,append=TRUE,sep="")
  aux=colMeans(xb); print(aux)
  write(t(aux), ncolumns = p, file=fname,append=TRUE)
  
  cat("\nBACON estimate of the covariance matrix is: \n")
  cat("\nBACON estimate of the covariance matrix is: \n", file=fname,append=TRUE,sep="")
  aux=var(xb);  print(aux)
  write(t(aux), ncolumns = p, file=fname,append=TRUE)
  #}
  
  if((sum(known.out==0)+sum(known.out==1))==n){
    list(dist=d,outliers=bacon.out,center=cb,cov=sb,iterations=itr,singsub=singsub,n.out=n.out,fn=fn,fp=fp)
  }
  else{list(dist=d,outliers=bacon.out,center=cb,cov=sb,cutoff=sqrt(cutoff),iterations=itr,singsub=singsub,n.out=length(bacon.out))
  }
  #----------------
  
}

df=x[,-13]
head(df)
d=as.matrix(df)
out=bacon(d, v=2, bss=1, known.out=0, c=5, alpha=.05, pr.opt=0,fname="OutputFile.txt")

output= mvBACON(df); 
y = cbind(1:nrow(df),output$dis)
colnames(y) <- c("Index","Distance");
plot(y, pch=19, main = "BACON Distances White Wine Group")
abline(h=output$limit,  col= "red",  lty=2)
points(y[ ! output$subset, ], pch = 4, col = 2, cex = 1.5)
identify(y)


#basic subset free outliers 

##Hotellings 
# Classical (non-robust) Hotelling T2 test
ht2=function(x, y) {
n=dim(x)[1]; m=dim(y)[1]; p=dim(x)[2]
xcov=cov(x); ycov=cov(y)
Sp=(n-1)*xcov+(m-1)*ycov; Sp=Sp/(n+m-2)
xcenter=colMeans(x); ycenter=colMeans(y)
d=xcenter-ycenter
T2=t(d)%*%solve(Sp)%*%d
T2=T2*n*m/(n+m)
F=T2*(n+m-p-1)/(p*(n+m-2))
pv=1-pf(F,p,n+m-p-1)
list(xcenter=xcenter,ycenter=ycenter,xcov=xcov,
ycov=ycov, Sp=Sp,T2=T2,F=F,df=c(p,n+m-p-1),pv=pv)
}

x=df1[1:300,]
head(x)
dim(x)

y=df1[301:600,]
head(y)
dim(y)

out1=ht2(x, y)
out1

#install.packages("rrcov")
library(rrcov)
out2=T2.test(x,y)
out2


