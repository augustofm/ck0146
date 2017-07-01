#required libraries
library(stringr) 
require(data.table)
require(ggplot2)
require(MASS)
require(mvtnorm)
require(ggthemes)
require(viridis)

dt<-fread("bigdata_xy_2017_1.txt")
names(dt)<-c("femur","height")

dmatrix<-matrix(nrow = 8, ncol = 2)

#Basis function
dmatrix[,1]<-1
dmatrix[,2]<-dt$femur

# Maximum likelihood estimate
w.ml <- solve(t(dmatrix)%*%dmatrix)%*%t(dmatrix)%*%as.matrix(dt$height)
w.ml <- as.vector(w.ml)
w.ml

ggplot()+
  geom_abline(color="blue",slope=w.ml[2],intercept = w.ml[1],size=1.5)+
  scale_x_continuous(limits = c(40,60))+
  scale_y_continuous(limits = c(150,190))+
  theme_bw()+
  labs(title=paste0("y = ",round(as.numeric(w.ml[1]),2)," + ",round(as.numeric(w.ml[2]),2),"*x"),
       x="Femur length [cm]",
       y="Height [cm]")


# Batch Bayesian estimate
noise<-0.5
precision<-(1/noise)^2
alpha<-0.01
mo <- c(65,2.6)
So <- (1/alpha)*diag(2)

bivn <-   mvrnorm(1000, mu = mo, Sigma = So)
bivn.kde <- kde2d(bivn[,1],bivn[,2],n=50)
image(bivn.kde,
xlab = "w0", ylab = "w1")
contour(bivn.kde,add=TRUE)

Sn.solved <- solve(So)+precision*t(dmatrix)%*%dmatrix  
Sn<-solve(Sn.solved)

mn <- Sn%*%(solve(So)%*%mo+precision*t(dmatrix)%*%as.matrix(dt$height))
mn <- as.vector(mn)

bivn <-   mvrnorm(1000, mu = mn, Sigma = Sn)
bivn.kde <- kde2d(bivn[,1],bivn[,2],n=50)
image(bivn.kde,
xlab = "w0", ylab = "w1")
contour(bivn.kde,add=TRUE)

#selection 8 lines fitted lines from the data space, as an example
x<-rmvnorm(n=8,mean=mn,sigma=Sn)
x<-as.data.table(x)
names(x)<-c("wo","w1")
x[,x:=dt[,1]]
x[,y:=wo+w1*x]

ggplot()+
  geom_abline(color="red",slope=x$w1,intercept = x$wo)+
  geom_point(data=dt,aes(x=femur,y=height),color="blue",size=1.5)+
  scale_x_continuous(limits = c(40,60))+
  scale_y_continuous(limits = c(150,190))+
  theme_bw()+
  labs(title=paste0("Data space after regression"),
       x="Femur length [cm]",
       y="Height [cm]")


# Sequential Bayesian estimate

mo <- c(65,2.6)
So <- (1/alpha)*diag(2)

mlist<-list()
Slist<-list()
ppost<-list()
pprior<-list()


mlist[[1]]<-mo
Slist[[1]]<-So

wo.scan <- seq(0,100,0.05)
w1.scan <- seq(1,5,0.01)

for(i in 1:nrow(dt)){

  Sn.solved <- solve(So)+precision*as.matrix(dmatrix[i,])%*%t(dmatrix[i,])  
  Sn<-solve(Sn.solved)
  mn <- Sn%*%(solve(So)%*%mo+precision*t(dmatrix)%*%as.matrix(dt$height))
  mn <- Sn%*%(solve(So)%*%as.matrix(mo)+precision*as.matrix(dmatrix[i,])%*%as.matrix(dt[i,2]))
  mn <- as.vector(mn)
  w.post.grid<-expand.grid(wo=wo.scan,w1=w1.scan)
  w.post.grid<-as.data.table(w.post.grid)
  w.post.grid[,p:=dmvnorm(w.post.grid,mean=mn,sigma=Sn)]
  
  ppost[[i]]<-ggplot(w.post.grid,aes(wo,w1))+geom_point(aes(colour=p),alpha=0.25)+
    scale_colour_viridis()+labs(title=paste0("# of observations: ",i))+
    scale_x_continuous(limits = c(40,100))+
    scale_y_continuous(limits = c(1.5,3.5))
  
  mlist[[i+1]]<-mn
  Slist[[i+1]]<-Sn
  mo<-mn
  So<-Sn
}


#OUTPUT EXAMPLES
#For first prior

i<-0
bivn <-   mvrnorm(1000, mu = mlist[[i+1]], Sigma = Slist[[i+1]])
bivn.kde <- kde2d(bivn[,1],bivn[,2],n=100)
image(bivn.kde,xlab = "w0", ylab = "w1")
contour(bivn.kde,add=TRUE)

x<-rmvnorm(n=32,mean=mlist[[i+1]],sigma=Slist[[i+1]])
x<-as.data.table(x)
names(x)<-c("wo","w1")
x[,x:=dt[,1]]
x[,y:=wo+w1*x]

ggplot()+
geom_abline(color="red",slope=x$w1,intercept = x$wo)+
scale_x_continuous(limits = c(40,60))+
scale_y_continuous(limits = c(150,190))+
theme_bw()+
labs(title=paste0("Data space with no observation"),
x="Femur length [cm]",
y="Height [cm]")

#After first iteration
i<-1
bivn <-   mvrnorm(1000, mu = mlist[[i+1]], Sigma = Slist[[i+1]])
bivn.kde <- kde2d(bivn[,1],bivn[,2],n=100)
image(bivn.kde,xlab = "w0", ylab = "w1")
contour(bivn.kde,add=TRUE)

x<-rmvnorm(n=8,mean=mlist[[i+1]],sigma=Slist[[i+1]])
x<-as.data.table(x)
names(x)<-c("wo","w1")
x[,x:=dt[,1]]
x[,y:=wo+w1*x]

ggplot()+
  geom_abline(color="red",slope=x$w1,intercept = x$wo)+
  geom_point(data=dt[1:i,],aes(x=femur,y=height),color="blue",size=1.5)+
  scale_x_continuous(limits = c(40,60))+
  scale_y_continuous(limits = c(150,190))+
  theme_bw()+
  labs(title=paste0("Data space with no observation"),
       x="Femur length [cm]",
       y="Height [cm]")


#After two iterations
i<-2
bivn <-   mvrnorm(1000, mu = mlist[[i+1]], Sigma = Slist[[i+1]])
bivn.kde <- kde2d(bivn[,1],bivn[,2],n=100)
image(bivn.kde,xlab = "w0", ylab = "w1")
contour(bivn.kde,add=TRUE)

x<-rmvnorm(n=8,mean=mlist[[i+1]],sigma=Slist[[i+1]])
x<-as.data.table(x)
names(x)<-c("wo","w1")
x[,x:=dt[,1]]
x[,y:=wo+w1*x]

ggplot()+
  geom_abline(color="red",slope=x$w1,intercept = x$wo)+
  geom_point(data=dt[1:i,],aes(x=femur,y=height),color="blue",size=1.5)+
  scale_x_continuous(limits = c(40,60))+
  scale_y_continuous(limits = c(150,190))+
  theme_bw()+
  labs(title=paste0("Data space with no observation"),
       x="Femur length [cm]",
       y="Height [cm]")


#PREDICTING THE HEIGHT OF A PERSON WHOSE FEMUR IS 48CM LONG
aux<-c(1,48)
var<-1/precision+t(as.matrix(aux))%*%Sn%*%as.matrix(aux)
var<-as.numeric(var)

aux.x<-seq(90,190,0.01)
aux.y<-dnorm(aux.x,mean=as.numeric(t(as.matrix(mn))%*%as.matrix(aux)), sd=var)
df<-data.table(height=aux.x,prob=aux.y)

ggplot(df,aes(x=height,y=prob))+geom_area(fill="#59afc6")+
  labs(x="height [cm]",
       y="P(height | femur = 48 [cm])")+
  theme_classic()+
  scale_x_continuous(limits = c(170,180))+
  theme(legend.position = "none") 
