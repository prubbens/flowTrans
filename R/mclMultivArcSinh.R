mclMultivArcSinh<-function(lambda,y)
{  
   a<-lambda[1]
   b<-lambda[2]
   cc<-0
   d<-dim(y)[2]         
   m<-nrow(y)
   n<-ncol(y)
   yt<- log((a+b*y)+sqrt((a+b*y)^2.0+1.0))+cc	         
   Gama<-cov(yt)          
   detcov<-det(Gama)
   partial<-eval(D(expression(log((a+b*y)+sqrt((a+b*y)^2.0+1.0))),"y"))
   partial<-apply(partial,1,function(x)sum(log(abs(x))))
   logG<-sum(2*partial)/length(partial)  
   return((log(detcov)-logG))
}
