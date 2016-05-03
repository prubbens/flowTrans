mclMultivLinLog<-function(lambda,y){
   y<-apply(as.matrix(y),2,function(x)x+abs(min(x,.Machine$double.eps))+.Machine$double.eps)
   yt<- flowTrans:::linlog(y,lambda)   
   detcov<-det(cov(yt));
   partial<-y;
   partial[y<=lambda]<-1/lambda;
   partial[y>lambda]<-log(y[y>lambda]);
   partial<-apply(abs(partial),1,prod)
   logG<-sum(log(partial^2))/length(partial)
   return(log(detcov)-logG)
}
