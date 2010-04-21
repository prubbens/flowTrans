
mclMultivBoxCox<-function(lambda,y){
    yp<-apply(y,2,function(x)x+abs(min(x,0))+.Machine$double.eps)
    yt<-flowClust::box(yp,lambda);
    Gama<-cov(yt)
    detcov<-det(Gama);
    if(lambda!=0){
    partial<-(abs(yp)^(lambda-1))
    partial<-apply(log(abs(partial)),1,sum)
    }else{
     partial<-apply(1/yp,1,function(x)sum(log(abs(x))))
    }
    logG<-sum(2*partial)/length(partial)
    return(log(detcov)-logG)
}

