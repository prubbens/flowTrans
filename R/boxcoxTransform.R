boxcoxTransform<-function(transformationId,pars){
    t<-new("transform",.Data=function(data)data<-flowClust::box(data,lambda=pars[1]));
    t@transformationId<-transformationId;
    t
}