
biexponentialTransformWrapper<-function(transformationId,pars=pars){
    t<-new("transform",.Data=function(x){
        a<-pars[1];
        b<-pars[2];
        c<-pars[3];
        d<-pars[4]
        w<-pars[6]
        f<-pars[5]
        x<-.Call("biexponential_transform",x,a=a,b=b,c=c,d=d,f=f,w=w,tol=.Machine$double.eps^0.25,maxit=as.integer(5000),PACKAGE="flowCore")
    });
    t@transformationId=transformationId;
    t   
}