mclMultivBiexp<-function(lambda,y)
{ 
   a<-lambda[1]
   b<-lambda[2]
   c<-lambda[3]
   d<-lambda[4]
   p.w<-lambda[5];
   p.ac<-c(a,c);
   p.bd<-c(b,d);
   p.ab<-c(a,b)
   p.cd<-c(c,d);
   convergence<-1;
   iter<-1;
   while((convergence==1)&(iter!=0)){
    repeat{
        p.ab<-try(optim(par=p.ab,mclMultivBiexpAB,y=y,c=p.cd[1],d=p.cd[2],w=p.w,method="L-BFGS-B",control=list(factr=1,pgtol=1e-8,maxit=10000),lower=c(1e-6,1e-6),upper=c(1,Inf)),silent=F)
        if(!inherits(p.bd,"try-error")){
            break;
        }else{
            p.ab<-runif(2,0,1);
        }
    }
   if(p.ab$convergence!=0){
    convergence<-0;
   }
   p.ab<-p.ab$par;
   repeat{
    p.cd<-try(optim(par=p.cd,mclMultivBiexpCD,y=y,a=p.ab[1],b=p.ab[1],w=p.w,method="L-BFGS-B",control=list(factr=1,pgtol=1e-8,maxit=10000),lower=c(1e-6,1e-6),upper=c(1,Inf)),silent=F)
    if(!inherits(p.cd,"try-error")){
        break;
    }else{
        p.cd<-runif(2,0,1);
    }
   }
    if(p.cd$convergence!=0){
        convergence<-0;
    }
    p.cd<-p.cd$par
    p.w<-optimize(mclMultivBiexpW,interval=c(-10,10),y=y,a=p.ab[1],b=p.ab[2],c=p.cd[1],d=p.cd[2])
    p.w<-p.w$minimum;
   iter<-iter-1;
   }
   if(convergence==1){
        return(list(par=c(p.ab[1],p.ab[2],p.cd[1],p.cd[2],0,p.w[1])))
   }else
   return(list(par=c(0.5,1,0.5,1,0,0))) #return default parameters if we haven't converged
}   

#mclMultivBiexpAC<-function(lambda,y,b,d,w)
#{ 
#    a<-lambda[1]
#    c<-lambda[2]
#    f<-0
#    w<-w;
#    xt<-apply(as.matrix(y),2,function(x)x+abs(min(x,0))+.Machine$double.eps)
#    xt<-.Call("biexponential_transform",y,a,b,c,d,f,w,.Machine$double.eps^0.25,5000,PACKAGE="flowCore");
#    ll<-log(det(cov(xt)))-sum(2*apply(log(abs(1/eval(D(expression(a*exp(b*(xt-w))-c*exp(-d*(xt-w))+f),"xt")))),1,sum)/dim(xt)[1])
#    return(ll);   
#
#}
mclMultivBiexpW<-function(lambda,y,a,c,b,d)
{ 
    w<-lambda[1];
    f<-0;
    xt<-apply(as.matrix(y),2,function(x)x+abs(min(x,0))+.Machine$double.eps)
    xt<-.Call("biexponential_transform",y,a,b,c,d,f,w,.Machine$double.eps^0.25,5000,PACKAGE="flowCore");
    ll<-log(det(cov(xt)))-sum(2*apply(log(abs(1/eval(D(expression(a*exp(b*(xt-w))-c*exp(-d*(xt-w))+f),"xt")))),1,sum)/dim(xt)[1])
    return(ll);   
}
#mclMultivBiexpBD<-function(lambda,y,a,c,w)
#{ 
#    b<-lambda[1]
#    d<-lambda[2]
#    f<-0
#    w<-w;
#    xt<-apply(as.matrix(y),2,function(x)x+abs(min(x,0))+.Machine$double.eps)
#    xt<-.Call("biexponential_transform",y,a,b,c,d,f,w,.Machine$double.eps^0.25,5000,PACKAGE="flowCore");
#    ll<-log(det(cov(xt)))-sum(2*apply(log(abs(1/eval(D(expression(a*exp(b*(xt-w))-c*exp(-d*(xt-w))+f),"xt")))),1,sum)/dim(xt)[1])
#    return(ll);   
#}
 mclMultivBiexpAB<-function(lambda,y,c,d,w){
a<-lambda[1]
b<-lambda[2]
f<-0
w<-w
xt <- .Call("biexponential_transform", y, a, b, c, d, f, 
        w, .Machine$double.eps^0.25, 5000,PACKAGE="flowCore")
ll <- log(det(cov(xt))) - sum(2 * apply(log(abs(1/eval(D(expression(a * 
        exp(b * (xt - w)) - c * exp(-d * (xt - w)) + f), "xt")))), 
        1, sum)/dim(xt)[1])
    return(ll)
}
mclMultivBiexpCD<-function(lambda,y,a,b,w){
f<-0
c<-lambda[1]
d<-lambda[2]
xt <- .Call("biexponential_transform", y, a, b, c, d, f, 
        w, .Machine$double.eps^0.25, 5000,PACKAGE="flowCore")
    ll <- log(det(cov(xt))) - sum(2 * apply(log(abs(1/eval(D(expression(a * 
        exp(b * (xt - w)) - c * exp(-d * (xt - w)) + f), "xt")))), 
        1, sum)/dim(xt)[1])
    return(ll)
}
