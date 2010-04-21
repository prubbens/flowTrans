InvArcSinh<-function(x,a,b){
    (sinh(x)-a)/b
}
InvBiexp<-function(lambda,x){
    a<-lambda[1]
    b<-lambda[2]
    c<-lambda[3]
    d<-lambda[4]
    f<-0
    w<-lambda[5];
    y<-a*exp(b*(x-w))-c*exp(-d*(x-w))+f;
    return(y);
}
InvLinLog<-function(lambda,x){
    x<-as.matrix(x)
    wh<-which(x>=log(lambda),T);
    wl<-which(x<log(lambda),T);
    x[wh]<-exp(x[wh])
    x[wl]<-lambda*(x[wl]-log(lambda))+lambda
    return(x)
}

setGeneric("flowTrans",function(dat,fun,dims,n2f,parameters.only){standardGeneric("flowTrans")});
setMethod("flowTrans",signature(dat="flowFrame",fun="character",dims="character",n2f="logical",parameters.only="logical"),function(dat,fun,dims,n2f,parameters.only=FALSE){
    d<-length(dims);
    fTable<-c("mclMultivArcSinh","mclMultivBiexp","mclMultivBoxCox","mclMultivLinLog");

    tTable<-c("arcsinhTransformWrapper","biexponentialTransformWrapper","boxcoxTransform","linlogTransform")
    npar<-c(3,6,1,1);

     default<-list(c(1,1),c(0.5,1,0.5,1,0),1,0)
     
 
   lower<-list(c(0,1e-6),c(1e-6,1e-6,1e-6,1e-6,-Inf),1e-10,".Machine$double.eps");
    upper<-list(c(1,1),c(Inf,Inf,Inf,Inf,Inf),Inf,"max(y)");
    nms<-list(c("a","b","c"),c("a","b","c","d","f","w"),"theta","theta");
    
    fun.c<-match.arg(fun,fTable);
    w<-match(fun.c,fTable);
    result<-NULL;
       
    fun<-try(getFunction(fun.c));
    if(class(fun)=="try-error"){
        stop(paste("function",fun.c,"not found",sep=""));
    }
if(w%in%(1:4)){
    if(d<2){
        stop("Multivariate transformations need more than two dimensions");
    }
    if(n2f){
            nf<-norm2Filter(filterId="n2f",x=dims[1:2],scale.factor=1);
            dat.t<-transformList(dims,biexponentialTransform())%on%dat
            fres<-filter(dat.t,nf);
            subs<-Subset(dat.t,fres)
            y<-exprs(subs)[,dims]
            y<-InvBiexp(lambda=c(0.5,1,0.5,1,0),as.matrix(y))
        }else{
            y<-exprs(dat)[,dims]
    }
    if(class(upper[[w]])=="character"){
        up <- as.vector(sapply(upper[[w]],function(x)eval(parse(text=x)))) 
    }else{
        up <- upper[[w]];
    }
    if(class(lower[[w]])=="character"){
        low <-as.vector(sapply(lower[[w]],function(x)eval(parse(text=x))))
    }else{
        low <- lower[[w]];
    }
    #Test whether we're doing biexponential transformation, and call a different optimization routine.
    if(class(all.equal(fun,mclMultivBiexp))=="logical"){
        result<-fun(default[[w]],y);
    }else if(class(all.equal(fun,mclMultivLinLog))=="logical"){
        o<-optimize(fun,interval=c(low,up),y=y)$minimum;
        result<-list(par=o);
    }else{ 
        result<-optim(par=default[[w]],fun,y=y,method="L-BFGS-B",lower=low,upper=up);
        if(result$convergence!=0){
            print("Failed to converge: applying default parameters");
            result<-list(par=default[[w]],convergence=0);
        }
    }
    tr<-sapply(1:length(dims),function(i){
           p1<-result$par;
           if(fun.c=="mclMultivArcSinh"){ #Add the fixed c=0 in case of arcsinh
            p1<-c(p1,0);
           }
           (eval(parse(text=paste(tTable[w],"(\"",fun.c,"\",","pars=p1",")",sep=""))));
    })
    trans<-transformList(dims,tr);
    result<-flowFrame(exprs(trans%on%dat));
    result@description<-dat@description;
    result@parameters@data$desc<-dat@parameters@data$desc;
    if(parameters.only){
        p<-summary(tr[[1]])$par #pull the first transformation, since the parameters are common across dimensions.
        names(p)<-nms[[w]];
        return(p);
    }else{
        #Re<-new("flowTransResult",result=result,trans=tr[[1]],dims=dims); #again, pull the first 
        Re<-list(result=result,trans=tr[[1]],dims=dims);
        class(Re)<-"flowTransResult";#transformation object since the parameters are common across all dimensions.
        return(Re);
    }
 }
})