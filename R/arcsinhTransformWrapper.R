arcsinhTransformWrapper<-function(transformationId,pars){
    t<-new("transform",.Data=function(x)x<-asinh(pars[1]+pars[2]*(x))+pars[3])
    t@transformationId<-transformationId
    t   
}