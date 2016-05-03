summary.flowTransResult<-function(object,...){
    if(class(object)!="flowTransResult"){
        stop("object must be of class flowTransResult");
    }
    cat("Transformation: \n");
    print(summary(object$trans)$transformationId);
    cat("Transformed Results: \n");
    print(object$result);
    cat("\n");
    cat("Transformation parameters:\n");
    pars<-extractParams(object)[[1]];
    print(pars);

}
extractParams<-function(x,dims=NULL){
    id<-c("mclMultivArcSinh","mclMultivBiexp","mclMultivBoxCox","mclMultivLinLog")
    pnames<-list(c("a","b","c"),c("a","b","c","d","f","w"),"theta","theta")
    if(class(x)!="flowTransResult"){
        stop("The input argument must be a 'flowTransResult' object")
    }   
    if(!is.null(dims)){
     idx <- which(dims %in% x$dims)
      if (identical(idx, integer(0))) {
      	   message("Channels do not match. Extract transformation parameters
for all the transformed channels.")
          dims <- NULL 
      }
    }
     if(is.null(dims)){
        idx<-1:length(x$dims)
        dims<-x$dims
     }
     p<-summary(x$trans)$pars
     names(p)<-pnames[[match(summary(x$trans)$transformationId,id)]]
     coef<-rep(list(p),length(idx))
     names(coef)<-dims
     coef;
}

#setMethod("summary",signature(object="flowTransResult"),function(object,...){
#cat("Transformation: \n");
#print(summary(object@trans)$transformationId);
#cat("Transformed Results: \n");
#print(object@result);
#cat("\n");
#cat("Transformation parameters:\n");
#pars<-extractParams(object)[[1]];
#print(pars);
#})

#extractParams<-function(x,dims=NULL){
#    id<-c("mclMultivArcSinh","mclMultivBiexp","mclMultivBoxCox","mclMultivLinLog")
#    pnames<-list(c("a","b","c"),c("a","b","c","d","f","w"),"theta","theta")
#    if(!is(x,"flowTransResult")){
#        stop("The input argument must be a 'flowTransResult' object")
#    }   
#    if(!is.null(dims)){
#     idx <- which(dims %in% x@dims)
#      if (identical(idx, integer(0))) {
#      	   message("Channels do not match. Extract transformation parameters
#for all the transformed channels.")
#          dims <- NULL 
#      }
#    }
#     if(is.null(dims)){
#        idx<-1:length(x@dims)
#        dims<-x@dims
#     }
#     p<-summary(x@trans)$pars
#     names(p)<-pnames[[match(summary(x@trans)$transformationId,id)]]
#     coef<-rep(list(p),length(idx))
#     names(coef)<-dims
#     coef;
#}