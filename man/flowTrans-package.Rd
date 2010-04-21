\name{flowTrans-package}
\Rdversion{1.1}
\alias{flowTrans-package}
\docType{package}
\title{
Parameter Optimization for Flow Cytometry Data Transformation
}
\description{
Maximum likelihood estimation of parameters for common flow cytometry data transformations.
}
\details{
\tabular{ll}{
Package: \tab flowTrans\cr
Type: \tab Package\cr
Version: \tab 0.6.0\cr
Date: \tab 2010-03-25\cr
License: \tab Artistic-2.0\cr
LazyLoad: \tab yes\cr
biocViews: \tab Bioinformatics, FlowCytometry\cr
Depends: \tab  flowCore, methods\cr
Imports: \tab flowCore, methods, stats, flowViz, flowClust\cr
Collate: \tab logicleTransformWrapper.R mclMultivArcSinh.R mclMultivBiexp.R mclMultivLinLog.R mclMultivBoxCox.R arcsinhTransformWrapper.R biexponentialTransformWrapper.R boxcoxTransform.R flowTrans.R linLog.R \cr
Packaged: \tab 2010-03-25 16:49:03 UTC; finak\cr
Built: \tab R 2.10.0; ; 2009-11-30 16:49:05 UTC; unix\cr
}

Index:
\preformatted{
flowTrans               Optimizing transformations for flow cytometry
                        data
}
}
\author{
Greg Finak <greg.finak@ircm.qc.ca>, Juan Manuel-Perez <jperez@ircm.qc.ca>, Raphael Gottardo <raphael.gottardo@ircm.qc.ca>

Maintainer: Greg Finak <greg.finak@ircm.qc.ca>
}
\references{
    Finak G, Perez J M, Weng A, Gottardo R. Optimizing Transformations for Flow Cytometry. (Submitted)
}
\keyword{ package }
\keyword{ models }
\keyword{ optimize }
\examples{
#Load some data
data(GvHD)
#transform the first sample, forward and side scatter.
result<-flowTrans(GvHD[[1]],"mclMultivArcSinh",colnames(GvHD[[1]])[1:2],n2f=FALSE,parameters.only=FALSE);
plot(result$result);
summary(result);
}
