\name{flowTrans}
\Rdversion{1.1}
\alias{flowTrans}

\title{
Optimizing transformations for flow cytometry data.
}
\description{
This is the principal function in the package. The function takes data,
a transformation function name, dimension names, and optional
preselection parameter as input, and runs the maximum likelihood
optimization on the data, performs the transformation, and returns the
transformed data together with a list of transformation paramters for
each dimension. The optimization is set to maximize the likelihood of
the parameters given the data, such that the transformed data is
normally distributed. An optional parameter \code{n2f=[T,F]} allows an
automated preselection of an approximately bivariate normal population
via the \code{norm2Filter} function. The transformation parameters are
then optimized for this preselected region, and finally applied globally
to the entire range of data. The optional argument
\code{parameters.only=[T,F]} specifies whether to return the parameters
or the transformed data. 
}
\usage{
flowTrans(dat, fun, dims, n2f,parameters.only)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
    The data to be transformed, should be a \code{flowFrame} object. 
}
  \item{fun}{
    A character string naming the transformation function to be
    applied. Can be one of: \code{mclMultivBiexp},
    \code{mclMultivBoxCox},\code{mclMultivArcSinh}, and \code{mclMultivLinLog}.
   The \code{Multiv} series of functions estimate common transformation
   parameters for all dimensions.} 
  \item{dims}{
    A character vector identifying the dimensions to be transformed.    
}
  \item{n2f}{
    An optional \code{TRUE, FALSE} value identifying whether the data
    should be prefiltered with \code{norm2Filter} to select a region of
    the data for which to optimize the transformation. Default is
    \code{FALSE}. 
}
    \item{parameters.only}{
        A logical specifying whether to return only the parameters, and
	not the transformed data. 
        \code{mclMultivArcSinh} a, b, c; c=0
        \code{mclMultivBiexp}   a, b, c, d, f, w; f=0
        \code{mclMultivBoxCox}  theta
        \code{mclMultivLinLog}  theta
    }
}
\details{
    The transformation functions are multivariate, common parameter
    transformations. The implementation utilizes a look up table to call
    optimization routines for the different transformations as well as
    optional \code{flowCore} and \code{flowClust} transformation
    implementations to optimize parameters and transform data from a
    single interface. Parameters are optimized such that the transformed
    data is multivariate-normal. 
}
\value{
    if \code{parameters.only = FALSE} then:
\item{flowTransResult}{A \code{flowTransResult} object}
If \code{parameters.only=TRUE} the returned value will be a vector of
common transformation parameters: 
        if \code{mclMultivArcSinh} is called, the returned vector will
	contain a, b, and c=0;
        if \code{mclMultivBiexp}   is called, the returned vector will
	contain a, b, c, d, w; 
        if \code{mclMultivBoxCox}  is called, the returned vector will
	contain theta;
        if \code{mclMultivLinLog}  is called, the returned vector will
	contain theta.
}

\references{
        Finak G, Perez J M, Weng A, Gottardo R. Optimizing
        Transformations for Flow Cytometry. (Submitted) 
}
\author{
Greg Finak <greg.finak@ircm.qc.ca>, Raphael Gottardo
        <raphael.gottardo@ircm.qc.ca> 
}
\examples{
data(GvHD)
result <- flowTrans(GvHD[[1]], "mclMultivArcSinh",
          colnames(GvHD[[1]])[1:2], n2f=FALSE, parameters.only=FALSE)
par(mfrow=c(1,2))
contour(GvHD[[1]][,1:2], main="Untransformed FSC vs SSC")
contour(result$result[,1:2], main="Transformed FSC vs SSC")
summary(result)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ optimize }
