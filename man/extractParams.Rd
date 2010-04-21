\name{extractParams}
\Rdversion{1.1}
\alias{extractParams}
\title{
Extract the transformation parameters from a flowTransResult
}
\description{
Extracts the transformation parameters from a \code{flowTransResult} object
}
\usage{
extractParams(x, dims = NULL)
}
\arguments{
  \item{x}{
An object of type \code{flowTransResult} returned by a call to \code{flowTrans}
}
  \item{dims}{
A character vector specifying the dimensions for which to extract transformation parameters.
}
}
\value{
    A list of length \code{length(dims)}. Each element contains a vector of parameters for transforming the particular dimension.The names of the list elements correspond to the names of the dimensions.
    
}
    \references{Finak G, Perez JM, Weng A, Gottardo R. Optimizing Data Transformation for Flow Cytometry. (Submitted)}
\author{
Greg Finak <greg.finak@ircm.qc.ca>, Raphael Gottardo <raphael.gottardo@ircm.qc.ca>
}
\note{
    \code{flowTrans} estimates common transformation parameters.
}


\seealso{
\code{\link{summary}},
\code{\link{flowTrans}},
\code{\link{flowTransResult}}
}

\examples{
m <- t(matrix(rnorm(10000),2))
colnames(m) <- c("A","B")
m <- flowFrame(m)
res <- flowTrans(m,"mclMultivBoxCox", c("A","B"), n2f=FALSE,
                 parameters.only=FALSE)
summary(res)
extractParams(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
