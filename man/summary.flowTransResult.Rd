\name{summary.flowTransResult}
\alias{summary.flowTransResult}
\title{
Summarize a flowTransResult object
}
\description{
Print a summary of the contents of a \code{flowTransResult} object, including the type of transformation, the transformed data frame, and the transformed dimensions and transformation parameters.
}
\usage{
summary.flowTransResult(object, ...)
}
\arguments{
  \item{object}{
An object of \code{flowTransResult} returned by a call to \code{flowTrans}
}
  \item{\dots}{
    Additional parameters passed through.
}
}
\value{
    No return value
}
\references{
Finak G, Perez JM, Weng A, Gottardo R. Optimizing Transformations for Flow Cytometry. (Submitted)
}
\author{
Greg Finak <greg.finak@ircm.qc.ca>, Raphael Gottardo <raphael.gottardo@ircm.qc.ca>
}

\seealso{
\code{\link{extractParams}}, \code{\link{flowTransResult}}
}
\examples{
m <- t(matrix(rnorm(10000),2))
colnames(m) <- c("A","B")
m <- flowFrame(m)
res <- flowTrans(m, "mclMultivBoxCox", c("A","B"),
                 n2f=FALSE, parameters.only=FALSE)
summary(res)

}
\keyword{ misc }
