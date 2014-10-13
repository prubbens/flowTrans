\name{flowTransResult-class}
\Rdversion{1.1}
\docType{class}
\alias{flowTransResult-class}
\alias{flowTransResult}
\alias{summary,flowTransResult-class}

\title{Class "flowTransResult" }
\description{Class that contains the results of a call to \code{flowTrans} on a \code{flowFrame}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("flowTransResult", ...)}.
	 The object contains three slots: the transformed flowFrame, a \code{transform} object from which the parameters can be extracted, and a \code{dims} slot containing a character vector of the dimensions that have been transformed.
}
\section{Slots}{
  \describe{
    \item{\code{result}:}{Object of class \code{"flowFrame"}. The transformed data.  }
    \item{\code{trans}:}{Object of class \code{"transform"}. The transform applied to the flowFrame. Parameters can be extracted via \code{summary(trans)}. }
    \item{\code{dims}:}{Object of class \code{"character"}. The names of the dimensions that have been transformed. }
  }
}
\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "flowTransResult")}: Summarize the transformed data. }
	 }
}
\references{Finak G, Perez JM, Weng A, Gottardo R. Optimizing Data Transformation for Flow Cytometry. }
\author{ Greg Finak <greg.finak@ircm.qc.ca>, Raphael Gottardo <raphael.gottardo@ircm.qc.ca>. }

 

\seealso{
	\code{\link{extractParams}}
	\code{\link{summary}}
}
\examples{
m <- t(matrix(rnorm(10000),2))
colnames(m) <- c("A","B")
m <- flowFrame(m)
res <- flowTrans(m,"mclMultivBoxCox", c("A","B"),
                 n2f=FALSE, parameters.only=FALSE)
summary(res)
extractParams(res)
}
\keyword{classes}
