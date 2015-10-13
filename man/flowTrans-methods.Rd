\name{flowTrans-methods}
\docType{methods}
\alias{flowTrans-methods}
\alias{flowTrans,flowFrame,character,character,logical,logical-method}
\title{ S4 method for optimizing parameters of a flow transformation and
  optionally transforming a flowFrame } 
\description{
Transforms a flowFrame and optimizes the parameters for the selected
transformation. 
} 
\section{Methods}{
\describe{

\item{dat = "flowFrame", fun = "character", dims = "character", n2f =
  "logical", parameters.only = "logical"}{ Transforms the flowFrame
  using the function fun over the dimensions dims, with optional
  bivaraite normal prefiltering. Parameters are optimized to make the
  transformed data multivariate normal depending on the selected
  transformation } 
}}
\keyword{methods}
\keyword{ models }
\keyword{ optimize }
