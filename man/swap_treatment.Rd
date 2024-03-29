\name{swap_treatment}
\alias{swap_treatment}

\title{
Swap Treatment
}
\description{
Scrambles the elements of the vector and swaps a predefined number of elements. Afterwards, the vector is ordered according to the original ordering and returned.
}
\usage{
swap_treatment(current_treatment,
               change,
               prev_index_list = NULL)
}

\arguments{
  \item{current_treatment}{a treatment vector to be changed.}
  \item{change}{number of elements that will be changed in the treatment vector.}
  \item{prev_index_list}{not used, needed for compatibility.}
}

\value{
Returns a new treatment vector.
}

\author{
Sebastian Schneider \email{sschneider@coll.mpg.de}; \email{sebastian@sebastianschneider.eu},
Giulia Baldini \email{giulia.baldini@uni-bonn.de}
}


\examples{
swap_treatment(current_treatment = c(0, 2, 0, 1, 1, 0, 2, 1, 1, 0),
               change = 2)
}

\keyword{treatment}
\keyword{assignment}
\keyword{optim}
\keyword{mse}

