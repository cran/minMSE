\name{swap_treatment_prev}
\alias{swap_treatment_prev}

\title{
Swap Treatment
}
\description{
Scrambles the elements of the vector and swaps a predefined number of elements. Afterwards, the vector is ordered according to the original ordering and returned.
}
\usage{
swap_treatment_prev(current_treatment,
                    change,
                    prev_index_list)
}

\arguments{
  \item{current_treatment}{a treatment vector to be changed.}
  \item{change}{number of elements that will be changed in the treatment vector.}
  \item{prev_index_list}{index list of the elements that can be changed. The current treatment vector may belong to a previous, unchangeable assignment.}
}

\value{
Returns a new treatment vector.
}

\author{
Sebastian Schneider \email{sschneider@coll.mpg.de}; \email{sebastian@sebastianschneider.eu},
Giulia Baldini \email{giulia.baldini@uni-bonn.de}
}


\examples{
swap_treatment_prev(current_treatment = c(0, 2, 0, 1, 1, 0, 2, 1, 1, 0),
                    change = 2,
                    prev_index_list = c(1, 2, 3, 4))
}

\keyword{treatment}
\keyword{assignment}
\keyword{optim}
\keyword{mse}

