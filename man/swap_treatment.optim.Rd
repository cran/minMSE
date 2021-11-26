\name{swap_treatment.optim}
\alias{swap_treatment.optim}

\title{
Swap Treatment (using \link{optim})
}
\description{
This function calls \link[minMSE]{swap_treatment}, but since \link{optim} requires fn and gr to have the same parameters, it has an additional one.
}
\usage{
swap_treatment.optim(current_treatment,
                     data = NULL,
                     evaluation_function = NULL,
                     swap_treatment_function = swap_treatment,
                     mse_weights = NULL,
                     change,
                     prev_index_list = NULL)
}

\arguments{
  \item{current_treatment}{a treatment vector to be changed.}
  \item{data}{the parameter is only needed for optim, it does not play any role.}
  \item{evaluation_function}{the parameter is only needed for optim, it does not play any role.}
  \item{swap_treatment_function}{the function used to create new treatments. Default is \link[minMSE]{swap_treatment}. Other options are \link[minMSE]{swap_treatment_prev} which, given a previous treatment, creates a new treatment assignment that takes the previous one into account.}
  \item{change}{number of elements that will be changed in the treatment vector.}
  \item{prev_index_list}{index list of the elements that can be changed. The current treatment vector may belong to a previous, unchangeable assignment.}
  \item{mse_weights}{the parameter is only needed for optim, it does not play any role.}
}

\value{
Returns a new treatment vector.
}

\author{
Sebastian Schneider \email{sschneider@coll.mpg.de}; \email{sebastian@sebastianschneider.eu},
Giulia Baldini \email{giulia.baldini@uni-bonn.de}
}

\seealso{
\code{\link{optim}}
}

\examples{
swap_treatment.optim(current_treatment = c(0, 2, 0, 1, 1, 0, 2, 1, 1, 0),
                     change = 2)

swap_treatment.optim(current_treatment = c(0, 2, 0, 1, 1, 0, 2, 1, 1, 0),
                     change = 2,
                     prev_index_list = c(1, 2, 3, 4))
}

\keyword{treatment}
\keyword{assignment}
\keyword{optim}
\keyword{mse}
