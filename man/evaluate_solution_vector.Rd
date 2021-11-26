\name{evaluate_solution_vector}
\alias{evaluate_solution_vector}

\title{
Evaluate MSE Equation
}
\description{
The function computes the mean squared error for a given treatment assignment.
More precisely: it computes the mean squared error of the treatment effect estimator resulting from the treatment groups as specified by the argument, the treatment assignment vector.
The function uses matrix multiplication and the Moore-Penrose generalized inverse.
}
\usage{
evaluate_solution_vector(treatment, data, mse_weights)
}

\arguments{
  \item{treatment}{a treatment assignment. The treatment and the data must have the same number of observations (rows).}
  \item{data}{a matrix containing the covariate vectors for each attribute.}
  \item{mse_weights}{a vector containing the mse_weights for each treatment.}
}
\value{
  Returns the mean squared error value for the current treatment assignment.
}

\references{
\href{https://www.econstor.eu/handle/10419/161931}{Schneider and Schlather (2017)},
}
\author{
Sebastian Schneider \email{sschneider@coll.mpg.de}; \email{sebastian@sebastianschneider.eu},
Giulia Baldini \email{giulia.baldini@uni-bonn.de}
}

\seealso{
\code{\link[MASS]{ginv}}
}

\examples{
input <- matrix(1:30, nrow = 10, ncol = 3)

evaluate_solution_vector(treatment = c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0),
                         input,
                         mse_weights = c(1, 2))
}

\keyword{treatment}
\keyword{assignment}
\keyword{optim}
\keyword{mse}

