\name{evaluate_solution}
\alias{evaluate_solution}

\title{
Evaluate MSE Equation
}
\description{
The function computes the mean squared error for a given treatment assignment.
More precisely: it computes the mean squared error of the treatment effect estimator resulting from the treatment groups as specified by the argument, the treatment assignment vector.
The function uses matrix multiplication and the Moore-Penrose generalized inverse.
}
\usage{
evaluate_solution(par, data)
}

\arguments{
  \item{par}{a treatment assignment. The treatment and the data must have the same number of observations (rows).}
  \item{data}{a dataframe containing the covariate vectors for each attribute.}
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
input <- data.frame(c(10, 20, 30, 40, 130, 40, 120, 5, 10, 80),
                    c(2, 6, 2, 8, 1, 10, 9, 8, 7, 5),
                    c(1, 0, 2, 1, 0, 1, 0, 2, 1, 0))
colnames(input) <- c("IQ", "grade_maths", "both_parents")

evaluate_solution(par = c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0),
                  input)
}

\keyword{treatment}
\keyword{assignment}
\keyword{optim}
\keyword{mse}

