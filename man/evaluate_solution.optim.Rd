\name{evaluate_solution.optim}
\alias{evaluate_solution.optim}

\title{
Evaluate MSE Equation (using \link{optim})
}
\description{
This function calls \link[minMSE]{evaluate_solution}, but since \link{optim} requires fn and gr to have the same parameters, it has two additional ones.
}
\usage{
evaluate_solution.optim(par,
                        data,
                        change = NULL,
                        prev_index_list = NULL)
}

\arguments{
  \item{par}{a treatment assignment. The treatment and the data must have the same number of observations (rows).}
  \item{data}{a dataframe containing the covariate vectors for each attribute.}
  \item{change}{the parameter is only needed for optim, it does not play any role.}
  \item{prev_index_list}{the parameter is only needed for optim, it does not play any role.}
}
\value{
  Returns the mean square error value for the current treatment assignment.
}

\references{
\href{https://www.econstor.eu/handle/10419/161931}{Schneider and Schlather (2017)},
}
\author{
Sebastian Schneider \email{sschneider@coll.mpg.de}; \email{sebastian@sebastianschneider.eu},
Giulia Baldini \email{giulia.baldini@uni-bonn.de}
}

\seealso{
\code{\link[MASS]{ginv}}, \code{\link{optim}}
}

\examples{
input <- data.frame(c(10, 20, 30, 40, 130, 40, 120, 5, 10, 80),
                    c(2, 6, 2, 8, 1, 10, 9, 8, 7, 5),
                    c(1, 0, 2, 1, 0, 1, 0, 2, 1, 0))
colnames(input) <- c("IQ", "grade_maths", "both_parents")

evaluate_solution.optim(par = c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0),
                        input)
}

\keyword{treatment}
\keyword{assignment}
\keyword{optim}
\keyword{mse}
