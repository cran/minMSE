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
                        evaluation_function = evaluate_solution,
                        swap_treatment_function = NULL,
                        mse_weights = NULL,
                        change = NULL,
                        prev_index_list = NULL)
}

\arguments{
  \item{par}{a treatment assignment. The treatment and the data must have the same number of observations (rows).}
  \item{data}{a matrix containing the covariate vectors for each attribute.}
  \item{evaluation_function}{the function used to evaluate the MSE treatment. Default is \link[minMSE]{evaluate_solution}, which does not take into account outcome or treatment weights. Other options are \link[minMSE]{evaluate_solution_vector} and \link[minMSE]{evaluate_solution_matrix}.}
  \item{swap_treatment_function}{the parameter is only needed for optim, it does not play any role.}
  \item{mse_weights}{a vector containing the mse_weights for each treatment, or a matrix containing the mse_weights for treatments and outcomes and scaling factors.}
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
input <- matrix(1:30, nrow = 10, ncol = 3)

evaluate_solution.optim(par = c(0, 1, 1, 1, 1, 0, 0, 0, 0, 0),
                        input)
}

\keyword{treatment}
\keyword{assignment}
\keyword{optim}
\keyword{mse}
