\name{scale_vars}
\alias{scale_vars}

\title{
Covariate Vectors Scaling
}
\description{
Scales the data such that the empty fiels (NA) are the mean of the column and all variables are scaled to have variance 1.
In case a variable has cero variance, the variable internally is treated as if it was 0, that way it is not taken into account for treatment assignment.
}
\usage{
scale_vars(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a dataframe containing the covariate vectors for each attribute.}
}

\value{
Returns a dataframe where the empty fields are filled with the mean of the column, and for all variables the variance is 1.
}

\author{
Sebastian Schneider \email{sschneider@coll.mpg.de}; \email{sebastian@sebastianschneider.eu},
Giulia Baldini \email{giulia.baldini@uni-bonn.de}
}

\examples{
input <- data.frame(c(10, 20, 30, 40, 130, 40, 120, 5, 10, 80),
                    c(2, 6, 2, 8, 1, 10, 9, 8, 7, 5),
                    c(1, 0, 2, 1, 0, 1, 0, 2, 1, 0))
colnames(input) <- c("IQ", "grade_maths", "both_parents")

scale_vars(input)
}

\keyword{treatment}
\keyword{assignment}
\keyword{optim}
\keyword{mse}
