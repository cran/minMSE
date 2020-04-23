\name{sample_with_prev_treatment}
\alias{sample_with_prev_treatment}

\title{
Sample Under Consideration of an Already Treated Subset of Units
}
\description{
Given a previous treatment assignment vector for a subset of all observations that treatment assigment is desired for, the function computes a treatment assignment vector for which the previously assigned units are not changed.
At a later step, the previously assigned units are also taken into consideration for computation of the score value, the min MSE function, to achieve balanced treatment groups.
}
\usage{
sample_with_prev_treatment(prev_treatment, n_treatments, n_per_group)
}

\arguments{
  \item{prev_treatment}{takes a numerical vector of partial treatment assignment as argument, and -- for a start -- assigns the missing units (where the value is NA) to a random treatment group, while maintaining the same proportions in the groups.}
  \item{n_treatments}{specifies the number of treatment groups desired (in addition to the control group). They might be more than the ones already defined in prev_treatment.}
  \item{n_per_group}{specifies the distribution of participants per experimental group. It is either an integer, which produces even-sized groups, or a vector which has the same length as the number of experimental groups.}
}

\value{
Returns a treatment assignment vector where the observations given by prev_treatment are unmodified, and the others are assigned to a group.
}

\author{
Sebastian Schneider \email{sschneider@coll.mpg.de}; \email{sebastian@sebastianschneider.eu},
Giulia Baldini \email{giulia.baldini@uni-bonn.de}
}

\examples{
sample_with_prev_treatment(prev_treatment = c(0, NA, NA, NA, 1, NA, NA, NA, NA, NA),
                           n_treatments = 2,
                           n_per_group = c(2, 4, 4))
}

\keyword{treatment}
\keyword{assignment}
\keyword{optim}
\keyword{mse}
