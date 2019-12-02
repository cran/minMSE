\name{count_occurrences}
\alias{count_occurrences}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Count of Equal Treatment Vectors
}
\description{
Checks if the treatment vector given as argument already exists in the dataframe, i.e., has been produced by one or more earlier iteration(s).
}
\usage{
count_occurrences(df_treatments, curr_treatment)
}

\arguments{
  \item{df_treatments}{dataframe containing all the discovered treatment vectors.}
  \item{curr_treatment}{treatment vector to be investigated. }
}

\value{
Returns the number of treatment assignment vectors which are equal to the one being investigated.
}

\author{
Sebastian Schneider \email{sschneider@coll.mpg.de}; \email{sebastian@sebastianschneider.eu},
Giulia Baldini \email{giulia.baldini@uni-bonn.de}
}

\examples{
df_treatments <- data.frame(c(0, 2, 0, 1, 1, 0, 2, 1, 1, 0),
                            c(0, 2, 0, 1, 2, 0, 1, 2, 1, 0),
                            c(0, 2, 1, 1, 2, 0, 0, 2, 0, 0))
colnames(df_treatments) <- c("treatment_iter_1", "treatment_iter_2", "treatment_iter_3")

count_occurrences(df_treatments,
                  c(0, 2, 1, 1, 2, 0, 0, 2, 0, 0))
}

\keyword{treatment}
\keyword{assignment}
\keyword{optim}
\keyword{mse}
