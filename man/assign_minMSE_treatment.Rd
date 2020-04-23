\name{assign_minMSE_treatment}
\alias{assign_minMSE_treatment}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{minMSE Treatment Assignment for One or Multiple Treatment Groups -- Wrapper}
\description{Calling \link{assign_treatment}, this user-friendly wrapper function computes a given number of treatment assignment vectors that will contain at most a certain percentage of duplicates, specified in percentage_equal_treatments. This is useful if non-parametric inference (randomization inference, sometimes called Fisher's exact test or permutation test) is desired (which is often advised), for analysis of significance of the treatment effect. The main function, assign_treatment, computes the treatment assignment vector according to available data (observable characteristics, covariate vectors) given about the units (individuals or clusters, such as schools, hospitals, ...)}
\usage{
assign_minMSE_treatment(data,
                        prev_treatment = NULL,
                        n_treatments = 1,
                        n_per_group = NULL,
                        mse_weights = NULL,
                        iterations = 50,
                        change = 3,
                        cooling = 1,
                        t0 = 10,
                        tmax = 10,
                        built_in = 0,
                        desired_test_vectors = 100,
                        percentage_equal_treatments = 1,
                        plot = 0,
                        trace_output = 1,
                        filename = NULL)
}

\arguments{
  \item{data}{a dataframe or a matrix containing the covariate vectors for each attribute. The values might be missing or on different scales as the software deals with missing values and scaling automatically.}
  \item{prev_treatment}{takes a numerical vector of partial treatment assignment as argument, and assigns the missing units (where the value is NA) to a treatment group while minimizing the objective function. Non-missing values are copied to the new vector, i.e., treatment group assignment of these observations is unaffected, but taken into consideration for achieving balanced treatment groups.}
  \item{n_treatments}{specifies the number of treatment groups desired (in addition to the control group); minimum and default value is n_treatments = 1.}
  \item{n_per_group}{specifies a vector containing uneven sizes for the treatment groups. Default value is NULL, which yields even sized groups. The sum of the elements in the vector should be equal to the total number of observations.}
  \item{mse_weights}{a vector containing the mse_weights for each treatment, or a matrix containing the mse_weights for treatments and outcomes and scaling factors.}
  \item{iterations}{specifies the number of iterations the algorithm performs; the default value is iterations = 50. Depending on the number of units and the number of covariates to consider for group assignment, a high value could result in a long run-time.}
  \item{change}{sets the number of units to exchange treatment in each iteration; the default value is change = 3. In case of big datasets (e.g. with more than 100 units), one might consider increasing the default value.}
  \item{cooling}{specifies the cooling scheme for the simulated annealing algorithm to use. cooling = 1, which is the default scheme, sets the temperature to
  \deqn{t0/log(floor((k - 1)/tmax ) * tmax + exp(1)),} whereas cooling = 2 sets the temperature to the faster decreasing sequence \deqn{t0 /(floor((k - 1)/tmax) * tmax + 1).} In praxis, cooling schemes are mostly of one of these forms. One might want to change the cooling scheme if the plot indicates a too slow decrease of objective values. For a theoretical discussion of cooling schemes Belisle (see 1992, p. 890).}
  \item{t0}{sets the starting temperature for the simulated annealing algorithm, see Belisle (1992) for theoretical convergence considerations. In praxis, a lower starting temperature t0 decreases the acceptance rate of a worse solution more rapidly. Specifying a negative number allows values proportional to the objective function, i.e. t0 = -5 sets the starting temperature to 1/5 of the objective function for the starting point, and thus - for the first tmax iterations of the algorithm - the difference of the old and the proposed solution is scaled by 1/5. When changing the default value, it should be considered that also worse solutions have to be accepted in order for the algorithm to escape a local minimum, so it should be chosen high enough. The default value is t0 = 10.}
  \item{tmax}{specifies the number of function evaluations at each temperature: For instance, tmax = 10 makes the algorithm evaluate 10 treatment assignments that are found based on the current solution, before the temperature is decreased and thus the probability of accepting a worse solution is decreased. The default value is tmax = 10.}
  \item{built_in}{if built_in = 1 the R built-in function \link{optim} with method 'SANN' (Simulated ANNealing) will be used to optimize the function. Otherwise, if built_in = 0, our implementation of the simulated annealing will be used. The function built_in = 0 uses our first cooling function and this cannot be changed. To used the second cooling function, set built_in = 0. All the other parameters, such as iterations, change, t0, tmax are taken into account.}
 \item{desired_test_vectors}{specifies the number of treatment assignment vectors that will be produced to perform Fischer's exact test (sometimes also called permutation test) for assessment of significance of the treatment effect (this, of course, will be done after treatment has been conducted and measurement of the outcome of interest has occurred). The number of possible treatment vectors will not exceed this number. The default value is desired_test_vectors = 100. For small datasets, one might consider increasing it without affecting performance. Note that this will affect your significance level: If desired_test_vectors = 100 and all of them are unique (see `percentage_equal_treatments` below), you can achieve a significance level of at most p < 0.01\%.If desired_test_vectors = 1, then the program returns a single vector that can be used for treatment assignment. Note that Fischer's exact test might still be possible and that alternative treatment vectors might also be produced after treatment has been conducted; yet, it is not sure how many *different* vectors can be produced with a given number of iterations. It is, therefore, good practice to produce the desired number of vectors with treatment assignment. For testing purposes, however, one might want to produce just one vector.}
 \item{percentage_equal_treatments}{the percentage of non-unique treatment vectors that we allow. The default value is percentage_equal_treatments = 1. Note that this will affect your significance level: If desired_test_vectors = 100 and percentage_equal_treatments = 1, you can achieve a significance level of at most p < 0.01\%.}
 \item{plot}{can be used to draw a plot showing the value of the objective function for the a percentage of the iterations by setting plot = 1. The default setting is plot = 0, which suppresses the plot.}
 \item{trace_output}{trace_output = 1 prints helpful output such as the current iteration. To avoid the program output to be too cumbersome, a more detailed output is saved in a txt file called program_output.txt.}
 \item{filename}{takes a string that represents the name of the csv file where the possible treatment assignments will be stored. If filename = NULL, then the file will not be saved.}
}

\value{
The program returns a dataframe containing all the unique treatments generated by the program. It also outputs the maximum number of iterations that were reached before finding a non-unique vector.
}
\references{
\href{https://www.econstor.eu/handle/10419/161931}{Schneider and Schlather (2017)},
\href{https://www.jstor.org/stable/3214721}{Belisle (1992)}
}
\author{
Sebastian Schneider \email{sschneider@coll.mpg.de}; \email{sebastian@sebastianschneider.eu},
Giulia Baldini \email{giulia.baldini@uni-bonn.de}
}
\note{
With the default setting of plotting and using the trace output, the program writes to different files. To avoid this, set plot = 0 and trace_output = 0. For the built-in function \link{optim}, the trace output is necessary for printing, because we pipe the output of the program to file to obtain the intermediate values of the optimization function.
}

\seealso{
\code{\link[MASS]{ginv}}, \code{\link{optim}}, \code{\link{assign_treatment}}
}

\examples{
input <- data.frame(c(10, 20, 30, 40, 130, 40, 120, 5, 10, 80),
                    c(2, 6, 2, 8, 1, 10, 9, 8, 7, 5),
                    c(1, 0, 2, 1, 0, 1, 0, 2, 1, 0))
colnames(input) <- c("IQ", "grade_maths", "both_parents")

assign_minMSE_treatment(input,
                        prev_treatment = c(0, NA, NA, NA, 1, NA, NA, NA, NA, NA),
                        n_treatments = 2,
                        mse_weights = c(1, 2),
                        iterations = 100,
                        trace_output = 1,
                        built_in = 0,
                        desired_test_vectors = 100,
                        plot = 0,
                        filename = NULL)
}

\keyword{treatment}
\keyword{assignment}
\keyword{optim}
\keyword{mse}
