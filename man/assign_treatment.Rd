\name{assign_treatment}
\alias{assign_treatment}

\title{
Min MSE Treatment Assignment
}
\description{
Computes the treatment assignment vector according to available data (observable characteristics, covariate vectors) given about the units (individuals or clusters, such as schools, hospitals, ...).
Consider using the user-friendly wrapper function \link[minMSE]{assignMinMSETreatment}.
}
\usage{
assign_treatment(current_data,
                 prev_treatment = NULL,
                 n_treatments = 1,
                 iterations = 50,
                 change = 3,
                 cooling = 1,
                 t0 = 10,
                 tmax = 10,
                 built_in = 0,
                 plot = 0,
                 create_plot_file = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{current_data}{a dataframe containing the covariate vectors for each attribute. If the values are missing or on different scales, please use \link[minMSE]{assignMinMSETreatment}, which automatically scales the data.}
  \item{prev_treatment}{takes a numerical vector of partial treatment assignment as argument, and assigns the missing units (where the value is NA) to a treatment group while minimizing the objective function. Non-missing values are copied to the new vector, i.e., treatment group assignment of these observations is unaffected, but taken into consideration for achieving balanced treatment groups.}
  \item{n_treatments}{specifies the number of treatment groups desired (in addition to the control group); minimum and default value is n_treatments = 1.}
  \item{iterations}{specifies the number of iterations the algorithm performs; the default value is iterations = 50. Depending on the number of units and the number of covariates to consider for group assignment, a high value could result in a long run-time.}
  \item{change}{sets the number of units to exchange treatment in each iteration; the default value is change = 3. In case of big datasets (e.g., with more than 100 units), one might consider increasing the default value.}
  \item{cooling}{specifies the cooling scheme for the simulated annealing algorithm to use. cooling = 1, which is the default scheme, sets the temperature to
  \deqn{t0/log(floor((k - 1)/tmax ) * tmax + exp(1)),} whereas cooling = 2 sets the temperature to the faster decreasing sequence \deqn{t0 /(floor((k - 1)/tmax) * tmax + 1).} In praxis, cooling schemes are mostly of one of these forms. One might want to change the cooling scheme if the plot indicates a too slow decrease of objective values. For a theoretical discussion of cooling schemes, see Belisle (see 1992, p. 890).}
  \item{t0}{sets the starting temperature for the simulated annealing algorithm, see Belisle (1992) for theoretical convergence considerations. In praxis, a lower starting temperature t0 decreases the acceptance rate of a worse solution more rapidly. Specifying a negative number allows values proportional to the objective function, i.e. t0 = -5 sets the starting temperature to 1/5 of the objective function for the starting point, and thus - for the first tmax iterations of the algorithm - the difference of the old and the proposed solution is scaled by 1/5. When changing the default value, it should be considered that also worse solutions have to be accepted in order for the algorithm to escape a local minimum, so it should be chosen high enough. The default value is t0 = 10.}
  \item{tmax}{specifies the number of function evaluations at each temperature: For instance, tmax = 10 makes the algorithm evaluate 10 treatment assignments that are found based on the current solution, before the temperature is decreased and thus the probability of accepting a worse solution is decreased. The default value is tmax = 10.}
  \item{built_in}{if built_in = 1 the R built-in function \link{optim} with method 'SANN' (Simulated ANNealing) will be used to optimize the function. Otherwise, if built_in = 0, our implementation of the simulated annealing will be used. The function built_in = 0 uses our first cooling function and this cannot be changed. To use the second cooling function, set built_in = 0. All the other parameters, such as iterations, change, t0, tmax are taken into account.}
 \item{plot}{can be used to draw a plot showing the value of the objective function for the a percentage of the iterations by setting plot = 1. The default setting is plot = 0, which suppresses the plot.}
 \item{create_plot_file}{Used to overwrite the plot file, in case there already exists one. It should only be 1 (true) when this method is called without the wrapper \link[minMSE]{assignMinMSETreatment}. This method alone is not capable of plotting, but it will create an auxiliary file that contains the information for plotting. To include plotting, use \link[minMSE]{assignMinMSETreatment} with desired_test_vectors = 1.}

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
With the default setting of plotting and using the trace output, the program writes to different files. To avoid this, set plot = 0 and trace_output = 0. For the built-in function \link{optim}, the trace output is necessary for printing, because we pipe the output of the program to a file to obtain the intermediate values of the optimization function.
}

\seealso{
\code{\link[MASS]{ginv}}, \code{\link{optim}}
}

\examples{
input <- data.frame(c(10, 20, 30, 40, 130, 40, 120, 5, 10, 80),
                    c(2, 6, 2, 8, 1, 10, 9, 8, 7, 5),
                    c(1, 0, 2, 1, 0, 1, 0, 2, 1, 0))
colnames(input) <- c("IQ", "grade_maths", "both_parents")

assign_treatment(input,
                 prev_treatment = c(0, NA, NA, NA, 1, NA, NA, NA, NA, NA),
                 n_treatments = 2,
                 iterations = 100,
                 built_in = 0,
                 plot = 0)
}

\keyword{treatment}
\keyword{assignment}
\keyword{optim}
\keyword{mse}
