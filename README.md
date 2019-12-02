# Min Mean Squared Error Treatment Assignment for One or Multiple Treatment Groups

## Description
Performs treatment assignment for (field) experiments considering available, possibly multivariate and continuous, information (covariates, observable characteristics), that is: forms balanced treatment groups, according to the minMSE-method as proposed by Schneider and Schlather (2017).

## Usage
```r
assignMinMSETreatment(data, 
                      prev_treatment = NULL, 
                      n_treatments = 1, 
                      iterations = 50, 
                      change = 3, 
                      cooling = 1, 
                      t0 = 10, 
                      tmax = 10, 
                      built_in = 0, 
                      desired_test_vectors = 100, 
                      percentage_equal_treatments = 1, 
                      plot = 1, 
                      trace_output = 1, 
                      filename = NULL)
```
The only mandatory argument is `data`. All other arguments have a default, as shown above.

## Arguments
`data`&nbsp;&nbsp;&nbsp;&nbsp;a dataframe containing the covariate vectors for each attribute. The values might be missing or on different scales as the software deals with missing values and scaling automatically.

`prev_treatment`&nbsp;&nbsp;&nbsp;&nbsp;takes a numerical vector of partial treatment assignment as argument, and assigns the missing units (where the value is NA) to a treatment group while minimizing the objective function. Non-missing values are copied to the new vector, i.e., treatment group assignment of these observations is unaffected, but taken into consideration for achieving balanced treatment groups.

`n_treatments`&nbsp;&nbsp;&nbsp;&nbsp;specifies the number of treatment groups desired (in addition to the control group); minimum and default value is `n_treatments = 1`.

`iterations`&nbsp;&nbsp;&nbsp;&nbsp;specifies the number of iterations the algorithm performs; the default value is `iterations = 50`. Depending on the number of units and the number of covariates to consider for group assignment, a high value could result in a long run-time.

`change`&nbsp;&nbsp;&nbsp;&nbsp;sets the number of units to exchange treatment in each iteration; the default value is `change = 3`. In case of big datasets (e.g. with more than 100 units), one might consider increasing the default value.

`cooling`&nbsp;&nbsp;&nbsp;&nbsp;specifies the cooling scheme for the simulated annealing algorithm to use. `cooling = 1`, which is the default scheme, sets the temperature to
$`t0 / \log(floor((k − 1)/tmax ) \cdot tmax + \exp(1))`$,
whereas `cooling = 2` sets the temperature to the faster decreasing sequence
$`t0 /(floor((k − 1)/tmax) \cdot tmax + 1)`$.

In praxis, cooling schemes are mostly of one of these forms. One might want to change the cooling scheme if the plot indicates a too slow decrease of objective values. For a theoretical discussion of cooling schemes [Belisle (see 1992, p. 890)](https://www.jstor.org/stable/3214721?seq=1#page_scan_tab_contents).

`t0`&nbsp;&nbsp;&nbsp;&nbsp;sets the starting temperature for the simulated annealing algorithm, see [Belisle (1992)](https://www.jstor.org/stable/3214721?seq=1#page_scan_tab_contents) for theoretical convergence considerations. In praxis, a lower starting temperature `t0` decreases the acceptance rate of a worse solution more rapidly. Specifying a negative number allows values proportional to the objective function, i.e. `t0 = -5` sets the starting temperature to 1/5 of the objective function for the starting point, and thus - for the first `tmax` iterations of the algorithm - the difference of the old and the proposed solution is scaled by 1/5. When changing the default value, it should be considered that also worse solutions have to be accepted in order for the algorithm to escape a local minimum, so it should be chosen high enough. The default value is `t0 = 10`.

`tmax`&nbsp;&nbsp;&nbsp;&nbsp;specifies the number of function evaluations at each temperature: For instance, `tmax = 10` makes the algorithm evaluate 10 treatment assignments that are found based on the current solution, before the temperature is decreased and thus the probability of accepting a worse solution is decreased. The default value is `tmax = 10`.

`built_in`&nbsp;&nbsp;&nbsp;&nbsp;if `built_in = 1` the R built-in function `optim` with method 'SANN' (Simulated ANNealing) will be used to optimize the function. Otherwise, if `built_in = 0`, our implementation of the simulated annealing will be used. The function `built_in = 0` uses our first cooling function and this cannot be changed. To used the second cooling function, set `built_in = 0`. All the other parameters, such as `iterations`, `change`, `t0`, `tmax` are taken into account.

`desired_test_vectors`&nbsp;&nbsp;&nbsp;&nbsp;specifies the number of treatment assignment vectors that will be produced to perform Fischer's exact test (sometimes also called permutation test) for assessment of significance of the treatment effect (this, of course, will be done after treatment has been conducted and measurement of the outcome of interest has occurred). The number of possible treatment vectors will not exceed this number. The default value is `desired_test_vectors = 100`. For small datasets, one might consider increasing it without affecting performance. Note that this will affect your significance level: If `desired_test_vectors = 100` and all of them are unique (see `percentage_equal_treatments` below), you can achieve a significance level of at most p < 0.01\%.

If desired_test_vectors = 1, then the program returns a single vector that can be used for treatment assignment. Note that Fischer's exact test might still be possible and that alternative treatment vectors might also be produced after treatment has been conducted; yet, it is not sure how many *different* vectors can be produced with a given number of iterations. It is, therefore, good practice to produce the desired number of vectors with treatment assignment. For testing purposes, however, one might want to produce just one vector.

`percentage_equal_treatments`&nbsp;&nbsp;&nbsp;&nbsp;the percentage of non-unique treatment vectors that we allow. The default value is `percentage_equal_treatments = 1`. Note that this will affect your significance level: If `desired_test_vectors = 100` and `percentage_equal_treatments = 1`, you can achieve a significance level of at most p < 0.01\%..

`plot`&nbsp;&nbsp;&nbsp;&nbsp;can be used to suppress drawing a plot showing the value of the objective function for the a percentage of the iterations by setting `plot = 0`. The default setting is `plot = 1`, which shows a plot.

`trace_output`&nbsp;&nbsp;&nbsp;&nbsp;`trace_output = 1` prints helpful output such as the current iteration. To avoid the program output to be too cumbersome, a more detailed output is saved in a txt file called `program_output.txt`.

`filename`&nbsp;&nbsp;&nbsp;&nbsp;takes a string that represents the name of the `csv` file where the possible treatment assignments will be stored. If `filename = NULL`, then the file will not be saved.

## Value
The program returns a dataframe containing all the unique treatment assignment vectors generated by the program, where each treatment assignment vector has been optimized according to the score function, the MSE, as proposed by Schneider and Schlather (2017). It also outputs the number of unique vectors that were obtained.

## References
[Schneider and Schlather, 2017](https://www.econstor.eu/handle/10419/161931)

[Belisle, 1992](https://www.jstor.org/stable/3214721?seq=1#page_scan_tab_contents)

## See also
[optim](https://www.rdocumentation.org/packages/stats/versions/3.6.1/topics/optim), [ginv](https://www.rdocumentation.org/packages/MASS/versions/7.3-51.4/topics/ginv)

## Notes
With the default setting of plotting and using the trace output, the program writes to different files. To avoid this, set `plot = 0` and `trace_output = 0`. To allow plotting when using the built-in function `optim` we save the output of optim in an external file, we read its content and we plot the values using a regex. Since this is quite cumbersome, we suggest to use our method when wanting to plot convergence curves. 

## Examples
```r
input <- data.frame(c(10, 20, 30, 40, 130, 40, 120, 5, 10, 80),
                    c(2, 6, 2, 8, 1, 10, 9, 8, 7, 5),
                    c(1, 0, 2, 1, 0, 1, 0, 2, 1, 0))
colnames(input) <- c("IQ", "grade_maths", "both_parents")

assignMinMSETreatment(input,
                      prev_treatment = c(0, NA, NA, NA, 1, NA, NA, NA, NA, NA),
                      n_treatments = 2,
                      iterations = 100,
                      trace_output = 1,
                      built_in = 0,
                      desired_test_vectors = 100,
                      plot = 1,
                      filename = 'possible_treatments.csv')
```

Running this will perform an optimization of the treatment assignment over the `input` dataset. It will generate three different groups (0, 1, 2) where 0 is the control group. The existing treatment assignment will be taken into consideration and will not be changed during the iterations.

There will be 100 different treatment assignment vectors produced, and for each of them, the simulated annealing algorithm optimizing the score function, the MSE, will be run 100 times. There will be a plot and some feedback output will appear both on the screen and some will be saved in the output file. The optimization will be done with the implemented simulated annealing procedure.
