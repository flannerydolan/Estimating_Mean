# Code and data for "The Consequences of Using the Sample Mean on Highly Skewed and Heavy Tailed Data" by Dolan, Lamontagne, and Vogel, 2022. 

All data is located in the `datasets` folder. 

1. Plot the fraction of the sample mean resulting from dropping the data above a certain percentile for each dataset (`PercentileFigure.R`)

2. Determine goodness-of-fit of theoretical distributions for the datasets with L-moment ratio diagrams (`Lmomdiagrams.R`).

3. Check robustness of distributional hypotheses obtained from Step 2 using probability-probability plots. Correlation coefficients and plots are generated in `PPplots.R`.

4. Perform controlled Monte Carlo experiments evaluating estimators of the mean with lognormal (`LN2.R`), generalized Pareto (`GPD.R`), and power law (`powerlaw.R`) data. Plot the Mean Absolute Error (MAE), L-scale, and bias of the intended distribution's MLE, the sample mean, and the median of means at varying degrees of skew (`alldist_rmse_bias_fig.R`).

5. Test estimators on datasets from Steps 1-3 (`SI_Table.R`).

6. Test estimators on USGS sites and plot relative difference between the sample mean and median of means estimators.
      i. Download and process USGS data (`process_USGS.R`)
      ii. Plot relative percent difference against Lskew differentiating by distribution and record length (`USGS_sites_fig.R`)
