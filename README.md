# TPDS-nonparametric
Nonparametric estimation of the welfare change under TPDS policy reform using Bihar and Rajasthan PoC data.



-------------
bihar_bybpl.m
-------------

Main code for Bihar. Read in price experiment data. Plot the observed cash take up, estimate and plot the CDF and the demand curve, and compute the welfare change. 

x: amount of cash offered in lieu of ration entitlement (deviation from the fiscal value of the ration)

y: 0 or 1 cash take up chocie.
 
Call "bSplineSieve.m" to Estimate the CDF using linear sieve with B-spline basis imposing monotonicity constraint. 

Call "CNS2015.m" to Find the 95% confidence intervals for the welfare change using the method proposed by Chernozhukov, Newey, and Santos (2015)

Note that bootstrapping cannot be used due to the presence of inequality constraints.


-------------
raj_pooled.m
-------------
Main code for Rajasthan. Read in  stated CEV data. Plot the observed cash take up, estimate and plot the CDF and the demand curve, and compute the welfare change. 

x: CEV

Estimate the CDF using kernel density estimation method

Bootstrap the 95% confidence intervals for the welfare change. 

