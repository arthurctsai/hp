# Hyperbolic Power transformation demonstration

This is a supplementary materials of the paper:

Arthur C. Tsai*, Michelle Liou*, Maria Simak, and Philip E. Cheng. On hyperbolic transformations to normality. Computational Statistics & Data Analysis, 115:250–266, 2017.

The Matlab scripts demonstrate the efficacy of hyperbolic power transformation in the four essential types of transformations: 
concave, convex, concave- to-convex, and convex-to-concave functions. Six families of probability distributions are 
used in the demonstration: lognormal, Gamma, beta, Laplacian, uniform, and bimodal. The closed-form solutions of 
both power and scale parameters by equi-percentile point estimation and maximum likelihood estimation are demonstrate 
by a demo program.

## The demo script:
-demo_hpt.m - A demo program of Hyperbolic power transformation with parameter estimation by (1) the method of 
              percentile and (2) Maximum Likelihood Estimation (MLE). The histogram of the transformed data and 
              estimated parameters as well as Jarque-Bera test are shown in demp_hpt.jpg.
## Usage: 
>>demo_hpt;


## Subroutines:
- hyperdistmop() - A program to implement equations (3)-(6) in the paper to estimate parameters 
                  by the method of equating percentile points in terms of the linear order statistic. 
Usage: 
 >>[alpha, betaminus, lambdaminus, betaplus, lambdaplus] = hyperdistmop(x); % x: original data to be transformed
 
- hyperdistfminsearch() - The parameter estimation by MLE.
 Usage: 
 >>[alpha, betaminus, lambdaminus, betaplus, lambdaplus] = hyperdistfminsearch(x, ...
   alpha0, betaminus0, lambdaminus0, betaplus0, lambdaplus0);
   

-pG() - A function to plot Gaussian curve.



------
When you use the code please cite:

Arthur C. Tsai*, Michelle Liou*, Maria Simak, and Philip E. Cheng. On hyperbolic transformations to normality. Computational Statistics & Data Analysis, 115:250–266, 2017.



