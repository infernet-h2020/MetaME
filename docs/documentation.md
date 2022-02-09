# MetaME: Metabolic Maximum Entropy modeling 
**M**aximum **E**ntropy modeling of **Meta**bolic flux distributions using Expectation Propagation

## Package features
+ Inference of metabolic flux distributions from a stoichiometric model and a set of experimental flux profiles;
+ Plot the marginal densities of non-measured fluxes.

## Overview
This package aims at determining the least possible constrained model (according to the Maximum Entropy principle) of metabolic fluxes under particular conditions. More precisely: 
+ fluxes satisfy a mass-balance relationship and they take value from a limited range of variability (i.e. they belong to the feasible space &Omega;);
+ their statistics mirror the experimental evidence encoded as observed means and variances of a sub-set of fluxes, i.e. the measured fluxes.

The resulting model is
<!-- $$\begin{align}
P\left( \boldsymbol{\nu}, \boldsymbol{\nu}^{e}\right) = \frac{1}{Z} \mathbb{I}\left[\boldsymbol{\nu} \in \Omega\right]e^{-\frac{\gamma}{2} \left(\boldsymbol{\nu} - \boldsymbol{\nu}^{e} \right)^2 + \boldsymbol{c}^{t}\boldsymbol{\nu}^{e} }
\end{align}$$ --> 

<div align="center"><img style="background: transparent;" src="https://render.githubusercontent.com/render/math?math=P%5Cleft(%20%5Cboldsymbol%7B%5Cnu%7D%2C%20%5Cboldsymbol%7B%5Cnu%7D%5E%7Be%7D%5Cright)%20%3D%20%5Cfrac%7B1%7D%7BZ%7D%20%5Cmathbb%7BI%7D%5Cleft%5B%5Cboldsymbol%7B%5Cnu%7D%20%5Cin%20%5COmega%5Cright%5De%5E%7B-%5Cfrac%7B%5Cgamma%7D%7B2%7D%20%5Cleft(%5Cboldsymbol%7B%5Cnu%7D%20-%20%5Cboldsymbol%7B%5Cnu%7D%5E%7Be%7D%20%5Cright)%5E2%20%2B%20%5Cboldsymbol%7Bc%7D%5E%7Bt%7D%5Cboldsymbol%7B%5Cnu%7D%5E%7Be%7D%20%7D"></div>

where 
+ &nu; is the set of fluxes in the feasible space;
+ &nu;<sup>e</sup> is a set of noisy fluxes distributed according to the experimental means and variances;
+ *c* and &gamma; are the Lagrange multipliers enforcing the constraints on &nu;<sup>e</sup>. 

Note that the two sets of fluxes are coupled by the presence of &gamma;;
Marginalizing over &nu;<sup>e</sup>, one gets
<!-- $$\begin{align}
P\left( \boldsymbol{\nu} \right) = \frac{1}{Z}\mathbb{I}\left[\boldsymbol{\nu} \in \Omega\right] e^{\boldsymbol{c}^{t}\boldsymbol{\nu}}
\end{align}$$ --> 

<div align="center"><img style="background: transparent;" src="https://render.githubusercontent.com/render/math?math=P%5Cleft(%20%5Cboldsymbol%7B%5Cnu%7D%20%5Cright)%20%3D%20%5Cfrac%7B1%7D%7BZ%7D%5Cmathbb%7BI%7D%5Cleft%5B%5Cboldsymbol%7B%5Cnu%7D%20%5Cin%20%5COmega%5Cright%5D%20e%5E%7B%5Cboldsymbol%7Bc%7D%5E%7Bt%7D%5Cboldsymbol%7B%5Cnu%7D%7D"></div>

The two distributions are hard to compute as the computation of the partition function *Z* is intractable. Therefore, we resort to an analytic approximation provided by Expectation Propagation which allows us to obtain a multivariate Gaussian density associated with the joint distribution of the fluxes, and a set of univariate truncated normal distributions approximating the marginal flux probabilities. We stress that within this framework it is also possible to estimate the unknown Lagrange multipliers *c* and &gamma; with no extra computation. See the works in the Reference section for a more detailed description of the model.

## Usage
The main script `MetaMEP.m` in `src` takes as input (in this order):
- `S` and `b`: the stoichiometric matrix and the vector constraining the fluxes &nu;

<!-- $$\begin{align}
\mathbf{S}\boldsymbol{v} = \boldsymbol{b}
\end{align}$$ --> 

<div align="center"><img style="background: transparent;" src="https://render.githubusercontent.com/render/math?math=%5Cmathbf%7BS%7D%5Cboldsymbol%7Bv%7D%20%3D%20%5Cboldsymbol%7Bb%7D"></div>

- `lb` and `ub` the lower and upper bounds of the range of variability of the fluxes

<!-- $$\begin{align}
lb\left(i\right) < v_{i} < ub\left(i\right)
\end{align}$$ --> 

<div align="center"><img style="background: transparent;" src="https://render.githubusercontent.com/render/math?math=lb%5Cleft(i%5Cright)%20%3C%20v_%7Bi%7D%20%3C%20ub%5Cleft(i%5Cright)"></div>

- `av_exp`: the experimental means of the measured fluxes;
- `va_exp`: the experimental variances (i.e. the square of the experimental errors) of the measured fluxes;
- `idx_exp`: the indices of the measured fluxes as they appear in the columns of the stoichiometric matrix;
- `type` specifies how to constrain the experimental variances; choose between 
    - `unique` (recommended) infers a unique Lagrange multiplier which ensures the fit of the average experimental variance;
    - `fixed` uses a fixed value for the Lagrange multiplier &gamma;  (no fit of the experimental variances is performed);
    - `adapt` infers as many Lagrange multipliers as the measured fluxes.
- `gamma0`: initial value(s) for the &gamma; Lagrange multiplier(s) (ex. 1e5).

EP parameters: 
- `precision`: precision required to the EP update equations to converge (ex. 1e-4);
- `precision` *c*: precision required to the fitting of the experimental means (ex. 1e-5);
- `precision` &gamma;: precision required to the fitting of the (average) experimental variances (ex. 1e-5);
- `damp`: damping used within the update of the EP equations (ex. 0.7);
- `max_iter`: the maxinum number of iterations (ex. 5e4);
- `minvar` and `maxvar`: the minimum and maximum allowed value for the approximate variances, respectively (ex. 1e-50 and 1e50).

The m-function returns:
+ &mu; and &sigma;: the parameters of the truncated Gaussian densities approximating each flux marginal probability:
<!-- $$\begin{align}
P\left(\nu_{i}\right)=\frac{1}{Z}e^{-\frac{1}{2\sigma_{i}}\left(\nu_{i} - \mu_{i}\right)^{2}}\mathbb{I}\left[lb\left(i\right) < \nu_{i} < ub\left(i\right)\right]
\end{align}$$ --> 

<div align="center"><img style="background: transparent;" src="https://render.githubusercontent.com/render/math?math=P%5Cleft(%5Cnu_%7Bi%7D%5Cright)%3D%5Cfrac%7B1%7D%7BZ%7De%5E%7B-%5Cfrac%7B1%7D%7B2%5Csigma_%7Bi%7D%7D%5Cleft(%5Cnu_%7Bi%7D%20-%20%5Cmu_%7Bi%7D%5Cright)%5E%7B2%7D%7D%5Cmathbb%7BI%7D%5Cleft%5Blb%5Cleft(i%5Cright)%20%3C%20%5Cnu_%7Bi%7D%20%3C%20ub%5Cleft(i%5Cright)%5Cright%5D"></div>

+ `av` and `va`: the mean and the variance of each flux &nu;<sub>i</sub>, that is
<!-- $$\begin{align}
<\nu_{i}>_{P}\qquad<\nu_{i}^{2}>_{P}-<\nu_{i}>_{P}^{2}
\end{align}$$ --> 

<div align="center"><img style="background: transparent;" src="https://render.githubusercontent.com/render/math?math=%3C%5Cnu_%7Bi%7D%3E_%7BP%7D%5Cqquad%3C%5Cnu_%7Bi%7D%5E%7B2%7D%3E_%7BP%7D-%3C%5Cnu_%7Bi%7D%3E_%7BP%7D%5E%7B2%7D"></div>

+ `a` and `d`: the means and variances of the univariate Gaussian of the approximation;
+ `c`: the coefficient or Lagrange multipliers **c** ensuring the fit of the experimental means;
+ `G`: the &gamma; Lagrange multipliers ensuring the fit of the (average) experimental variances;
+ &nu;<sup>e</sup> : the set of auxiliary fluxes at convergence (their values mirror *exactly* the experimental means);
+ &Sigma;: the covariance matrix of the fluxes &nu;;
+ H: the entropy of the multivariate Gaussian approximating the joint flux distribution.

### Example

The `data` folder contains two files, `nanchen_2006.mat` and `schuetz_2012.mat`, containing the datasets used in [arXiv:2104.02594](https://arxiv.org/abs/2104.02594), i.e. the metabolic models and the experimental data. The MATLAB structures contain the necessary input data and some extra information about the two models, i.e. the reaction and metabolite names. 

The main script `Run_inference.m` launches `MetaMEP.m` on the two datasets to get 
+ the Lagrange multipliers **c** and &gamma;;
+ the parameters of the approximate flux densities.
As an example, we report a function to plot the marginal probabilities of non-measured fluxes.

### Important pre-processing steps

While preparing the input data one should note that:
+ Any fixed flux should be removed from the columns of the stoichiometric matrix. For this data, the glucose uptake has been fixed to the experimental mean (the values are reported in the `.mat` files) and, therefore, the associated flux has been removed from **S** and encoded into the constant term **b**;
+ The lower and upper bounds of the fluxes should be those obtained from Flux Variability Analysis (those in `nanchen_2006.mat` and `schuetz_2012.mat` have been already pre-processed).


## Reference
All the mathematical and implementation details (as well as a deep analysis of the results obtained for these datasets) are present in:
+ AP Muntoni, A Braunstein, A Pagnani, D De Martino, and A De Martino. *Relationship between fitness and heterogeneity in exponentially growing microbial populations.* [arXiv:2104.02594](https://arxiv.org/abs/2104.02594), April 2021

Please, cite this work if you use (even partially) this code.

More about Expectation Propagation for metabolic modelling in:
+ A Braunstein, A Muntoni, A Pagnani. *An analytic approximation of the feasible space of metabolic networks.* Nature Communications 8, Article number: 14915 (2017) [doi:10.1038/ncomms14915](https://www.nature.com/articles/ncomms14915)
+ https://github.com/anna-pa-m/Metabolic-EP
