# MetaME: Metabolic Maximum Entropy modelling 
**M**aximum **E**ntropy modelling of **Meta**bolic flux distributions using Expectation Propagation

## Package features
+ Inference of metabolic flux distributions from a stoichiometric model and a set of experimental flux profiles;
+ Plot the marginal densities of non-measured fluxes.


## Usage
The main script `MetaMEP.m` in `src` takes as input (in this order):
- `S` and `b`: the stoichiometric matrix and the vector constraining the fluxes $\boldsymbol{\nu}$
$$\begin{align}
\mathbf{S}\boldsymbol{\nu} = \boldsymbol{b}
\end{align}$$

- `lb` and `ub` the lower and upper bounds of the range of variability of the fluxes
$$\begin{align}
lb\left(i\right) < \nu_{i} < ub\left(i\right)
\end{align}$$

- `av_exp`: the experimental means of the measured fluxes;
- `va_exp`: the experimental variances (i.e. the square of the experimental errors) of the measured fluxes;
- `idx_exp`: the indices of the measured fluxes as they appear in the columns of the stoichiometric matrix;
- `type` specifies how to constrain the experimental variances; choose between 
    - `unique` (recommended) infers a unique Lagrange multiplier which ensures the fit of the average experimental variance;
    - `fixed` uses a fixed value for the Lagrange multiplier $\gamma$  (no fit of the experimental variances is performed);
    - `adapt` infers as many Lagrange multipliers as the measured fluxes.
- `gamma0`: initial value(s) for the $\gamma$ Lagrange multiplier(s) (ex. $10^5$).

EP parameters: 
- `precision`: precision required to the EP update equations to converge (ex. $10^{-4}$);
- `precision` $c$: precision required to the fitting of the experimental means (ex. $10^{-5}$);
- `precision` $\gamma$: precision required to the fitting of the (average) experimental variances (ex. $10^{-5}$);
- `damp`: damping used within the update of the EP equations (ex. $0.7$);
- `max_iter`: the maxinum number of iterations (ex. $5\cdot 10^{4}$);
- `minvar` and `maxvar`: the minimum and maximum allowed value for the approximate variances, respectively (ex. $10^{-50}$ and $10^{50}$).

The routine returns:
+ $\mu$ and $\sigma$: the parameters of the truncated Gaussian densities approximating each $\nu_{i}$ flux marginal probability:
$$\begin{align}
P\left(\nu_{i}\right)=\frac{1}{Z}e^{-\frac{1}{2\sigma_{i}}\left(\nu_{i} - \mu_{i}\right)^{2}}\mathbb{I}\left[lb\left(i\right) < \nu_{i} < ub\left(i\right)\right]
\end{align}$$
+ `av` and `va`: the mean and the variance of each flux $\nu_{i}$, that is
$$\begin{align}
<\nu_{i}>_{P}\qquad<\nu_{i}^{2}>_{P}-<\nu_{i}>_{P}^{2}
\end{align}$$
+ `a` and `d`: the means and variances of the univariate Gaussian of the approximation;
+ `c`: the coefficient or Lagrange multipliers $c$ ensuring the fit of the experimental means;
+ `G`: the $\gamma$ Lagrange multipliers ensuring the fit of the (average) experimental variances;
+ $\boldsymbol{\nu}^{e}$: the set of auxiliary fluxes at convergence (their values are *exactly* the experimental means);
+ $\Sigma$: the covariance matrix of the fluxes $\boldsymbol{\nu}$;
+ H: the entropy of the multivariate Gaussian approximating the joint flux distribution.

## Example

The `data` folder contains two files, `nanchen_2006.mat` and `schuetz_2012.mat`, containing the datasets used in [arXiv:2104.02594](https://arxiv.org/abs/2104.02594), i.e. the metabolic models and the experimental data. 

The main script `Run_inference.m` launches `MetaMEP.m` on the two datasets to get 
+ the Lagrange multipliers $c$ and $\gamma$;
+ the parameters of the approximate flux densities.
As an example, we report a function to plot the marginal probabilities of non-measured fluxes.

While preparing the input data one should note that:
+ Any fixed flux should be removed from the columns of the stoichiometric matrix. For this data, the glucose uptake has been fixed to the experimental mean (the values are reported in the `.mat` files) and, therefore, the associated flux has been removed from $\mathbf{S}$ and encoded into the constant term $\boldsymbol{b}$;
+ The lower and upper bounds of the fluxes should be those obtained from Flux Variability Analysis (those in `nanchen_2006.mat` and `schuetz_2012.mat` have been already pre-processed).




## Reference
This code is associated with the work in:
+ AP Muntoni, A Braunstein, A Pagnani, D De Martino, and A De Martino. *Relationship between fitness and heterogeneity in exponentially growing microbial populations.* [arXiv:2104.02594](https://arxiv.org/abs/2104.02594), April 2021

More about the use of Expectation Propagation for metabolic modelling in:
+ A Braunstein, A Muntoni, A Pagnani. *An analytic approximation of the feasible space of metabolic networks.* Nature Communications 8, Article number: 14915 (2017) [doi:10.1038/ncomms14915](https://www.nature.com/articles/ncomms14915)
+ https://github.com/anna-pa-m/Metabolic-EP