# README

## Aptitude-Achievement Paradox Simulation and Data Analysis Code

This repository contains the simulation and data analysis code used for the manuscript:

**The Aptitude-Achievement Paradox: A Theory of Rank Reversals under Selective Observation**

The code reproduces the numerical simulations and empirical Major League Baseball example reported in the main manuscript and Supplementary Materials.

## Contents

The zip file contains the following main files:

```text
simulation.R
luis-arraez-2023.R
raw_batter.csv
baseball_simu.RData
README.md
```

## File descriptions

### `simulation.R`

This file contains the simulation code for the theoretical and sensitivity analyses.

It includes three main components:

1. **Simulation for Theorem 3.4**

   This part studies the aptitude-achievement paradox probability under identical achievement distributions. It generates achievement values $Y_1$ and $Y_2$, maps them to latent aptitude values $X_1$ and $X_2$ using the Full House order-statistic transformation, and estimates

   $$
   P(Y_1 < Y_2,\ X_1 > X_2).
   $$

   The simulations consider both Normal achievement distributions,

   $$
   Y \sim \mathcal{N}(0,1),
   $$

   and Gamma achievement distributions,

   $$
   Y \sim \mathrm{Gamma}(2,1),
   $$

   to check robustness to skewness.

2. **Sensitivity analysis for Theorem 3.4**

   This part creates the four-panel sensitivity surfaces used in the manuscript and supplement. The four panels correspond to:

   - **Panel (A):** fixed $(n_1,n_2)$, varying $(N_1,N_2)$;
   - **Panel (B):** fixed $(N_1,N_2)$, varying $(n_1,n_2)$;
   - **Panel (C):** fixed $(n_1,N_2)$, varying $(n_2,N_1)$;
   - **Panel (D):** fixed $(n_2,N_1)$, varying $(n_1,N_2)$.

   The code runs the analysis for the first through fourth largest achievers.

3. **Simulation for Theorem 4.4**

   This part verifies the Gaussian representation in Theorem 4.4. It compares Monte Carlo estimates of the paradox probability with the theoretical probability computed using the multivariate normal representation:

   $$
   n_2 F_W(w) - \frac{N_2}{N_1+N_2}.
   $$

   The multivariate normal probabilities are computed using the `mvtnorm` package.

### `luis-arraez-2023.R`

This file contains the empirical Major League Baseball analysis used in the manuscript.

The analysis treats the 2023 season as the reference system and compares it with seasons from 1871 to 2022. The goal is to estimate paradox probabilities for Luis Arraez's 2023 batting average relative to historical batting environments.

The file performs the following steps:

1. Loads adjusted batter data from the `fullhouse` package and additional historical baseball data.
2. Constructs season-level full-time batter samples.
3. Computes season-specific estimates of:
   - $\mu_t$: mean batting average or adjusted batting average;
   - $\sigma_t$: standard deviation;
   - $n_t$: number of full-time batters;
   - $N_t$: estimated eligible talent-pool size.
4. Computes paradox probability estimates using the Gaussian representation from Theorem 4.4.
5. Produces plots comparing estimated paradox probabilities with empirical exceedance patterns for Luis Arraez's 2023 batting average.

### `raw_batter.csv`

This file contains raw batter-level Major League Baseball data used in the empirical analysis.

It is used to compute raw batting averages and empirical exceedance indicators for historical seasons.

### `baseball_simu.RData`

This file contains preprocessed baseball data objects used by `luis-arraez-2023.R`.

In particular, it contains objects required to construct the season-level batter data and talent-pool estimates used in the Full House Modeling analysis.

## Required R packages

The code uses the following R packages:

```r
library(tidyverse)
library(parallel)
library(doParallel)
library(Pareto)
library(ggplot2)
library(reshape2)
library(akima)
library(patchwork)
library(mvtnorm)
library(fullhouse)
```

Some scripts may also use:

```r
library(dplyr)
library(viridis)
```

Please install any missing packages before running the scripts.

For example:

```r
install.packages(c(
  "tidyverse",
  "doParallel",
  "Pareto",
  "ggplot2",
  "reshape2",
  "akima",
  "patchwork",
  "mvtnorm",
  "viridis"
))
```

The `fullhouse` package can be installed from GitHub if it is not already available:

```r
# install.packages("remotes")
remotes::install_github("DEck13/fullhouse")
```

## Reproducibility

The simulations use the random seed

```r
set.seed(13)
```

to improve reproducibility.

Because several simulations use parallel computation through `mclapply`, exact results may vary slightly across operating systems, R versions, and hardware. The qualitative conclusions are stable across repeated runs.

The default number of Monte Carlo repetitions in several sections is:

```r
size = 1000
nrep = 10
```

For more precise estimates, increase `size` or `nrep`, but this will increase computation time.

## Running the simulations

To reproduce the simulation results, run:

```r
source("simulation.R")
```

The script generates simulation results and saves figures such as:

```text
paradox_probability_1st_to_4th_achievers.png
paradox_probability_1st_to_4th_achievers_gamma.png
sensitivity_normal_rank_1.png
sensitivity_normal_rank_2.png
sensitivity_normal_rank_3.png
sensitivity_normal_rank_4.png
sensitivity_gamma_rank_1.png
sensitivity_gamma_rank_2.png
sensitivity_gamma_rank_3.png
sensitivity_gamma_rank_4.png
```

The exact list of generated files depends on which sections of the script are run.

## Running the MLB analysis

To reproduce the Major League Baseball analysis, place the following files in the working directory:

```text
luis-arraez-2023.R
raw_batter.csv
baseball_simu.RData
```

Then run:

```r
source("luis-arraez-2023.R")
```

The script computes season-level quantities and paradox probability estimates for Luis Arraez's 2023 season.

Please check the file paths inside `luis-arraez-2023.R`. Some paths may need to be changed depending on where the files are stored. For example, replace paths such as

```r
"~/research/full_house_model/baseball_simu.RData"
"~/website/era-adjustment-website/src/data/raw_batter.csv"
```

with local relative paths:

```r
"baseball_simu.RData"
"raw_batter.csv"
```

## Main simulation settings

### Ratio analysis for $N_2/N_1$

The ratio analysis fixes

$$
n_1 = 60, \qquad n_2 = 40,
$$

varies

$$
N_1 \in \{100,200,\ldots,2000\},
$$

and sets

$$
N_2 = cN_1,
\qquad
c \in \{0.1,0.2,\ldots,1.0\}.
$$

The analysis is performed for the first through fourth largest achievers.

For the largest achiever case, Corollary 3.4.1 implies that the paradox probability is positive only when

$$
\frac{n_1}{N_1} \leq \frac{n_2}{N_2}.
$$

Under the above design, this condition becomes

$$
c \leq \frac{n_2}{n_1} = \frac{2}{3}.
$$

Therefore, curves for $c > 2/3$ collapse to zero and may overlap with the horizontal axis in the plots.

### Sensitivity analysis for Theorem 3.4

The sensitivity analysis uses the following four designs:

1. Fixed $n_1=40$, $n_2=60$, varying $N_1,N_2$.
2. Fixed $N_1=N_2=1000$, varying $n_1,n_2$.
3. Fixed $n_1=40$, $N_2=1000$, varying $n_2,N_1$.
4. Fixed $n_2=60$, $N_1=1000$, varying $n_1,N_2$.

The latent aptitude distribution is taken to be Pareto with parameter

$$
\alpha = 1.16.
$$

The achievement distribution is either Normal or Gamma.

### Simulation for Theorem 4.4

The Theorem 4.4 simulation compares two quantities:

1. A Monte Carlo estimate of

   $$
   P(Y_1 < Y_2,\ X_1 > X_2),
   $$

   based on simulated achievement and latent aptitude values.

2. The theoretical Gaussian approximation

   $$
   n_2 F_W(w) - \frac{N_2}{N_1+N_2},
   $$

   computed using `mvtnorm::pmvnorm`.

The parameter grid includes combinations of:

```r
n1_vec     <- c(10, 20)
n2_vec     <- c(20, 40)
N1_vec     <- c(200, 400)
N2_vec     <- c(300, 600)
mu1_vec    <- c(-1, 0)
mu2_vec    <- c(0, 1)
sigma1_vec <- c(1, 1.5)
sigma2_vec <- c(2, 3)
```

Only configurations satisfying

$$
\frac{n_1}{N_1} \leq \frac{n_2}{N_2},
\qquad
\sigma_1 \leq \sigma_2,
\qquad
\mu_1 \leq \mu_2
$$

are retained.

## Notes on computation time

Some simulations can be computationally intensive because they use nested loops and repeated Monte Carlo sampling. Runtime depends on:

- the number of Monte Carlo repetitions;
- the number of grid points;
- the number of achiever ranks;
- the number of CPU cores available.

To reduce runtime, decrease `size`, `nrep`, or the number of grid points.

To increase precision, increase `size` or `nrep`.

## Output interpretation

The main output of each simulation is the estimated paradox probability:

$$
\widehat P
=
\frac{1}{M}
\sum_{m=1}^M
\mathbf{1}
\{Y_{1,m}<Y_{2,m},\ X_{1,m}>X_{2,m}\},
$$

where $M$ is the number of Monte Carlo repetitions.

A larger value indicates a higher probability that the individual from system 1 has larger latent aptitude but lower observed achievement than the corresponding individual from system 2.

## Contact

For questions about the code or data, please contact the Shen Yan: shenyan3@illinois.edu. 
