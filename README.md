# ivamse: IV-LASSO with approximate-MSE selection of the dictionary and penalty

## Overview

**ivamse** implements the approximate mean squared error (AMSE) criterion of
Ma, Navjeevan and Salahub, *Choosing the Dictionary and Penalty for IV-LASSO*.
We consider linear instrumental variables models with one endogenous regressor
in which the first stage is estimated using the LASSO. Estimation in these
models requires choosing both a dictionary of technical instruments and a
penalty level. The package selects these jointly to minimize a feasible AMSE
criterion for the coefficient on the endogenous regressor.

When candidate first stages consistently estimate the optimal instrument, the
corresponding structural estimators have the same first-order limiting
distribution. This limits the usefulness of first-order asymptotics for choosing
among them. The proposed criterion uses higher-order terms that account for
first-stage approximation error and many-instrument bias. For each candidate
`c`, consisting of a dictionary and a penalty level, we evaluate

```
S_c = ( sigma_eps^2 * E_n[Pi_hat_c^2]  +  sigma_epsv^2 * (d_c^2 + d_c) / n ) / h_c^2
```

where `Pi_hat_c` is the fitted first stage, `h_c = E_n[Pi_hat_c * x]` is the IV
denominator, `d_c` is the effective dimension of the LASSO fit, and `n` is the
effective sample size. After division by `h_c^2`, the first term measures
first-stage approximation error up to a common offset. The second captures
many-instrument bias as studied by Donald and Newey (2001). The weight on this
bias depends on the squared covariance `sigma_epsv^2` between the structural
and first-stage errors. In particular, greater endogeneity increases the cost
of a larger effective dimension in the criterion.

The package also implements selection by cross-validation and by the plug-in
penalty of Belloni, Chen, Chernozhukov and Hansen (2012). Cross-validation
minimizes first-stage prediction error, and the plug-in penalty is calibrated
to dominate the first-stage score. These procedures do not use the outcome and
therefore do not account for the estimated endogeneity of the regressor.

The main function, `ivamse()`, accepts a model formula. The matrix interface is
provided by `ivamse_fit()`. Both return an `"ivamse"` object with methods such
as `summary()`, `coef()`, `confint()` and `predict()`. The `selected()` method
reports the chosen dictionary, penalty and first-stage coefficients. Methods
are also provided for
[**sandwich**](https://CRAN.R-project.org/package=sandwich) covariance estimators
and [**broom**](https://CRAN.R-project.org/package=broom) tidiers. These methods
are described below.

## Installation

Install from GitHub with
[**remotes**](https://CRAN.R-project.org/package=remotes):

```r
# install.packages("remotes")
remotes::install_github("mnavjeev/ivamse")
```

It requires R (>= 3.6.0) and imports
[**Formula**](https://CRAN.R-project.org/package=Formula) and
[**glmnet**](https://CRAN.R-project.org/package=glmnet). The **sandwich**,
[**lmtest**](https://CRAN.R-project.org/package=lmtest), **broom** and
[**ivreg**](https://CRAN.R-project.org/package=ivreg) packages are suggested
rather than required. The illustration below takes its data from **ivreg** and
its robust covariances from **sandwich** and **lmtest**.

```r
library("ivamse")
library("sandwich")
library("lmtest")
data("SchoolingReturns", package = "ivreg")
```

## IV-LASSO with an approximate-MSE criterion

The model is

```
y_i = beta * x_i + w_i'gamma + eps_i,        x_i = Pi_i + v_i,
```

with one endogenous regressor `x`, a vector of included exogenous controls `w`,
and a dictionary of technical instruments `Z` used in the first stage. Before
fitting, the package partials the controls out of the outcome, the endogenous
regressor and every dictionary. The effective sample size is `n = N - q`, where
`q` is the rank of the controls. Every average written `E_n[.]` divides by `n`.
Controls are never penalized and are the same across dictionaries.

After residualization, dictionary columns are rescaled to unit empirical second
moment, as required by Assumption 2(i) of Ma et al. Columns that become zero
after projecting out the controls are dropped with a message.

For each candidate the first stage solves

```
minimize  (1/(2n)) * ||x - Z pi||^2  +  lambda * ||pi||_1
```

and the resulting fitted values serve as the single excluded instrument in a
just-identified IV regression of `y` on `(x, w)`. The first stage uses the LASSO
coefficients without a least-squares refit on the selected support.

The effective dimension `d_c` is the rank of the dictionary columns in the
equicorrelation set, defined by the columns whose score against the LASSO
residual attains the penalty. With collinear columns, the LASSO coefficients
need not be unique, so the number of nonzero coefficients can depend on the
solver. The equicorrelation set depends only on the fitted values, which are
unique. When the dictionary has full column rank, the solution is unique and
the effective dimension is generically equal to the number of selected columns.
The package uses this equality to simplify the calculation in that case.

### Model specification

The formula syntax follows `ivreg::ivreg()`. In the two-part form, list the
regressors before `|` and the instruments after it:

```r
ivamse(y ~ x + w1 + w2 | z1 + z2 + z3 + w1 + w2, data = d)
```

Exogenous regressors appear in both parts as their own instruments. Exactly one
regressor must be absent from the instrument list; it is treated as endogenous.
Listing it in both parts would make it exogenous.

The three-part form specifies the controls, endogenous regressor and excluded
instruments separately:

```r
ivamse(y ~ w1 + w2 | x | z1 + z2 + z3, data = d)
```

The two calls fit the same model. The three-part form avoids repeating the
controls and is required when using the `dictionaries` argument. An intercept
is included unless it is removed with `-1` in the first part.

All parts of the formula, together with `cluster` and `cv_group`, are assembled
into a single model frame, so `subset`, `na.action` and `contrasts` apply
uniformly to every variable in the fit. The arguments `x_name`, `intercept`,
`cluster`, `w` and `z` of `ivamse_fit()` are determined by the formula and cannot
be passed as additional matrix-interface arguments through `ivamse()`.
For direct specification of the design matrices, use `ivamse_fit()`.

### Dictionaries

The criterion selects over the dictionary as well as over the penalty. Candidate
dictionaries are supplied to `dictionaries` as a named list of one-sided
formulas, evaluated in `data`:

```r
ivamse(
  y ~ w1 + w2 | x | z1,
  dictionaries = list(
    small  = ~ z1 + z2 + z3,
    medium = ~ z1 + z2 + z3 + rain + snow,
    wide   = ~ (z1 + z2 + z3) * (rain + snow)
  ),
  data = d
)
```

With `dictionaries`, the third part of the model formula is a placeholder for
the excluded instruments; the named list supplies the dictionaries used for
estimation. The formula must still have three parts, as in the example above.

All dictionaries use the same controls, outcome and criterion, so their scores
are comparable. Candidates follow the order of the dictionary list, then
increasing penalty within each dictionary. Ties go to the earlier candidate,
so list simpler dictionaries first. For a matrix specification, pass a named
list of matrices to the `z` argument of `ivamse_fit()`.

The package rejects dictionary columns that are numerically collinear with the
endogenous regressor, since such a column would reproduce the regressor in the
first stage. This can occur when `.` is used on the instrument side of a formula.

## Illustration: Returns to schooling

We illustrate the procedure by estimating returns to schooling using the data
of Card (1995), an extract from the U.S. National Longitudinal Survey of Young
Men available as `SchoolingReturns` in **ivreg**. The outcome is log wage and
the endogenous regressor is years of `education`. We control for `ethnicity`,
residence in a metropolitan area (`smsa`), residence in the `south`, and a
quadratic in `age`.

We use a quadratic in age in place of the usual quadratic in labor-market experience.
Experience is defined as age minus education minus six and would introduce
additional endogenous regressors, which the criterion does not cover. As a
comparison, the ordinary least-squares estimate is

```r
m_ols <- lm(log(wage) ~ education + ethnicity + smsa + south + poly(age, 2),
            data = SchoolingReturns)
round(coef(summary(m_ols))["education", ], 4)
```
```
  Estimate Std. Error    t value   Pr(>|t|) 
    0.0341     0.0027    12.4898     0.0000 
```

To account for the possible endogeneity of education, we use indicators for
having grown up near a two-year or a four-year college as instruments. We first
take these two indicators as the dictionary and use the default penalty grid:

```r
m_iv <- ivamse(
  log(wage) ~ ethnicity + smsa + south + poly(age, 2) | education |
    nearcollege2 + nearcollege4,
  data = SchoolingReturns)
summary(m_iv)
```
```

Call:
ivamse(formula = log(wage) ~ ethnicity + smsa + south + poly(age, 
    2) | education | nearcollege2 + nearcollege4, data = SchoolingReturns)

Coefficients:
              Estimate Std. Error z value Pr(>|z|)    
education      0.12420    0.03840   3.235  0.00122 ** 
(Intercept)    4.60518    0.50466   9.125  < 2e-16 ***
ethnicityafam -0.05618    0.05998  -0.937  0.34900    
smsayes        0.07896    0.04076   1.937  0.05272 .  
southyes      -0.08377    0.02645  -3.168  0.00154 ** 
poly(age, 2)1  7.01046    0.45307  15.473  < 2e-16 ***
poly(age, 2)2 -0.08250    0.55565  -0.148  0.88196    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 0.4401 on 3003 degrees of freedom

First-stage selection by minimum AMSE criterion:
  dictionary:          instruments
  penalty:             0.01425 (kappa = 0.0986)
  effective dimension: 1
  IV denominator h:    0.04075
  first-stage F:       7.014 on the reference dictionary; 20.82 on the fitted
                       instrument, which is chosen using x and so is biased upward
  chosen from 24 eligible of 25 candidates

Error moments from pilot candidate 1 (dictionary "instruments"):
  sigma_eps^2: 0.1904,  sigma_epsv: -0.5539,  implied correlation: -0.505

3010 observations, 6 control columns, effective sample size 3004.
The selected penalty does not satisfy Assumption 2(iii) (implied alpha -0.842); 22 of 25 candidates are outside it.

```

The estimated coefficient on education is 0.12420, compared with 0.0341 under
least squares. The IV estimate therefore implies a larger return to an
additional year of schooling in this specification.

### The coefficient table and the selection block

The coefficient table reports the endogenous coefficient followed by the
controls. By default, standard errors assume homoscedasticity and inference uses
the normal distribution. Passing `df` to `summary()` specifies a t distribution
with the given degrees of freedom. The residual standard error uses `N` minus
the number of estimated coefficients. The criterion uses the effective sample
size after partialling out the controls: in this example, 3010 observations and
6 control columns give `n = 3004`.

The selection block reports the chosen dictionary and penalty. A single
dictionary is named `instruments` by default. The penalty is reported on the
scale of the objective and as a multiple `kappa` of the dictionary-specific
plug-in penalty of Belloni et al. This normalization permits comparison across
dictionaries of different widths and is reported under each grid specification.
Here `kappa = 0.0986`, so the selected penalty is substantially below the plug-in
level. The selected first stage has effective dimension one. Of the 25 grid
points, 24 are eligible; the remaining penalty sets all coefficients to zero
and does not define an IV estimator.

The two reported first-stage F statistics differ in how they use the data. The
first is the least-squares F statistic for `education` on the reference
dictionary after partialling out the controls. It does not depend on the
selected candidate. The second uses the fitted instrument, which was selected
using `x`. This statistic is biased upward and may be misleading as a measure
of identification strength, as indicated in the printed summary.

The pilot structural estimate supplies the error moments `sigma_eps^2` and
`sigma_epsv`. These moments determine the relative weight on approximation error
and many-instrument bias. The implied error correlation is -0.505 in this
example. When the estimated covariance is close to zero, the contribution of
the bias term is small.

### Diagnostics

`summary()` reports several diagnostics:

- **Penalty condition.** Assumption 2(iii) of Ma et al. requires the penalty to
  dominate the first-stage score by a margin that grows with the sample size.
  A candidate's `implied_alpha` is the largest constant it supports in that
  condition and is positive exactly when the condition holds. Candidates that
  fail it are still fitted and scored. The summary reports the selected
  candidate's constant and how many candidates fail the condition.
- **Effective dimension.** A dimension above `sqrt(n)` falls outside the regime
  covered by the consistency result.
- **IV denominator.** A denominator within a factor of two of `1/sqrt(n)`
  indicates a weak first stage. The reported standard errors are unreliable
  in this case.
- **Close scores.** A runner-up within 0.1% of the selected candidate indicates
  that the criterion offers little separation between them.
- **Fallback.** If the selected candidate gives no usable estimate, estimation
  uses the pilot. The reason is recorded and the original selected index is
  retained.

### Selecting over dictionaries

To illustrate joint selection of the dictionary and penalty, we consider three
dictionaries. The first, `proximity`, contains the two college-proximity
indicators. The second, `background`, adds residence in 1966 and family
circumstances at age 14. The third, `expanded`, includes interactions between
the two groups:

```r
dictionaries <- list(
  proximity  = ~ nearcollege2 + nearcollege4,
  background = ~ nearcollege2 + nearcollege4 + smsa66 + south66 + parents14 + library14,
  expanded   = ~ (nearcollege2 + nearcollege4) * (smsa66 + south66 + parents14 + library14)
)
m_dict <- ivamse(
  log(wage) ~ ethnicity + smsa + south + poly(age, 2) | education | nearcollege,
  dictionaries = dictionaries, pilot = "background", data = SchoolingReturns)
m_dict
```
```

Call:
ivamse(formula = log(wage) ~ ethnicity + smsa + south + poly(age, 
    2) | education | nearcollege, data = SchoolingReturns, dictionaries = dictionaries, 
    pilot = "background")

First stage: dictionary "expanded", lambda = 0.001464 (kappa = 0.0086), effective dimension 25

Coefficients:
    education    (Intercept)  ethnicityafam        smsayes       southyes  
      0.05337        5.53515       -0.15942        0.14708       -0.12053  
poly(age, 2)1  poly(age, 2)2  
      6.79958       -0.74011  

```

The criterion selects the interacted dictionary with a penalty of 0.001464,
giving an estimated return of 0.05337. The third formula part, `nearcollege`,
serves as a placeholder because `dictionaries` specifies the candidate
instruments. Setting `pilot = "background"` restricts the pilot to the second
dictionary. Its residuals supply the error moments used in the criterion.

The fitted object retains the scores for all candidates in `$candidates`,
ordered by dictionary and then by increasing penalty:

```r
head(m_dict$candidates[, c("dictionary", "lambda", "kappa", "h", "d",
                           "criterion", "implied_alpha")])
```
```
  dictionary       lambda        kappa          h d criterion implied_alpha
1  proximity 2.103026e-05 0.0001516927 0.04473814 3  3.296366    -0.8614800
2  proximity 3.086820e-05 0.0002226545 0.04473564 3  3.296371    -0.8614799
3  proximity 4.530833e-05 0.0003268121 0.04473199 3  3.296378    -0.8614798
4  proximity 6.650353e-05 0.0004796945 0.04472662 3  3.296390    -0.8614796
5  proximity 9.761383e-05 0.0007040952 0.04471874 3  3.296406    -0.8614791
6  proximity 1.432775e-04 0.0010334705 0.04470718 3  3.296432    -0.8614779
```

The `kappa` column expresses each penalty relative to the plug-in level. The
full table also reports the dictionary width `p`, the candidate's position
`column` within its dictionary, the fitted instrument's second moment
`fitted_moment`, and the number of nonzero coefficients `nonzero`. The fields
`empty`, `eligible` and `used` indicate whether the fit is zero, whether it can
be selected, and whether it was used for the reported estimate. Cross-validation
scores appear in `cvm` and `cvsd` when computed. Scores are retained for
ineligible candidates so that the criterion can be inspected over the full grid.

The `selected()` method reports the chosen candidate, the best candidate within
each dictionary, and the nonzero first-stage coefficients on the normalized
dictionary, ordered by decreasing absolute value. The `plot()` method shows the criterion
and effective dimension against `log(lambda)`, with a line for each dictionary
and a marker for the selected candidate:

```r
selected(m_dict)
plot(m_dict)
```

## The penalty grid

The default grid follows the construction used by **glmnet**, with
`nlambda = 25` geometrically spaced points between `lambda_max` and
`lambda.min.ratio * lambda_max`. Here `lambda_max` is the smallest penalty
that sets every coefficient to zero.
The ratio defaults to `1e-2` when the dictionary is wider than the effective
sample size and to `1e-4` otherwise. The default `nlambda` is smaller than
glmnet's 100 because the theory considers selection over a fixed, moderate
number of candidates. The grids of Ma et al. have 13 to 29 points per dictionary.
The package prints a message when the total exceeds 100 candidates.

| argument | effect |
| --- | --- |
| `penalty` | `"path"`, the default, for the geometric grid; `"bcch"` for fixed multiples of the plug-in penalty |
| `nlambda` | number of points per dictionary, used only when `penalty = "path"` |
| `lambda.min.ratio` | smallest penalty as a fraction of the largest, used only when `penalty = "path"` |
| `kappa` | the multipliers themselves, used only when `penalty = "bcch"`; the default is `2^seq(-4, 2, by = 0.5)` |
| `lambda` | penalties given directly, on the scale of the objective above |
| `c_lambda` | inflation constant in the plug-in penalty, `1.1` in Ma et al. and in Belloni et al. |
| `sigma_v` | first-stage error scale behind a `"bcch"` grid; `NULL` estimates it by least squares of `x` on the reference dictionary |

```r
ivamse(..., nlambda = 40, lambda.min.ratio = 1e-3)               # a finer path
ivamse(..., lambda = c(0.001, 0.005, 0.01, 0.05))                # penalties given directly
ivamse(..., penalty = "bcch", kappa = 2^seq(-4, 2, by = 0.5))    # the grid of Ma et al.
```

Under `penalty = "path"` the grid includes `lambda_max` itself, whose fit is
identically zero, so one candidate per dictionary is always empty and therefore
ineligible. Under `penalty = "bcch"` the grid is instead `kappa` times

```
lambda_plugin = c_lambda * sigma_v * qnorm(1 - gamma / (2p)) / sqrt(n),
        gamma = 0.1 / log(max(p, n)),
```

the homoscedastic plug-in penalty of Belloni et al. Thus `kappa = 1` gives the
plug-in level. By default, `sigma_v` is the degrees-of-freedom corrected
residual scale from least squares of `x` on the reference dictionary. A supplied
`sigma_v` changes the grid but does not enter the reported error correlation.

Supplying `lambda` overrides `penalty`, `nlambda`, `lambda.min.ratio` and
`kappa`. The package warns if you explicitly set an argument that has no effect.

For `lambda`, `nlambda`, `lambda.min.ratio` and `kappa`, a single specification
applies to all dictionaries. A named list gives each dictionary its own settings,
as in this example with grids of different lengths:

```r
m_bcch <- ivamse(
  log(wage) ~ ethnicity + smsa + south + poly(age, 2) | education | nearcollege,
  dictionaries = dictionaries, pilot = "background", penalty = "bcch",
  kappa = list(proximity  = 2^seq(-12, 2, by = 0.5),
               background = 2^seq(-12, 2, by = 0.5),
               expanded   = 2^seq(-4, 2, by = 0.5)),
  data = SchoolingReturns)
cat("candidates:", nrow(m_bcch$candidates), "\n")
table(m_bcch$candidates$dictionary)
```
```
candidates: 71 

background   expanded  proximity 
        29         13         29 
```

The list must contain an entry for every dictionary; missing names cause an error.

Penalties are on the scale of the objective written above, in which the loss is
divided by the effective sample size `n`. glmnet normalizes by the number of
rows, and the package handles that conversion internally. Use `bcch_lambda()`
to compute a plug-in level on the package's scale when building a grid by hand.

## Selection rules

The `select` argument specifies the selection rule.

`select = "amse"`, the default, minimizes the criterion. `select = "cv"` ranks
candidates by K-fold cross-validated prediction error of the first stage, with
`cv_s = "lambda.min"` taking the largest penalty attaining the minimum error and
`cv_s = "lambda.1se"` the largest penalty within one standard error of it, as in
`cv.glmnet()`. Folds are set by `nfolds`, which defaults to 10, or given
directly through `foldid`, and `cv_group` keeps rows sharing a value in the same
fold. The scores are stored in the `cvm` and `cvsd` columns of `$candidates` when
they are computed.

`select = "bcch"` chooses the candidate in the reference dictionary whose penalty
is closest on the log scale to the plug-in level, with
the first-stage scale refined iteratively on the LASSO residual as in the
original programs of Belloni et al. The `penalty = "bcch"` grid uses the
least-squares scale, so its plug-in level can differ from the selection target.

Use the full argument name `cv_s`: R would partially match `s` to `subset`.
The package catches values such as `s = "lambda.1se"` and reports the mistake.
It also warns if `cv_s` is supplied when `select` is not `"cv"`.

We compare the three selection rules using the same candidate list:

```r
rules <- lapply(c("amse", "cv", "bcch"), function(rule)
  ivamse(log(wage) ~ ethnicity + smsa + south + poly(age, 2) | education | nearcollege,
         dictionaries = dictionaries, pilot = "background", select = rule,
         data = SchoolingReturns))
data.frame(
  rule       = c("amse", "cv", "bcch"),
  dictionary = sapply(rules, function(f) f$selected$dictionary),
  dimension  = sapply(rules, function(f) f$selected$d),
  education  = round(sapply(rules, function(f) unname(coef(f)[1])), 4))
```
```
  rule dictionary dimension education
1 amse   expanded        25    0.0534
2   cv background         5    0.0579
3 bcch background         3    0.0556
```

In this example, AMSE selects the expanded dictionary with effective dimension
25, while CV and BCCH select the background dictionary with dimensions 5 and 3.
The resulting estimates of the schooling coefficient are similar. The larger
dimension selected by AMSE reflects the balance between approximation error and
bias: the improvement in first-stage approximation is sufficient to offset the
additional bias cost in the estimated criterion.

A candidate is eligible if its first stage is nonzero and, when `screen = TRUE`,
its IV denominator satisfies the threshold described below. AMSE and
cross-validation minimize over eligible candidates. The plug-in rule chooses
the nearest grid point in the reference dictionary without imposing eligibility.
If that candidate cannot be used for estimation, the package uses the pilot
and reports this in `summary()`.

### The pilot and the reference dictionary

The error moments `sigma_eps2` and `sigma_epsv` require a preliminary structural
estimate. With `pilot = "auto"`, the pilot is the candidate with the strongest
observed first stage in the reference dictionary. Fixing that dictionary in
advance keeps the choice of pilot candidate a function of `(x, Z)` alone and
allows a narrower dictionary to be used for this preliminary fit.

The moments are computed from the pilot's structural residual and `x`, without
a separate estimate of the first-stage error. The reported error correlation
uses the least-squares first-stage scale, even when `sigma_v` was supplied to
build the penalty grid.

The reference dictionary also supplies the default `sigma_v` for a `"bcch"`
grid and the candidates considered by `select = "bcch"`. It defaults to the
first dictionary in the list. A different dictionary can be specified by name,
as in `pilot = "background"`, and should be fixed before examining the outcome.
For replication, an integer value of `pilot` specifies a candidate index directly.

`summary()` reports which candidate served as pilot and the moments it produced.
The implied correlation is attenuated because it uses a preliminary estimate
of the structural coefficient. The paper quantifies the effect of this
attenuation on the criterion.

### Screening on the IV denominator

`screen = TRUE` restricts selection to candidates with `abs(h) >= 1/sqrt(n)`,
the root-`n` scale at which the IV denominator itself fluctuates. The default is
`FALSE`. The pilot is screened under either setting because a near-zero
denominator would affect the estimated error moments used for all candidates.
If no candidate in the reference dictionary satisfies the threshold, the pilot
is chosen from the unscreened set. This is recorded in `$pilot$screened`.

## Standard errors, clustering and cross-validation groups

After selection, the package computes the usual just-identified IV variance,
treating the fitted first stage as fixed. `vcov()` returns the clustered
variance when the fit was made with `cluster` and the homoscedastic variance
otherwise; `vcov(fit, cluster = NULL)` forces the homoscedastic form. A fit made
with `cluster` also has the clustering and the number of clusters reported in
its `summary()`.

The `cluster` argument accepts a one-sided formula, a bare column name, or a
vector. A formula or a name is looked up in the model frame, so it is subset and
filtered for missingness along with everything else. The same three forms are
accepted by `cv_group`, which is used only when `select = "cv"`. Clustering on
`age` below illustrates the interface rather than being a substantive choice.

```r
m_cl <- ivamse(
  log(wage) ~ ethnicity + smsa + south + poly(age, 2) | education | nearcollege,
  dictionaries = dictionaries, pilot = "background", cluster = ~ age,
  data = SchoolingReturns)
lmtest::coeftest(m_cl)[1:2, ]
```
```
              Estimate Std. Error  t value      Pr(>|t|)
education   0.05337377 0.01137054  4.69404  2.799127e-06
(Intercept) 5.53515221 0.14630974 37.83174 2.809290e-256
```
```r
confint(m_cl)["education", ]
```
```
     2.5 %     97.5 % 
0.03108792 0.07565963 
```

The point estimate remains 0.05337377 because clustering affects the reported
variance but does not change the selected candidate or first-stage fit.
The `coeftest()` function uses the fitted object's residual degrees of freedom
and reports a t statistic, while `summary()` defaults to normal inference.

The package provides `bread()`, `estfun()` and `vcovHC()` methods for covariance
estimation with **sandwich**. For example, HC1 standard errors are obtained by

```r
round(sqrt(diag(sandwich::vcovHC(m_dict, type = "HC1")))[1:2], 5)
```
```
  education (Intercept) 
    0.00866     0.11470 
```

The same applies to the cluster-robust estimators of **sandwich** and
`lmtest::coeftest()`, given a covariance function:

```r
sandwich::vcovCL(m_dict, cluster = SchoolingReturns$age)
lmtest::coeftest(m_dict, vcov = sandwich::vcovHC, type = "HC1")
```

`summary()`, `confint()` and `predict()` each take a `vcov.` argument that is
either a matrix or a function of the fitted object, as in `ivreg`, and further
arguments in `...` are passed on to that function:

```r
summary(m_dict, vcov. = sandwich::vcovHC, type = "HC1")
confint(m_dict, vcov. = sandwich::vcovHC, type = "HC1")
```

`summary()` and `confint()` also take `df`: the default `Inf` uses the normal
distribution, which is why the summary above reports `z value`, and passing
`df = m_iv$df.residual` switches to the t distribution.

A `vcovHC()` method is supplied because **sandwich** recovers residuals by
dividing the estimating functions by the model matrix, and the estimating
functions here use the projected regressors rather than the regressors
themselves. `hatvalues()` returns the leverage from the projected regressors, so
HC2 and HC3 are also available. Residuals are retained as `N`-vectors in the
original observation coordinates. This matters for robust and clustered
covariances, which can depend on how the controls are partialled out even when
the homoscedastic covariance is unchanged.

Clustering and heteroscedasticity adjustments apply to the reported standard
errors. Selection still uses the homoscedastic criterion covered by the theory.

## The matrix interface

The matrix interface `ivamse_fit()` accepts the outcome `y`, endogenous
regressor `x`, controls `w`, and instruments `z`. The instruments may be supplied
as a single matrix or a named list of dictionary matrices. An intercept is added
to the controls unless `intercept = FALSE`, and `x_name` labels the endogenous
coefficient.

The grid, selection, pilot and screening arguments described above are accepted
by `ivamse_fit()` and passed through `...` when using `ivamse()`. A matrix fit
does not retain formula terms, levels or a model frame, so `predict()` is not
available for new data.

```r
fit <- ivamse_fit(y, x, z = list(small = z1, wide = z2), w = W, select = "amse")
```

## The fitted object

An `"ivamse"` object is a list with the following main components:

| component | contents |
| --- | --- |
| `coefficients` | the endogenous coefficient first, then the intercept and the controls |
| `residuals`, `fitted.values` | `N`-vectors in the original observation coordinates |
| `candidates` | one row per candidate, the columns described above |
| `selected` | the chosen candidate: `dictionary`, `lambda`, `kappa`, `d`, `h`, `fitted_moment`, `implied_alpha`, `criterion`, the index selected, the index actually used, and whether estimation fell back to the pilot together with the reason |
| `beta` | the structural estimate implied by each candidate |
| `instrument` | the fitted first stage used as the instrument, an `N`-vector |
| `first_stage`, `scale` | coefficients on the normalized dictionary, and the column scales that normalized it |
| `sigma_eps2`, `sigma_epsv`, `rho_hat` | the error moments from the pilot, and the implied correlation |
| `sigma_v` | the first-stage scale used to build a `"bcch"` grid |
| `first_stage_F`, `instrument_F` | the F on the reference dictionary and the F on the fitted instrument |
| `pilot` | the pilot's index, its dictionary, its estimate, and whether it was screened |
| `n`, `N`, `q` | effective sample size, observations, and the rank of the controls |
| `regressors`, `projected` | the regressor matrix and its projection on the instruments, which the variance estimators use |
| `cov.unscaled` | `(X_hat' X_hat)^{-1}` |

The fields `select`, `screen`, `cluster`, `endogenous`, `controls`, `nobs`,
`df.residual` and `foldid` record how the fit was made, and `ivamse()` adds
`call`, `formula`, `terms`, `levels`, `contrasts`, `na.action` and, when
`model = TRUE`, the model frame in `model`.

## Methods

`print()` gives the call, the selected dictionary, penalty and effective
dimension, and the coefficients. `summary()` adds the coefficient table, the
residual standard error, the selection block shown above, the best candidate in
each dictionary when there is more than one, the pilot's error moments, and the
diagnostic notes. `selected()` returns the selection alone, with the nonzero
first-stage coefficients.

`coef()` takes `complete`, which governs whether aliased regressors are kept.
`fitted()`, `residuals()`, `nobs()` and `sigma()` behave as for `lm` objects,
with `sigma()` the residual standard error of the structural equation.
`formula()` and `terms()` return the model formula and its terms.

`model.matrix()` takes a `component` argument: `"regressors"` returns `(x, w)`,
`"projected"` returns the projected regressors that the variance estimator uses,
and `"instrument"` returns the fitted first stage. `hatvalues()` computes
leverage from the projected regressors, which is what HC2 and HC3 need.

`predict()` supports `type = "response"` and `se.fit`. Prediction from new data
requires a model fitted with `ivamse()` rather than `ivamse_fit()`, since it
rebuilds the regressors from the stored terms.

`plot()` draws the criterion and the effective dimension against `log(lambda)`,
one line per dictionary, with the selected candidate marked and the panel titles
taken from `main`.

Three **broom** methods are available:

- `tidy()` returns the coefficient table and accepts `conf.int` and `conf.level`.
- `glance()` reports `nobs`, `n.effective`, `n.controls`, `n.candidates`, `rule`,
  `dictionary`, `lambda`, `kappa`, `dimension`, `criterion`, `rho.hat` and `sigma`.
- `augment()` adds `.fitted` and `.resid` to the stored model frame, or `.fitted`
  to `newdata`. Using the stored model frame requires `model = TRUE` at fitting.

## Exported building blocks

Four functions are exported for separate use. The function `bcch_lambda()`
computes the plug-in penalty, and `implied_alpha()` computes the largest
score-domination constant supported by a penalty. The function
`effective_dimension()` computes the rank of the equicorrelation set, and
`feasible_criterion()` evaluates the AMSE score from the error moments and
candidate quantities supplied to it.

```r
bcch_lambda(p = 50, n = 500, sigma_v = 1)
implied_alpha(2^seq(-4, 2, by = 0.5) * bcch_lambda(50, 500), p = 50, n = 500)
```

## References

Ma, Y., Navjeevan, M., and Salahub, B. *Choosing the Dictionary and Penalty for
IV-LASSO.*

Belloni, A., Chen, D., Chernozhukov, V., and Hansen, C. (2012). Sparse models and
methods for optimal instruments with an application to eminent domain.
*Econometrica* **80**, 2369-2429.

Card, D. (1995). Using geographic variation in college proximity to estimate the
return to schooling. In L. N. Christofides, E. K. Grant and R. Swidinsky (eds.),
*Aspects of Labour Market Behaviour: Essays in Honour of John Vanderkamp*.
University of Toronto Press, Toronto.

Donald, S. G., and Newey, W. K. (2001). Choosing the number of instruments.
*Econometrica* **69**, 1161-1191.
