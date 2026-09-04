# ivamse: IV-LASSO with approximate-MSE selection of the dictionary and penalty

## Overview

The **ivamse** package estimates a linear instrumental-variables model whose
first stage is fitted by the LASSO, and chooses both the dictionary of technical
instruments and the penalty level by minimizing a feasible estimate of the
approximate mean squared error (AMSE) of the resulting estimator of the
structural coefficient. It implements the criterion of Ma, Navjeevan and
Salahub, *Choosing the Dictionary and Penalty for IV-LASSO*.

The two choices are not settled by first-order asymptotics. Any first stage that
is consistent for the optimal instrument leaves the structural estimator with the
same limiting distribution, so candidates cannot be separated by their
first-order behavior. The criterion separates them by the terms that the
first-order argument discards. Writing a candidate `c` for one pairing of a
dictionary with a penalty, it evaluates

```
S_c = ( sigma_eps^2 * E_n[Pi_hat_c^2]  +  sigma_epsv^2 * (d_c^2 + d_c) / n ) / h_c^2
```

where `Pi_hat_c` is the fitted first stage, `h_c = E_n[Pi_hat_c * x]` is the IV
denominator, `d_c` is the effective dimension of the LASSO fit, and `n` is the
effective sample size. The first term in the numerator is an approximation error
and rewards first-stage fit. The second is a many-instrument bias of the kind
studied by Donald and Newey (2001), weighted by the squared covariance
`sigma_epsv^2` of the structural and first-stage errors. Cross-validation ranks
candidates by the first term alone; the plug-in penalty of Belloni, Chen,
Chernozhukov and Hansen (2012) is constructed to dominate the first-stage score
and uses neither term. Because the weight on the bias rises with the endogeneity
of the regressor, the criterion regularizes a badly endogenous regressor more
heavily than a nearly exogenous one.

The main fitting function is `ivamse()`, a formula interface to the work-horse
`ivamse_fit()`, which takes design matrices. Both return an object of class
`"ivamse"`, for which the usual methods are available, including `print()`,
`summary()`, `coef()`, `vcov()`, `confint()`, `predict()`, `fitted()`,
`residuals()`, `nobs()`, `sigma()`, `hatvalues()`, `model.matrix()`,
`formula()`, `terms()` and `plot()`, together with `selected()`, which reports
what the selection rule chose. Additionally, `bread()`, `estfun()` and
`vcovHC()` methods make the variance estimators of the
[**sandwich**](https://CRAN.R-project.org/package=sandwich) package available
directly, and `tidy()`, `glance()` and `augment()` methods connect the fit to
[**broom**](https://CRAN.R-project.org/package=broom).

## Installation

The package is hosted on GitHub and can be installed with the
`install_github()` function of the
[**remotes**](https://CRAN.R-project.org/package=remotes) package:

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
and a dictionary of technical instruments `Z` that enters only through the first
stage. The included controls are partialled out of the outcome, the endogenous
regressor and every dictionary before anything else happens, so the effective
sample size is `n = N - q` with `q` the rank of the controls, and every average
written `E_n[.]` divides by `n` rather than by `N`. Controls are never
penalized, and the same controls apply to every dictionary. Dictionary columns
are rescaled to unit empirical second moment after residualization, which is the
normalization Assumption 2(i) of Ma et al. requires; columns that the projection
annihilates are dropped with a message.

For each candidate the first stage solves

```
minimize  (1/(2n)) * ||x - Z pi||^2  +  lambda * ||pi||_1
```

and the LASSO fitted values themselves, rather than a refit on the selected
support, serve as the single excluded instrument in a just-identified IV
regression of `y` on `(x, w)`.
The effective dimension `d_c` is the rank of the dictionary columns in the
equicorrelation set, the columns whose score against the LASSO residual attains
the penalty. Counting nonzero coefficients would not do, since with collinear
columns the coefficient vector is not unique and two solvers can report different
supports for the same fit; the equicorrelation set is a function of the fitted
values, which are unique. Note that when the dictionary has full column rank the
solution is unique and the effective dimension is generically the number of
selected columns, which the code uses as a shortcut.

### Model specification

Regressors and instruments are specified in a formula of the same two shapes
that `ivreg::ivreg()` accepts. The two-part form lists regressors and then
instruments,

```r
ivamse(y ~ x + w1 + w2 | z1 + z2 + z3 + w1 + w2, data = d)
```

where exogenous regressors appear on both sides and are their own instruments.
The regressors that do not also appear among the instruments are taken to be
endogenous, and exactly one of them is allowed. Note that repeating a variable in
both parts by accident therefore makes it exogenous silently. The three-part form
instead names the included exogenous variables, the endogenous regressor, and the
excluded instruments only,

```r
ivamse(y ~ w1 + w2 | x | z1 + z2 + z3, data = d)
```

The two calls fit the same model. The three-part form is the more convenient one
here, since a dictionary of technical instruments can run to hundreds of columns
and the two-part form requires writing all of them twice; it is also the form
required by the `dictionaries` argument. Both parts carry an implied intercept
unless it is removed with `-1` in the first part.

All parts of the formula, together with `cluster` and `cv_group`, are assembled
into a single model frame, so `subset`, `na.action` and `contrasts` apply
uniformly to every variable in the fit. The arguments `x_name`, `intercept`,
`cluster`, `w` and `z` of `ivamse_fit()` are determined by the formula and cannot
be set through `ivamse()`; use `ivamse_fit()` when direct control over the
matrices is wanted.

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

Supplying `dictionaries` requires the three-part formula, whose third part is
then ignored; it still has to be there, and something short such as a single
instrument is the usual filler. Every dictionary faces the same controls, the
same outcome and the same criterion, so candidates from different dictionaries
are directly comparable. Candidates are ordered by dictionary in the order they
are listed and, within a dictionary, by increasing penalty, and ties in the
criterion are broken toward the earlier candidate. Simpler dictionaries should
therefore be listed first. To supply dictionaries as matrices rather than as
formulas, pass a named list of matrices as the `z` argument of `ivamse_fit()`.

A dictionary column that is numerically collinear with the endogenous regressor
makes the first stage trivially perfect and the IV estimate meaningless, so it is
refused. This happens most easily by writing a `.` on the instrument side of a
formula.

## Illustration: Returns to schooling

The data are the extract from the U.S. National Longitudinal Survey of Young Men
used by Card (1995), supplied as `SchoolingReturns` by the **ivreg** package.
The outcome is the log wage, the endogenous regressor is years of `education`,
and the controls are `ethnicity`, residence in a metropolitan area (`smsa`),
residence in the `south`, and a quadratic in `age`. The quadratic in labor-market
experience of the usual wage equation is replaced by one in age, since experience
is age minus education minus six and would be a second endogenous regressor,
which the criterion does not cover. Estimating the wage equation by ordinary
least squares gives

```r
m_ols <- lm(log(wage) ~ education + ethnicity + smsa + south + poly(age, 2),
            data = SchoolingReturns)
round(coef(summary(m_ols))["education", ], 4)
```
```
  Estimate Std. Error    t value   Pr(>|t|) 
    0.0341     0.0027    12.4898     0.0000 
```

Education is plausibly endogenous, and Card's instruments are indicators for
having grown up near a two-year and near a four-year college. Taking those two
indicators as the dictionary and leaving everything else at its default:

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

The IV estimate of the return to a year of schooling, 0.12420, is well above the
least-squares estimate of 0.0341, which is the familiar finding in these data.

### The coefficient table and the selection block

The coefficient table is ordered with the endogenous coefficient first and the
controls after it. Standard errors are the homoscedastic ones by default and the
reference distribution is the normal, so the columns are labeled `z value` and
`Pr(>|z|)`; passing `df` to `summary()` switches to a t reference with that many
degrees of freedom. The residual standard error is reported on `N` minus the
number of estimated coefficients. Controls are partialled out rather than
penalized, at the cost of the degrees of freedom they use: with 3010 observations
and 6 control columns the criterion is evaluated at an effective sample size of
3004.

The block below the coefficient table describes the first stage that produced the
estimate. Only one dictionary was supplied, so it is named `instruments`; the
chosen penalty is reported both on the scale of the objective and as the multiple
`kappa` of the plug-in penalty of Belloni et al. for that dictionary, which makes
it comparable across dictionaries of different widths and is reported whichever
rule built the grid. Here it is well below one, so the criterion selected a
penalty lighter than the plug-in level. The effective dimension is one, so the
selected first stage is a single direction in the two instruments. Of the 25
candidates on the default grid, 24 are eligible: the grid includes the smallest
penalty at which every coefficient is zero, whose fit is empty and therefore
defines no IV estimator.

The summary reports two first-stage F statistics, and they answer different
questions. The first is the conventional least-squares F of `education` on the
reference dictionary, computed after the controls are partialled out; it does not
use the selected candidate and so is not a post-selection quantity. The second is
the F on the fitted instrument actually used, which was chosen with `x` in hand
and is biased upward, and it is printed with that qualification attached.

The error moments come from the pilot candidate, whose structural estimate
supplies `sigma_eps^2` and `sigma_epsv`. The implied correlation, here -0.505,
determines how much weight the criterion places on the many-instrument term; a
regressor whose endogeneity is mild leaves that term small and the criterion
close to a fit criterion.

### Diagnostics

`summary()` prints a line for each condition of the theory that the fit does not
meet. The selected penalty may fail Assumption 2(iii) of Ma et al., which
requires the penalty to dominate the first-stage score by a margin that grows in
the sample size; the `implied_alpha` of a candidate is the largest constant it
supports in that condition, positive exactly when the assumption holds.
Candidates outside it are still fitted and still scored, and both the implied
constant and the number of candidates in the same position are reported, as in
the illustration above. The effective dimension may exceed `sqrt(n)`, outside the
regime the consistency result covers. The IV denominator may lie within a factor
of two of the `1/sqrt(n)` scale at which it fluctuates, in which case the first
stage is weak and the reported standard errors should not be relied on. The
runner-up may be within a tenth of a percent of the selected candidate, which
says that the ranking is close and that little rides on it. Finally, if the
selected candidate gives no usable estimate, estimation falls back to the pilot
and the reason is recorded; the selected index is not overwritten.

### Selecting over dictionaries

These data carry further background variables, and the design of the paper puts
several dictionaries built from them in competition. Here `proximity` is the
pair of college-proximity indicators, `background` adds residence in 1966 and
family circumstances at age 14, and `expanded` interacts the two groups:

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

The criterion prefers the interacted dictionary at a light penalty, and the
estimated return falls to 0.05337. Note that the third part of the formula,
`nearcollege`, plays no role; it is displaced by `dictionaries` and is present
only because the three-part form requires a third part. Setting
`pilot = "background"` makes the middle dictionary the reference, so the pilot
candidate is chosen within it and the error moments that weight the bias term
come from there rather than from the first dictionary in the list.

Every candidate that was scored is kept in the fit, one per row of
`$candidates`, in the order dictionaries were listed and, within each, by
increasing penalty:

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

`kappa` expresses each penalty as a multiple of the plug-in level whichever rule
built the grid, so the plug-in penalty can always be located among the
candidates. The full table adds `p`, the number of columns in the dictionary;
`column`, the position within it; `fitted_moment`, the second moment of the
fitted instrument; `nonzero`, the naive dimension; `empty`, marking fits the
LASSO shrank to zero; `eligible`; `cvm` and `cvsd` when cross-validation was
computed; and `used`, marking the candidate the reported estimate came from.
Scores for ineligible candidates are kept in the table for inspection rather than
overwritten, so the criterion can be examined along the whole grid.

Two further displays summarize the same information. `selected()` returns the
chosen candidate, the best candidate within each dictionary, and the nonzero
first-stage coefficients on the normalized dictionary, largest in absolute value
first. `plot()` draws two panels in the style of `plot.cv.glmnet()`, the
criterion and the effective dimension against `log(lambda)`, with one line per
dictionary and the selected candidate marked:

```r
selected(m_dict)
plot(m_dict)
```

## The penalty grid

By default the grid for each dictionary is built as **glmnet** builds one:
`nlambda = 25` points running geometrically down from `lambda_max`, the smallest
penalty at which every coefficient is zero, to `lambda.min.ratio * lambda_max`.
The ratio defaults to `1e-2` when the dictionary is wider than the effective
sample size and to `1e-4` otherwise. The default `nlambda` is smaller than
glmnet's 100 because the theory ranks a fixed and moderate list of candidates
rather than tracing a path; the grids of Ma et al. have 13 to 29 points per
dictionary. Scoring more than 100 candidates in total draws a message saying so.

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

the homoscedastic plug-in penalty of Belloni et al., so a multiplier of
`kappa = 1` reproduces that penalty exactly. The scale `sigma_v` behind the
anchor defaults to the degrees-of-freedom corrected residual scale from least
squares of `x` on the reference dictionary, which uses `(x, Z)` alone; supplying
it shapes the grid and never the reported error correlation. Supplying `lambda`
makes `penalty`, `nlambda`, `lambda.min.ratio` and `kappa` inert, and an argument
that cannot have any effect given the others draws a warning rather than being
ignored silently.

Each of `lambda`, `nlambda`, `lambda.min.ratio` and `kappa` may be given once
for every dictionary or as a named list with one entry per dictionary, which is
how dictionaries of very different widths are given grids of different lengths:

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

Such a list must have an entry for every dictionary; a missing name is an error
rather than a silent default.

Penalties are on the scale of the objective written above, in which the loss is
divided by the effective sample size `n`. Note that glmnet instead normalizes by
the number of rows; the conversion is handled internally. The
exported `bcch_lambda()` computes the plug-in level on the same scale, so a grid
can be built by hand if wanted.

## Selection rules

`select` names the rule that picks the candidate used for estimation.

`select = "amse"`, the default, minimizes the criterion. `select = "cv"` ranks
candidates by K-fold cross-validated prediction error of the first stage, with
`cv_s = "lambda.min"` taking the largest penalty attaining the minimum error and
`cv_s = "lambda.1se"` the largest penalty within one standard error of it, as in
`cv.glmnet()`. Folds are set by `nfolds`, which defaults to 10, or given
directly through `foldid`, and `cv_group` keeps rows sharing a value in the same
fold. The scores are stored in the `cvm` and `cvsd` columns of `$candidates` when
they are computed. `select = "bcch"` takes the candidate in the reference
dictionary whose penalty is closest on the log scale to the plug-in level, with
the first-stage scale refined iteratively on the LASSO residual as in the
original programs of Belloni et al. Note that the `penalty = "bcch"` grid instead
anchors on the least-squares scale, so the two levels differ.

The cross-validation rule is named `cv_s` rather than `s` because `s` would
partially match the `subset` argument of `ivamse()`, so that `s = "lambda.1se"`
would silently subset the data instead; a value passed to `subset` that looks
like an `s` argument is caught and reported. Setting `cv_s` under a rule other
than `"cv"` draws a warning, as do the other arguments that cannot take effect
given the rest.

Both comparators look at the first stage alone and neither uses the outcome, so
neither responds to the endogeneity of the regressor. Running the three rules
over the same candidate list gives

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

The three rules disagree about the dictionary and about the complexity of the
first stage, and the estimated return to schooling differs accordingly. Note that
the criterion here selects the richest dictionary and the largest effective
dimension. The bias term does not push mechanically toward simplicity; it is
traded against first-stage approximation error, and which of the two dominates
depends on the estimated endogeneity of the regressor and on how much fit the
additional columns buy.

A candidate is eligible when its first stage is not empty, so that it defines an
IV estimator at all, and when `screen = TRUE` also when its denominator clears
the screen described below. The AMSE rule and cross-validation both minimize
over that set. The plug-in rule instead takes the nearest grid point in the
reference dictionary whether or not it is eligible, so it can land on a
candidate the other two exclude; estimation then falls back to the pilot and
`summary()` reports that it did.

### The pilot and the reference dictionary

The criterion needs the two error moments `sigma_eps2` and `sigma_epsv`, and
those need a preliminary structural estimate. That estimate comes from the pilot
candidate, which `pilot = "auto"` takes to be the candidate with the strongest
observed first stage within a reference dictionary. Restricting the pilot to one
dictionary keeps a wide interacted dictionary from contaminating the preliminary
estimate and makes the pilot a function of `(x, Z)` alone. The moments are
computed from the pilot's structural residual and from `x` itself, so no separate
estimate of the first-stage error is needed, and the reported error correlation
always uses the least-squares first-stage scale rather than any `sigma_v`
supplied to shape a grid.

The reference dictionary is the dictionary the pilot is drawn from. It also
supplies the default `sigma_v` behind a `"bcch"` grid, and `select = "bcch"`
searches within it. By default it is the first dictionary in the list, and
`pilot` may name a different one, as `pilot = "background"` does above. It should
be a dictionary fixed before the outcome is examined. Passing an integer to
`pilot` instead names a candidate index directly, which is useful for replication
but bypasses the rule.

`summary()` reports which candidate served as pilot and the moments it produced.
Note that the implied correlation is attenuated relative to the truth, since it
is computed from a preliminary estimate rather than from the true coefficient;
the paper quantifies the cost of that attenuation.

### Screening on the IV denominator

`screen = TRUE` restricts selection to candidates with `abs(h) >= 1/sqrt(n)`,
the root-`n` scale at which the IV denominator itself fluctuates. The default is
`FALSE`, so that the full candidate list is scored. The pilot is screened in
either case, since a pilot with a near-zero denominator would distort both error
moments and hence every candidate's score; if nothing in the reference dictionary
survives, the unscreened set is used and that fact is recorded in
`$pilot$screened`.

## Standard errors, clustering and cross-validation groups

Once a candidate is chosen the estimator is an ordinary just-identified IV
regression and its variance is the ordinary one. `vcov()` returns the clustered
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

The point estimate 0.05337377 is the one `m_dict` already reported, since
clustering enters the reported variance and nothing else: the same candidate is
selected and the same first stage is fitted. Note that `coeftest()` takes the
residual degrees of freedom off the fitted object and reports a t statistic,
where `summary()` defaults to the normal.

Because `bread()` and `estfun()` methods are provided, any variance estimator in
**sandwich** applies unchanged:

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
HC2 and HC3 work as well. Residuals are kept as `N`-vectors in the original
observation coordinates so that robust and clustered variances are available at
all: the homoscedastic variance is invariant to how the controls are partialled
out, but a robust one is not.

Note that clustering and heteroscedasticity affect the reported standard errors
only. Selection always evaluates the homoscedastic criterion, since that is what
is proved.

## The matrix interface

`ivamse_fit()` is the work-horse behind `ivamse()` and takes design matrices
directly: the outcome `y`, the endogenous regressor `x`, the dictionaries `z` as
one matrix or a named list of matrices, and the controls `w`. It adds an
intercept column to the controls unless `intercept = FALSE`, and labels the
endogenous coefficient with `x_name`. Every grid, selection, pilot and screening
argument described above is an argument of `ivamse_fit()` and reaches it through
the `...` of `ivamse()`. It does not parse a formula, so it returns no `terms`,
`levels` or model frame, and `predict()` on new data is unavailable.

```r
fit <- ivamse_fit(y, x, z = list(small = z1, wide = z2), w = W, select = "amse")
```

## The fitted object

An object of class `"ivamse"` is a list. The components most often wanted are

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

Additionally, `select`, `screen`, `cluster`, `endogenous`, `controls`, `nobs`,
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

The **broom** methods are `tidy()`, with `conf.int` and `conf.level`;
`glance()`, which reports the selection as well as the fit, since which candidate
was chosen is usually the quantity of interest, and which carries `nobs`,
`n.effective`, `n.controls`, `n.candidates`, `rule`, `dictionary`, `lambda`,
`kappa`, `dimension`, `criterion`, `rho.hat` and `sigma`; and `augment()`, which
adds `.fitted` and `.resid` to the model frame, or `.fitted` to `newdata`.
`augment()` on a stored model frame requires `model = TRUE` at fitting.

## Exported building blocks

Four functions from the criterion are exported for use on their own:
`bcch_lambda()` for the plug-in penalty level, `implied_alpha()` for the largest
score-domination constant a penalty supports, `effective_dimension()` for the
rank of the equicorrelation set of a LASSO fit, and `feasible_criterion()` for
the criterion itself given the pieces it is built from.

```r
bcch_lambda(p = 50, n = 500, sigma_v = 1)
implied_alpha(2^seq(-4, 2, by = 0.5) * bcch_lambda(50, 500), p = 50, n = 500)
```

## Scope of the criterion

The theory covers one endogenous regressor. A formula implying more than one is
refused.

The criterion is derived under homoscedastic Gaussian errors and is always
evaluated in that form. Simulations with Laplace errors leave the ranking of
candidates essentially unchanged, so departures from Gaussianity are not fatal,
but heteroscedastic selection is covered by no result. Clustering and
heteroscedasticity affect the reported standard errors only and never enter
selection, so a clustered fit selects the same candidate as an unclustered one on
the same data.

Standard errors treat the selected first stage as fixed. They reflect neither
first-stage estimation error nor the fact that the candidate was chosen on the
same data, and they are not robust to weak identification. This is first-order
valid, since the limiting distribution of the estimator does not involve the
first-stage estimate.

The level of the criterion is not an estimate of mean squared error. It carries
an offset that is common to every candidate and drops out of comparisons, so
only differences and rankings across candidates are meaningful, and a criterion
value on its own says nothing about the precision of the estimate.

## Verification

The package reproduces the authors' replication code to numerical precision.
Point the environment variable `IVAMSE_PAPER_REPO` at a checkout of the paper
repository and the test suite compares penalties, IV denominators, fitted second
moments, effective dimensions and criterion values candidate by candidate. The
comparison runs on a simulated design, since the data used in the paper is not
redistributable, and the test is skipped when the variable is unset.

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
